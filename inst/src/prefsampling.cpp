/* By Xiangjie Xue, updated on 01/05/2021. */
#include <TMB.hpp>
#include <numeric>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density; // this where the structure for GMRF is defined
  using namespace Eigen;  // probably for sparseness clacs
  // assume the right number of dimention of parameter and data is passed.
  DATA_MATRIX(yresp); // A matrix of response, each column should have the same distribution
  DATA_VECTOR(ypp) // A vector for point process response (poisson).
  DATA_SPARSE_MATRIX(lmat); // A sparse matrix mapping mesh points to the observations
  DATA_STRUCT(spde,spde_t); // this as structure of spde object is defined in r-inla hence using that namespace
  DATA_VECTOR(w); // weight of mesh
  DATA_IMATRIX(idx); 
  /*
    a matrix of size (no of resp) * (no of resp + 1). 
    The first column is for sharing random field with point process. 
    For the first column, the values are indicated as >0 if sharing random field with point process.
    For the i-th column (idx.col(i + 1)), the values are also indicated as >0 if sharing random field defined in i-th random field.
    If idx(i, i + 1) = 0, then the entire column is ignored.
  */ 
  DATA_MATRIX(strfixed);
  /*
    a matrix of fixed structural parameters. For some distributions, non-NaN values will overwrites the corresponding strparam (consider it fixed),
    and NaN values indicate that it should be estimated (there is no NA in C++). Non-NaN values can be set partially.
    The entire column shares one estimate if it needs estimating (and hence strparam is a vector of length no. of resp).
    Details:
    Normal distribution, this is the log of element-wise standard deviation. If non-NaN, it is taken as a fixed value 
      of the log of standard deviation.
    Poisson distribution, this is the effort. The default input passed from R should be 1s.
    Binomial distribution, this is the number of trials. The default input from R should be 1s, i.e., the bernoulli distribution.
    Gamma distribution, this is the log of scale. If non-NaN, this is taken as a fixed value of the log of the scale.
  */
  DATA_IVECTOR(methods);
  /*
    About methods:
    0 - normal distribution, strparam/strfixed as log_sigma.
    1 - poisson distribution, strparam is not referenced, strfixed as effort.
    2 - binomial distribution, strparam is not referenced, strfixed as number of trials.
    3 - gamma distribution, the implementation in TMB is shape-scale. strparam/strfixed as log_scale;
  */
  PARAMETER_VECTOR(betaresp); // intercept term of each response
  PARAMETER(betapp); // intercept of point process
  PARAMETER_MATRIX(beta); 
  /*
    A matrix of the same size as idx.
    This matrix contains the estimates for the sharing effects.
    Note that for the right no. of resp-by-no. of resp matrix, the diagonal terms WILL always be 0.
    This might be resulting from the features of TMB. Consider a manual fixed from R side if needed.
    All the values in this matrix is constrained to be a normal variable.

    NB: Tried to remove this constraint before, the estimates will be unbounded. The diagonal terms can not
      be constrained as a normal variable with very small sd, this will cause erroneous results as well.
  */ 
  PARAMETER_VECTOR(log_kappa); 
  PARAMETER_VECTOR(log_tau);
  /*
    log of kappas/taus for the random field.
    The lengths of these vectors are no. of resp + 1.
    The first element in both vectors are for the random field of the point process.
  */
  PARAMETER_VECTOR(strparam); // see strfixed.
  PARAMETER_MATRIX(x); 
  /*
    the random field/effect. First column is always the random field for point process.
    The number of columns is sum(diag(idx[, -1]) > 0) + 1, and must have the right size before passing.

    NB: Tried to use the matrix with the column (1 + no. of resp). This will give an error when using nlminb.
  */
  vector<Type> kappa = exp(log_kappa);
  vector<Type> tau = exp(log_tau);
  // spde part
  SparseMatrix<Type> Q = Q_spde(spde, kappa[0]) * tau[0]; // create the precision matrix from the spde model for the GMRF
  // Type nll = 0.0;
  vector<Type> tempx = x.col(0);
  Type nll = GMRF(Q)(tempx); // field the random effect is a GMRF with precision Q

  // point process
  vector<Type> lambdapp = exp(tempx + betapp) * w.cwiseEqual(0).select(vector<Type>::Ones(w.size()), w); // Eexp(eta);
  nll -= w.cwiseEqual(0).select(vector<Type>::Zero(ypp.size()), dpois(ypp, lambdapp, true)).sum(); // construct likelihood
  
  matrix<Type> lambdaresp(x.rows(), yresp.cols());
  // for responses interact with point process and set intercept
  for (int i = 0; i < yresp.cols(); ++i){
    lambdaresp.col(i).setConstant(betaresp[i]); // set intercept of responses
    // if ith response has sharing effect with point process
    if (idx(i, 0) > 0){
      vector<Type> temp = lambdaresp.col(i);
      temp += beta(i, 0) * tempx; 
      lambdaresp.col(i) = temp;
    }
  }
  // for other random fields
  int RFcount = 1;
  for (int i = 1; i < idx.cols(); ++i)
    if (idx(i - 1, i) > 0){
      // inspect if each of diag(idx[, -1]) is positive. If 0, the whole column is ignored.
      Q = Q_spde(spde, kappa[i]) * tau[i];
      tempx = x.col(RFcount);
      nll += GMRF(Q)(tempx); // field the random effect is a GMRF with precision Q
      for (int j = 0; j < yresp.cols(); ++j)
        if (idx(j, i) > 0){
          vector<Type> temp = lambdaresp.col(j);
          if (j == i - 1) {
            temp += tempx; // For i-th variable itself.
            /*
              An IMPORTANT NOTE, the corresponding term of beta WILL NOT be 1.
              Assigning beta(j, i)=1 will produce the same result, this issue might be explained by TMB evaluating its derivative.
              As a result, a manual fix is needed in R side.
            */
          }else{
            temp += beta(j, i) * tempx;  // for other variables.
          }
          lambdaresp.col(j) = temp;
        }
      RFcount++;
    }

  matrix<Type> matchresp = lmat * lambdaresp;
  for (int i = 0; i < yresp.cols(); ++i){
    // construct likelihood
    vector<Type> resp = yresp.col(i); // response 
    vector<Type> pred = matchresp.col(i); // predictor
    vector<Type> fixedstr = strfixed.col(i); // fixed structural parameters
    switch(methods[i]){
      case 0 : {
        // normal
        vector<Type> sigma = fixedstr.isNaN().select(vector<Type>::Constant(fixedstr.size(), strparam[i]), fixedstr); // replace NaN terms with estimates
        sigma = exp(sigma);
        nll -= sum(dnorm(resp, pred, sigma, true));
        break;
      }
      case 1 : {
        // poisson
        pred = exp(pred) * fixedstr; // incorporate efforts.
        nll -= sum(dpois(resp, pred, true));
        break;
      }
      case 2 : {
        // binomial
        nll -= sum(dbinom_robust(resp, fixedstr, pred, true)); // this function is linking density of binomial to predictor directly.
        break;
      }
      case 3 : {
        // gamma.
        vector<Type> scale = fixedstr.isNaN().select(vector<Type>::Constant(fixedstr.size(), strparam[i]), fixedstr); // replace NaN with estimates
        scale = exp(scale);
        pred = exp(pred) / scale;
        nll -= sum(dgamma(resp, pred, scale, true));
      }
    }
  }

  // ? constrain beta? (see PARAMETER_MATRIX(beta)).
  nll -= sum(dnorm(beta.vec(), Type(0.), Type(1. / sqrt(2.)), true));
  ADREPORT(betaresp);
  ADREPORT(betapp);
  ADREPORT(beta);
  ADREPORT(strparam);
  // ADREPORT parameters for RF.
  ADREPORT(tau);
  ADREPORT(kappa);
  ADREPORT(x);
  return nll;
}
