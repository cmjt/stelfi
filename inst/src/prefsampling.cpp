/* By Xiangjie Xue, updated on 01/05/2021. */
#include <TMB.hpp>
#include <numeric>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density; // this where the structure for GMRF is defined
  using namespace Eigen;  // probably for sparseness clacs
  // assume the right number of dimention of parameter and data is passed.
  DATA_MATRIX(yresp); // A matrix of marks, each column contain data from the same distribution
  //DATA_MATRIX(yresp_ind); // A matrix of marks with additional fields... could not work out how to index here
  DATA_VECTOR(ypp) // A vector for point process response (poisson).
  DATA_SPARSE_MATRIX(lmat); // A sparse matrix mapping mesh points to the observations
  DATA_STRUCT(spde,spde_t); // this as structure of spde object is defined in r-inla hence using that namespace
  DATA_VECTOR(w); // weights of mesh
  DATA_MATRIX(strfixed);
  /*
   A matrix of fixed structural parameters.
   Normal distribution: this is the log of element-wise standard deviation.
   Poisson distribution: not used
   Binomial distribution, this is the number of trials.
   The default input from R should be 1s, i.e., the Bernoulli distribution.
   Gamma distribution, this is the log of scale.
  */
  DATA_IVECTOR(methods);
  /*
    About methods:
    0 - normal distribution, strparam/strfixed as log_sigma.
    1 - poisson distribution, strparam is not referenced, strfixed as effort.
    2 - binomial distribution, strparam is not referenced, strfixed as number of trials.
    3 - gamma distribution, the implementation in TMB is shape-scale. strparam/strfixed as log_scale;
  */
  PARAMETER_VECTOR(betaresp); // intercept term of each mark (n_mark)
  PARAMETER_VECTOR(betapp); // intercept of point process
  //std::cout << betaresp << "betaresp\n";
  PARAMETER_VECTOR(beta_coefs_pp); // coefficients of shared field
  DATA_IVECTOR(mark_field); // indicator if mark should have its own field
  int n_mark = 0;
  for (int i = 0; i < mark_field.size(); ++i){
    if (mark_field(i) == 1){
      n_mark = n_mark + 1;
	  }
  }
  //std::cout << n_mark << "n_mark\n";
  DATA_MATRIX(designmat);  // the design matrix for the fixed effects of the point process.
  PARAMETER_VECTOR(log_kappa); //same length as number of fields (i.e., max shared + n marks)
  PARAMETER_VECTOR(log_tau);
  /*
    log of kappas/taus for the random field.
    The lengths of these vectors are no. of resp + 1 if all marks have own field.
    The first element in both vectors are for the random field of the point process.
  */
  //PARAMETER_VECTOR(strparam); // see strfixed.
  PARAMETER_MATRIX(x); 
  /*
    the random field/effect. First column is always the random field for point process.
    The number of columns is + 1 + n(mark__field == 1), and must have the right size before passing.
  */
  vector<Type> kappa = exp(log_kappa);
  vector<Type> tau = exp(log_tau);
  // spde part
  SparseMatrix<Type> Q = Q_spde(spde, kappa[0]) * pow(tau[0],2);
  // create the precision matrix from the spde model for the GMRF
  // Type nll = 0.0;
  vector<Type> tempx = x.col(0);
  Type nll = GMRF(Q)(tempx); // field the random effect is a GMRF with precision Q

  // point process
  vector<Type> lambdapp = exp(tempx + designmat*betapp) * w.cwiseEqual(0).select(vector<Type>::Ones(w.size()), w); // Eexp(eta);
  nll -= w.cwiseEqual(0).select(vector<Type>::Zero(ypp.size()), dpois(ypp, lambdapp, true)).sum(); // construct likelihood
  
  matrix<Type> lambdaresp(x.rows(), yresp.cols());
  // for responses interact with point process and set intercept
  int RFcount = 1;
  for (int i = 0; i < yresp.cols(); ++i){
    lambdaresp.col(i).setConstant(betaresp[i]); // set intercept of responses
    vector<Type> temp = lambdaresp.col(i);
    temp += beta_coefs_pp(i) * log(lambdapp); 
    lambdaresp.col(i) = temp;
    if(mark_field(i) == 1){
	     Q = Q_spde(spde, kappa[i + 1]) * pow(tau[i + 1],2);
	     tempx = x.col(RFcount);
	     nll += GMRF(Q)(tempx); // field the random effect is a GMRF with precision Q
	     //vector<Type> temp = lambdaresp.col(i);
	     temp += tempx; // For i-th variable itself.
	     lambdaresp.col(i) = temp;
       RFcount++;
    }
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
        //vector<Type> sigma = fixedstr.isNaN().select(vector<Type>::Constant(fixedstr.size(), strparam[i]), fixedstr); // replace NaN terms with estimates
        vector<Type> sigma = fixedstr;
        sigma = exp(sigma);
        nll -= sum(dnorm(resp, pred, sigma, true));
        break;
      }
      case 1 : {
        // poisson
        pred = exp(pred);
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
        //vector<Type> scale = fixedstr.isNaN().select(vector<Type>::Constant(fixedstr.size(), strparam[i]), fixedstr); // replace NaN with estimates
        vector<Type> scale = fixedstr;
        scale = exp(scale);
        pred = exp(pred) / scale;
        nll -= sum(dgamma(resp, pred, scale, true));
      }
    }
  }

  // ? constrain beta? (see PARAMETER_MATRIX(beta)).
  //nll -= sum(dnorm(beta_coef.vec(), Type(0.), Type(1. / sqrt(2.)), true));
  ADREPORT(betaresp);
  ADREPORT(betapp);
  ADREPORT(beta_coefs_pp);
  //ADREPORT(strparam);
  // ADREPORT parameters for RF.
  ADREPORT(tau);
  ADREPORT(kappa);
  REPORT(x);
  return nll;
}
