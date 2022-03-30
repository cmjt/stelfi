/* written by Xiangjie Xue, updated on 01/05/2021. */
/* Linearisation of a matrix can be directly done using .vec() */
/* Added simulation code on 20/05/2021. */
#include <TMB.hpp>
#include <numeric>
#include <math.h> // for pi

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla; // Where Q_spde is defined.
  using namespace density; // this where the structure for GMRF and AR_t is defined
  using namespace Eigen;  // probably for sparseness class
  DATA_VECTOR(y); // the vector of observations. 
  /*
    For spatial only, this is the vector of counts.
    For spatial + AR1 temporal, this vector needs to be 
    arranged in the following way:
    y                 t
    ---------------------
    y_1               1
    y_2               1
    ...
    y_n               1
    y_{n+1}           2
    y_{n+2}           2
    ... 
    y_{n * t_n - 1}   t_n
    y_{n * t_n}       t_n
  */
  DATA_SPARSE_MATRIX(A); // the predicator matrix A.
  /*
    A matrix of size y.size() by x.rows() * x.cols() (or xref2.size()).
    This sparse matrix A maps ALL of the random effects (linearised 
    version of x, i.e., xref2) to the observation vector y. 
  */
  DATA_MATRIX(designmat);  // the design matrix for the fixed effects.
  /*
    A matrix of size y.size() by no. of fixed effects including intercepts.
    Set it as cbind(1, matrix of covariates).
    If only the intercept needs to be estimated, set it as a MATRIX of 1s
    with ncol = 1 and nrow = y.size().
  */
  DATA_STRUCT(spde,spde_t); // this as structure of spde object.
  /*
   This is defined in r-inla. By inspection of the template, the minimal
   required components are M0, M1, M2.
  */
  DATA_VECTOR(w); 
  /*
    This vector corresponds to the E term for poisson models in INLA;
    see INLA::inla.doc("poisson") for more detail.
    If not required, set it to a vector of 1s or 0s (by definition in INLA docs).
  */
  DATA_IVECTOR(idx);
  /*
    This vector should be a binary vector of the same size as the observation vector y.
    With this vector, the log-likelihood can be computed using a subset
    of the observations. 1 for contributing to the log-likelihood, and 0 otherwise.
  */
  PARAMETER_VECTOR(beta); // A vector of fixed effects needs to be estimated.
  /*
    This vector is of length designmat.cols(). The appropriate length should be
    set before passing from R.
  */
  PARAMETER_ARRAY(x); // the random field/effects
  /*
    This is the random field/effects. Set this variable random in MadeADFun().
    This array is of size no. of random effect for each time knots (x.rows()) by no. of temporal 
    knots (x.cols()), and hence the total number of random effects is 
    x.rows() * x.cols().
    Each column represents a time knot:
              t
              1           2           3           ...           t_n
  ----------------------------------------------------------------------
  random  
  effects x.col(0)     x.col(1)    x.col(2)       ...         x.col(t_n - 1)
  ----------------------------------------------------------------------
  */
  PARAMETER(log_tau); // tau parameter for the RF
  PARAMETER(log_kappa); // kappa parameter for the RF
  Type kappa = exp(log_kappa);
  Type tau = exp(log_tau);
  
  // Section 5.1.1 of https://becarioprecario.bitbucket.io/spde-gitbook/ch-nonstationarity.html.
  SparseMatrix<Type> Q = Q_spde(spde, kappa) * pow(tau,2); 
  Type nll = 0.0;
  /*
    The following chunk is for the log-likelihood of the random effects.
    x.cols() == 1 corresponds to spatial only. In this case, the AR_1 
    parameter atanh_rho is not required from R side, and the ADREPORT for
    rho would also not be shown in the exported object.
    For spatial-temporal models, atanh_rho is required from R side and 
    ADREPORT for rho will be exported. 
    x.col(i) = rho * x.col(i - 1) + \epsilon_i, where \epsilon_i ~ N(0, Q).
  */
  if (x.cols() == 1){
    nll = GMRF(Q)(x);

    SIMULATE{
      array<Type> xmat(x.size(), 2);
      AR1_t<GMRF_t<Type> >(Type(0.), GMRF(Q)).simulate(xmat);
      x = xmat.col(0);
    }
  }else{
    PARAMETER(atanh_rho); // AR1 parameter
    Type rho = tanh(atanh_rho); // translate from R to [-1, 1].
    nll += AR1_t<GMRF_t<Type> >(rho, GMRF(Q))(x); 
    /*
      rho: AR1 parameter.
      Here x has to be an array in order for computation in one go.
      See https://kaskr.github.io/adcomp/classdensity_1_1AR1__t.html.
    */

    SIMULATE{
      AR1_t<GMRF_t<Type> >(rho, GMRF(Q)).simulate(x);
    }

    ADREPORT(rho);
  }
  
  /*
    construct the linear predictor eta* = Eexp(eta), with random effects
    transformed by A and add the fixed effects.
  */
  vector<Type> lambda = exp(A * x.vec() + designmat * beta) * w.cwiseEqual(0).select(vector<Type>::Ones(w.size()), w);
  /*
    construct the likelihood with binary variable idx.
  */
  nll -= idx.cwiseEqual(0).select(vector<Type>::Zero(y.size()), dpois(y, lambda, true)).sum(); // construct likelihood

  // simulation chunk
  SIMULATE {

    lambda = exp(A * x.vec() + designmat * beta) * w.cwiseEqual(0).select(vector<Type>::Ones(w.size()), w);

    y = rpois(lambda);

    REPORT(y);
  }
  ADREPORT(beta);
  // ADREPORT parameters for RF.
  ADREPORT(log_tau);
  ADREPORT(log_kappa);
  // ADREPORT transformed params
  Type range = sqrt(8)/exp(log_kappa); //as reported by INLA
  Type stdev = 1/(4*M_PI*pow(exp(log_kappa),2)*pow(exp(log_tau),2)) ; // as reported by INLA
  ADREPORT(range);
  ADREPORT(stdev);
  REPORT(x);
  return nll;
}
