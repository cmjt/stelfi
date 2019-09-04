#include <TMB.hpp>
#include <vector>
#include <iostream>
using namespace R_inla;
using namespace density; // this where the structure for GMRF is defined
using namespace Eigen;  // probably for sparseness calcs
template<class Type> // dpois to deal with low expectations
Type dpois_stable (const Type &x, const Type &lambda, const int &give_log){
  Type out;
  out = pow(lambda, x)*exp(-lambda)/exp(lgamma(x + 1));
  if (give_log){
    out = log(out + DBL_MIN);
  }
  return out;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(resp);
  DATA_MATRIX(covariates); // matrix of covariates, one column per covariate
  DATA_VECTOR(area);
  DATA_VECTOR(times);
  DATA_STRUCT(spde,spde_t); // this as structure of spde object is defined in r-inla hence using that namespce
  PARAMETER_VECTOR(gammas); // regression coefficients
  PARAMETER(log_kappa); // kappa of random field
  PARAMETER_VECTOR(x); //the random field/effect each matrix row is a time step
  PARAMETER(mu); // baseline parameter of the hawkes temporal process
  PARAMETER(alpha); // jump parameter of the hawkes temporal process
  PARAMETER(beta); // decay parameter of the hawkes temporal process
  Type kappa = exp(log_kappa);
  // return the kappa parameter of the random field
  SparseMatrix<Type> Q = Q_spde(spde,kappa); // create the precision matrix from the spde model for the GMRF
  // t_n
  Type last = times[times.size() - 1];
  Type term_1 = -mu*last;
  vector<Type> term_2vec = exp(-beta * (last - times)) - 1;
  Type term_2 = alpha/beta*sum(term_2vec);
  vector<Type> A(times.size());
  // initialize negative log-likelihood
  Type nll = 0;
  // hawkes
  for(int i = 1; i <= times.size(); i++){
    vector<Type> sub = times(i - 1) - times.head(i - 1);
    vector<Type> temp = -beta*sub;
    A(i - 1) = sum(exp(temp));
  }
  vector<Type> As = log(mu + alpha*A);
  Type term_3 =  sum(As);
  vector<Type> design = covariates*gammas;
  for(int i = 0; i <resp.size(); i++){
    Type eta = design(i) + log(area(i)) + x(i);
    Type lambda = exp(eta); // intensity
    nll -= dpois_stable(resp(i),lambda,true) + term_1 + term_2 + term_3; 
  }
  //nll += -term_1 - term_2 - term_3;
  ADREPORT(kappa);
  Type sigma2 = 1/(4*M_PI*pow(kappa,2));   // sigma2 = 1/(4*pi*k^2) with tau = 1
  ADREPORT(sigma2);
  Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011)
  ADREPORT(range);
  return nll;
}
