#include <TMB.hpp>

template<class Type> // covariate list of matrices for each time step
struct covariate_list : vector<matrix <Type> > {
  covariate_list(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

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
  using namespace R_inla;
  using namespace density; // this where the structure for GMRF is defined
  using namespace Eigen;  // probably for sparseness calcs
  // Poisson response
  DATA_MATRIX(y);
  // covariates as a list of matricies
  DATA_STRUCT(covariates, covariate_list); // list of design matrcies for each time step
  // weight of each mesh node
  DATA_VECTOR(area);
  // number of time steps
  DATA_INTEGER(tsteps);
  // this as structure of spde object is defined in r-inla hence using that namespce
  DATA_STRUCT(spde,spde_t);
  // vector of regression coefficients
  PARAMETER_VECTOR(beta);
  // Parameters for the AR(1) part of the spatiotemporal field
  PARAMETER(log_rho);
  PARAMETER(log_sigma);
  // Parameter for the spatial GMRF of the spatiotemporal field
  PARAMETER(log_kappa);
  // array (n.y * tsteps) of spatiotemporal field random variables.
  PARAMETER_ARRAY(field);
  // transform, rho of AR1 process
  Type rho = exp(log_rho);
  // transform, sigma of AR1 process
  Type sigma = exp(log_sigma);
  // transform, kappa parameter of the field
  Type kappa = exp(log_kappa);
  // initialise neg log-likelihood
  Type nll = 0;
  // create the precision matrix from the spde model for the GMRF
  SparseMatrix<Type> Q = Q_spde(spde,kappa);
  // component due to spatiotemporal field with precision Q
  nll = SEPARABLE(SCALE(AR1(rho),sigma), GMRF(Q))(field);
  // data contribution
  for(int i = 0; i< (tsteps - 1); i++){
    // design matrix and regression coeffs. (fixed effects)
    vector<Type> eta = covariates(i)*beta; 
    vector<Type> gmrf = (vector<Type> (field.col(i)));
    vector<Type> respi = y.col(i);
    for(int j = 0; j <respi.size(); j++){
      Type mu;
      mu = eta(j) +  gmrf(j);
      Type lambda = area(j)*exp(mu); // intensity
      nll -= dpois_stable(respi(j), lambda, true);
    }
   
  }
  // ADREPORT(kappa);
  // Type sigma2 = 1/(4*M_PI*pow(kappa,2));   // sigma2 = 1/(4*pi*k^2) with tau = 1
  // ADREPORT(sigma2);
  // Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011)
  // ADREPORT(range);
  // SIMULATE {
  //   SEPARABLE(AR1(rho), GMRF(Q)).simulate(x);
  //   for(int i = 0; i< (tsteps - 1); i++){
  //     vector<Type> eta_s = covariates(i)*beta; // design matrix and regression coeffs. (fixed effects)
  //     vector<Type> gmrf_s = (vector<Type> (x.col(i)));
  //     //vector<Type> respi_s = resp(i);
  //     for(int j = 0; j <resp(i).size(); j++){
  // 	Type mu_s;
  // 	mu_s = eta_s(j) +  gmrf_s(j);
  // 	Type lam_s = area(j)*exp(mu_s); // intensity
  // 	resp(i)(j)= rpois(lam_s);
  //     }
  // REPORT(resp);          // Report the simulation
  // }
  // }
  return nll;
}
