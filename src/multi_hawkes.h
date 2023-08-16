#ifndef multi_hawkes_hpp
#define multi_hawkes_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/* Estimating multivariavte Hawkes process(es) 15/08/2023 */
template<class Type>
Type multi_hawkes(objective_function<Type>* obj) {
  using namespace Eigen;
  // vector of time
  DATA_INTEGER(N); //Number of dimensions/streams
  DATA_VECTOR(times);
  DATA_IVECTOR(events);
  DATA_IVECTOR(EPS); //Events per stream (excluding the very last event)
  int T = times.size();
  // parameters of the Hawkes process
  PARAMETER_VECTOR(log_mu);
  PARAMETER_MATRIX(logit_abratio);
  PARAMETER_VECTOR(log_beta);
  
  //Indexing scheme:
  //i is for (self-excited) streams/dimensions, j is for time/events, k for self-exciting streams
  
  vector<Type> mu = exp(log_mu);
  vector<Type> beta = exp(log_beta);
  

  // t_n
  Type last = times[(T-1)];
  matrix<Type> A(N,T); //Self-exciting effects
  A.setZero();
  
  for(int j = 1; j < T; ++j){
    int M = events[(j-1)]; // Most recent stream to have had an event
    for(int k = 0; k < N; ++k){
      if (k == M){
        A(k,j) = exp(-beta[k] * (times[j] - times[(j - 1)])) * (Type(1.) + A(k,(j-1)));
      } else {
        A(k,j) = exp(-beta[k] * (times[j] - times[(j - 1)])) * A(k,(j-1));
      }
    }
  }
  
  // Calculate Alphas
  matrix<Type> alpha(N,N);
  alpha.setZero();
  for (int i = 0; i < N; ++i){
    for (int k = 0; k < N; ++k){
      alpha(i,k) = exp(logit_abratio(i,k) + log_beta(k)) / (Type(1.) + exp(logit_abratio(i,k)));
    }
  }
  

  // Third term
  Type nll = 0;
  for(int j = 0; j < T; ++j){
    int M = events[j];
    //calculate lambda of stream M at event j
    nll -= log(mu[M] + ((alpha.row(M) * A.col(j)).sum()));
  }
  
  
  for(int i = 0; i < N; ++i){
    for(int k = 0; k < N; ++k){
      nll += (alpha(i,k)/ beta[k]) * (EPS[k] - A(k, (T-1)));
    }
  }
  nll += sum(mu) * last;
  ADREPORT(mu);
  ADREPORT(alpha);
  ADREPORT(beta);
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
