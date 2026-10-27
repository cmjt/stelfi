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
  //j is for the excited stream, k for the exciting stream, n for time/events.
  //The decay rate of the excitement is that of the excited stream (beta_j), i.e.,
  //lambda_j(t) = mu_j + sum_k alpha_jk sum_{i: stream_i = k, t_i < t} exp(-beta_j (t - t_i))
  
  vector<Type> mu = exp(log_mu);
  vector<Type> beta = exp(log_beta);
  
  // Calculate alphas, enforcing 0 <= alpha_jk <= beta_j
  matrix<Type> alpha(N,N);
  alpha.setZero();
  for (int j = 0; j < N; ++j){
    for (int k = 0; k < N; ++k){
      alpha(j,k) = exp(logit_abratio(j,k) + log_beta(j)) / (Type(1.) + exp(logit_abratio(j,k)));
    }
  }

  // t_n
  Type last = times[(T-1)];
  // A(j,k) = sum_{i: stream_i = k, t_i < t_n} exp(-beta_j (t_n - t_i)), updated recursively
  matrix<Type> A(N,N);
  A.setZero();
  Type nll = 0;
  for(int n = 0; n < T; ++n){
    if (n > 0){
      int M = events[(n-1)]; // stream of the previous event
      Type dt = times[n] - times[(n - 1)];
      for(int j = 0; j < N; ++j){
        for(int k = 0; k < N; ++k){
          Type add = (k == M) ? Type(1.) : Type(0.);
          A(j,k) = exp(-beta[j] * dt) * (A(j,k) + add);
        }
      }
    }
    // log of the intensity of the stream of event n at time t_n
    int S = events[n];
    Type lambda = mu[S];
    for(int k = 0; k < N; ++k){
      lambda += alpha(S,k) * A(S,k);
    }
    nll -= log(lambda);
  }
  // compensator up to the last event; A excludes the last event and EPS
  // is the number of events in each stream excluding the last event
  for(int j = 0; j < N; ++j){
    for(int k = 0; k < N; ++k){
      nll += (alpha(j,k) / beta[j]) * (EPS[k] - A(j,k));
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
