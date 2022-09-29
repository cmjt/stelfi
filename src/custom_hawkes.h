#ifndef custom_hawkes_hpp
#define custom_hawkes_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/* Modified version of estimating Hawkes process */
/* The simulation code is added on 20/05/2021. */
/* Generic Inhomogenous Self-exciting Hawke's processes 23/05/2022 */
template<class Type>
Type custom_hawkes(objective_function<Type>* obj) {
  using namespace Eigen;
  // vector of time
  DATA_VECTOR(times);
  DATA_VECTOR(lambda);
  DATA_VECTOR(marks);
  Type marks_mean = marks.sum()/marks.size(); // Average mark
  DATA_SCALAR(lambda_integral);
  
  // parameters of the hawkes process
  PARAMETER(logit_abratio);
  PARAMETER(log_beta);

  Type beta = exp(log_beta);
  Type alpha = exp(logit_abratio) / (Type(1.) + exp(logit_abratio)) * (beta/marks_mean);
  // enforcing 0<=alpha<=beta
  vector<Type> A = vector<Type>::Zero(times.size());
  Type nll = 0;
  for(int i = 1; i < times.size(); ++i){
    // Page 28 of https://pat-laub.github.io/pdfs/honours_thesis.pdf
    A[i] = exp(-beta * (times[i] - times[i - 1])) * (marks[i-1] + A[i - 1]);
  }
  vector<Type> term_3vec = log(lambda + alpha * A);
  nll = lambda_integral + ((alpha/beta) * Type(sum(marks) - marks.template tail<1>()[0] - A.template tail<1>()[0])) - sum(term_3vec);
  ADREPORT(alpha);
  ADREPORT(beta);
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
