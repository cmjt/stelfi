#ifndef neg_alpha_custom_hawkes_hpp
#define neg_alpha_custom_hawkes_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/* Modified version of estimating Hawkes process */
/* The simulation code is added on 20/05/2021. */
/* Generic Inhomogenous Hawke's processes(either self-exciting or self-inhibiting) 24/05/2022 */
template<class Type>
Type neg_alpha_custom_hawkes(objective_function<Type>* obj) {
  using namespace Eigen;
  // vector of time
  DATA_VECTOR(times);
  DATA_VECTOR(marks);
  Type marks_mean = marks.sum()/marks.size(); // Average mark
  DATA_VECTOR(lambda);
  DATA_VECTOR(lambda_min); // min value of lambda between tk and tk+1
  DATA_SCALAR(lambda_integral);
  // parameters of the hawkes process
  PARAMETER(a_par);
  PARAMETER(log_beta);
  Type beta = exp(log_beta);
  vector<Type> A = vector<Type>::Zero(times.size());
  Type nll = 0;
  // A[i]*alpha is self-exciting component immediately before event i
  for(int i = 1; i < times.size(); ++i){
    // Page 28 of https://pat-laub.github.io/pdfs/honours_thesis.pdf
    A[i] = exp(-beta * (times[i] - times[i - 1])) * (marks[i-1] + A[i - 1]);
  }
  // B[i]*alpha is self-exciting component immediately after event i
  vector<Type> B = vector<Type>::Zero(times.size());
  for(int i = 0; i < times.size(); ++i){
    B[i] = A[i] + marks[i];
  }
  // -min(lambda_min/B) <= alpha <= beta/mean(marks)
  // A conservative formula for ensuring the intensity is never negative
  vector<Type> C = vector<Type>::Zero(times.size());
  for(int i = 0; i < times.size(); ++i){
    C[i] = lambda_min[i]/B[i];
  }
  Type Max = beta/marks_mean;
  Type alpha = ((exp(a_par) / (Type(1.) + exp(a_par)))*(Max+min(C))) - min(C);
  vector<Type> term_3vec = log(lambda + alpha * A);
  nll = lambda_integral + ((alpha/beta) * Type(sum(marks) - marks.template tail<1>()[0] - A.template tail<1>()[0])) - sum(term_3vec);
  ADREPORT(alpha);
  ADREPORT(beta);
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
