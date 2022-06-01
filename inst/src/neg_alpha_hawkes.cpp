/* Modified version of estimating Hawkes process */
/* The simulation code is added on 20/05/2021. */
#include <TMB.hpp>
#include <vector>
#include <iostream>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace Eigen;
  // vector of time
  DATA_VECTOR(times);
  DATA_VECTOR(marks);
  Type marks_mean = marks.sum()/marks.size(); // Average mark
  // parameters of the hawkes process
  PARAMETER(log_mu);
  PARAMETER(a_par);
  PARAMETER(log_beta);
  Type mu = exp(log_mu);
  Type beta = exp(log_beta);
 
  // t_n
  Type last = times.template tail<1>()[0];
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
  
  // -min(lambda/B) <= alpha <= beta
  Type Max = beta/marks_mean;
  Type Min = mu/max(B);
  Type alpha = (exp(a_par) / (Type(1.) + exp(a_par)))*(Max+Min) - Min;

  
  vector<Type> term_3vec = log(mu + alpha * A);
  // part of part 2 is computed in A[times.size() - 1].
  nll = (mu * last) - ((alpha/beta)*A.template tail<1>()[0])+ ((alpha / beta) * Type(sum(marks)-marks.template tail<1>()[0])) - sum(term_3vec);

  SIMULATE {
    Type eps = 1e-10, t = 0, M = mu, U;
    int index = 0;
    while (index < times.size()){
      M = mu + alpha * (-beta * (t + eps - times.array().head(index))).exp().sum();
      t += rexp(Type(1.) / M); U = runif(Type(0.), M); // There is currently a bug as at TMB-1.7.20, 14/05/2021.
      if (U <= mu + alpha * (-beta * (t - times.array().head(index))).exp().sum()){
        times[index] = t;
        index++;
      }
    }
    REPORT(times);
  }

  ADREPORT(mu);
  ADREPORT(alpha);
  ADREPORT(beta);

  return nll;
}
