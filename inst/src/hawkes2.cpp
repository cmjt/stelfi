/* Modified version of estimating Hawkes process */
/* The simulation code is added on 20/05/2021. */
/* Expanded to self-inhibiting Hawke's processes 19/5/2022 */
#include <TMB.hpp>
#include <vector>
#include <iostream>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace Eigen;
  // vector of time
  DATA_VECTOR(times);
  // parameters of the hawkes process
  PARAMETER(log_mu);
  PARAMETER(logit_abratio);
  PARAMETER(log_beta);
  PARAMETER(logit_amuratio);
  PARAMETER(log_b);
  PARAMETER(atanh_c);
  Type b = exp(log_b); //enforcing b > 0
  Type c = tanh(atanh_c) * M_PI; //enforcing -pi<c<pi
  Type mu = exp(log_mu);
  Type a = exp(logit_amuratio) / (Type(1.) + exp(logit_amuratio)) * mu; // enforcing 0<=a<=mu
  Type beta = exp(log_beta);
  Type alpha = exp(logit_abratio) / (Type(1.) + exp(logit_abratio)) * beta; // enforcing 0<=alpha<=beta

  // t_n
  Type last = times.template tail<1>()[0];
  vector<Type> A = vector<Type>::Zero(times.size());
  
  // Background lambda
  vector<Type> lambda = vector<Type>::Zero(times.size());
  for(int i = 0; i < times.size(); ++i){
    lambda[i] = mu + a * sin((b*times[i])+c);
  }
  
  Type nll = 0;
  for(int i = 1; i < times.size(); ++i){
    // Page 28 of https://pat-laub.github.io/pdfs/honours_thesis.pdf
    A[i] = exp(-beta * (times[i] - times[i - 1])) * (Type(1.0) + A[i - 1]);
  }
  vector<Type> term_3vec = log(lambda + alpha * A);
  Type lambda_integral = ((a * cos(c) - a * cos((b*last)+c)) / b) + mu * last;
  // part of part 2 is computed in A[times.size() - 1].
  nll = lambda_integral - ((alpha/beta)*A.template tail<1>()[0])+ ((alpha / beta) * Type(times.size() - 1)) - sum(term_3vec);

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
  ADREPORT(a);
  ADREPORT(b);
  ADREPORT(c);

  return nll;
}