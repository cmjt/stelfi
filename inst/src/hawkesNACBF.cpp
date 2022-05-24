/* Modified version of estimating Hawkes process */
/* The simulation code is added on 20/05/2021. */
/* Negative Alpha and Custom Background Function 24/05/2022 */
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
  DATA_VECTOR(lambda);
  DATA_VECTOR(lambda_min); // min value of lambda between tk and tk+1
  DATA_SCALAR(lambda_integral);
  // parameters of the hawkes process
  PARAMETER(a_par);
  PARAMETER(log_beta);
  Type beta = exp(log_beta);
  //Type alpha = exp(logit_abratio) / (Type(1.) + exp(logit_abratio)) * beta; // enforcing 0<=alpha<=beta

  // t_n
  //Type last = times.template tail<1>()[0];
  vector<Type> A = vector<Type>::Zero(times.size());
  
  Type nll = 0;
  for(int i = 1; i < times.size(); ++i){
    // Page 28 of https://pat-laub.github.io/pdfs/honours_thesis.pdf
    A[i] = exp(-beta * (times[i] - times[i - 1])) * (marks[i-1] + A[i - 1]);
  }
  vector<Type> B = vector<Type>::Zero(times.size());
  for(int i = 0; i < times.size(); ++i){
    B[i] = A[i] + marks[i];
  }
  
  // min(lambda_min/B) <= alpha <= beta
  // A conservative formula for ensuring the intensity is never negative
  vector<Type> C = vector<Type>::Zero(times.size());
  for(int i = 0; i < times.size(); ++i){
    C[i] = lambda_min[i]/B[i];
  }
  
  Type negative = ((a_par)/(1-a_par)) * min(C); // alpha if a_par negative
  Type positive = ((a_par)/(1+a_par)) * beta; // alpha if a_par positive
  Type m = fabs(a_par)/(a_par+1.e-35); // returns 1 if a_par>=0, -1 otherwise
  Type alpha = (0.5*(Type(1.0)+m)*positive) - (0.5*(Type(-1.0)+m)*negative);

  
  vector<Type> term_3vec = log(lambda + alpha * A);
  // part of part 2 is computed in A[times.size() - 1].
  nll = lambda_integral - ((alpha/beta)*A.template tail<1>()[0])+ ((alpha / beta) * Type(times.size() - 1)) - sum(term_3vec);

  //SIMULATE {
    //Type eps = 1e-10, t = 0, M = mu, U;
    //int index = 0;
    //while (index < times.size()){
      //M = mu + alpha * (-beta * (t + eps - times.array().head(index))).exp().sum();
      //t += rexp(Type(1.) / M); U = runif(Type(0.), M); // There is currently a bug as at TMB-1.7.20, 14/05/2021.
      //if (U <= mu + alpha * (-beta * (t - times.array().head(index))).exp().sum()){
        //times[index] = t;
        //index++;
      //}
    //}
    //REPORT(times);
  //}

  ADREPORT(alpha);
  ADREPORT(beta);

  return nll;
}
