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
  PARAMETER(mu);
  PARAMETER(alpha);
  PARAMETER(beta);
  // t_n
  Type last = times[times.size() - 1];
  Type term_1 = -mu*last;
  vector<Type> term_2vec = exp(-beta * (last - times)) - 1;
  Type term_2 = alpha/beta*sum(term_2vec);
  vector<Type> A(times.size());
  // initialize negative log-likelihood
  Type nll = 0;
  for(int i = 1; i <= times.size(); i++){
    vector<Type> sub = times(i - 1) - times.head(i - 1);
    vector<Type> temp = -beta*sub;
    A(i - 1) = sum(exp(temp));
  }
  vector<Type> As = log(mu + alpha*A);
  Type term_3 =  sum(As);
  nll = -term_1 - term_2 - term_3;
  return nll;
}

