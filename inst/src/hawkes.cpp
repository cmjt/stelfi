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
  Type last = times[times.size() - 1];
  Type term_1 = -mu*last;
  vector<Type> term_2vec = alpha/beta*(exp(-beta * (last - times)) - 1);
  Type term_2 = sum(term_2vec);
  vector<Type> A;
  Type nll = 0;
  for(int i = 1; i <=times.size(); i++){
    vector<Type> sub = times.head(i-1);
    vector<Type> temp = -beta*sub;
    Type etemp = sum(exp(temp));
    Type mualphaA = mu + alpha*etemp;
    Type term_3 = log(mualphaA);
    nll += -term_1 - term_2 - term_3;
  }
  std::cout << nll << "\n";
  return nll;
}

