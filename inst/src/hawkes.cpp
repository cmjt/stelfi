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
  std::cout << last << "\n";
  Type term_1 = -mu*last;
  std::cout << term_1 << "\n";
  vector<Type> term_2vec = alpha/beta*(exp(-beta * (last - times)) - 1);
  std::cout << term_2vec << "\n";
  Type term_2 = sum(term_2vec);
  vector<Type> A;
  Type nll = 0;
  for(int i = 1; i <times.size(); i++){
    Type sub = times[0];
    std::cout << sub << "\n";
    // vector<Type> temp = -beta*sub;
    // A(i) += sum(exp(temp));
  }
  // vector<Type> mualphaA = mu + alpha*A;
  // Type term_3 = sum(log(mualphaA));
  // nll = -term_1 - term_2 - term_3;
  std::cout << nll << "\n";
  return nll;
}

