#include <TMB.hpp>
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
  std::cout << beta << "\n";
  Type term_1 = -mu*times(times.size());
  std::cout << term_1 << "\n";
  vector<Type> term_2vec = alpha/beta*(exp(-beta * (times(times.size()) - times)) - 1);
  std::cout << term_2vec << "\n";
  Type term_2 = sum(term_2vec);
  vector<Type> A;
  Type nll = 0;
  for(int i = 1; i <times.size(); i++){
    vector<Type> sub(0,(i - 1));
    std::cout << sub << "\n";
    vector<Type> temp = -beta*sub;
    A(i) += sum(exp(temp));
  }
  vector<Type> mualphaA = mu + alpha*A;
  Type term_3 = sum(log(mualphaA));
  nll = -term_1 - term_2 - term_3;
  std::cout << nll << "\n";
  return nll;
}

