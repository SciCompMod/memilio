#include "ad/ad.hpp"
#include <iostream>

template<typename FP> FP square(FP x){
  return x*x;
}

int main() {
  ad::gt1s<double>::type x;
  ad::value(x) = 3.0;
  ad::derivative(x) = 1.0;
  ad::gt1s<double>::type y = square(x);
  std::cout << "value of square(" << ad::value(x) << ") is " << ad::value(y) << std::endl;
  std::cout << "derivative of square(" << ad::value(x) << ") is " << ad::derivative(y) << std::endl;
  return 0;
}
