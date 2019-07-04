#include "LinearAlgebra.hpp"
#include <iostream>

int main() {
  int n = 20;
  double a = 1.0;
  std::vector<std::vector<double>> A(n, std::vector<double>(n));
  for (double i = 1; i <= n; i++) {
    for (double j = 1; j <= n; j++) {
      A.at(i-1).at(j-1) = 1/(i+j-1);
    }
  }
  std::vector<double> b(n, 1.0);
  Gauss_elimination(A, b);
}
