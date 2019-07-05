#include "LinearAlgebra.hpp"
#include <iostream>

int main() {
  int n = 5;
  std::vector<std::vector<double>> A(n, std::vector<double>(n));
  for (double i = 1; i <= n; i++) {
    for (double j = 1; j <= n; j++) {
      A.at(i-1).at(j-1) = 1/(i+j-1);
    }
  }
  // b
  std::vector<double> x_tilde(n, 1);
  std::vector<double> b(n);
  b = MatrixVector(A, x_tilde);
  std::cout << "b = \n";
  printVector(b);
  std::cout << "----------\n\n";

  // LU_compostion
  std::vector<double> x(n, 0.0);
  std::vector<double> y(n, 0.0);
  std::vector<std::vector<double>> A_LU = LU_decomposition(A);
  std::vector<std::vector<double>> L = A_LU;
  std::vector<std::vector<double>> U = A_LU;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      if (i == j) {
        L.at(i).at(j) = 1.0;
      } else if (i < j) {
        L.at(i).at(j) = 0.0;
      } else if (i > j) {
        U.at(i).at(j) = 0.0;
      }
    }
  }

  y = Forward_sub(L, b);
  x = Backward_sub(U, y);
  std::cout << "x = \n";
  printVector_more_detail(x);
  std::cout << "----------\n\n";

  double b_Ax = VectorNormInfty(VectorSubstract(b, MatrixVector(A, x)));
  printf("|b - Ax|_inf = %.2e\n", b_Ax);
  std::cout << "----------\n\n";

  double x_tilde_x = VectorNormInfty(VectorSubstract(x_tilde, x));
  printf("|x* - x|_inf = %.2e\n", x_tilde_x);
  std::cout << "----------\n\n";


  return 0;
}
