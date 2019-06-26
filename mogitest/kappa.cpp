#include "LinearAlgebra.hpp"
#include <iostream>

int main() {
  int n = 5;
  std::vector<std::vector<double>> Hilbert(n, std::vector<double>(n));
  for (double i = 1; i <= n; i++) {
    for (double j = 1; j <= n; j++) {
      Hilbert.at(i-1).at(j-1) = 1/(i+j-1);
    }
  }

  // Hilbert
  std::cout << "Hilbert = \n";
  printMatrix(Hilbert);
  std::cout << "----------\n\n";

  // Inverse_matrix
  std::vector<std::vector<double>> Hilbert_Inverse;
  Hilbert_Inverse = Inverse_matrix(Hilbert);
  std::vector<std::vector<double>> solve = MatrixMatrix(Hilbert, Hilbert_Inverse);
  std::cout << "Hilbert_Inverse = \n";
  printMatrix(Hilbert_Inverse);
  std::cout << "----------\n\n";

  /* std::cout << "Hilbert*Hilbert_Inverse = \n"; */
  /* printMatrix(solve); */
  /* std::cout << "----------\n\n"; */

  // NormInf
  double Hilbert_Norm_inf = MatrixNormInfty(Hilbert);
  printf("Hilbert_Norm_inf = %.2e\n\n",Hilbert_Norm_inf);
  std::cout << "----------\n\n";
  double Hilbert_Inverse_Norm_inf = MatrixNormInfty(Hilbert_Inverse);
  printf("Hilbert_Inverse_Norm_inf = %.2e\n\n",Hilbert_Inverse_Norm_inf);
  std::cout << "----------\n\n";

  // Kappa
  std::cout << "Kappa = \n";
  double kappa = MatrixNormInfty(Hilbert) * MatrixNormInfty(Hilbert_Inverse);
  printf("%.2e\n\n",kappa);
  std::cout << "----------\n\n";


  // b
  std::vector<double> x_tilde(n, 1);
  std::vector<double> b(n);
  b = MatrixVector(Hilbert, x_tilde);
  std::cout << "b = \n";
  printVector(b);
  std::cout << "----------\n\n";

  // LU_compostion
  std::vector<double> x(n, 0.0);
  std::vector<double> y(n, 0.0);
  std::vector<std::vector<double>> Hilbert_LU = LU_decomposition(Hilbert);
  std::vector<std::vector<double>> L = Hilbert_LU;
  std::vector<std::vector<double>> U = Hilbert_LU;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
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

  double b_Ax = VectorNormInfty(VectorSubstract(b, MatrixVector(Hilbert, x)));
  printf("|b - Ax|_inf = %.2e\n", b_Ax);
  std::cout << "----------\n\n";

  double x_tilde_x = VectorNormInfty(VectorSubstract(x_tilde, x));
  printf("|x* - x|_inf = %.2e\n", x_tilde_x);
  std::cout << "----------\n\n";


  return 0;
}
