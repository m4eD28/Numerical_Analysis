#include "../src/LinearAlgebra.hpp"
#include "../src/Algo.hpp"
#include <iostream>

int main() {
  int n = 20;
  double a = 1.0;
  std::vector<std::vector<double>> A(n, std::vector<double>(n));
  for (double i = 1; i <= n; i++) {
    for (double j = 1; j <= n; j++) {
      if(i == j) A.at(i-1).at(i-1) = a;
      else if (0 < fabs(i-j) and fabs(i-j) < 4) A.at(i-1).at(j-1) = 1 / (fabs(i*i - j*j + 15) + 1);
      else A.at(i-1).at(j-1) = 0;
    }
  }
  std::vector<double> b(n, 1.0);
  Gauss_elimination(A, b);
  std::vector<std::vector<double>> A_Inverse;
  A_Inverse = Inverse_matrix(A);
  std::vector<std::vector<double>> solve = MatrixMatrix(A, A_Inverse);
  printMatrix(solve);
  std::cout << "A_Inverse = \n";
  printMatrix(A_Inverse);
  std::cout << "----------\n\n";

  /* std::cout << "A*Hilbert_Inverse = \n"; */
  /* printMatrix(solve); */
  /* std::cout << "----------\n\n"; */

  // NormInf
  double A_Norm_inf = MatrixNormInfty(A);
  printf("A_Norm_inf = %.2e\n\n",A_Norm_inf);
  std::cout << "----------\n\n";
  double A_Inverse_Norm_inf = MatrixNormInfty(A_Inverse);
  printf("A_Inverse_Norm_inf = %.2e\n\n",A_Inverse_Norm_inf);
  std::cout << "----------\n\n";

  // Kappa
  std::cout << "Kappa = \n";
  double kappa = MatrixNormInfty(A) * MatrixNormInfty(A_Inverse);
  printf("%.2e\n\n",kappa);
  std::cout << "----------\n\n";
}
