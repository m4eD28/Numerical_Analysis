#include "../src/LinearAlgebra.hpp"
#include "../src/Algo.hpp"
#include <iostream>
#include <random>

std::random_device rnd_dev;
std::mt19937 mt(rnd_dev());
std::uniform_real_distribution<> dist(0, 1);
inline double rnd() {
  return dist(mt);
}

int main() {
  int n = 4;
  std::vector<std::vector<double>> A(n, std::vector<double>(n));
  A = {{1, 2, 1, 2}, {2, 4+pow(10, -10), 2, 3}, {1, 2, 3, 4}, {4, 3, 8, 1}};
  std::vector<double> b(n, 7);
  printMatrix(A);
  printVector(b);
  Gauss_elimination(A, b, n);
  std::printf("-------------------------\n");
  printMatrix(A);
  printVector(b);
  Gauss_elimination_pivot(A, b, n);
  printMatrix(A);
  printVector(b);



  /* std::vector<std::vector<double>> A_Inverse; */
  /* A_Inverse = Inverse_matrix(A); */
  /* std::vector<std::vector<double>> solve = MatrixMatrix(A, A_Inverse); */
  /* printMatrix(solve); */
  /* std::cout << "A_Inverse = \n"; */
  /* printMatrix(A_Inverse); */
  /* std::cout << "----------\n\n"; */

  /* std::cout << "A*Hilbert_Inverse = \n"; */
  /* printMatrix(solve); */
  /* std::cout << "----------\n\n"; */

  // NormInf
  /* double A_Norm_inf = MatrixNormInfty(A); */
  /* printf("A_Norm_inf = %.2e\n\n",A_Norm_inf); */
  /* std::cout << "----------\n\n"; */
  /* double A_Inverse_Norm_inf = MatrixNormInfty(A_Inverse); */
  /* printf("A_Inverse_Norm_inf = %.2e\n\n",A_Inverse_Norm_inf); */
  /* std::cout << "----------\n\n"; */

  /* // Kappa */
  /* std::cout << "Kappa = \n"; */
  /* double kappa = MatrixNormInfty(A) * MatrixNormInfty(A_Inverse); */
  /* printf("%.2e\n\n",kappa); */
  /* std::cout << "----------\n\n"; */

  return 0;
}
