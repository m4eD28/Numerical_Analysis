#include "../src/LinearAlgebra.hpp"
#include "../src/Algo.hpp"
#include <iostream>

int main() {
  int n = 20;
  double a = 0.1;
  std::vector<std::vector<double>> Hilbert(n, std::vector<double>(n));
  for (double i = 1; i <= n; i++) {
    for (double j = 1; j <= n; j++) {
      if(i == j) Hilbert.at(i-1).at(i-1) = a;
      else if (0 < fabs(i-j) and fabs(i-j) < 4) 
        Hilbert.at(i-1).at(j-1) = 1 / (fabs(i*i - j*j + 15) + 1);
      else {
        Hilbert.at(i-1).at(j-1) = 0;
      }
    }
  }

  std::vector<double> gaus_b(n, 1.0);
  std::vector<double> x_tilde(n, 1.0);
  std::vector<double>x = Gauss_elimination(Hilbert, gaus_b);
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

  /* double b_Ax = VectorNormInfty(VectorSubstract(gaus_b, MatrixVector(Hilbert, gaus_b)))/VectorNormInfty(gaus_b); */
  double b_Ax = VectorNormInfty(VectorSubstract(MatrixVector(Hilbert, gaus_b), gaus_b))/VectorNormInfty(gaus_b);
  printf("|b - Ax|_inf = %.2e\n", b_Ax);
  std::cout << "----------\n\n";

  double x_tilde_x = VectorNormInfty(VectorSubstract(x_tilde, x))/VectorNormInfty(x_tilde);
  printf("|x* - x|_inf = %.2e\n", x_tilde_x);
  std::cout << "----------\n\n";

  /* std::vector<double> x_tilde(n, 1); */
  /* std::vector<double> b(n); */
  /* b = MatrixVector(Hilbert, x_tilde); */
  /* std::cout << "b = \n"; */
  /* printVector(b); */
  /* std::cout << "----------\n\n"; */

  /* std::vector<double> x(n, 0.0); */
  /* std::vector<double> y(n, 0.0); */
  /* std::vector<std::vector<double>> Hilbert_LU = LU_decomposition(Hilbert); */
  /* std::vector<std::vector<double>> L = Hilbert_LU; */
  /* std::vector<std::vector<double>> U = Hilbert_LU; */
  /* for (int i = 0; i < n; i++) { */
  /*   for (int j = 0; j < n; j++) { */
  /*     if (i == j) { */
  /*       L.at(i).at(j) = 1.0; */
  /*     } else if (i < j) { */
  /*       L.at(i).at(j) = 0.0; */
  /*     } else if (i > j) { */
  /*       U.at(i).at(j) = 0.0; */
  /*     } */
  /*   } */
  /* } */

  /* y = Forward_sub(L, b); */
  /* x = Backward_sub(U, y); */
  /* std::cout << "x = \n"; */
  /* printVector_more_detail(x); */
  /* std::cout << "----------\n\n"; */

  /* double b_Ax = VectorNormInfty(VectorSubstract(b, MatrixVector(Hilbert, x))); */
  /* printf("|b - Ax|_inf = %.2e\n", b_Ax); */
  /* std::cout << "----------\n\n"; */

  /* double x_tilde_x = VectorNormInfty(VectorSubstract(x_tilde, x)); */
  /* printf("|x* - x|_inf = %.2e\n", x_tilde_x); */
  /* std::cout << "----------\n\n"; */
/* 解の誤差と右辺項の摂動は次の関係式が成り立つ. */
/* ∥x − x ̃∥∞ ≤ κ∞(A)∥Ax ̃ − b∥∞ = κ∞(A)∥∆b∥∞ ∥x∥∞ ∥b∥∞ ∥b∥∞ */
/* 右辺項に混入した誤差が微小であっても，解に含まれる誤差は条件数に応じて大きくなる 可能性があることが分かる.以上の結果より，実際の解の誤差ノルムはほぼ「条件数 × 右 辺項の摂動」になった. */

  return 0;
}
