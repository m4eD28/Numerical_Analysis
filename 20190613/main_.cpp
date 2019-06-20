#include "LinearAlgebra.hpp"
#include <iostream>

int main() {
  int n = 11;
  vector<vector<double>> Hilbert(n, vector<double>(n));
  for (double i = 1; i <= n; i++) {
    for (double j = 1; j <= n; j++) {
      Hilbert.at(i-1).at(j-1) = 1/(i+j-1);
    }
  }

  // Hilbert
  cout << "Hilbert = \n";
  printMatrix(Hilbert);
  cout << "----------\n\n";

  // Inverse_matrix
  vector<vector<double>> Hilbert_Inverse;
  Hilbert_Inverse = Inverse_matrix(Hilbert);
  vector<vector<double>> solve = MatrixMatrix(Hilbert, Hilbert_Inverse);
  cout << "Hilbert_Inverse = \n";
  printMatrix(Hilbert_Inverse);
  cout << "----------\n\n";

  cout << "Hilbert*Hilbert_Inverse = \n";
  printMatrix(solve);
  cout << "----------\n\n";

  // NormInf
  double Hilbert_Norm_inf = MatrixNormInfty(Hilbert);
  printf("Hilbert_Norm_inf = %.2e\n\n",Hilbert_Norm_inf);
  cout << "----------\n\n";
  double Hilbert_Inverse_Norm_inf = MatrixNormInfty(Hilbert_Inverse);
  printf("Hilbert_Inverse_Norm_inf = %.2e\n\n",Hilbert_Inverse_Norm_inf);
  cout << "----------\n\n";

  // Kappa
  cout << "Kappa = \n";
  double kappa = MatrixNormInfty(Hilbert) * MatrixNormInfty(Hilbert_Inverse);
  printf("%.2e\n\n",kappa);
  cout << "----------\n\n";


  // b
  vector<double> x_tilde(n, 1);
  vector<double> b(n);
  b = MatrixVector(Hilbert, x_tilde);
  cout << "b = \n";
  printVector(b);
  cout << "----------\n\n";

  // LU_compostion
  vector<double> x(n, 0);
  vector<double> y(n, 0);
  vector<vector<double>> Hilbert_LU = LU_decomposition(Hilbert);
  vector<vector<double>> L = Hilbert_LU;
  vector<vector<double>> U = Hilbert_LU;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) {
        L.at(i).at(j) = 1;
      } else if (i < j) {
        L.at(i).at(j) = 0;
      } else if (i > j) {
        U.at(i).at(j) = 0;
      }
    }
  }

  y = Forward_sub(Hilbert_LU, b);
  x = Backward_sub(Hilbert_LU, y);
  cout << "x = \n";
  printVector_more_detail(x);
  cout << "----------\n\n";

  double b_Ax = VectorNormInfty(VectorSubstract(b, MatrixVector(Hilbert, x)));
  printf("|b - Ax|_inf = %.2e\n", b_Ax);
  cout << "----------\n\n";

  double x_tilde_x = VectorNormInfty(VectorSubstract(x_tilde, x));
  printf("|x* - x|_inf = %.2e\n", x_tilde_x);
  cout << "----------\n\n";


  return 0;
}
