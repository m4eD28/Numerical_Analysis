#include "LinearAlgebra.hpp"
#include <iostream>

int main() {
  int n = 7;
  std::vector<std::vector<double>> Hilbert(n, std::vector<double>(n));
  std::vector<std::vector<double>> Marix_5(n, std::vector<double>(n));
  std::vector<double> b_1;
  std::vector<double> b_2;
  std::vector<double> x_tilde;
  for (double i = 1; i <= n; i++) {
    for (double j = 1; j <= n; j++) {
      Hilbert.at(i-1).at(j-1) = 1/(i+j);
    }
  }
  for (double i = 1; i <= n; i++) {
    for (double j = 1; j <= n; j++) {
      if (i == j) Marix_5.at(i-1).at(j-1) = 5.0;
      else if (0 < fabs(i-j) and fabs(i-j) <= 2) Marix_5.at(i-1).at(j-1) = -1.0;
      else Marix_5.at(i-1).at(j-1) = 0;
    }
  }
  for (int i = 0; i < n; i++) {
    x_tilde.emplace_back(1);
  }

  std::cout << "Hilbert = " << std::endl;
  printMatrix(Hilbert);
  std::cout << "------------------\n b = " << std::endl;
  b_1 = MatrixVector(Hilbert, x_tilde);
  std::vector<double> b_Hilbert(b_1);
  printVector(b_1); //Hilbert*x_tilde
  std::cout << "------------------\n Forward_erase" << std::endl;
  std::vector<std::vector<double>> Hilbert_Forward = Forward_erase(Hilbert, b_1);
  printMatrix(Hilbert_Forward);
  std::cout << "------------------\n Backward_sub" << std::endl;
  std::vector<double> x_1 = Backward_sub(Hilbert_Forward, b_1);
  printVector_more_detail(x_1);
  std::cout << "------------------\n ResidualError" << std::endl;
  std::vector<double> Hilbert_Ax_b = ResidualError(Hilbert, x_1, b_Hilbert);
  std::cout << "Norm1 = " << VectorNorm1(Hilbert_Ax_b) << std::endl;
  std::cout << "NormInf = " << VectorNormInfty(Hilbert_Ax_b) << std::endl;
  std::cout << "------------------\n Error" << std::endl;
  std::vector<double> Hilbert_Error= VectorSubstract(x_1, x_tilde);
  std::cout << "Norm1 = " << VectorNorm1(Hilbert_Error) << std::endl;
  std::cout << "NormInf = " << VectorNormInfty(Hilbert_Error) << std::endl;
  std::cout << "------------------\n" << std::endl;



  std::cout << "Marix_5 = " << std::endl;
  printMatrix(Marix_5);
  std::cout << "------------------\n b = " << std::endl;
  b_2 = MatrixVector(Marix_5, x_tilde);
  std::vector<double> b_Marix5(b_2);
  printVector(b_2); //Marix_5*x_tilde
  std::cout << "------------------\n Forward_erase" << std::endl;
  std::vector<std::vector<double>> Marix_5_Forward = Forward_erase(Marix_5, b_2);
  printMatrix(Marix_5_Forward);
  std::cout << "------------------\n Backward_sub" << std::endl;
  std::vector<double> x_2 = Backward_sub(Marix_5_Forward, b_2);
  printVector_more_detail(x_2);
  std::cout << "------------------\n ResidualError = " << std::endl;
  std::vector<double> Marix_5_Ax_b = ResidualError(Marix_5, x_2, b_Marix5);
  std::cout << "Norm1 = " << VectorNorm1(Marix_5_Ax_b) << std::endl;
  std::cout << "NormInf = " << VectorNormInfty(Marix_5_Ax_b) << std::endl;
  std::cout << "------------------\n Error" << std::endl;
  std::vector<double> Marix_5_Error= VectorSubstract(x_2, x_tilde);
  std::cout << "Norm1 = " << VectorNorm1(Marix_5_Error) << std::endl;
  std::cout << "NormInf = " << VectorNormInfty(Marix_5_Error) << std::endl;

  /* std::cout << "------------------\n Forward_erase" << std::endl; */
  /* vector<vector<double>> Forward  = Forward_erase(A, b); */
  /* printMatrix(Forward); */
  /* std::cout << "------------------\n Backward_sub" << std::endl; */
  /* vector<double> x_x = Backward_sub(Forward, b); */
  /* printVector(x_x); */
  /* std::cout << "------------------\n ResidualError" << std::endl; */
  /* vector<double> Hilbert_Ax_b = ResidualError(A, x_1, b_1); */
  /* std::cout << "Norm1 = \n" << VectorNorm1(Hilbert_Ax_b) << std::endl; */
  /* std::cout << "NormInf = \n" << VectorNormInfty(Hilbert_Ax_b) << std::endl; */
  /* std::cout << "------------------\n Error" << std::endl; */
  /* vector<double> Hilbert_Error= VectorSubstract(x_1, x_tilde); */
  /* std::cout << "Norm1 = \n" << VectorNorm1(Hilbert_Error) << std::endl; */
  /* std::cout << "NormInf = \n" << VectorNormInfty(Hilbert_Error) << std::endl; */
  /* std::cout << "------------------\n" << std::endl; */

  /* std::cout << "Hilbert = " << std::endl; */
  /* printMatrix(Hilbert); */
  /* std::cout << "------------------\n b = " << std::endl; */
  /* b_1 = MatrixVector(Hilbert, x_tilde); */
  /* printVector(b_1); //Hilbert*x_tilde */
  /* Gauss_elimination(Hilbert, b_1); */

  std::vector<std::vector<double>> A_Inverse;

  std::vector<std::vector<double> > A = {{2, 3, 3}, {3, 2, -1}, {5, 4, 2}};
  std::vector<double> b = {5, -4, 3};
  A_Inverse = Inverse_matrix(A);
  std::vector<std::vector<double>> solve = MatrixMatrix(A, A_Inverse);
  printMatrix(solve);

  return 0;
}
