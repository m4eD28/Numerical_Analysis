#include "LinearAlgebra.hpp"
#include <iostream>

int main() {
  int n = 7;
  vector<vector<double>> Hilbert(n, vector<double>(n));
  vector<vector<double>> Marix_5(n, vector<double>(n));
  vector<double> b_1;
  vector<double> b_2;
  vector<double> x_tilde;
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

  cout << "Hilbert = " << endl;
  printMatrix(Hilbert);
  cout << "------------------\n b = " << endl;
  b_1 = MatrixVector(Hilbert, x_tilde);
  vector<double> b_Hilbert(b_1);
  printVector(b_1); //Hilbert*x_tilde
  cout << "------------------\n Forward_erase" << endl;
  vector<vector<double>> Hilbert_Forward = Forward_erase(Hilbert, b_1);
  printMatrix(Hilbert_Forward);
  cout << "------------------\n Backward_sub" << endl;
  vector<double> x_1 = Backward_sub(Hilbert_Forward, b_1);
  printVector_more_detail(x_1);
  cout << "------------------\n ResidualError" << endl;
  vector<double> Hilbert_Ax_b = ResidualError(Hilbert, x_1, b_Hilbert);
  cout << "Norm1 = " << VectorNorm1(Hilbert_Ax_b) << endl;
  cout << "NormInf = " << VectorNormInfty(Hilbert_Ax_b) << endl;
  cout << "------------------\n Error" << endl;
  vector<double> Hilbert_Error= VectorSubstract(x_1, x_tilde);
  cout << "Norm1 = " << VectorNorm1(Hilbert_Error) << endl;
  cout << "NormInf = " << VectorNormInfty(Hilbert_Error) << endl;
  cout << "------------------\n" << endl;



  cout << "Marix_5 = " << endl;
  printMatrix(Marix_5);
  cout << "------------------\n b = " << endl;
  b_2 = MatrixVector(Marix_5, x_tilde);
  vector<double> b_Marix5(b_2);
  printVector(b_2); //Marix_5*x_tilde
  cout << "------------------\n Forward_erase" << endl;
  vector<vector<double>> Marix_5_Forward = Forward_erase(Marix_5, b_2);
  printMatrix(Marix_5_Forward);
  cout << "------------------\n Backward_sub" << endl;
  vector<double> x_2 = Backward_sub(Marix_5_Forward, b_2);
  printVector_more_detail(x_2);
  cout << "------------------\n ResidualError = " << endl;
  vector<double> Marix_5_Ax_b = ResidualError(Marix_5, x_2, b_Marix5);
  cout << "Norm1 = " << VectorNorm1(Marix_5_Ax_b) << endl;
  cout << "NormInf = " << VectorNormInfty(Marix_5_Ax_b) << endl;
  cout << "------------------\n Error" << endl;
  vector<double> Marix_5_Error= VectorSubstract(x_2, x_tilde);
  cout << "Norm1 = " << VectorNorm1(Marix_5_Error) << endl;
  cout << "NormInf = " << VectorNormInfty(Marix_5_Error) << endl;

  /* cout << "------------------\n Forward_erase" << endl; */
  /* vector<vector<double>> Forward  = Forward_erase(A, b); */
  /* printMatrix(Forward); */
  /* cout << "------------------\n Backward_sub" << endl; */
  /* vector<double> x_x = Backward_sub(Forward, b); */
  /* printVector(x_x); */
  /* cout << "------------------\n ResidualError" << endl; */
  /* vector<double> Hilbert_Ax_b = ResidualError(A, x_1, b_1); */
  /* cout << "Norm1 = \n" << VectorNorm1(Hilbert_Ax_b) << endl; */
  /* cout << "NormInf = \n" << VectorNormInfty(Hilbert_Ax_b) << endl; */
  /* cout << "------------------\n Error" << endl; */
  /* vector<double> Hilbert_Error= VectorSubstract(x_1, x_tilde); */
  /* cout << "Norm1 = \n" << VectorNorm1(Hilbert_Error) << endl; */
  /* cout << "NormInf = \n" << VectorNormInfty(Hilbert_Error) << endl; */
  /* cout << "------------------\n" << endl; */

  /* cout << "Hilbert = " << endl; */
  /* printMatrix(Hilbert); */
  /* cout << "------------------\n b = " << endl; */
  /* b_1 = MatrixVector(Hilbert, x_tilde); */
  /* printVector(b_1); //Hilbert*x_tilde */
  /* Gauss_elimination(Hilbert, b_1); */

  vector<vector<double>> A_Inverse;

  vector< vector<double> > A = {{2, 3, 3}, {3, 2, -1}, {5, 4, 2}};
  vector<double> b = {5, -4, 3};
  A_Inverse = Inverse_matrix(A);
  vector<vector<double>> solve = MatrixMatrix(A, A_Inverse);
  printMatrix(solve);

  return 0;
}
