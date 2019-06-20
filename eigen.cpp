#include "LinearAlgebra.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include "Eigen/Core"
#include "Eigen/LU"
using namespace Eigen;


int main() {
  int n = 5;
  vector<vector<double>> Hilbert(n, vector<double>(n));
  /* MatrixXd Hilbert = MarixXd::(n, n); */
  for (double i = 1; i <= n; i++) {
    for (double j = 1; j <= n; j++) {
      Hilbert.at(i-1).at(j-1) = 1/(i+j-1);
    }
  }

  // b
  vector<double> x_tilde(5, 1);
  vector<double> b(5);
  b = MatrixVector(Hilbert, x_tilde);
  cout << "b = \n";
  printVector(b);
  cout << "----------\n\n";
  vector<double> x;
  FullPivLU<MarixXd> lu(Hilbert);
  VectorXD x;

  // LU_compostion
  /* vector<double> x(n); */
  /* vector<double> y(n); */
  /* /1* vector<vector<double>> Hilbert_Forward; *1/ */
  /* /1* Hilbert_Forward = Forward_erase(Hilbert); *1/ */
  /* vector<vector<double>> Hilbert_LU = LU_decomposition(Hilbert); */
  /* vector<vector<double>> L = Hilbert_LU; */
  /* vector<vector<double>> U = Hilbert_LU; */
  /* for (int i = 0; i < n; i++) { */
  /*   for (int j = 0; j < n; j++) { */
  /*     if (i == j) { */
  /*       L.at(i).at(j) = 1; */
  /*     } else if (i < j) { */
  /*       L.at(i).at(j) = 0; */
  /*     } else if (i > j) { */
  /*       U.at(i).at(j) = 0; */
  /*     } */
  /*   } */
  /* } */

  /* y = Forward_sub(Hilbert_LU, b); */
  /* x = Backward_sub(Hilbert_LU, y); */
  /*   // printVector(x); */
  cout << "x = \n";
  printVector_more_detail(x);
  cout << "----------\n\n";

  double b_Ax = VectorNormInfty(VectorSubstract(b, MatrixVector(Hilbert, x)));
  printf("|b - Ax|_inf = %.2e\n", b_Ax);
  cout << "----------\n\n";

  double x_tilde_x = VectorNormInfty(VectorSubstract(x_tilde, x));
  printf("|x* - x|_inf = %.2e\n", x_tilde_x);
  cout << "----------\n\n";


  /* // Hilbert */
  /* cout << "Hilbert = \n"; */
  /* printMatrix(Hilbert); */
  /* cout << "----------\n\n"; */

  /* // Inverse_matrix */
  /* vector<vector<double>> Hilbert_Inverse; */
  /* Hilbert_Inverse = Inverse_matrix(Hilbert); */
  /* vector<vector<double>> solve = MatrixMatrix(Hilbert, Hilbert_Inverse); */
  /* cout << "Hilbert_Inverse = \n"; */
  /* printMatrix(Hilbert_Inverse); */
  /* cout << "----------\n\n"; */

  /* cout << "Hilbert*Hilbert_Inverse = \n"; */
  /* printMatrix(solve); */
  /* cout << "----------\n\n"; */

  /* // NormInf */
  /* double Hilbert_Norm_inf = MatrixNormInfty(Hilbert); */
  /* printf("Hilbert_Norm_inf = %.2e\n\n",Hilbert_Norm_inf); */
  /* cout << "----------\n\n"; */
  /* double Hilbert_Inverse_Norm_inf = MatrixNormInfty(Hilbert_Inverse); */
  /* printf("Hilbert_Inverse_Norm_inf = %.2e\n\n",Hilbert_Inverse_Norm_inf); */
  /* cout << "----------\n\n"; */

  /* // Kappa */
  /* cout << "Kappa = \n"; */
  /* double kappa = MatrixNormInfty(Hilbert) * MatrixNormInfty(Hilbert_Inverse); */
  /* printf("%.2e\n\n",kappa); */
  /* cout << "----------\n\n"; */






  /* cout << "------------------\n Forward_erase" << endl; */
  /* vector<vector<double>> Hilbert_Forward = Forward_erase(Hilbert, b_1); */
  /* printMatrix(Hilbert_Forward); */
  /* cout << "------------------\n Backward_sub" << endl; */
  /* vector<double> x_1 = Backward_sub(Hilbert_Forward, b_1); */
  /* printVector_more_detail(x_1); */
  /* cout << "------------------\n ResidualError" << endl; */
  /* vector<double> Hilbert_Ax_b = ResidualError(Hilbert, x_1, b_Hilbert); */
  /* cout << "Norm1 = " << VectorNorm1(Hilbert_Ax_b) << endl; */
  /* cout << "NormInf = " << VectorNormInfty(Hilbert_Ax_b) << endl; */
  /* cout << "------------------\n Error" << endl; */
  /* vector<double> Hilbert_Error= VectorSubstract(x_1, x_tilde); */
  /* cout << "Norm1 = " << VectorNorm1(Hilbert_Error) << endl; */
  /* cout << "NormInf = " << VectorNormInfty(Hilbert_Error) << endl; */
  /* cout << "------------------\n" << endl; */

  return 0;
}
