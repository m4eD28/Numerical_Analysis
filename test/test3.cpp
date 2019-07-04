#include "LinearAlgebra.hpp"
#include <iostream>

int main() {
  int n = 5;
  std::vector<std::vector<double>> Hilbert(n, std::vector<double>(n));
  std::vector<std::vector<double>> Marix_5(n, std::vector<double>(n));
  std::vector<std::vector<double>> test(n, std::vector<double>(n));
  std::vector<double> b_1;
  std::vector<double> b_2;
  std::vector<double> b_3;
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
  for (double i = 1; i <= n; i++) {
    for (double j = 1; j <= n; j++) {
      test.at(i-1).at(j-1) = pow(j, i-1);
    }
  }
  for (int i = 0; i < n; i++) {
    x_tilde.emplace_back(1);
  }

  std::vector<std::vector<double>> test_Inverse;
  test_Inverse = Inverse_matrix(test);
  std::vector<std::vector<double>> solve = MatrixMatrix(test, test_Inverse);
  printMatrix(solve);

  std::cout << "test = " << std::endl;
  printMatrix(test);
  std::cout << "------------------\n b = " << std::endl;
  b_3 = MatrixVector(test , x_tilde);
  std::vector<double> b_test(b_3);
  printVector(b_3); //Hilbert*x_tilde
  b_3.at(4) = b_3.at(4) + 0.001*b_3.at(4);
  printVector(b_3); //Hilbert*x_tilde
  std::cout << "------------------\n Forward_erase" << std::endl;
  std::vector<std::vector<double>> test_forward = Forward_erase(test, b_3);
  printMatrix(test_forward);
  std::cout << "------------------\n Backward_sub" << std::endl;
  std::vector<double> x_3 = Backward_sub(test_forward, b_3);
  printVector(x_3);
  std::cout << "------------------\n ResidualError" << std::endl;
  std::vector<double> test_Ax_b = ResidualError(test, x_3, b_test);
  std::cout << "Norm1 = " << VectorNorm1(test_Ax_b) << std::endl;
  std::cout << "NormInf = " << VectorNormInfty(test_Ax_b) << std::endl;
  std::cout << "------------------\n Error" << std::endl;
  std::vector<double> test_error = VectorSubstract(x_3, x_tilde);
  std::cout << "Norm1 = " << VectorNorm1(test_error) << std::endl;
  std::cout << "NormInf = " << VectorNormInfty(test_error) << std::endl;
  std::cout << "------------------\n" << std::endl;

  double kappa = MatrixNormInfty(test) * MatrixNormInfty(test_Inverse);
  printf("kappa = %.2e\n", kappa);

  double relation_error = VectorNormInfty(test_error) / VectorNormInfty(x_tilde);
  printf("相対誤差 = %.2e\n", relation_error);

  std::vector<double> delta_b(5, 0);
  delta_b.at(4) = 0.001*b_test.at(4);
  double setudou = VectorNormInfty(delta_b)/VectorNormInfty(b_test);
  printf("摂動 = %.2e\n", setudou);



  /* std::cout << "Hilbert = " << std::endl; */
  /* printMatrix(Hilbert); */
  /* std::cout << "------------------\n b = " << std::endl; */
  /* b_1 = MatrixVector(Hilbert, x_tilde); */
  /* std::vector<double> b_Hilbert(b_1); */
  /* printVector(b_1); //Hilbert*x_tilde */
  /* std::cout << "------------------\n Forward_erase" << std::endl; */
  /* std::vector<std::vector<double>> Hilbert_Forward = Forward_erase(Hilbert, b_1); */
  /* printMatrix(Hilbert_Forward); */
  /* std::cout << "------------------\n Backward_sub" << std::endl; */
  /* std::vector<double> x_1 = Backward_sub(Hilbert_Forward, b_1); */
  /* printVector_more_detail(x_1); */
  /* std::cout << "------------------\n ResidualError" << std::endl; */
  /* std::vector<double> Hilbert_Ax_b = ResidualError(Hilbert, x_1, b_Hilbert); */
  /* std::cout << "Norm1 = " << VectorNorm1(Hilbert_Ax_b) << std::endl; */
  /* std::cout << "NormInf = " << VectorNormInfty(Hilbert_Ax_b) << std::endl; */
  /* std::cout << "------------------\n Error" << std::endl; */
  /* vector<double> Hilbert_Error= VectorSubstract(x_1, x_tilde); */
  /* std::cout << "Norm1 = " << VectorNorm1(Hilbert_Error) << std::endl; */
  /* std::cout << "NormInf = " << VectorNormInfty(Hilbert_Error) << std::endl; */
  /* std::cout << "------------------\n" << std::endl; */



  /* std::cout << "Marix_5 = " << std::endl; */
  /* printMatrix(Marix_5); */
  /* std::cout << "------------------\n b = " << std::endl; */
  /* b_2 = MatrixVector(Marix_5, x_tilde); */
  /* std::vector<double> b_Marix5(b_2); */
  /* printVector(b_2); //Marix_5*x_tilde */
  /* std::cout << "------------------\n Forward_erase" << std::endl; */
  /* std::vector<std::vector<double>> Marix_5_Forward = Forward_erase(Marix_5, b_2); */
  /* printMatrix(Marix_5_Forward); */
  /* std::cout << "------------------\n Backward_sub" << std::endl; */
  /* std::vector<double> x_2 = Backward_sub(Marix_5_Forward, b_2); */
  /* printVector_more_detail(x_2); */
  /* std::cout << "------------------\n ResidualError = " << std::endl; */
  /* std::vector<double> Marix_5_Ax_b = ResidualError(Marix_5, x_2, b_Marix5); */
  /* std::cout << "Norm1 = " << VectorNorm1(Marix_5_Ax_b) << std::endl; */
  /* std::cout << "NormInf = " << VectorNormInfty(Marix_5_Ax_b) << std::endl; */
  /* std::cout << "------------------\n Error" << std::endl; */
  /* std::vector<double> Marix_5_Error= VectorSubstract(x_2, x_tilde); */
  /* std::cout << "Norm1 = " << VectorNorm1(Marix_5_Error) << std::endl; */
  /* std::cout << "NormInf = " << VectorNormInfty(Marix_5_Error) << std::endl; */

  /* std::cout << "------------------\n Forward_erase" << std::endl; */
  /* std::vector<std::vector<double>> Forward  = Forward_erase(A, b); */
  /* printMatrix(Forward); */
  /* std::cout << "------------------\n Backward_sub" << std::endl; */
  /* std::vector<double> x_x = Backward_sub(Forward, b); */
  /* printVector(x_x); */
  /* std::cout << "------------------\n ResidualError" << std::endl; */
  /* std::vector<double> Hilbert_Ax_b = ResidualError(A, x_1, b_1); */
  /* std::cout << "Norm1 = \n" << VectorNorm1(Hilbert_Ax_b) << std::endl; */
  /* std::cout << "NormInf = \n" << VectorNormInfty(Hilbert_Ax_b) << std::endl; */
  /* std::cout << "------------------\n Error" << std::endl; */
  /* std::vector<double> Hilbert_Error= VectorSubstract(x_1, x_tilde); */
  /* std::cout << "Norm1 = \n" << VectorNorm1(Hilbert_Error) << std::endl; */
  /* std::cout << "NormInf = \n" << VectorNormInfty(Hilbert_Error) << std::endl; */
  /* std::cout << "------------------\n" << std::endl; */

  /* std::cout << "Hilbert = " << std::endl; */
  /* printMatrix(Hilbert); */
  /* std::cout << "------------------\n b = " << std::endl; */
  /* b_1 = MatrixVector(Hilbert, x_tilde); */
  /* printVector(b_1); //Hilbert*x_tilde */
  /* Gauss_elimination(Hilbert, b_1); */


  return 0;
}
