#include "LinearAlgebra.hpp"
#include <iostream>

int Jacobi_law(double a) {
  int n = 20;
  printf("a = %f\n", a);
  std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
  for (int i = 1; i <= A.size(); i++) {
    for (int j = 1; j <= A.at(0).size(); j++) {
      if (0 < fabs(i - j) and fabs(i - j) <= 2) {
        A.at(i-1).at(j-1) = -1.0;
      }
      if (i == j) {
        A.at(i-1).at(j-1) = a;
      }
    }
  }

  // A
  /* cout << "A = \n"; */
  /* printMatrix(A); */
  /* cout << "----------\n\n"; */

  // b
  std::vector<double> b(n, 1.0);

  Jacobi_law(A, b);
  /* Gauss_Seidel_law(A, b); */
  std::cout << "-----------------------" << std::endl;
  return 0;
}

int main() {
  /* int n = 20; */
  /* double a = 6.0; */
  /* vector<vector<double>> A(n, vector<double>(n, 0)); */
  /* for (int i = 1; i <= A.size(); i++) { */
  /*   for (int j = 1; j <= A.at(0).size(); j++) { */
  /*     if (0 < fabs(i - j) and fabs(i - j) <= 2) { */
  /*       A.at(i-1).at(j-1) = -1.0; */
  /*     } */
  /*     if (i == j) { */
  /*       A.at(i-1).at(j-1) = a; */
  /*     } */
  /*   } */
  /* } */

  Jacobi_law(28.0);



  return 0;
}

