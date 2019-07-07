#include "../src/LinearAlgebra.hpp"
#include "../src/Algo.hpp"
#include <iostream>

int Jacobi_law(double a) {
  int n = 20;
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
  // A
  /* cout << "A = \n"; */
  /* printMatrix(A); */
  /* cout << "----------\n\n"; */

  // b
  std::vector<double> b(n, 0);

  Jacobi_law(Hilbert, b);
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

  Jacobi_law(1.0);
  Jacobi_law(0.1);



  return 0;
}

