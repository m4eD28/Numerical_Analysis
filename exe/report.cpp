#include "../src/LinearAlgebra.hpp"
#include "../src/Algo.hpp"
#include <iostream>


std::vector<std::vector<double> > Generate_Matrix(int N, double d) {
  std::vector<std::vector<double> > A(N, std::vector<double>(N, 0));
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      if (i == j) {
        A.at(i-1).at(j-1) = (100 * 100 + d * (i - 1)) * i;
      } else {
        A.at(i-1).at(j-1) = (j + d * (i - 1)) * i;
      }
    }
  }

  return A;
}

void Jacobi(int N, double d) {
  std::vector<std::vector<double>> A = Generate_Matrix(N, d);
  // A
  /* cout << "A = \n"; */
  /* printMatrix(A); */
  /* cout << "----------\n\n"; */

  // b
  std::vector<double> b(N, 1.0);

  std::vector<double> x = Jacobi_law(A, b);
  printVector(x);
  std::cout << "-----------------------" << std::endl;
}

void Gauss_Seidel(int N, double d) {
  std::vector<std::vector<double>> A = Generate_Matrix(N, d);
  // A
  /* cout << "A = \n"; */
  /* printMatrix(A); */
  /* cout << "----------\n\n"; */

  // b
  std::vector<double> b(N, 1.0);

  std::vector<double> x = Gauss_Seidel_law(A, b);
  printVector(x);
  std::cout << "-----------------------" << std::endl;
}

void kadai(double a) {
  std::cout << "a = " << a << std::endl;
  int N = 20;
  std::vector<std::vector<double> > A = Generate_Matrix(N, a);
  std::vector<double> b(N, 1.0);
  std::vector<double> x = Jacobi_law(A, b);
  printf("x_1 = %.2e\n", x.at(0));
  std::cout << "----------" << std::endl;
}

int main() {
  Jacobi(10, 1);
  Gauss_Seidel(10, 1);

  return 0;
}

