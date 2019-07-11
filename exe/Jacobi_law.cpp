#include "../src/LinearAlgebra.hpp"
#include "../src/Algo.hpp"
#include <iostream>

/* int Jacobi_law(double a) { */
/*   int n = 20; */
/*   std::vector<std::vector<double>> Hilbert(n, std::vector<double>(n)); */
/*   for (double i = 1; i <= n; i++) { */
/*     for (double j = 1; j <= n; j++) { */
/*       if(i == j) Hilbert.at(i-1).at(i-1) = a; */
/*       else if (0 < fabs(i-j) and fabs(i-j) < 4) */
/*         Hilbert.at(i-1).at(j-1) = 1 / (fabs(i*i - j*j + 15) + 1); */
/*       else { */
/*         Hilbert.at(i-1).at(j-1) = 0; */
/*       } */
/*     } */
/*   } */
/*   // A */
/*   /1* cout << "A = \n"; *1/ */
/*   /1* printMatrix(A); *1/ */
/*   /1* cout << "----------\n\n"; *1/ */

/*   // b */
/*   std::vector<double> b(n, 0); */

/*   Jacobi_law(Hilbert, b); */
/*   /1* Gauss_Seidel_law(A, b); *1/ */
/*   std::cout << "-----------------------" << std::endl; */
/*   return 0; */
/* } */
std::vector<std::vector<double> > Generate_Matrix(int N, double a) {
  std::vector<std::vector<double> > A(N, std::vector<double>(N, 0));
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      if (i == j) {
        A.at(i-1).at(j-1) = a;
      } else if (0 < fabs(i-j) and fabs(i-j) <= 2) {
        A.at(i-1).at(j-1) = -1;
      }
    }
  }

  return A;
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
  kadai(3.0);
  kadai(4.0);
  kadai(5.0);
  kadai(6.0);
  kadai(7.0);
  kadai(14.0);
  kadai(28.0);

  return 0;
}

