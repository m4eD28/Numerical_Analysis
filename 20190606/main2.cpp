#include "LinearAlgebra.hpp"
#include <iostream>

/* void Gauss_elimination(const vector<vector<double>>& _A, const vector<double>& _b) { */
/*   vector<double> alpha; */
/*   vector<vector<double>> A(_A); */
/*   vector<double> b(_b); */
/*   vector<double> x(A.size()); */
/*   for (int k = 0; k < A.size()-1; k++) { */
/*     for (int i = k+1; i < A.size(); i++) { */
/*       double alpha = A[i][k] / A[k][k]; */
/*       for (int j = k+1; j < A.size(); j++) { */
/*         A[i][j] = A[i][j] - alpha*A[k][j]; */
/*       } */
/*       b[i] = b[i] - alpha*b[k]; */
/*     } */
/*   } */

/*   for (int k = A.size()-1; k >= 0; k--) { */
/*     double sum = 0; */
/*     for (int j = k+1; j < A.size(); j++) { */
/*       sum += A[k][j] * x[j]; */
/*     } */
/*     x[k] = 1.0/A[k][k] * (b[k] - sum); */
/*   } */

/*   cout << "A = " << endl; */
/*   printMatrix(A); */
/*   cout << "b = " << endl; */
/*   printVector(b); */
/*   cout << "x = " << endl; */
/*   printVector(x); */
/* } */

int main() {
  vector<vector<double>> A;
  vector<double> b;
  A = {{3, -1, -1}, {-1, 3, -1}, {-1, -1, 3}};
  b = {1, 1, 1};
  /* printMatrix(A); */
  Gauss_elimination(A, b);


}
