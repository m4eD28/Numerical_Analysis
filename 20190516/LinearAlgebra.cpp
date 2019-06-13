#include "LinearAlgebra.hpp"
#include <vector>
#include <cmath>
using namespace std;

int main() {
  int n = 3;
  vector<double> a(n);
  vector<vector<double>> A(n, vector<double>(n));
  for (int i = 0; i < n; i++) {
    a[i] = i;
    for (int j = 0; j < n; j++) {
      A[i][j] = i+j;
    }
  }

  cout << "a = " << endl;
  printVector(a);
  cout << "A = " << endl;
  printMatrix(A);

  vector<double> b = a;
  b[0] = 1000;

  vector<vector<double>> B = A;
  B[0][0] = 1000;

  cout << "b = " << endl;
  printVector(b);
  cout << "B = " << endl;
  printMatrix(B);

  cout << "a = " << endl;
  printVector(a);
  cout << "A = " << endl;
  printMatrix(A);

  vector<double> c = VectorSubstract(a, b);
  cout << "c = a-b = " << endl;
  printVector(c);

  c = MatrixVector(A, a);
  cout << "c = A*a = " << endl;
  printVector(c);

  vector<double> d = ResidualError(A, a, b);
  cout << "d = A*a - b = " << endl;
  printVector(d);

  printf("||a||_1 = %.2e\n", VectorNorm1(a));
  printf("||a||_2 = %.2e\n", VectorNorm2(a));
  printf("||a||_infty = %.2e\n", VectorNormInfty(a));

  A[0][1] = 100;
  printf("||A||_1 = %.2e\n", MatrixNorm1(A));
  printf("||A||_infty = %.2e\n", MatrixNormInfty(A));

  return 0;
}
