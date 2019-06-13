#include "LinearAlgebra.hpp"
#include <iostream>

int main() {
  vector<vector<double>> A;
  vector<double> b;
  A = {{3, -1, -1}, {-1, 3, -1}, {-1, -1, 3}};
  b = {1, 1, 1};
  Gauss_elimination(A, b);
}
