#include <iostream>
#include <cmath>
using namespace std;

const double heps = 1e-10;

double f(double x) {
  return ;
}

double df(double x) {
  return (f(x + heps) - f(x)) / heps;
}

int main() {
  const int Max = 100;
  int N = 0;
  double x0 = ;
  double eps = ;
  double X_n[Max + 1];
  X_n[0] = x0;
  for (int k = 0; k < Max; k++) {
    X_n[k+1] = X_n[k] - f(X_n[k]) / df(X_n[k]);
    N++;
    if(fabs((X_n[k+1] - X_n[k]) / X_n[k+1])){
      cout << X_n[k+1] << endl;
      break;
    }
  }
}
