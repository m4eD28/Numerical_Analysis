#include <iostream>
#include <cmath>
using namespace std;

const double heps = 1e-15;
const int Max = 50;
double x0 = 1.0;
double x1 = 1.5;
double eps =1e-12;
double X_n[Max + 1];
int N = 0;

double f(double x) {
  return 10 * sin(x) + pow(M_E, x);
}

double df(double x) {
  /* return (f(x + heps) - f(x)) / heps; */
  return 10 * cos(x) + pow(M_E, x);
}

int error(double *X_n) {
  for (int k = 1; k < N; k++) {
    printf("%d \t %.12e \t %.2e\n",k,X_n[k],fabs(X_n[k] - X_n[N]));
  }
  return 1;
}

int Newton(double *X_n) {
  N = 0;
  for (int k = 0; k < Max; k++) {
    X_n[k+1] = X_n[k] - f(X_n[k]) / df(X_n[k]);
    N++;
    if(fabs((X_n[k+1] - X_n[k]) / X_n[k+1]) <= eps){
      printf("Newton法\nx = %.12e\nCount = %d\n",X_n[k+1],N-1);
      error(X_n);
      printf("----------\n");
      return 0;
    }
  }
  printf("Newton法\n収束しない\n---------------\n");
  return 1;
}

int Secant(double *X_n) {
  N = 0;
  for (int k = 1; k < Max; k++) {
    X_n[k+1] = X_n[k] - f(X_n[k]) * ((X_n[k] -X_n[k-1]) / (f(X_n[k]) - f(X_n[k-1])));
    N++;
    if(fabs((X_n[k+1] - X_n[k]) / X_n[k+1]) < eps){
      printf("Secant法\nx = %.12e\nCount = %d\n",X_n[k+1],N);
      error(X_n);
      printf("----------\n");
      return 1;
    }
  }
  printf("Secant法\n収束しない\n---------------\n");
  return 0;
}

int main() {
  X_n[0] = x0;
  Newton(X_n);
  for (int i = 0; i < Max; i++) {
    X_n[i] = 0;
  }
  X_n[0] = x0;
  X_n[1] = x1;
  Secant(X_n);
}
