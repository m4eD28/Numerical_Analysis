#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include "../src/LinearAlgebra.hpp"
#include "../src/Algo.hpp"

#define dim 50

const int Max = 100;
int N;
const double error_eps = 1e-8;
std::vector<std::vector<double> > x(Max, std::vector<double>(dim));
std::vector<double> delta(dim);
std::vector<std::vector<double> > jacobi(dim, std::vector<double>(dim, 0));
std::vector<double> F(dim);
std::vector<double> x_alpha(dim, 1);

std::vector<std::vector<double> > Generate_Jacobi(std::vector<double>& x) {
  for (int i = 1; i < dim-1; ++i) {
    for (int j = 0; j < dim; ++j) {
      if (i == j) {
        jacobi.at(i).at(j) = 1200*pow(x.at(j), 2) - 400*x.at(j+1) + 202;
      } else if (i - j == 1) {
        jacobi.at(i).at(j) = -400 * x.at(j);
      } else if (j - i == 1) {
        jacobi.at(i).at(j) = -400 * x.at(j-1);
      }
    }
  }
  jacobi.front().front() = 1200*pow(x.at(0), 2) - 400*x.at(1) + 2;
  jacobi.front().at(1) = -400*x.at(0);
  jacobi.back().back() = 200;
  jacobi.back().at(dim-2) = -400*x.at(dim-2);
  return jacobi;
}

std::vector<double> Generate_negative_F(std::vector<double>& x) {
  for (int i = 1; i < dim-1; ++i) {
    F.at(i) = -(400*pow(x.at(i), 3) + (202 - 400*x.at(i+1))*x.at(i) - 200*pow(x.at(i-1), 2) - 2);
  }
  F.front() = -(400*pow(x.at(0), 3) + (2 - 400*x.at(1))*x.at(0) - 2);
  F.back() = -(200*(x.back() - pow(x.at(dim-2), 2)));
  return F;
}

void Newton_convergence_order(std::vector<double>& initial, std::vector<double>& x_alpha) {
  double error = std::numeric_limits<double>::infinity();
  double error_true;
  double convergence_order;
  printf("Newton法 (収束判定 残差) 収束次数\n");
  x.at(0) = initial;
  N = 0;
  for(int k = 0; k < Max; ++k){
    delta = Gauss_elimination_pivot(Generate_Jacobi(x.at(k)), Generate_negative_F(x.at(k)));
    x.at(k+1) = VectorAdd(x.at(k), delta);
    N++;
    error = VectorNormInfty(Generate_negative_F(x.at(k)));
    if(error < error_eps ){
      for (int i = k-2; i <= k; ++i) {
        printf("N = %d\n", i+1);
        printf("近似解 x[n]\n");
        printVector_more_detail(x.at(i));
        error = VectorNormInfty(Generate_negative_F(x.at(i)));
        printf("error = %.10e\n", error);
        error_true = VectorNormInfty(VectorSubstract(x.at(i), x_alpha));
        printf("誤差 ||x[n] - x|| = %.2e\n", error_true);
        convergence_order = log(VectorNormInfty(VectorSubstract(x.at(i), x_alpha))) / log(VectorNormInfty((VectorSubstract(x.at(i-1), x_alpha))));
        printf("収束次数 : %.2e\n", convergence_order);
      }
      return;
    }
  }
  printf("収束しない\n");
}

void Newton(std::vector<double>& initial) {
  double error = std::numeric_limits<double>::infinity();
  printf("Newton法 (収束判定 残差)\n");
  x.at(0) = initial;
  N = 0;
  for(int k = 0; k < Max; ++k){
    if(error < error_eps ){
      printf("反復回数 n = %d\n", N);
      printf("近似解 x[n]\n");
      printVector_more_detail(x.at(k));
      return;
    }
    std::vector<std::vector<double> > test;
    test = Generate_Jacobi(x.at(k));
    printMatrix(test);
    delta = Gauss_elimination_pivot(Generate_Jacobi(x.at(k)), Generate_negative_F(x.at(k)));
    x.at(k+1) = VectorAdd(x.at(k), delta);
    N++;
    error = VectorNormInfty(VectorSubstract(x.at(k+1), x.at(k)));
  }

  printf("収束しない\n");
}

int main(){
  std::vector<double> initial(dim, 1);
  for (int i = 0; i < dim; ++i) {
    if (i % 2 == 0) {
      initial.at(i) = -1.2;
    }
  }
  Newton_convergence_order(initial, x_alpha);
  Newton(initial);

  return 0;
}
