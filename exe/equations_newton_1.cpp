#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include "../src/LinearAlgebra.hpp"
#include "../src/Algo.hpp"

#define dim 2

const int Max = 100;
int N;
const double error_eps = 1e-10;
std::vector<std::vector<double> > x(Max, std::vector<double>(dim));
std::vector<double> delta(dim);
std::vector<std::vector<double> > jacobi(dim, std::vector<double>(dim));
std::vector<double> F(dim);
std::vector<double> x_alpha(dim);

std::vector<std::vector<double> > Generate_Jacobi(std::vector<double>& x) {
  jacobi = {{3 * pow(x.at(0), 2) - 2 * x.at(0) * x.at(1), - pow(x.at(0), 2) + 2 * x.at(1)},
            {2 * x.at(0) + x.at(1), x.at(0) - 2 * x.at(1)}};
  return jacobi;
}

std::vector<double> Generate_negative_F(std::vector<double>& x) {
  F = {-(pow(x.at(0), 3) - pow(x.at(0), 2) * x.at(1) + pow(x.at(1), 2) - 2),
       -(pow(x.at(0), 2) + x.at(0) * x.at(1) - pow(x.at(1), 2) - 2)};
  return F;
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
    delta = Gauss_elimination_pivot(Generate_Jacobi(x.at(k)), Generate_negative_F(x.at(k)));
    x.at(k+1) = VectorAdd(x.at(k), delta);
    N++;
    error = VectorNormInfty(VectorSubstract(x.at(k+1), x.at(k)));
  }

  printf("収束しない\n");
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
    error = VectorNormInfty(VectorSubstract(x.at(k+1), x.at(k)));
    if(error < error_eps ){
      for (int i = k-2; i <= k; ++i) {
        printf("N = %d\n", i+1);
        printf("近似解 x[n]\n");
        printVector_more_detail(x.at(i));
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


int main(){
  printf("-----(iii)-----\n");
  std::vector<double> initial1 = {1, 0};
  Newton(initial1);
  std::vector<double> initial2 = {1, 1};
  Newton(initial2);
  std::vector<double> initial3 = {0, -1};
  Newton(initial3);

  printf("-----(iv)-----\n");
  std::vector<double> initial = {1, 1};
  x_alpha = {sqrt(2), sqrt(2)};
  Newton_convergence_order(initial, x_alpha);

  return 0;
}
