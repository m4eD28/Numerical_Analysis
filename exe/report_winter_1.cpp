#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include "../src/LinearAlgebra.hpp"
#include "../src/Algo.hpp"

#define dim 2

const int Max = 100;
int N;
const double error_eps = 1e-8;
std::vector<std::vector<double> > x(Max, std::vector<double>(dim));
std::vector<double> delta(dim);
std::vector<std::vector<double> > jacobi(dim, std::vector<double>(dim, 0));
std::vector<double> F(dim);
std::vector<double> x_alpha(dim, 0);

double func(std::vector<double>& x) {
  return pow((pow(x.at(0), 2) + x.at(1) -11), 2) + pow((x.at(0) + pow(x.at(1), 2) - 7), 2);
      }

std::vector<std::vector<double> > Generate_Jacobi(std::vector<double>& x) {
  jacobi = {{2*(-21 + 6*pow(x.at(0), 2) + 2*x.at(1)), 4*(x.at(0) + x.at(1))},
            {4*(x.at(0) + x.at(1)), 2*(-13 + 2*x.at(0) + 6*pow(x.at(1), 2))}};
  return jacobi;
}

std::vector<double> Generate_negative_F(std::vector<double>& x) {
  F = {-2*(-7 - 21*x.at(0) + 2*pow(x.at(0), 3) + 2*x.at(0)*x.at(1) + pow(x.at(1), 2)), -2*(-11 + pow(x.at(0), 2) - 13*x.at(1) + 2*x.at(0)*x.at(1) + 2*pow(x.at(1), 3))};
  return F;
}

void Newton(std::vector<double>& initial) {
  double error = std::numeric_limits<double>::infinity();
  printf("-----Newton法-----\n");
  printf("初期値:\n");
  printVector(initial);
  x.at(0) = initial;
  N = 0;
  for(int k = 0; k < Max; ++k){
    delta = Gauss_elimination_pivot(Generate_Jacobi(x.at(k)), Generate_negative_F(x.at(k)));
    x.at(k+1) = VectorAdd(x.at(k), delta);
    N++;
    error = VectorNormInfty(Generate_negative_F(x.at(k)));
    if(error < error_eps ){
      printf("反復回数 n = %d\n", N);
      printf("関数値 f(x) = %.10e\n", func(x.at(k)));
      printf("近似解 x[n]\n");
      printVector_more_detail(x.at(k));
      printf("勾配 g(x)\n");
      printVector_more_detail(Generate_negative_F(x.at(k)));
      printf("ヘッセ行列  ∇2f(x)\n");
      printMatrix_more_detail(Generate_Jacobi(x.at(k)));
      return;
    }
  }

  printf("収束しない\n");
}

void BFGS(std::vector<double>& initial) {
  double error = std::numeric_limits<double>::infinity();
  printf("-----BFGS公式-----\n");
  printf("初期値:\n");
  printVector(initial);
  x.at(0) = initial;
  N = 1;
  std::vector<std::vector<std::vector<double> > > H(Max, std::vector<std::vector<double> >(dim, std::vector<double>(dim, 0)));
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      if (i == j) {
        H.at(0).at(i).at(j) = 1;
      }
    }
  }
  delta = MatrixVector(H.at(0), Generate_negative_F(x.at(0)));
  x.at(1) = VectorAdd(x.at(0), delta);

  // H


  // delta
  std::vector<double> Hf(dim);
  std::vector<double> Hy(dim);
  std::vector<double> y(dim);
  std::vector<double> s(dim);
  double sy;
  double sf;
  double scalar;

  for(int k = 1; k < Max; ++k){
    Hf = MatrixVector(H.at(k-1), Generate_negative_F(x.at(k)));
    y = VectorScalar(-1, VectorSubstract(Generate_negative_F(x.at(k)), Generate_negative_F(x.at(k-1))));
    Hy = MatrixVector(H.at(k-1), y);
    s = VectorSubstract(x.at(k), x.at(k-1));
    sy = VectorToScalar(s, y);
    sf = VectorToScalar(s, VectorScalar(-1, Generate_negative_F(x.at(k))));
    scalar = VectorToScalar(Hy, VectorScalar(-1, Generate_negative_F(x.at(k)))) / sy - (1 + VectorToScalar(y, Hy) / sy) * sf / sy;
    delta = VectorAdd(VectorAdd(Hf, VectorScalar(scalar, s)), VectorScalar(sf/sy, Hy));

    x.at(k+1) = VectorAdd(x.at(k), delta);

    /* H.at(k) = */ 
    N++;
    error = VectorNormInfty(Generate_negative_F(x.at(k)));

    if(error < error_eps ){
      printf("反復回数 n = %d\n", N);
      printf("関数値 f(x) = %.10e\n", func(x.at(k)));
      printf("近似解 x[n]\n");
      printVector_more_detail(x.at(k));
      printf("勾配 g(x)\n");
      printVector_more_detail(Generate_negative_F(x.at(k)));
      printf("ヘッセ行列  ∇2f(x)\n");
      printMatrix_more_detail(Generate_Jacobi(x.at(k)));
      return;
    }
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

int main(){
  std::vector<double> initial(dim);
  std::vector<double> alpha(dim);

  /* // 初期値(0, 0) */
  /* initial = {0, 0}; */
  /* Newton(initial); */

  /* // 初期値(-3, -3) */
  /* initial = {-3, -3}; */
  /* Newton(initial); */

  /* // 初期値(-3, 3) */
  /* initial = {-3, 3}; */
  /* Newton(initial); */

  /* // 初期値(5, 5) */
  /* initial = {5, 5}; */
  /* Newton(initial); */

  initial = {-3, 3};
  alpha = {-2.805118, 3.131312};
  Newton_convergence_order(initial, alpha);


  //BFGS
  /* BFGS(initial); */
  return 0;
}


