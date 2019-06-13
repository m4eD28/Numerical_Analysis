#ifndef NEWTON_H
#define NEWTON_H

#include <iostream>
#include <cmath>

using namespace std;

double f(double x){}

double df(double x){}

const int Max = 60;
double eps = 1e-7;
double X[Max + 1];


void Newton_zansa(){
  int N = 0;
  printf("Newton法 (収束判定 残差)\n");
  for(int k = 0; k < Max; k++){
    if( fabs( f(X[k]) ) < eps ){
      printf("反復回数 N = %d\n", N);
      printf("近似解 X[N] = %.16lf\n", X[k]);

      return ;
    }else{
      X[k + 1] = X[k] - (f(X[k]) / df(X[k]));
      N++;
    }
  }

  printf("収束しない\n");
}


void Newton(double x0, double alpha){
  X[0] = x0;
  int N;
  Newton_zansa();

  printf("#--------------------------------\n");
  printf("反復回数 k\t 近似解 X[k]\t\t 誤差 |X[k] - X[N]|\t\t 収束次数\n");
  for(int k = 1; k <= N; k++){
    printf("%d \t\t %.16lf \t %.16e \t %.16lf\n",
    k, X[k], fabs(X[k] - X[N]), log10(fabs(X[k] - alpha)) / log10(fabs(X[k-1] - alpha)) );
  }

  printf("誤差の履歴 |X[i] - alpha|\n");
  for(int i = 1; i <= N; i++){
    printf("%2d\t%.16e\n", i, fabs(X[i] - alpha));
  }

  printf("収束次数 %.16lf\n", log10(fabs(X[N] - alpha)) / log10(fabs(X[N-1] - alpha)));

}
#endif
