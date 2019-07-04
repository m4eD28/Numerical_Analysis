#include <iostream>
#include <cmath>

using namespace std;

double f(double x){
  /* return pow(M_E, (pow(x, 4) - 2 * pow(x, 2) + 1)) - 1; */
  /* return pow(x, 4) - 2*pow(x, 2) + 1; */
  return pow(x, 2) - 1;
  /* return 1.0 - 1.0/pow(x, 2); */
}

double df(double x){
  /* return (4*pow(x, 3) - 4*x)*pow(M_E, (pow(x, 4) - 2 * pow(x, 2) + 1)); */
  /* return 4*pow(x, 3) - 4*x; */
  return 2*x;
  /* return -2.0/pow(x, 3); */
}

const int Max = 100;
int N = 0;
double x0 = 3.0;
double eps = 1e-8;
double X[Max + 1];
double alpha = 1.0;


void Newton_zansa(){
    printf("Newton法 (収束判定 残差)\n");
    for(int k = 0; k < Max; k++){
        if( fabs( f(X[k]) ) < eps ){
            printf("反復回数 N = %d\n", N);
            printf("近似解 X[N] = %.12e\n", X[k]);

            return ;
        }else{
            X[k + 1] = X[k] - (f(X[k]) / df(X[k]));
            N++;
        }
    }

    printf("収束しない\n");
}


int main(){
    X[0] = x0;
    Newton_zansa();

    printf("#--------------------------------\n");
    printf("反復回数 k\t 近似解 X[k]\t\t 誤差 |X[k] - X[N]|\t\t 収束次数\n");
    for(int k = 1; k <= N; k++){
        printf("%d \t\t %.12e \t %.2e \t %.2e\n",
        k, X[k], fabs(X[k] - X[N]), log10(fabs(X[k] - alpha)) / log10(fabs(X[k-1] - alpha)) );
    }

    printf("誤差の履歴 |X[i] - alpha|\n");
    for(int i = 1; i <= N; i++){
        printf("%2d\t%.2e\n", i, fabs(X[i] - alpha));
    }

    printf("収束次数 %.2lf\n", log10(fabs(X[N] - alpha)) / log10(fabs(X[N-1] - alpha)));

    return 0;
}
