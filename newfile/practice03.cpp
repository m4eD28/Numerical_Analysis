/* // (2) Newton法 */

#include <iostream>
#include <cmath>

using namespace std;

double f(double x){
    /* return x*x-8/x; */
    /* return pow(x,3) - 8; */
  return pow(x,4) + 5*pow(x,3) + 6*pow(x,2) - 4*x -8;
}

double df(double x){
    /* return 2*x+8/(x*x); */
    /* return 3*pow(x,2); */
  return 4*pow(x,3) + 15*pow(x,2) + 12*x - 4;
}

const int Max = 60;
int N = 0;
double x0 = 30.0;
double eps = 1e-7;
double X[Max + 1];


void Newton_zansa(){
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


int main(){
    X[0] = x0;
    Newton_zansa();

    printf("#--------------------------------\n");
    printf("反復回数 k\t 近似解 X[k]\t\t 誤差 |X[k] - X[N]|\t\t 収束次数\n");
    for(int k = 1; k <= N; k++){
        printf("%d \t\t %.16lf \t %.16e \t %.16lf\n",
        k, X[k], fabs(X[k] - X[N]), log10(fabs(X[k] - 1.0)) / log10(fabs(X[k-1] - 1.0)) );
    }

    printf("誤差の履歴 |X[i] - 2.0|\n");
    for(int i = 1; i <= N; i++){
        printf("%2d\t%.16e\n", i, fabs(X[i] - -.0));
    }

    printf("収束次数 %.16lf\n", log10(fabs(X[N] - 1.0)) / log10(fabs(X[N-1] - 1.0)));

    return 0;
}