// (2) Secant法

#include <cstdio>
#include <cmath>

using namespace std;

double f(double x){
  /* return pow(x,3) - 8; */
  /* return pow(x,2) -8/x; */
  return pow(x,4) + 5*pow(x,3) + 6*pow(x,2) - 4*x -8;
}

const int Max = 60;
int N = 0;
double x0 = -30.0, x1 = -20.0;
double eps = 1e-7;
double X[Max + 2];
double alpha = -2.0;


void Secant_zansa(){
    printf("Secant法 (収束判定 残差)\n");
    for(int k = 1; k < Max; k++){
        if( fabs( f(X[k]) ) < eps ){
            printf("反復回数 N = %d\n", N);
            printf("近似解 X[N+1] = %.12e\n", X[k]);

            return ;
        }else{
            X[k+1] = X[k] - f(X[k]) * ( (X[k]-X[k-1]) / (f(X[k])-f(X[k-1])) );
            N++;
        }
    }

    printf("収束しない\n");

}


int main(){
    X[0] = x0;
    X[1] = x1;
    Secant_zansa();

    printf("#--------------------------------\n");
    printf("反復回数 k\t 近似解 X[k+1]\t\t 誤差 |X[k+1] - X[N+1]|\t\t 収束次数\n");
    for(int k = 1; k <= N; k++){
        printf("%d \t\t %.12e\t %.2e\t %.2e\n",
        k, X[k+1], fabs(X[k+1] - X[N+1]), log10(fabs(X[k+1] - alpha)) / log10(fabs(X[k] - alpha)) );
    }

    printf("誤差の履歴 |X[i] - alpha|\n");
    for(int i = 1; i <= N; i++){
        printf("%2d\t%.2e\n", i, fabs(X[i+1] - alpha));
    }

    printf("収束次数 %.2e\n", log10(fabs(X[N+1] - alpha)) / log10(fabs(X[N] - alpha)));

    return 0;
}
