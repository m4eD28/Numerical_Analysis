#include <cstdio>
#include <cmath>

using namespace std;

double f(double x){
  return pow(x, 3) + pow(x, 2)*11.0/2 + 6*x - 9.0/2;
}

double df(double x){
  return 3*pow(x, 2) + 11.0*x + 6;
}

const int Max = 100;
int N = 0;
double x0 = -10;
double eps = 1e-8;
double X_n[Max + 1];
double alpha = -3;


void Newton_zansa(){
    printf("Newton法 (収束判定 残差)\n");
    for(int k = 0; k < Max; k++){
        if( fabs( f(X_n[k]) ) < eps ){
            printf("反復回数 : %d\n", N);
            printf("近似解 : %.10e\n", X_n[k]);

            return ;
        }else{
            X_n[k + 1] = X_n[k] - (f(X_n[k]) / df(X_n[k]));
            N++;
        }
    }

    printf("収束しない\n");

    return ;
}


int main(){
    X_n[0] = x0;
    Newton_zansa();
    printf("反復回数 \t 解 \t 絶対誤差 \t 対数誤差\n");
    for(int k = 1; k <= N; k++){
        printf("%d \t %.10e \t %.2e \t %.2e\n", k, X_n[k], fabs(X_n[k] - alpha), log(fabs(X_n[k]-alpha)) / log(fabs(X_n[k-1]-alpha)));
    }

    return 0;
}
