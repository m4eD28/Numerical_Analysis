#include <cstdio>
#include <cmath>

using namespace std;

double f(double x){
    return pow(x, 4) + 5*pow(x, 3) + 6*pow(x, 2) - 4*x - 8;
    /* return pow(M_E, x-1) - x; */
}

double df(double x){
    return 4*pow(x, 3) + 15*pow(x, 2) + 12*x - 4;
    /* return (x-1)*pow(M_E, x-1) - 1; */
}

const int Max = 60;
int N = 0;
double x0 = 30;
double eps = 1e-7;
double X_n[Max + 1];
double alpha = 1.0;


void Newton_zansa(){
    printf("Newton法 (収束判定 残差)\n");
    for(int k = 0; k < Max; k++){
        if( fabs( f(X_n[k]) ) < eps ){
            printf("反復回数 : %d\n", N);
            printf("近似解 : %.12e\n", X_n[k]);

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
        printf("%d \t %.12e \t %.2e \t %.2e\n", k, X_n[k], fabs(X_n[k] - alpha), log(fabs(X_n[k]-alpha)) / log(fabs(X_n[k-1]-alpha)));
    }

    return 0;
}
