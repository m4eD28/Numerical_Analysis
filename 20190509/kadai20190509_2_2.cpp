#include <cstdio>
#include <cmath>

using namespace std;

double f(double x){
    return pow(x, 2) - 8/x;
}

double df(double x){
    return 2*x + 8*pow(x,2);
}

const int Max = 60;
int N = 0;
double x0 = 10, x1 = 5;
double eps = 1e-7;
double X_n[Max + 1];


void Secant_zansa(){
    printf("Secant法 (収束判定 残差)\n");
    for(int k = 1; k < Max; k++){
        if( fabs( f(X_n[k]) ) < eps ){
            printf("反復回数 : %d\n", N);
            printf("近似解 : %.12e\n", X_n[k]);

            return ;
        }else{
            X_n[k+1] = X_n[k] - f(X_n[k]) * ( (X_n[k] - X_n[k-1]) / (f(X_n[k]) - f(X_n[k-1])) );
            N++;
        }
    }

    printf("収束しない\n");

    return ;
}


int main(){
    X_n[0] = x0;
    X_n[1] = x1;
    Secant_zansa();
    for(int k = 1; k <= N; k++){
        printf("%d \t %.12e \t %.2e\n", k, X_n[k+1], fabs(X_n[k+1] - X_n[N+1]));
    }

    return 0;
}
