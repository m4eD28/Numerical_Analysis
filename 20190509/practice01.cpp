//(1) Newton法

#include <cstdio>
#include <cmath>

using namespace std;

double f(double x){
    return pow(x, 3) - 1;
}

double df(double x){
    return 3 * pow(x, 2);
}

const int Max = 50;
int N = 0;
double x0 = 7.0;
double eps = 1e-7;
double X_n[Max + 1];


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
    for(int k = 1; k <= N; k++){
        printf("%d \t %.12e \t %.12e \t %.12e\n", k, X_n[k], fabs(X_n[k] - 1), log(fabs(X_n[k]-1)) / log(fabs(X_n[k-1]-1)));
    }

    return 0;
}
