#include <iostream>
#include <cmath>

using namespace std;

double f(double x){
    return pow(x, 3) - 8;
}

double df(double x){
    return 3*pow(x, 2);
}

const int Max = 60;
int N = 0;
double x0 = 10;
double eps = 1e-7;
double X_n[Max + 1];


void Newton_zansa(){
    printf("Newton法 (収束判定 残差)\n");
    for(int k = 0; k < Max; k++){
        if( fabs( f(X_n[k]) ) < eps ){
            printf("反復回数 : %d\n", N);
            printf("近似解 : %.16e\n", X_n[k]);

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
    cout << "------------" <<
      "反復回数 k \t 近似解 X[k] \t\t 誤差 |X[k] - X[N]| \t\t 収束次数\n" ;
    for(int k = 1; k <= N; k++){
        printf("%d \t %.16e \t %.16e \t %.16e\n", k, X_n[k], fabs(X_n[k] - X_n[N]), log(fabs(X_n[k]-2)) / log(fabs(X_n[k-1]-2)) );
    }
    cout << "誤差 |X[i] - 2.0|\n";
    for (int i = 0; i <= N; i++) {
      printf("%2d \t %.16e\n", i, fabs(X_n[i] - 2.0));
    }

    printf("収束次数 %.16e\n", log(fabs(X_n[N] - 2.0) / log10(fabs(X_n[N-1] - 2.0))));

    return 0;
}
