
#include "LinearAlgebra.hpp"
#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;


//ベクトルaを出力する
void printVector(vector<double> a){
    unsigned long n = a.size();
    for(int i = 0; i < n; i++){
        printf("%.10e\n",a[i]);
    }
    printf("\n");
}

//行列Aを出力する
void printMatrix(vector< vector<double> > A){
    unsigned long m = A.size();
    unsigned long n = A[0].size();
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            printf("%.2e\t",A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

//ベクトルa-bを計算し，結果を返す
vector<double> VectorSubtract(vector<double> a,vector<double> b){
    unsigned long n = a.size();
    vector<double> c(n);
    for(int i = 0; i < n; i++){
        c[i] = a[i] - b[i];
    }

    return c;
}

//行列Aとベクトルbの乗算を行い，結果を返す
vector<double> MatrixVector(vector< vector<double> > A, vector<double> b){
    unsigned long m = A.size();
    unsigned long n = A[0].size();
    vector<double> c(m);
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            c[i] += A[i][j] * b[j];
        }
    }

    return c;
}

//行列Aと行列Bの積を返す
vector< vector<double> > MatrixMatrix(vector< vector<double> > A, vector< vector<double> > B){
    unsigned long m = A.size();
    unsigned long n = A[0].size(); // == B.size();
    unsigned long l = B[0].size();
    vector< vector<double> > C(m, vector<double>(l));

    double sum;
    for(int i = 0; i < m; i++){
        for(int j = 0; j < l; j++){
            sum = 0;
            for(int k = 0; k < n; k++){
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }

    return C;
}

//A*x-bを返す
vector<double> ResidualError(vector< vector<double> > A, vector<double> x, vector<double> b){
    unsigned long m = A.size();
    vector<double> c(m);
    for(int i = 0; i < m; i++){
        c[i] = MatrixVector(A,x)[i] - b[i];
    }

    return c;
}



//ベクトルaの1ノルムを返す
double VectorNorm1(vector<double> a){
    double norm_1 = 0;
    unsigned long n = a.size();
    for(int i = 0; i < n; i++){
        norm_1 += fabs(a[i]);
    }

    return norm_1;
}


//ベクトルaの2ノルムを返す
double VectorNorm2(vector<double> a){
    double norm_2 = 0;
    unsigned long n = a.size();
    for(int i = 0; i < n; i++){
        norm_2 += fabs(a[i] * a[i]);
    }
    norm_2 = sqrt(norm_2);

    return norm_2;
}


//ベクトルaのInftyノルムを返す
double VectorNormInfty(vector<double> a){
    double norm_inf = 0;
    unsigned long n = a.size();
    for(int i = 0; i < n; i++){
        norm_inf = max(norm_inf, fabs(a[i]));
    }

    return norm_inf;
}


//行列Aの1ノルムを返す
double MatrixNorm1(vector< vector<double> > A){
    double m_norm = 0;
    double num;
    unsigned long m = A.size();
    unsigned long n = A[0].size();
    for(int i = 0; i < n; i++){
        num = 0;
        for(int j = 0; j < m; j++){
                num += fabs(A[j][i]);
        }
        m_norm = max(m_norm, num);
    }

    return m_norm;
}


//行列AのInftyノルムを返す
double MatrixNormInfty(vector< vector<double> > A){
    double n_norm_inf = 0;
    double num;
    unsigned long m = A.size();
    unsigned long n = A[0].size();
    for(int i = 0; i < m; i++){
        num = 0;
        for(int j = 0; j < n; j++){
            num += fabs(A[i][j]);
        }
        n_norm_inf = max(n_norm_inf, num);
    }

    return n_norm_inf;
}


////////////////////////

vector<double> GaussianElimination(vector< vector<double> > a, vector<double> c){ // a:m*n c:m*1
    unsigned long m = a.size(); // = c.size();
    unsigned long n = a[0].size();

    vector< vector<double> > A(m, vector<double>(n));
    vector<double> b(m);
    A = a; b = c;

    // 前進消去
    double alpha;
    for(int k = 0; k < n - 1; k++){
        for(int i = k + 1; i < n; i++){
            alpha = A[i][k] / A[k][k];
            for(int j = k + 1; j < n; j++){
                A[i][j] = A[i][j] - alpha * A[k][j];
            }
            b[i] = b[i] - alpha * b[k];
        }
    }

    // cout << "前進消去" << endl;
    // cout << "A = " << endl;
    // printMatrix(A);
    // cout << "b = " << endl;
    // printVector(b);

    // cout << "後退代入" << endl;
    vector<double>  x(m);
    x = Backward_sub(A, b);
    // cout << "x = " << endl;
    // printVector(x);
    
    return x;
}


// 後退代入
vector<double> Backward_sub(vector< vector<double> > A, vector<double> b){
    int m = int(b.size());
    vector<double> x(m);

    double sum;
    for(int k = m - 1; k >= 0; k--){
        sum = 0;
        for(int j = k + 1; j < m; j++){
            sum += A[k][j] * x[j];
        }
        x[k] = (1.0 / A[k][k]) * (b[k] - sum);
    }

    return x;
}


// LU分解
vector< vector<double> > LU_decomposition(vector< vector<double> > A){
    unsigned long m = A.size();
    unsigned long n = A[0].size(); // = A.size();

    double alpha;
    for(int k = 0; k < n - 1; k++){
        for(int i = k + 1; i < n; i++){
            alpha = A[i][k] / A[k][k];
            A[i][k] = alpha;
            for(int j = k + 1; j < n; j++){
                A[i][j] = A[i][j] - alpha * A[k][j];
            }
        }
    }

    return A;
}


// 前進代入
vector<double> Forward_sub(vector< vector<double> > A, vector<double> b){
    unsigned long n = b.size();
    vector<double> y(n);

    double sum;
    for(int k = 0; k < n; k++){
        sum = 0;
        for(int j = 0; j <= k - 1; j++){
            sum += A[k][j] * y[j];
        }
        y[k] = b[k] - sum;
    }

    return y;
}


// 逆行列
vector< vector<double> > Inverse_matrix(vector< vector<double> > A){
    vector< vector<double> > A_LU = LU_decomposition(A);
    unsigned long n = A.size();
    vector<double> e(n);
    vector<double> x(n);
    vector<double> y(n);
    vector< vector<double> > A_Inverse(n, vector<double>(n));

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            e[j] = 0.0;
        }
        e[i] = 1.0;
        y = Forward_sub(A_LU, e);
        x = Backward_sub(A_LU, y);
        // printVector(x);

        for(int j = 0; j < n; j++){
            A_Inverse[j][i] = x[j];
        }
    }

    return A_Inverse;
}