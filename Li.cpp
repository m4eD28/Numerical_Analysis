#ifndef MATRIX_H
#define MATRIX_H
#include <cassert>
#include <cmath>
#include <numeric>
#include <limits>
#include <vector>
#include <algorithm>
#include <iostream>

class Matrix {
public:
  std::vector<std::vector<double> > gyous;
  int letu_num;
  int gyou_num;
  Matrix(std::vector<std::vector<double>> _gyous) {
      gyous = _gyous;
      gyou_num = _gyous.size();
      letu_num = _gyous.at(0).size();
  }
};

std::ostream& operator<<(std::ostream& stream, const std::vector<double>& vec) {
    stream << "[";
    for(auto& p : vec) {
        stream << p << " ";
    }
    stream << "]" << std::endl;
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const Matrix& A) {
    for(const auto& p : A.gyous) {
        for(const auto& q : p) {
            stream << q << " ";
        }
        stream << std::endl;
    }
    return stream;
}

std::vector<double> operator*(const Matrix& A, const std::vector<double>& b) {
    assert(b.size() == A.letu_num);

    int m = A.gyou_num;
    int n = b.size();

    std::vector<double> result(m);

    for(int i = 0; i < m; i++) {
        double sum = 0;
        for(int j = 0; j < n; j++) {
            sum += A.gyous.at(i).at(j) * b.at(j);
        }
        result.at(i) = sum;
    }
    return result;
}

std::vector<double> operator+(const std::vector<double>& v1, const std::vector<double>& v2){
    assert(v1.size() == v2.size());
    int n = v1.size();
    std::vector<double> result(n);
    for(int i = 0; i < n; i++) {
        result.at(i) = v1.at(i) + v2.at(i);
    }
    return result;
}

std::vector<double> operator*(const double& t, const std::vector<double>& v) {
    std::vector<double> result;
    std::for_each(v.begin(), v.end(), [&](double x){result.push_back(t * x);});
    return result;
}
std::vector<double> operator*(const std::vector<double>& v, const double& t) {
    return t * v;
}

std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double> v2) {
    return v1 + (-1) * v2;
}

double L_inf_norm(const std::vector<double> vec) {
    double biggest = std::abs(vec.at(0));
    for(int i = 1; i < vec.size(); i++) {
        biggest = std::max(biggest, std::abs(vec.at(i)));
    }
    return biggest;
}

double L_inf_norm(const Matrix& A) {
    double biggest = std::numeric_limits<double>::min();
    for(int i = 0; i < A.letu_num; i++) {
        double now_letu_sum = 0;
        for(int j = 0; j < A.gyou_num; j++) {
            now_letu_sum += std::abs(A.gyous.at(j).at(i));
        }
        biggest = std::max(biggest, now_letu_sum);
    }
    return biggest;
}

double L_1_norm(const std::vector<double> vec) {
    double sum = 0;
    for(auto& p : vec) {
        sum += std::abs(p);
    }
    return sum;
}

double L_1_norm(const Matrix& A) {
    double biggest = std::numeric_limits<double>::min();
    for(const auto& p : A.gyous) {
        double now_gyou_sum = 0;
        for(const auto& q : p) {
            now_gyou_sum += std::abs(q);
        }
        biggest = std::max(biggest, now_gyou_sum);
    }
    return biggest;
}

double L_2_norm(const std::vector<double>& vec) { // L2ノルムの行列verは難しいので作らない
    double sum = 0;
    for(auto& p : vec) {
        sum += p*p;
    }
    return std::sqrt(sum);
}

void solve_with_gauss(const Matrix& _A, const std::vector<double>& _b, std::vector<double>& x,bool display_content = false) {
    assert(_A.gyou_num == _A.letu_num && _A.letu_num == _b.size());

    const int n =  _A.gyou_num;
    x.resize(n);
    Matrix A = _A;
    std::vector<double> b(_b);
    for(int k = 0; k < n-1; k++) { //forward
        for(int i = k + 1; i < n; i++) {
            double alpha = A.gyous.at(i).at(k)/A.gyous.at(k).at(k);
            for(int j = k + 1; j < n; j++) {
                A.gyous.at(i).at(j) = A.gyous.at(i).at(j) - alpha * A.gyous.at(k).at(j);  
            } 
            b.at(i) = b.at(i) - alpha * b.at(k);
        }
    }

    for(int k = n-1; k >= 0; k--) { //back
        double sum = 0;
        for(int j = k + 1; j < n; j++) {
            sum += A.gyous.at(k).at(j) * x.at(j);
        }
        x.at(k) = 1.0/A.gyous.at(k).at(k) * (b.at(k) - sum);
    }
    
    if(display_content) {
        std::cout << "チルダAは" <<std::endl << A << std::endl;
        std::cout << "チルダbは" << b << std::endl;
        std::cout << "解xは" << x << std::endl;
    }

}

void solve_with_LU(const Matrix& _A, const std::vector<double>& _b, std::vector<double>& x, bool display_content = false) {
    Matrix A = _A;
    std::vector<double> b = _b;
    int n = A.gyou_num;

    ///1. Lの殻
    std::vector<std::vector<double> > l;
    l.resize(n);
    for(auto& p : l) {
        p.resize(n);
    }
    Matrix L(l);

    for(int i = 0; i <= n-1; i++) { //Lの対角成分は1である
        L.gyous.at(i).at(i) = 1.0;
    }

    ///2. Uの殻
    std::vector<std::vector<double> > u;
    u.resize(n);
    for(auto& p : u) {
        p.resize(n);
    }
    Matrix U(u);
    U.gyous.at(0) = A.gyous.at(0); //Uの0行目はAの0行目と同じ

    ///3. LU分解
    ///3.1 前進消去しつつLとUの更新
    

    for(int k = 0; k <= n - 2; k++) { // k行目をかなめにして
        for(int i = k+1; i <= n - 1; i++) { // k+1行目 ~ n-1行目中のi行目において
            double alpha_i_k = A.gyous.at(i).at(k) / A.gyous.at(k).at(k);
            L.gyous.at(i).at(k) = alpha_i_k;
            for(int j = k+1; j <= n - 1; j++) { // (i,k+1) ~ (i,j) ~ (i,n-1)に操作を加える。つまりjは列
                A.gyous.at(i).at(j) -= alpha_i_k * A.gyous.at(k).at(j);
                if(i <= j) {
                    U.gyous.at(i).at(j) = A.gyous.at(i).at(j);
                }
            } 
        }
    }

    if(display_content) {
        std::cout << "Lは" << std::endl;
        std::cout << L << std::endl;
        
        std::cout << "Uは" << std::endl;
        std::cout << U << std::endl;
    }
    
    //4. 前進代入でLy = bのyを求める
    std::vector<double> y(n);
    for(int k = 0; k <= n-1; k++) { //0~y~n-1のk番目は
        double sum = b.at(k);
        for(int j = 0; j <= k-1; j++) { //0~j~k-1において足し算を行い
            double alpha_k_j = L.gyous.at(k).at(j);
            sum -= alpha_k_j * y.at(j);
        }
        y.at(k) = sum;
    }

    if(display_content) {
        std::cout << "yは" << std::endl;
        std::cout << y << std::endl;
    }

    //5. 後退代入で Ux = y の x を求める
    x.resize(n);
    for(int k = n-1; k >= 0; k--) {
        double sum = 0;
        for(int j = n-1; j > k; j--) {
            sum += U.gyous.at(k).at(j) * x.at(j);
        }
        x.at(k) = (y.at(k) - sum)/U.gyous.at(k).at(k) ;
    }

    if(display_content) {
        std::cout << "xは" << std::endl;
        std::cout << x << std::endl;
    }
}


void Inverse(const Matrix& A, Matrix& Inv_A) {
    const int n = A.gyou_num;
    std::vector<std::vector<double> > vecs(n); //n個のベクトルの組
    for(auto& p : vecs) {
        p.resize(n);
    }
    std::vector<std::vector<double> > e(n); //標準基底たち
    for(auto& p : e) {
        p.resize(n);
    }
    for(int i = 0; i < n; i++) {
        e.at(i).at(i) = 1.0;
    }
    for(int i = 0; i < n; i++) {
        solve_with_LU(A, e.at(i), vecs.at(i)); // Todo: 無駄が多いのでここでLU分解を再実装すべき
                                               // ---> solve_with_LUをLUを作る部分と前進/後進代入をする部分で分ける?
    }
    for(int i = 0; i < n ; i++) {
        for(int j = 0; j < n; j++) {
            Inv_A.gyous.at(j).at(i) = vecs.at(i).at(j);
        }
    }
}

#endif
