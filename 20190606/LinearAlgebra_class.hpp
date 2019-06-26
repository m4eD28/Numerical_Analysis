#ifndef LinearAlgebra_hpp
#define LinearAlgebra_hpp

#include <iostream>
#include <vector>
#include <cmath>

class Matrix {
  public:
    std::vector<std::vector<double>> line;
    int Column_num;
    int Line_num;
    Matrix(std::vector<std::vector<double>> _line) {
      line = _line;
      Column_num = _line.size();
      Line_num = _line.at(0).size();
    }
};

std::ostream& operator<<(std::ostream& stream, const std::vector<double>& vec) {
  stream << "[";
  for (auto& p : vec) {
    stream << p << " ";
  }
  stream << "]" << std::endl;
  return stream;
}

std::ostream& operator<<(std::ostream& stream, const Matrix& A) {
    for(const auto& p : A.line) {
        for(const auto& q : p) {
            stream << q << " ";
        }
        stream << std::endl;
    }
    return stream;
}

std::vector<double> operator*(const Matrix& A, const std::vector<double>& b) {
    assert(b.size() == A.Line_num);

    int m = A.Column_num;
    int n = b.size();

    std::vector<double> result(m);

    for(int i = 0; i < m; i++) {
        double sum = 0;
        for(int j = 0; j < n; j++) {
            sum += A.line.at(i).at(j) * b.at(j);
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
    std::for_each(v.begin(), v.end(), [&](double x){result.emplace_back(t * x);});
    return result;
}
std::vector<double> operator*(const std::vector<double>& v, const double& t) {
    return t * v;
}

std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double> v2) {
    return v1 - v2;
}

double VectorNormInfty(const std::vector<double> vec) {
    double biggest = std::abs(vec.at(0));
    for(int i = 1; i < vec.size(); i++) {
        biggest = std::max(biggest, std::abs(vec.at(i)));
    }
    return biggest;
}

double MatrixNormInfty(const Matrix& A) {
    double biggest = std::numeric_limits<double>::min();
    for(int i = 0; i < A.Line_num; i++) {
        double now_line_sum = 0;
        for(int j = 0; j < A.Line_num; j++) {
            now_line_sum += std::abs(A.line.at(j).at(i));
        }
        biggest = std::max(biggest, now_line_sum);
    }
    return biggest;
}

double VectorNorm1(const std::vector<double> vec) {
    double sum = 0;
    for(auto& p : vec) {
        sum += std::abs(p);
    }
    return sum;
}

double MatrixNorm1(const Matrix& A) {
    double biggest = std::numeric_limits<double>::min();
    for(const auto& p : A.line) {
        double now_gyou_sum = 0;
        for(const auto& q : p) {
            now_gyou_sum += std::abs(q);
        }
        biggest = std::max(biggest, now_gyou_sum);
    }
    return biggest;
}

double VectorNorm2(const std::vector<double>& vec) {
    double sum = 0;
    for(auto& p : vec) {
        sum += p*p;
    }
    return std::sqrt(sum);
}

void Gauss_elimination(const Matrix& _A, const std::vector<double>& _b, std::vector<double>& x,bool display_content = false) {
    assert(_A.Column_num == _A.Line_num && _A.Line_num == _b.size());

    const int n =  _A.Column_num;
    x.resize(n);
    Matrix A = _A;
    std::vector<double> b(_b);
    for(int k = 0; k < n-1; k++) {
        for(int i = k + 1; i < n; i++) {
            double alpha = A.line.at(i).at(k)/A.line.at(k).at(k);
            for(int j = k + 1; j < n; j++) {
                A.line.at(i).at(j) = A.line.at(i).at(j) - alpha * A.line.at(k).at(j);
            }
            b.at(i) = b.at(i) - alpha * b.at(k);
        }
    }

    for(int k = n-1; k >= 0; k--) {
        double sum = 0;
        for(int j = k + 1; j < n; j++) {
            sum += A.line.at(k).at(j) * x.at(j);
        }
        x.at(k) = 1.0/A.line.at(k).at(k) * (b.at(k) - sum);
    }

    if(display_content) {
        std::cout << "チルダAは" <<std::endl << A << std::endl;
        std::cout << "チルダbは" << b << std::endl;
        std::cout << "解xは" << x << std::endl;
    }

}

void LU_solve(const Matrix& _A, const std::vector<double>& _b, std::vector<double>& x, bool display_content = false) {
    Matrix A = _A;
    std::vector<double> b = _b;
    int n = A.Column_num;

    std::vector<std::vector<double> > l;
    l.resize(n);
    for(auto& p : l) {
        p.resize(n);
    }
    Matrix L(l);

    for(int i = 0; i <= n-1; i++) {
        L.line.at(i).at(i) = 1.0;
    }

    std::vector<std::vector<double> > u;
    u.resize(n);
    for(auto& p : u) {
        p.resize(n);
    }
    Matrix U(u);
    U.line.at(0) = A.line.at(0);


    for(int k = 0; k <= n - 2; k++) {
        for(int i = k+1; i <= n - 1; i++) {
            double alpha_i_k = A.line.at(i).at(k) / A.line.at(k).at(k);
            L.line.at(i).at(k) = alpha_i_k;
            for(int j = k+1; j <= n - 1; j++) {
                A.line.at(i).at(j) -= alpha_i_k * A.line.at(k).at(j);
                if(i <= j) {
                    U.line.at(i).at(j) = A.line.at(i).at(j);
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

    std::vector<double> y(n);
    for(int k = 0; k <= n-1; k++) {
        double sum = b.at(k);
        for(int j = 0; j <= k-1; j++) {
            double alpha_k_j = L.line.at(k).at(j);
            sum -= alpha_k_j * y.at(j);
        }
        y.at(k) = sum;
    }

    if(display_content) {
        std::cout << "yは" << std::endl;
        std::cout << y << std::endl;
    }

    x.resize(n);
    for(int k = n-1; k >= 0; k--) {
        double sum = 0;
        for(int j = n-1; j > k; j--) {
            sum += U.line.at(k).at(j) * x.at(j);
        }
        x.at(k) = (y.at(k) - sum)/U.line.at(k).at(k) ;
    }

    if(display_content) {
        std::cout << "xは" << std::endl;
        std::cout << x << std::endl;
    }
}


void Inverse(const Matrix& A, Matrix& Inv_A) {
    const int n = A.Column_num;
    std::vector<std::vector<double> > vecs(n);
    for(auto& p : vecs) {
        p.resize(n);
    }
    std::vector<std::vector<double> > e(n);
    for(auto& p : e) {
        p.resize(n);
    }
    for(int i = 0; i < n; i++) {
        e.at(i).at(i) = 1.0;
    }
    for(int i = 0; i < n; i++) {
        LU_solve(A, e.at(i), vecs.at(i));
    }
    for(int i = 0; i < n ; i++) {
        for(int j = 0; j < n; j++) {
            Inv_A.line.at(j).at(i) = vecs.at(i).at(j);
        }
    }
}

#endif
