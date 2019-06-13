#ifndef MATRIX_H
#define MATRIX_H
#include <cassert>
#include <vector>
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
  std::vector<double> operator*(std::vector<double> vec) {
      assert(vec.size() == letu_num);
      std::vector<double> res(letu_num);
      for(auto& p : gyous) {
          for(int i = 0; i < letu_num; i++) {
              res.at(i) += vec.at(i) * p.at(i);
          }
      }
      return res;
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


void solve(const Matrix& _A, const std::vector<double>& _b, std::vector<double>& x) {
    assert(_A.gyou_num == _A.letu_num && _A.letu_num == _b.size());

    const int n =  _A.gyou_num;
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

    std::cout << "チルダA: " << std::endl << A << std::endl;
    std::cout << "bは" << b << std::endl;
    std::cout << "xは" << x << std::endl;
}


#endif