#ifndef LinearAlgebra_hpp
#define LinearAlgebra_hpp

#include <iostream>
#include <vector>
#include <cmath>

class Matrix{
  public:
    std::vector<std::vector<double> > Line;
    int Column_num;
    int Row_num;
    Matrix(std::vector<std::vector<double> > _Line) {
      Line = _Line;
      Column_num = _Line.size();
      Row_num = _Line.at(0).size();
    }
};

std::ostream& operator<<(std::ostream& stream, const std::vector<double>& vec) {
  stream << "|";
  for (auto& p : vec) {
    stream << p << " ";
  }
  stream << "|" << std::endl;
  return stream;
}

std::ostream& operator<<(std::ostream& stream, const Matrix& A) {
    for(const auto& p : A.Line) {
        for(const auto& q : p) {
            stream << q << " ";
        }
        stream << std::endl;
    }
    return stream;
}

std::vector<double> operator*(const Matrix& A, const std::vector<double>& b) {
    assert(b.size() == A.Row_num);

    int m = A.Column_num;
    int n = b.size();

    std::vector<double> result(m);

    for(int i = 0; i < m; i++) {
        double sum = 0;
        for(int j = 0; j < n; j++) {
            sum += A.Line.at(i).at(j) * b.at(j);
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

void printVector(const std::vector<double>& a) {
  for (int i = 0; i < a.size(); i++) {
    printf("|%.2e|\n",a[i]);
  }
  std::cout << std::endl;

}

void printVector_more_detail(const std::vector<double>& a) {
  for (int i = 0; i < a.size(); i++) {
    printf("|%.14e|\n",a[i]);
  }
  std::cout << std::endl;
}

void printMatrix(const std::vector<std::vector<double> >& A) {
  for (int i = 0; i < A.size(); i++) {
    std::cout << "|";
    for (int j = 0; j < A[0].size(); j++) {
      if(j == A[0].size()-1) {
        printf("%.2e|\n",A.at(i).at(j));
        continue;
      }
      printf("%.2e\t",A.at(i).at(j));
    }
  }
  std::cout << std::endl;
}

std::vector<double> VectorSubstract(const std::vector<double>& a, const std::vector<double>& b) {
  unsigned long n = a.size();
  std::vector<double> c(n);
  for (int i = 0; i < n; i++) {
    c.at(i) = a.at(i) - b.at(i);
  }

  return c;
}

std::vector<double> MatrixVector(const std::vector<std::vector<double> >& A, const std::vector<double>& b) {
  std::vector<double> c;
  double sum;
  for (int i = 0; i < A.size(); i++) {
    sum = 0;
    for (int j = 0; j < A[0].size(); j++) {
      sum += A.at(i).at(j) * b.at(j);
    }
    c.emplace_back(sum);
  }

  return c;
}

std::vector<std::vector<double> > MatrixMatrix(const std::vector<std::vector<double> >& A, const std::vector<std::vector<double> >& B) {
  std::vector<std::vector<double> > C(A.size(), std::vector<double>(B.at(0).size()));

  double sum;
  for (int i = 0; i < A.size(); i++) {
    for (int j = 0; j < A.at(0).size(); j++) {
      sum = 0;
      for (int k = 0; k < A.at(0).size(); k++) {
        sum += A.at(i).at(k) * B.at(k).at(j);
      }
      C.at(i).at(j) = sum;
    }
  }

  return C;
}

std::vector<double> ResidualError(const std::vector<std::vector<double> >& A, const std::vector<double>& x, const std::vector<double>& b) {
  std::vector<double> Ax;
  Ax = MatrixVector(A, x);
  std::vector<double> Ax_b;
  Ax_b = VectorSubstract(Ax, b);

  return Ax_b;
}

double VectorNorm1(const std::vector<double>& a) {
  double norm = 0;
  for (int i = 0; i < a.size(); i++) {
    norm += fabs(a.at(i));
  }
  return norm;
}

double VectorNorm2(const std::vector<double>& a) {
  double norm = 0;
  for (int i = 0; i < a.size(); i++) {
    norm += sqrt(pow(a.at(i), 2));
  }
  return norm;
}

double VectorNormInfty(const std::vector<double>& a) {
  double max = 0;
  for (int i = 0; i < a.size(); i++) {
    if(max <= fabs(a.at(i))) max = fabs(a.at(i));
  }
  return max;
}

double MatrixNorm1(const std::vector<std::vector<double> >& A) {
  double max = 0.0;
  double sum;
  for (int j = 0; j < A.at(0).size(); j++) {
    sum = 0.0;
    for (int i = 0; i < A.size(); i++) {
      sum += fabs(A.at(i).at(j));
    }
    if(max <= sum) max = sum;
  }
  return max;
}

double MatrixNormInfty(const std::vector<std::vector<double> >& A) {
  double max = 0, sum;
  for (int i = 0; i < A.size(); i++) {
    sum = 0;
    for (int j = 0; j < A.at(0).size(); j++) {
      sum += fabs(A.at(i).at(j));
    }
    if(max <= sum) max = sum;
  }
  return max;
}

#endif
