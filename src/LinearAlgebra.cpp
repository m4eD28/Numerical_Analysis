#ifndef LinearAlgebra_hpp
#define LinearAlgebra_hpp

#include <iostream>
#include <vector>
#include <cmath>

void printVector(const std::vector<double>& a) {
  for (int i = 0; i < a.size(); i++) {
    printf("|%.2e|\n",a[i]);
  }
  std::cout << std::endl;

}

void printVector_more_detail(const std::vector<double>& a) {
  for (int i = 0; i < a.size(); i++) {
    printf("|%.6e|\n",a[i]);
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

void printMatrix_more_detail(const std::vector<std::vector<double> >& A) {
  for (int i = 0; i < A.size(); i++) {
    std::cout << "|";
    for (int j = 0; j < A[0].size(); j++) {
      if(j == A[0].size()-1) {
        printf("%.6e|\n",A.at(i).at(j));
        continue;
      }
      printf("%.6e\t",A.at(i).at(j));
    }
  }
  std::cout << std::endl;
}

std::vector<double> VectorAdd(const std::vector<double>& a, const std::vector<double>& b) {
  unsigned long n = a.size();
  std::vector<double> c(n);
  for (int i = 0; i < n; i++) {
    c.at(i) = a.at(i) + b.at(i);
  }

  return c;
}

std::vector<double> VectorSubstract(const std::vector<double>& a, const std::vector<double>& b) {
  unsigned long n = a.size();
  std::vector<double> c(n);
  for (int i = 0; i < n; i++) {
    c.at(i) = a.at(i) - b.at(i);
  }

  return c;
}

double VectorToScalar(const std::vector<double>& a, const std::vector<double>& b) {
  double val = 0;
  for (int i = 0; i < a.size(); ++i) {
    val += a.at(i) + b.at(i);
  }
  return val;
}

std::vector<double> VectorScalar(const double a, const std::vector<double>& b) {
  std::vector<double> vec(b);
  for (double& elem : vec) {
    elem = a * elem;
  }
  return vec;
}

std::vector<std::vector<double> > MatrixSclar(const double a, const std::vector<std::vector<double> >& b) {
  std::vector<std::vector<double> > mat(b);
  for (std::vector<double>& vec : mat) {
    for (double& elem : vec) {
      elem = a * elem;
    }
  }
  return mat;
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
