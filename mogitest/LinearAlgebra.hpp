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
    printf("|%.12e|\n",a[i]);
  }
  std::cout << std::endl;
}

void printMatrix(const std::vector<std::vector<double>>& A) {
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

std::vector<double> MatrixVector(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
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

std::vector<std::vector<double>> MatrixMatrix(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
  std::vector<std::vector<double>> C(A.size(), std::vector<double>(B.at(0).size()));

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

std::vector<double> ResidualError(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b) {
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

double MatrixNorm1(const std::vector<std::vector<double>>& A) {
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

double MatrixNormInfty(const std::vector<std::vector<double>>& A) {
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

std::vector<std::vector<double>> Forward_erase(const std::vector<std::vector<double>>& _A, std::vector<double>& b) {
  std::vector<std::vector<double>> A(_A);
  double alpha;
  for (int k = 0; k < A.size()-1; k++) {
    for (int i = k+1; i < A.size(); i++) {
      alpha = A.at(i).at(k) / A.at(k).at(k);
      for (int j = k+1; j < A.size(); j++) {
        A.at(i).at(j) = A.at(i).at(j) - alpha*A.at(k).at(j);
      }
      b.at(i) = b.at(i) - alpha*b.at(k);
    }
  }

  return A;
}

std::vector<std::vector<double>> Forward_erase(const std::vector<std::vector<double>>& _A) {
  std::vector<std::vector<double>> A(_A);
  double alpha;
  for (int k = 0; k < A.size()-1; k++) {
    for (int i = k+1; i < A.size(); i++) {
      alpha = A.at(i).at(k) / A.at(k).at(k);
      for (int j = k+1; j < A.size(); j++) {
        A.at(i).at(j) = A.at(i).at(j) - alpha*A.at(k).at(j);
      }
      A.at(i).at(k) = alpha;
    }
  }

  return A;
}

std::vector<double> Backward_sub(const std::vector<std::vector<double>>& _A,const std::vector<double>& _b) {
  std::vector<std::vector<double>> A(_A);
  std::vector<double> b(_b);
  std::vector<double> x(A.size());
  double sum = 0;
  for (int k = A.size()-1; k >= 0; k--) {
    sum = 0;
    for (int j = k+1; j < A.size(); j++) {
      sum += A.at(k).at(j) * x.at(j);
      /* sum += A.at(k).at(j) * x.at(j) / A.at(k).at(k); */
    }
    x.at(k) = (b.at(k) - sum) / A.at(k).at(k);
    /* x.at(k) = b.at(k) / A.at(k).at(k) - sum; */
  }

  return x;
}

std::vector<double> Forward_sub(const std::vector<std::vector<double>>& _A,const std::vector<double>& _b) {
  std::vector<std::vector<double>> A(_A);
  std::vector<double> b(_b);
  std::vector<double> y(A.size());
  double sum;
  for (int k = 0; k < b.size(); k++) {
    sum = 0;
    for (int j = 0; j < k; j++) {
      sum += A.at(k).at(j) * y.at(j);
    }
    y.at(k) = b.at(k) - sum;
  }

  return y;
}

void Gauss_elimination(const std::vector<std::vector<double>>& _A, const std::vector<double>& _b) {
  std::vector<std::vector<double>> A(_A);
  std::vector<double> b(_b);
  double alpha;
  for (int k = 0; k < A.size()-1; k++) {
    for (int i = k+1; i < A.size(); i++) {
      alpha = A.at(i).at(k) / A.at(k).at(k);
      for (int j = k+1; j < A.size(); j++) {
        A.at(i).at(j) = A.at(i).at(j) - alpha*A.at(k).at(j);
      }
      b.at(i) = b.at(i) - alpha*b.at(k);
    }
  }

  std::vector<double> x = Backward_sub(A, b);

  std::cout << "A = " << std::endl;
  printMatrix(A);
  std::cout << "b = " << std::endl;
  printVector(b);
  std::cout << "x = " << std::endl;
  printVector_more_detail(x);
}

std::vector<std::vector<double>> LU_decomposition(const std::vector<std::vector<double>>& _A) {
  std::vector<std::vector<double>> A(_A);
  double alpha;
  for (int k = 0; k < A.size()-1; k++) {
    for (int i = k+1; i < A.size(); i++) {
      alpha = A.at(i).at(k) / A.at(k).at(k);
      A.at(i).at(k) = alpha;
      for (int j = k+1; j < A.size(); j++) {
        A.at(i).at(j) = A.at(i).at(j) - alpha*A.at(k).at(j);
      }
    }
  }

  return A;
}

std::vector<std::vector<double>> Inverse_matrix(const std::vector<std::vector<double>>& A) {
  std::vector<std::vector<double>> A_LU = LU_decomposition(A);
  std::vector<double> e(A.size());
  std::vector<double> x(A.size());
  std::vector<double> y(A.size());
  std::vector<std::vector<double>> A_Inverse(A.size(), std::vector<double>(A.size()));

  for (int i = 0; i < A.size(); i++) {
    for (int j = 0; j < A.size(); j++) {
      e.at(j) = 0.0;
    }
    e.at(i) = 1.0;
    y = Forward_sub(A_LU, e);
    x = Backward_sub(A_LU, y);
    for (int j = 0; j < A.size(); j++) {
      A_Inverse.at(j).at(i) = x.at(j);
    }
  }
  return A_Inverse;
}

void Jacobi_law(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
  int M = 200;
  double eps = 1e-8;
  std::vector<std::vector<double>> x(M, std::vector<double>(A.size(), 1.0));
  double sum;
  for (int m = 1; m < M; m++) {
    for (int i = 0; i < A.size(); i++) {
      sum = 0;
      for (int j = 0; j < A.at(0).size(); j++) {
        if (j == i) continue;
        /* sum += A.at(i).at(j) * x.at(m-1).at(j) / A.at(i).at(i); */
        sum += A.at(i).at(j) * x.at(m-1).at(j);
      }
      /* x.at(m).at(i) = b.at(i)/A.at(i).at(i) - sum ; */
      x.at(m).at(i) = (b.at(i) - sum) / A.at(i).at(i);
    }
    if ((VectorNormInfty(VectorSubstract(x.at(m), x.at(m-1)))/VectorNormInfty(x.at(m-1))) <= eps) {
      std::cout << "x = " << std::endl;
      printVector(x.at(m));
      std::cout << "m = " << m << std::endl;
      return;
    }
  }
  std::cout << "収束しない" << std::endl;
}

void Gauss_Seidel_law(const std::vector<std::vector<double>>& _A, const std::vector<double>& _b) {
  std::vector<std::vector<double>> A(_A);
  std::vector<double> b(_b);
  int M = 200;
  double eps = 1e-8;
  std::vector<std::vector<double>> x(M, std::vector<double>(A.size(), 1.0));
  double sum1;
  double sum2;
  for (int m = 1; m < M; m++) {
    for (int i = 0; i < A.size(); i++) {
      sum1 = 0;
      for (int j = 0; j < i; j++) {
        sum1 += A.at(i).at(j) * x.at(m).at(j);
      }
      sum2 = 0;
      for (int j = i+1; j < A.size(); j++) {
        sum2 += A.at(i).at(j) * x.at(m-1).at(j);
      }
      x.at(m).at(i) = (b.at(i) - sum1 - sum2)/A.at(i).at(i);
    }
    if ((VectorNormInfty(VectorSubstract(x.at(m-1), x.at(m)))/VectorNormInfty(x.at(m))) <= eps) {
      std::cout << "x = " << std::endl;
      printVector(x.at(m));
      std::cout << "m = " << m << std::endl;
      return;
    }
  }
  std::cout << "収束しない" << std::endl;
}
#endif
