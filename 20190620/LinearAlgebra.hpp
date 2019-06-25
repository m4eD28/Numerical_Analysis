#ifndef LinearAlgebra_hpp
#define LinearAlgebra_hpp

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void printVector(const vector<double>& a) {
  for (int i = 0; i < a.size(); i++) {
    printf("|%.2e|\n",a[i]);
  }
  cout << endl;

}

void printVector_more_detail(const vector<double>& a) {
  for (int i = 0; i < a.size(); i++) {
    printf("|%.12e|\n",a[i]);
  }
  cout << endl;
}

void printMatrix(const vector<vector<double>>& A) {
  for (int i = 0; i < A.size(); i++) {
    cout << "|";
    for (int j = 0; j < A[0].size(); j++) {
      if(j == A[0].size()-1) {
        printf("%.2e|\n",A.at(i).at(j));
        continue;
      }
      printf("%.2e\t",A.at(i).at(j));
    }
  }
  cout << endl;
}

vector<double> VectorSubstract(const vector<double>& a, const vector<double>& b) {
  unsigned long n = a.size();
  vector<double> c(n);
  for (int i = 0; i < n; i++) {
    c.at(i) = a.at(i) - b.at(i);
  }

  return c;
}

vector<double> MatrixVector(const vector<vector<double>>& A, const vector<double>& b) {
  vector<double> c;
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

vector<vector<double>> MatrixMatrix(const vector<vector<double>>& A, const vector<vector<double>>& B) {
  vector<vector<double>> C(A.size(), vector<double>(B.at(0).size()));

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

vector<double> ResidualError(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b) {
  vector<double> Ax;
  Ax = MatrixVector(A, x);
  vector<double> Ax_b;
  Ax_b = VectorSubstract(Ax, b);

  return Ax_b;
}

double VectorNorm1(const vector<double>& a) {
  double norm = 0;
  for (int i = 0; i < a.size(); i++) {
    norm += fabs(a.at(i));
  }
  return norm;
}

double VectorNorm2(const vector<double>& a) {
  double norm = 0;
  for (int i = 0; i < a.size(); i++) {
    norm += sqrt(pow(a.at(i), 2));
  }
  return norm;
}

double VectorNormInfty(const vector<double>& a) {
  double max = 0;
  for (int i = 0; i < a.size(); i++) {
    if(max <= fabs(a.at(i))) max = fabs(a.at(i));
  }
  return max;
}

double MatrixNorm1(const vector<vector<double>>& A) {
  double max = 0, sum;
  for (int j = 0; j < A[0].size(); j++) {
    sum = 0;
    for (int i = 0; i < A.size(); i++) {
      sum += fabs(A.at(i).at(j));
    }
    if(max <= sum) max = sum;
  }
  return sum;
}

double MatrixNormInfty(const vector<vector<double>>& A) {
  double max = 0, sum;
  for (int i = 0; i < A[0].size(); i++) {
    sum = 0;
    for (int j = 0; j < A.size(); j++) {
      sum += fabs(A.at(i).at(j));
    }
    if(max <= sum) max = sum;
  }
  return sum;
}

vector<vector<double>> Forward_erase(const vector<vector<double>>& _A, vector<double>& b) {
  vector<vector<double>> A(_A);
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

vector<vector<double>> Forward_erase(const vector<vector<double>>& _A) {
  vector<vector<double>> A(_A);
  double alpha;
  for (int k = 0; k < A.size()-1; k++) {
    for (int i = k+1; i < A.size(); i++) {
      alpha = A.at(i).at(k) / A.at(k).at(k);
      for (int j = k+1; j < A.size(); j++) {
        A.at(i).at(j) = A.at(i).at(j) - alpha*A.at(k).at(j);
      }
    }
  }

  return A;
}

vector<double> Backward_sub(const vector<vector<double>>& _A,const vector<double>& _b) {
  vector<vector<double>> A(_A);
  vector<double> b(_b);
  vector<double> x(A.size());
  double sum = 0;
  for (int k = A.size()-1; k >= 0; k--) {
    sum = 0;
    for (int j = k+1; j < A.size(); j++) {
      sum += A.at(k).at(j) * x.at(j);
    }
    x.at(k) = 1.0/A.at(k).at(k) * (b.at(k) - sum);
  }

  return x;
}

vector<double> Forward_sub(const vector<vector<double>>& _A,const vector<double>& _b) {
  vector<vector<double>> A(_A);
  vector<double> b(_b);
  vector<double> y(A.size());
  double sum;
  for (int k = 0; k < b.size(); k++) {
    sum = 0;
    for (int j = 0; j <= k-1; j++) {
      sum += A.at(k).at(j) * y.at(j);
    }
    y.at(k) = b.at(k) - sum;
  }

  return y;
}

void Gauss_elimination(const vector<vector<double>>& _A, const vector<double>& _b) {
  vector<vector<double>> A(_A);
  vector<double> b(_b);
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

  vector<double> x = Backward_sub(A, b);

  cout << "A = " << endl;
  printMatrix(A);
  cout << "b = " << endl;
  printVector(b);
  cout << "x = " << endl;
  printVector_more_detail(x);
}

vector<vector<double>> LU_decomposition(const vector<vector<double>>& _A) {
  vector<vector<double>> A(_A);
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

vector<vector<double>> Inverse_matrix(const vector<vector<double>>& A) {
  vector<vector<double>> A_LU = LU_decomposition(A);
  vector<double> e(A.size());
  vector<double> x(A.size());
  vector<double> y(A.size());
  vector<vector<double>> A_Inverse(A.size(), vector<double>(A.size()));

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

void Jacobi_law(const vector<vector<double>>& A, const vector<double>& b) {
  int M = 200;
  double eps = 1e-8;
  vector<vector<double>> x(M, vector<double>(A.size(), 1.0));
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
    if ((VectorNorm1(VectorSubstract(x.at(m), x.at(m-1)))/VectorNorm1(x.at(m-1))) <= eps) {
      cout << "x = " << endl;
      printVector(x.at(m));
      cout << "m = " << m << endl;
      return;
    }
  }
  cout << "収束しない" << endl;
  /* for (int i = 0; i < A.size(); i++) { */
  /*   x.at(0).at(i) = 1.0; */
  /* } */
  /* vector<double> x_before(A.size(), 1); */
  /* vector<double> x_after(A.size()); */
  /* for (int m = 1; m <= M; m++) { */
  /*   for (int i = 0; i < A.size(); i++) { */
  /*     sum = 0; */
  /*     for (int j = 0; j < A.size(); j++) { */
  /*       sum += A.at(i).at(j) * x_before.at(j); */
  /*     } */
  /*     x_after.at(i) = (b.at(i) - sum); */
  /*   } */
  /*   if ((VectorNormInfty(VectorSubstract(x_before, x_after))/VectorNormInfty(x_after)) <= eps) { */
  /*     cout << "x = " << endl; */
  /*     printVector(x_after); */
  /*     cout << "m = " << m << endl; */
  /*     break; */
  /*   } */
  /*   x_before = x_after; */
  /* } */

  /* cout << "A = " << endl; */
  /* printMatrix(A); */
  /* cout << "b = " << endl; */
  /* printVector(b); */
  /* cout << "x = " << endl; */
  /* printVector_more_detail(x.at()); */
}

void Gauss_Seidel_law(const vector<vector<double>>& _A, const vector<double>& _b) {
  vector<vector<double>> A(_A);
  vector<double> b(_b);
  int M = 200;
  double eps = 1e-8;
  vector<vector<double>> x(M, vector<double>(A.size(), 1.0));
  double sum1;
  double sum2;
  for (int m = 1; m < M; m++) {
    for (int i = 0; i < A.size(); i++) {
      sum1 = 0;
      for (int j = 0; j < i-1; j++) {
        sum1 += A.at(i).at(j) * x.at(m).at(j);
      }
      sum2 = 0;
      for (int j = i+1; j < A.size(); j++) {
        sum2 += A.at(i).at(j) * x.at(m-1).at(j);
      }
      x.at(m).at(i) = (b.at(i) - sum1 - sum2)/A.at(i).at(i);
    }
    if ((VectorNorm1(VectorSubstract(x.at(m-1), x.at(m)))/VectorNorm1(x.at(m))) <= eps) {
      cout << "x = " << endl;
      printVector(x.at(m));
      cout << "m = " << m << endl;
      return;
    }
  }
  cout << "収束しない" << endl;
}
#endif
