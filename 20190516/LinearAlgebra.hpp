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

void printMatrix(const vector<vector<double>>& A) {
  for (int i = 0; i < A.size(); i++) {
    cout << "|";
    for (int j = 0; j < A[0].size(); j++) {
      if(j == A[0].size()-1) {
        printf("%.2e|\n",A[i][j]);
        continue;
      }
      printf("%.2e\t",A[i][j]);
    }
  }
}

vector<double> VectorSubstract(const vector<double>& a, const vector<double>& b) {
  unsigned long n = a.size();
  vector<double> c(n);
  for (int i = 0; i < n; i++) {
    c[i] = a[i] - b[i];
  }

  return c;
}

vector<double> MatrixVector(const vector<vector<double>>& A, const vector<double>& b) {
  vector<double> c;
  double sum;
  for (int i = 0; i < A.size(); i++) {
    sum = 0;
    for (int j = 0; j < A[0].size(); j++) {
      sum += A[i][j] * b[j];
    }
    c.emplace_back(sum);
  }

  return c;
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
    norm += fabs(a[i]);
  }
  return norm;
}

double VectorNorm2(const vector<double>& a) {
  double norm = 0;
  for (int i = 0; i < a.size(); i++) {
    norm += sqrt(pow(a[i], 2));
  }
  return norm;
}

double VectorNormInfty(const vector<double>& a) {
  double max = 0;
  for (int i = 0; i < a.size(); i++) {
    if(max <= fabs(a[i])) max = fabs(a[i]);
  }
  return max;
}

double MatrixNorm1(const vector<vector<double>>& A) {
  double max = 0, sum;
  for (int j = 0; j < A[0].size(); j++) {
    sum = 0;
    for (int i = 0; i < A.size(); i++) {
      sum += fabs(A[i][j]);
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
      sum += fabs(A[i][j]);
    }
    if(max <= sum) max = sum;
  }
  return sum;
}

void Gauss_elimination(const vector<vector<double>>& _A, const vector<double>& _b) {
  vector<double> alpha;
  vector<vector<double>> A(_A);
  vector<double> b(_b);
  vector<double> x(A.size());
  for (int k = 0; k < A.size()-1; k++) {
    for (int i = k+1; i < A.size(); i++) {
      double alpha = A[i][k] / A[k][k];
      for (int j = k+1; j < A.size(); j++) {
        A[i][j] = A[i][j] - alpha*A[k][j];
      }
      b[i] = b[i] - alpha*b[k];
    }
  }

  for (int k = A.size()-1; k >= 0; k--) {
    double sum = 0;
    for (int j = k+1; j < A.size(); j++) {
      sum += A[k][j] * x[j];
    }
    x[k] = 1.0/A[k][k] * (b[k] - sum);
  }

  cout << "A = " << endl;
  printMatrix(A);
  cout << "b = " << endl;
  printVector(b);
  cout << "x = " << endl;
  printVector(x);
}
#endif
