#ifndef ALOG_H
#define ALOG_H


#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include "LinearAlgebra.hpp"

std::vector<std::vector<double> > Forward_erase(const std::vector<std::vector<double> >& _A, std::vector<double>& b) {
  std::vector<std::vector<double> > A(_A);
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

std::vector<std::vector<double> > Forward_erase(const std::vector<std::vector<double> >& _A) {
  std::vector<std::vector<double> > A(_A);
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

/* std::vector<double> Backward_sub(const std::vector<std::vector<double> >& _A,const std::vector<double>& _b) { */
/*   std::vector<std::vector<double> > A(_A); */
/*   std::vector<double> b(_b); */
/*   std::vector<double> x(A.size()); */
/*   double sum = 0; */
/*   for (int k = A.size()-1; k >= 0; k--) { */
/*     sum = 0; */
/*     for (int j = k+1; j < A.size(); j++) { */
/*       sum += A.at(k).at(j) * x.at(j); */
/*       /1* sum += A.at(k).at(j) * x.at(j) / A.at(k).at(k); *1/ */
/*     } */
/*     x.at(k) = (b.at(k) - sum) / A.at(k).at(k); */
/*     /1* x.at(k) = b.at(k) / A.at(k).at(k) - sum; *1/ */
/*   } */

/*   return x; */
/* } */

std::vector<double> Backward_sub(const std::vector<std::vector<double> >& A,const std::vector<double>& _b) {
  std::vector<double> b(_b);
  std::vector<double> x(A.size());
  for (int k = A.size()-1; k >= 0; k--) {
    for (int j = k+1; j < A.size(); j++) {
      b.at(k) -= A.at(k).at(j) * x.at(j);
    }
    x.at(k) = b.at(k) / A.at(k).at(k);
  }

  return x;
}

/* std::vector<double> Forward_sub(const std::vector<std::vector<double> >& _A,const std::vector<double>& _b) { */
/*   std::vector<std::vector<double> > A(_A); */
/*   std::vector<double> b(_b); */
/*   std::vector<double> y(A.size()); */
/*   double sum; */
/*   for (int k = 0; k < b.size(); k++) { */
/*     sum = 0; */
/*     for (int j = 0; j < k; j++) { */
/*       sum += A.at(k).at(j) * y.at(j); */
/*     } */
/*     y.at(k) = b.at(k) - sum; */
/*   } */

/*   return y; */
/* } */

std::vector<double> Forward_sub(const std::vector<std::vector<double> >& _A,const std::vector<double>& _b) {
  std::vector<std::vector<double> > A(_A);
  std::vector<double> b(_b);
  std::vector<double> y(b);
  for (int k = 0; k < b.size(); k++) {
    for (int j = 0; j < k; j++) {
      y.at(k) -= A.at(k).at(j) * y.at(j);
    }
  }

  return y;
}


std::vector<std::vector<double> > LU_decomposition(const std::vector<std::vector<double> >& _A) {
  std::vector<std::vector<double> > A(_A);
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

std::vector<std::vector<double> > Inverse_matrix(const std::vector<std::vector<double> >& A) {
  std::vector<std::vector<double> > A_LU = LU_decomposition(A);
  std::vector<double> e(A.size());
  std::vector<double> x(A.size());
  std::vector<double> y(A.size());
  std::vector<std::vector<double> > A_Inverse(A.size(), std::vector<double>(A.size()));

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

std::vector<double> Jacobi_law(const std::vector<std::vector<double> >& A, const std::vector<double>& b) {
  int M = 200;
  double eps = 1e-8;
  std::vector<std::vector<double> > x(M, std::vector<double>(A.size(), 1.0));
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
      return x.at(m);
    }
  }
  std::cout << "収束しない" << std::endl;
  return x.at(0);
}

std::vector<double> Gauss_Seidel_law(const std::vector<std::vector<double> >& _A, const std::vector<double>& _b) {
  std::vector<std::vector<double> > A(_A);
  std::vector<double> b(_b);
  int M = 200;
  double eps = 1e-8;
  std::vector<std::vector<double> > x(M, std::vector<double>(A.size(), 1.0));
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
      return x.at(m);
    }
  }
  std::cout << "収束しない" << std::endl;
  return x.at(0);
}

std::vector<double> Gauss_elimination(const std::vector<std::vector<double> >& _A, const std::vector<double>& _b, int n) {
  std::vector<std::vector<double> > A(_A);
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

  std::vector<double> x1 = Backward_sub(A, b);

  /* std::vector<double> gaus_b(n, 1.0); */
  /* std::vector<std::vector<double> > A_Inverse; */
  /* A_Inverse = Inverse_matrix(A); */


  // prints
  std::cout << "A = " << std::endl;
  printMatrix(A);
  std::cout << "b = " << std::endl;
  printVector(b);
  std::cout << "x = " << std::endl;
  printVector_more_detail(x1);


  /* double A_Norm_inf = MatrixNormInfty(A); */
  /* printf("Hilbert_Norm_inf = %.2e\n\n",A_Norm_inf); */
  /* std::cout << "----------\n\n"; */
  /* double A_Inverse_Norm_inf = MatrixNormInfty(A_Inverse); */
  /* printf("Hilbert_Inverse_Norm_inf = %.2e\n\n",A_Inverse_Norm_inf); */
  /* std::cout << "----------\n\n"; */

  /* // Kappa */
  /* std::cout << "Kappa = \n"; */
  /* double kappa = MatrixNormInfty(A) * MatrixNormInfty(A_Inverse); */
  /* printf("%.2e\n\n",kappa); */
  /* std::cout << "----------\n\n"; */

  double b_Ax = VectorNormInfty(VectorSubstract(_b, MatrixVector(_A, x1)));
  printf("|b - Ax|_inf = %.8e\n", b_Ax);
  std::cout << "----------\n\n";

  /* double x_tilde_x = VectorNormInfty(VectorSubstract(x_tilde, x)); */
  /* printf("|x* - x|_inf = %.2e\n", x_tilde_x); */
  /* std::cout << "----------\n\n"; */
  return x1;
}

std::vector<double> Gauss_elimination_pivot(const std::vector<std::vector<double> >& _A, const std::vector<double>& _b, int n) {
  std::vector<std::vector<double> > A(_A);
  std::vector<double> b(_b);
  std::vector<std::vector<double> > A_true(_A);
  std::vector<double> b_true(_b);
  double alpha;
  for (int k = 0; k < A.size()-1; k++) {

    /* pivot選択 */
    double max = abs(A.at(k).at(k));
    int index = k;
    for (int i = k; i < n; i++) {
      if(max < abs(A.at(i).at(k))){
        max = abs(A.at(i).at(k));
        index = i;
      }
    }
    /* int index = std::max_element(A.at(k).at(k), A.at(k).at(k+1)); */
    /* std::vector<double>::iterator index = *std::max_element(A.at(k).at(k), A.at(n-1).at(k)); */

    for (int i = k; i < n; i++) {
      std::swap(A.at(index).at(i), A.at(k).at(i));
    }
    for (int i = 0; i < n; i++) {
      std::swap(A_true.at(index).at(i), A_true.at(k).at(i));
    }

    std::swap(b.at(index), b.at(k));
    std::swap(b_true.at(index), b_true.at(k));

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

  double b_Ax = VectorNorm2(VectorSubstract(b_true, MatrixVector(A_true, x)));
  printf("|b - Ax|_inf = %.8e\n", b_Ax);
  std::cout << "----------\n\n";
  return x;
}

double Gauss_elimination_norm(const std::vector<std::vector<double> >& _A, const std::vector<double>& _b, int n) {
  std::vector<std::vector<double> > A(_A);
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

  double b_Ax = VectorNorm2(VectorSubstract(_b, MatrixVector(_A, x)));

  return b_Ax;
}

double Gauss_elimination_pivot_norm(const std::vector<std::vector<double> >& _A, const std::vector<double>& _b, int n) {
  std::vector<std::vector<double> > A(_A);
  std::vector<std::vector<double> > A_true(_A);
  std::vector<double> b(_b);
  std::vector<double> b_true(_b);
  double alpha;
  for (int k = 0; k < A.size()-1; k++) {

    /* pivot選択 */
    double max = abs(A.at(k).at(k));
    int index = k;
    for (int i = k; i < n; i++) {
      if(max < abs(A.at(i).at(k))){
        max = abs(A.at(i).at(k));
        index = i;
      }
    }

    for (int i = k; i < n; i++) {
      std::swap(A.at(index).at(i), A.at(k).at(i));
    }
    for (int i = 0; i < n; i++) {
      std::swap(A_true.at(index).at(i), A_true.at(k).at(i));
    }
    std::swap(b.at(index), b.at(k));
    std::swap(b_true.at(index), b_true.at(k));

    for (int i = k+1; i < A.size(); i++) {
      alpha = A.at(i).at(k) / A.at(k).at(k);
      for (int j = k+1; j < A.size(); j++) {
        A.at(i).at(j) = A.at(i).at(j) - alpha*A.at(k).at(j);
      }
      b.at(i) = b.at(i) - alpha*b.at(k);
    }
  }

  std::vector<double> x = Backward_sub(A, b);

  double b_Ax = VectorNorm2(VectorSubstract(b_true, MatrixVector(A_true, x)));
  return b_Ax;
}

std::vector<double> Euler_method(double x_0, double a, double T, int N) {
  std::vector<double> x(N);
  x.at(0) = x_0;
  double h = T/N;
  double t = 0;
  for (int i = 1; i < N; i++) {
    t += h;
    x.at(i) = x.at(i-1) + h * a * x.at(i-1);
  }
  return x;
}
#endif
