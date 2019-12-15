#ifndef ALOG_H
#define ALOG_H


#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include "LinearAlgebra.hpp"

std::vector<std::vector<double> > Forward_erase(const std::vector<std::vector<double> >& _A, std::vector<double>& b);

std::vector<std::vector<double> > Forward_erase(const std::vector<std::vector<double> >& _A);

std::vector<double> Backward_sub(const std::vector<std::vector<double> >& _A,const std::vector<double>& _b);

std::vector<double> Forward_sub(const std::vector<std::vector<double> >& _A,const std::vector<double>& _b);

std::vector<std::vector<double> > LU_decomposition(const std::vector<std::vector<double> >& _A);

std::vector<std::vector<double> > Inverse_matrix(const std::vector<std::vector<double> >& A);

std::vector<double> Jacobi_law(const std::vector<std::vector<double> >& A, const std::vector<double>& b);

std::vector<double> Gauss_Seidel_law(const std::vector<std::vector<double> >& _A, const std::vector<double>& _b);

std::vector<double> Gauss_elimination(const std::vector<std::vector<double> >& _A, const std::vector<double>& _b);

std::vector<double> Gauss_elimination_pivot(const std::vector<std::vector<double> >& _A, const std::vector<double>& _b);

double Gauss_elimination_norm(const std::vector<std::vector<double> >& _A, const std::vector<double>& _b);

double Gauss_elimination_pivot_norm(const std::vector<std::vector<double> >& _A, const std::vector<double>& _b);

std::vector<double> Euler_method(double x_0, double a, double T, int N);
#endif
