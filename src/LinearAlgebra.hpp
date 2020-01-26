#ifndef LinearAlgebra_hpp
#define LinearAlgebra_hpp

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

constexpr double eps = std::numeric_limits<double>::epsilon();

void printVector(const std::vector<double>& a);

void printVector_more_detail(const std::vector<double>& a);

void printMatrix(const std::vector<std::vector<double> >& A);

void printMatrix_more_detail(const std::vector<std::vector<double> >& A);

std::vector<double> VectorAdd(const std::vector<double>& a, const std::vector<double>& b);

std::vector<double> VectorSubstract(const std::vector<double>& a, const std::vector<double>& b);

double VectorToScalar(const std::vector<double>& a, const std::vector<double>& b);

std::vector<double> VectorScalar(const double a, const std::vector<double>& b);

std::vector<std::vector<double> > MatrixSclar(const double a, const std::vector<std::vector<double> >& b);

std::vector<double> MatrixVector(const std::vector<std::vector<double> >& A, const std::vector<double>& b);

std::vector<std::vector<double> > MatrixMatrix(const std::vector<std::vector<double> >& A, const std::vector<std::vector<double> >& B);

std::vector<double> ResidualError(const std::vector<std::vector<double> >& A, const std::vector<double>& x, const std::vector<double>& b);

double VectorNorm1(const std::vector<double>& a);

double VectorNorm2(const std::vector<double>& a);

double VectorNormInfty(const std::vector<double>& a);

double MatrixNorm1(const std::vector<std::vector<double> >& A);

double MatrixNormInfty(const std::vector<std::vector<double> >& A);

#endif
