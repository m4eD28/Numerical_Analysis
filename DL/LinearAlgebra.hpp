
#ifndef LinearAlgebra_hpp
#define LinearAlgebra_hpp

#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

void printVector(vector<double>);
void printMatrix(vector< vector<double> >);
vector<double> VectorSubtract(vector<double>, vector<double>);
vector<double> MatrixVector(vector< vector<double> >, vector<double>);
vector< vector<double> > MatrixMatrix(vector< vector<double> >, vector< vector<double> >);
vector<double> ResidualError(vector< vector<double> >, vector<double>, vector<double>);
double VectorNorm1(vector<double>);
double VectorNorm2(vector<double>);
double VectorNormInfty(vector<double>);
double MatrixNorm1(vector< vector<double> >);
double MatrixNormInfty(vector< vector<double> >);

////////////////////

vector<double> GaussianElimination(vector< vector<double> >, vector<double>);
vector<double> Backward_sub(vector< vector<double> >, vector<double>);
vector< vector<double> > LU_decomposition(vector< vector<double> >);
vector<double> Forward_sub(vector< vector<double> >, vector<double>);
vector< vector<double> > Inverse_matrix(vector< vector<double> >);

/*
 void printVector(double*, int);
 void printMatrix(double**, int, int);
 void copyVector(double*, double**, int);
 void copyMatrix(double**, double**, int , int);
*/

#endif /* LinearAlgebra_hpp */
