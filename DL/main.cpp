
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "LinearAlgebra.hpp"

using namespace std;

int main(int argc, const char * argv[]){

    // vector< vector<double> > B = {{3, -1, -1}, {-1, 3, -1}, {-1, -1, 3}};
    vector< vector<double> > A = {{2, 3, 3}, {3, 2, -1}, {5, 4, 2}};
    // vector<double> b = {1, 1, 1};
    vector<double> b = {5, -4, 3};

    // cout << "A" << endl;
    // printMatrix(A);
    // cout << "b" << endl;
    // printVector(b);

    // vector<double> x_1 = GaussianElimination(A, b);
    // cout << "x_1" << endl;
    // printVector(x_1);
    

    // vector< vector<double> > LU = LU_decomposition(A);
    // cout << "LU" << endl;
    // printMatrix(LU);

    // vector<double> y = Forward_sub(LU, b);
    // cout << "y" << endl;
    // printVector(y);


    // vector<double> x_2 = Backward_sub(LU, y);
    // cout << "x_2" << endl;
    // printVector(x_2);


    vector< vector<double> > A_Inverse;
    A_Inverse = Inverse_matrix(A);
    cout << "A_Inverse" << endl;
    printMatrix(A_Inverse);

    // A * A_Inverse = I を確認
    printMatrix(MatrixMatrix(A, A_Inverse));



    return 0;
}
