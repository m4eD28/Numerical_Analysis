#include "../src/LinearAlgebra.hpp"
#include "../src/Algo.hpp"
#include <cstdlib>
#include <iostream>
#include <random>
#include <fstream>
#include <string>

std::random_device rnd_dev;
std::mt19937 mt(rnd_dev());
std::uniform_real_distribution<> dist(0, 1);
inline double rnd() {
  return dist(mt);
}

int main() {
  int n = 500;
  int count = 100;
  std::vector<std::vector<double>> A(n, std::vector<double>(n));
  std::vector<double> b(n);
  std::vector<double> norm_no_pivot(count);
  std::vector<double> norm_pivot(count);

  for (double i = 0 ; i < count; ++i) {
    for (std::vector<double>& vec : A) {
      for (double& elem : vec) {
        elem = rnd();
        /* elem = (double)rand()/RAND_MAX; */
      }
    }
    for (double& elem : b) {
      elem = rnd();
      /* elem = (double)rand()/RAND_MAX; */
    }

      /* std::cout << Gauss_elimination_norm(A, b, n) << std::endl; */
      /* std::cout << Gauss_elimination_pivot_norm(A, b, n) << std::endl; */
      /* ofs << Gauss_elimination_norm(A, b) << ","  << Gauss_elimination_pivot_norm(A, b) << std::endl; */

    norm_no_pivot.at(i) = Gauss_elimination_norm(A, b);
    norm_pivot.at(i) = Gauss_elimination_pivot_norm(A, b);
    fprintf(stderr, "Running : %2.2f%%\r", i/count*100);
  }
  /* std::cout << Gauss_elimination_norm(A, b) << ","  << Gauss_elimination_pivot_norm(A, b) << std::endl; */

  std::ofstream ofs("output.csv");
  ofs << "1" << "," << "2" << std::endl;
  for (int i = 0;i < count;i++) {
    ofs << norm_no_pivot.at(i) << ","  << norm_pivot.at(i) << std::endl;
  }
  return 0;
}
