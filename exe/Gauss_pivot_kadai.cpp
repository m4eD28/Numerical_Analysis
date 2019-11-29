#include "../src/LinearAlgebra.hpp"
#include "../src/Algo.hpp"
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
  std::ofstream ofs("test.csv");
  ofs << "1" << "," << "2" << std::endl;
  int n = 500;
  int count = 100;
  std::vector<std::vector<double>> A(n, std::vector<double>(n, rnd()));
  std::vector<double> b(n);
  for (double i = 0 ; i < count; i++) {
    for (std::vector<double>& vec : A) {
      for (double& elem : vec) {
        elem = rnd();
      }
    }
    for (double& elem : b) {
      elem = rnd();
    }

    /* std::cout << Gauss_elimination_norm(A, b, n) << std::endl; */
    /* std::cout << Gauss_elimination_pivot_norm(A, b, n) << std::endl; */

    ofs << Gauss_elimination_norm(A, b, n) << ","  << Gauss_elimination_pivot_norm(A, b, n) << std::endl;
    fprintf(stderr, "Running : %2.2f%%\r", i/count*100);
  }

  return 0;
}
