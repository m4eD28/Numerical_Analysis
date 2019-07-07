#include<iostream>
#include "src/matplotlibcpp.h"
#include "src/Algo.hpp"

namespace plt = matplotlibcpp;

int main(){
  std::cout<<"Plot start"<<std::endl;

  double a = 3.0, T = 4.0, x0 = 1.0;
  int N = 200;
  double h = T/N, t = 0;
  std::vector<double> xx(N), x_ex(N);
  xx.at(0) = x0;
  x_ex.at(0) = x0;
  for(int i=0; i<N-1; ++i) {
    xx.at(i+1) = xx.at(i) + h*a*xx.at(i);
    t += h;
    x_ex.at(i+1) = x0*pow(M_E, a*t);
  }

  plt::plot(xx, "b");
  plt::plot(x_ex, "r");
	plt::title("Plot");
	plt::xlabel("t");
	plt::ylabel("x(t)");
  plt::show();

  return 0;
}
