#include<iostream>
#include "../src/matplotlibcpp.h"
#include "../src/Algo.hpp"

namespace plt = matplotlibcpp;

int main(){
  std::cout<<"Plot another window"<<std::endl;

  double r = 3.0, T = 8.0, u0 = 1.0;
  int N = 400;
  std::vector<double> u_400;
  u_400.emplace_back(u0);
  for(int i=0; i<N; ++i) {
    u_400.emplace_back(u0*i);
  }

  plt::plot(u_400, "b");
	plt::title("Plot");
	plt::xlabel("t");
	plt::ylabel("x(t)");
  plt::show();

  return 0;
}
