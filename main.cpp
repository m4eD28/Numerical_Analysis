#include <iostream>
#include "Newton.h"
using namespace std;

double f(double x){
  return pow(x,4) + 5*pow(x,3) + 6*pow(x,2) - 4*x - 8;
}

double df(double x) {
  return 4*pow(x,3) + 15*pow(x,2) + 12*x - 4;
}

int main() {
  Newton_zansa();
  Newton(-30, -2);

  return 0;
}
