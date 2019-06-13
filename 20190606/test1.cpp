#include <iostream>
using namespace std;

int main() {
  double a1 = 1.0 + 1e-16;
  double b1 = 1.0 + 1.2*1e-16;
  double c1 = 1.0 + 1.1*1e-16;
  printf("a1 = %.20e\n",a1);
  printf("b1 = %.20e\n",b1);
  printf("c1 = %.20e\n",c1);
}
