#include <iostream>
#include <cmath>
using namespace std;

int main() {
  double kakko1;
  double kakko2;
  cout << log(M_E) << endl;
  kakko1 = pow(1-(log(5)/10000), 10000);
  printf("(1) = %.1e\n",kakko1);
  kakko2 = pow(pow(3, -147), 100);
  printf("(2) = %.1e", kakko2);
  kakko2 = pow(3, -14700);
  printf("(2) = %.1e", kakko2);
}
