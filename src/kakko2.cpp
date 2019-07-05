#include<iostream>
using namespace std;

int main() {
  double sum1 = 0, sum2 = 0, sum3 = 0;
  int n = 5000000, N = 1000000;
  for (int i = 1; i <= n ; i++) {
    sum1 = sum1 + 1.0/i/i/i;
    sum2 = sum2 + 1.0/(n+1-i)/(n+1-i)/(n+1-i);
  }
  for (int i = 1; i <= N; i++) {
    sum3 = sum3 + 1.0/i/i/i;
  }
  printf("sum1 = %.16e\n", sum1);
  printf("sum2 = %.16e\n", sum2);
  printf("sum3 = %.16e\n", sum3);
  bool flag = (sum1 == sum3);
  cout << flag << endl;
}
