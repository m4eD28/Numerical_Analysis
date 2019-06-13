#include<iostream>
#include<cmath>
using namespace std;

int toi1(){
  cout << "問1" << endl;
  double a = 0.1, b = 0.2, c = 0.3;
  printf("(a+b)+c = %.16e\n",(a+b)+c);
  printf("a+(b+c) = %.16e\n",a+(b+c));

  return 0;
}

double toi2_fuc_1(double x){
  return (1 - cos(x)) / x*x;
}

double toi2_fuc_2(double x){
  return 2 * sin(x/2) * sin(x/2) / x*x;
}

int toi2(){
  cout << "問2" << endl;
  double x = pow(10,-7);
  printf("f(x) = %.16e\n",toi2_fuc_1(x));
  printf("f(x) = %.16e\n",toi2_fuc_2(x));

  return 0;
}

double toi3_fuc_1(long n){
  double sum = 0;
  for (double i = 0; i <= n; i++) {
    sum += 2 / (4 * i + 1) / (4 * i + 3);
  }
  return 4 * sum;
}

double toi3_fuc_2(long n){
  double sum = 0;
  for (double i = n; i >= 0; i--) {
    sum += 2 / (4 * i + 1) / (4 * i + 3);
  }

  return 4 * sum;
}

double toi3(){
  cout << "問3" << endl;
  long n[3];
  n[0] = 5000000;
  n[1] = 50000000;
  n[2] = 500000000;

  cout << "(1)" << endl;

  for (int i = 0; i < 3; i++) {
    printf("4Sn_%9ld = %.16e\n", n[i], toi3_fuc_1(n[i]));
  }
  cout << "(2)" << endl;

  for (int i = 0; i < 3; i++) {
    printf("4Sn_%9ld = %.16e\n", n[i], toi3_fuc_2(n[i]));
  }
  printf("Pi = %.16e\n",M_PI);

  if ((toi3_fuc_1(n[2]) - M_PI) >= (toi3_fuc_2(n[2])-M_PI)) {
    cout << "(1)の方が精度高め" << endl;
  } else {
    cout << "(2)の方が精度高め" << endl;
  }
  return 0;
}

int main(){
  toi1();
  toi2();
  toi3();
}
