#include<iostream>
#include<cmath>
using namespace std;
double machine_eps(){
  double eps = 1.0;
  while(!(eps + 1.0 == 1.0))eps /= 2;
  return eps * 2;
}

int main(){
  double eps= 1.0;
  while(!(eps + 1.0 == 1.0)){
    eps /= 2;
  }
  printf("Epsilon = %.20e\n",machine_eps());
  printf("2 ^ -52= %.20e\n",pow(2,-52));
  printf("Epsilon = %.20e\n",2 * eps);
  printf("Epsilon = %g\n",2 * eps);
  printf("Epsilon = %g\n",pow(2,-52));

  printf("1.0 + 1.0 /2 * Epsilon = %.20e\n",1.0 + 1.0 /2.0 * machine_eps());
  printf("1.0 + 1.0 /2 * Epsilon = %.20e\n",1.0 /2.0 * machine_eps()+ 1);
  printf("1.0 + 3.0 /4 * Epsilon = %.20e\n",1.0 + 3.0 /4.0 * machine_eps());
  printf("1.0 + Epsilon = %.20e\n",1.0 + machine_eps());
  printf("2.0 - 1.0 / 2 * Epsilon = %.20e\n",2.0 - 1.0 /2.0 * machine_eps());
  printf("2.0 - 3.0 / 4.0 * Epsilon = %.20e\n",2.0 - 3.0 /4.0 * machine_eps());
  printf("2.0 + Epsilon = %.20e\n",2.0 + machine_eps());

  /* cout << pow(10,45) - pow(10,45) - pow(10,35)+ pow(10,35)+ 123.4 + 432.1 << endl; */
  /* cout << 123.4 + 432.1 + pow(10,45) - pow(10,45) - pow(10,35)+ pow(10,35)  << endl; */
  /* cout << pow(10,45) - pow(10,35) - pow(10,45) + pow(10,35)+ 123.4 + 432.1 << endl; */
}
