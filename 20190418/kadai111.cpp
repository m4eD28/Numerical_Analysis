#include<iostream>
#include<cmath>
#define MAX_SIZE 1000000

double func(int n){
  if(n == 0) return 1.00000;
  return double(1 / (1 + double(func(n-1))) + 1);
}

int main(){
  int n;
  printf("Number? : ");
  std::cin >> n;
  printf(" a[%2d]  = %.18e \nsqrt(2) = %.18e \n",n,func(n), sqrt(2));
  for (int i = 0; i <= n; i++){
    printf(" a[%2d]  = %.18e \n",i,func(i));
  }
  for (int i = 0; i <= n; i++){
    printf(" |a[%2d] - a| = %.2e \n",i,fabs(func(i)-sqrt(2)));
  }

  return 0;
}
