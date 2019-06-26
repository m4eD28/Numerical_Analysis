#include<iostream>
using namespace std;
#define MAX_SIZE 1000000

int main(){
  double a[MAX_SIZE];
  a[0] = 1;
  for(int i = 0;i <= MAX_SIZE-1;i++){
    a[i+1] = 1 + 1 / double(1 + a[i]);
  }

  for (int i = 0; i <= 23; i++) {
    printf("a[%d] = %.18e\n",i,a[i]);
  }

  printf("a[%d] = %.18e\n",MAX_SIZE,a[MAX_SIZE-1]);

  for (int i = 0; i <= 23; i++) {
    printf("|a[%d] - a| = %.18e\n",i,double(a[i] - a[MAX_SIZE-1]));
  }

  return 0;
}
