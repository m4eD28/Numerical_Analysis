#include <stdio.h>

int main(void){
  double eps = 1.0;
  for (;;) {
    double nextEps = eps / 2;
    if (1.0 + nextEps == 1.0)
      break;
    eps = nextEps;
  }
  printf("%g\n", eps);
}
