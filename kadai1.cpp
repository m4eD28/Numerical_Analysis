#include<iostream>
#include<vector>
using namespace std;

int main(){
  vector<double> a[23];
  a[0] = 1;
  for(int i = 0;i < a.size();i++){
    a[i+1] = 1 + 1 / (1 + a[i]);
    cout << a[i+1] << endl;
  }

  return 0;
}
