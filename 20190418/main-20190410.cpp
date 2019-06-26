//
//  main.cpp
//  NumericalEx1
//
//  Created by Guanyu Zhou on 2019/04/08.
//  Copyright © 2019年 Guanyu Zhou. All rights reserved.
//

#include <iostream>
#include<vector>
#include<cmath>

using namespace::std;


double sum(vector<double> a){
    double s = 0;
    for(int i = 0;i<a.size();i++){
        s = s + a[i];
    }
    return s;
}

int main(int argc, const char * argv[]) {
    cout << "Hello, World!\n";
    int aa = 12345;
    double bb = 12345.6789;
    printf("%8d, %9.3f, %.10e\n", aa, bb, bb);

    cout<<"sin(aa)="<<sin(aa)<<endl;


    int x = 3, m = 6, n = 7;
    double y = sin(M_PI/2);
    cout<<"x+y="<<exp(x)+y<<endl;
    cout<<"m/n="<<m/n<<endl;
    printf("%7d, %9.3f, %.15e\n",12345,12345.6789,12345.6789);

    double d[3] = {1,11,111};
    double B[3][2]={
        {1,2},
        {3,4},
        {5,6}
    };
    d[2] = 111.111; B[0][0] = 1.1;
    cout<<"d=\n";
    for(int i = 0;i<3;i++){
        printf("d[%d] = %.6f\n",i,d[i]);
    }
    cout<<"B=\n";
    for(int i = 0;i<3;i++){
        for(int j = 0;j<2;j++){
            printf("B[%d][%d] = %.6f ",i,j,B[i][j]);
        }
        cout<<"\n";
    }
    double *pd = new double[m];
    pd[2] = 8.8;
    for(int i=0;i<m;i++){
        cout<<"pd[i]="<<pd[i]<<endl;
    }
    delete[] pd;

    vector<double> b(n,log(1000));
    cout<<"b = \n";
    for(int i=0;i<b.size();i++){
        printf("b[%d] = %.3e\n",i,b[i]);
    }
    for(int i=0;i<b.size();i++){
        b[i] = i*i;
    }
    cout<<"b = \n";
    for(int i=0;i<b.size();i++){
        printf("b[%d] = %.3e\n",i,b[i]);
    }

    vector<double> a(5), c(5);
    vector<vector<double>> A(5,vector<double>(3));
    for(int i=0;i<A.size();i++){
        a[i] = 10; c[i] = 20;
        for(int j=0;j<A[0].size();j++){
            A[i][j] = i+j;
        }
    }
    cout<<"A = \n";
    A[2][2] = pow(10,10);
    for(int i=0;i<A.size();i++){
        for(int j=0;j<A[0].size();j++){
            printf("%.2e\t",A[i][j]);
        }
        printf("\n");
    }
    cout<<"sum(a)="<<sum(a)<<", sum(c)="<<sum(c)<<endl;

    return 0;
}


