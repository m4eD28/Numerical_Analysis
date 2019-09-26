//
//  main.cpp
//  FuDou
//
//  Created by Guanyu Zhou on 2019/09/11.
//  Copyright © 2019年 Guanyu Zhou. All rights reserved.
//

#include <iostream>
#include <cmath>

using namespace::std;

int main(int argc, const char * argv[]) {
    //単精度
    float af = 1.0f, bf = 501.023f, cf = 1.0f;
    float xf1 = (-bf - sqrt(bf*bf-4.0f*af*cf))/(2.0f*af);
    float xf2 = (-bf + sqrt(bf*bf-4.0f*af*cf))/(2.0f*af);
    cout<<"単精度\n";
    printf("xf1 = %.12f\n",xf1);
    printf("xf2 = %.12f\n",xf2);

    //倍精度
    double aa = 1.0, bb = 501.023, cc = 1.0;
    double xx1 = (-bb - sqrt(bb*bb-4.0*aa*cc))/(2.0*aa);
    double xx2 = (-bb + sqrt(bb*bb-4.0*aa*cc))/(2.0*aa);
    cout<<"倍精度\n";
    printf("xx1 = %.12f\n",xx1);
    printf("xx2 = %.12f\n",xx2);


    //Logistic Map
    cout<<"Logistic Mapと初期値誤差\n";
    double x0 = 0.4, x0e = x0 + 0.001;
    double xn = 0, xne = 0;
    double a = 3.5;
    for(int i = 0;i<50;i++){
        xn = a*x0*(1-x0);
        x0 = xn;
        xne = a*x0e*(1-x0e);
        x0e = xne;
    }
    printf("xn = %.12f, \t xne = %.12f, \n |xn-xne|= %.12f \n",xn,xne,fabs(xn-xne));

    return 0;
}
