#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 16:07:38 2019

@author: takaseyuusuke
"""

import math
import numpy
from matplotlib import pyplot as plt

def f(x):
    # return pow(x, 2) - 8/x
    return pow(x, 2) - 8/x

def df(x):
    # return pow(x,2) * 2
    return 2*x + 8/pow(x, 2)

if __name__ == '__main__':
    print('Newton')
    x0=10; xs=2; eps=1e-7;Max=60;
    x=x0; N=0; x_N=0;
    Error=[];
    print('step \t x_k \t error \t order')
    for k in range(Max+1):
        if(abs(f(x))<eps):
            N=k+1; x_N=x;
            print('N={},x_N={:.12e}'.format(N,x_N))
            break
        x=x0-f(x0)/df(x0)
        print('{}\t{:.12e}\t{:.2e}\t{:.2e}'.format(k+1,x,abs(x-xs),math.log(abs(x-xs))))
        x0=x;
        Error.append(abs(x-xs));
    plt.plot(range(1,N),Error,marker='o')
    plt.yscale('log')
    plt.show()
