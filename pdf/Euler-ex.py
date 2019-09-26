#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 15:40:18 2019

@author: zhou
"""

import matplotlib.pyplot as plt
import math
import numpy as np

r = 3; K = 20;
T = 8;
u0 = 30; t = 0;
u = u0;
N = 400; h = T/N;

tt = np.linspace(0, T, N+1)
uu = [];
uu.append(u0);

for i in range(N):
    u_new = u + h*r*u*(1 - u/K)
    uu.append(u_new);
    u = u_new;

print("u_N = {:.10e}".format(u))
print("h = {:.2e}".format(h))
print("|u_N - u(T)| = {:.2e}".format(abs(u - K/(1+ (K/u0 -1)*math.exp(-r*T))) ))

plt.plot(tt,uu,color='b',label = 'u')
plt.title('u0={0}, r={1}, K={2}'.format(u0,r,K),fontsize=20)
plt.legend()
plt.xlabel('t')
plt.ylabel('u')
plt.show()
