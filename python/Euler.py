import matplotlib.pyplot as plt
import math
import numpy as np

a, b, tau, T = -0.1, 0.5, 0.1, 100.0
u0, w0, t = 0.5, 0.3, 0
u, w = u0, w0
N = 20000
h = T/N
tt = np.linspace(0, T, N+1)
uu = []
ww = []
uu.append(u0), ww.append(w0)

for i in range(N):
    u_new = u + h*(u*(u-a)*(1-u) - w)
    w_new = w + h*tau*(u - b*w)
    uu.append(u_new); ww.append(w_new)
    u = u_new; w = w_new;

plt.plot(tt, uu, color='b', label='u')
plt.plot(tt, ww, color='r', label='w')
plt.legend()
plt.xlabel('t')
plt.ylabel('u, w')
plt.show()
