import matplotlib.pyplot as plt
import math
import numpy as np

a = 3; T = 4.0; x0 = 1;
t = 0; x = x0;
N = 200; tau = T/N;
tt = np.linspace(0, T, N+1)
xx = []; x_ex = [];
xx.append(x0); x_ex.append(x0);

for i in range(N):
    x_new = x + tau*a*x
    t = t + tau
    x_ex.append(x0*math.exp(a*t))
    xx.append(x_new)
    x = x_new

plt.xkcd()
plt.plot(tt, xx, color='b')
plt.plot(tt, x_ex, color='r')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.show()
