import matplotlib.pyplot as plt
import math
import numpy as np

T, N = 8, 400
u0, r, K = 1, 1, 10
tau = T/N
tt = np.linspace(0, T, N+1)
uu = []
uu.append(u0)

for i in range(N):
    u_new = u + tau*
    t = t + tau
    x_ex.append(x0*math.exp(a*t))
    xx.append(x_new)
    x = x_new
    x_1.append(np.sin(2*i))

fig = plt.figure(figsize=(10, 8))
ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=2)


ax1.plot(xx, linewidth=2, color='r')
ax1.grid(True)
ax2.plot(x_ex, linewidth=2, color='r')
ax2.grid(True)
ax3.plot(x_1, linewidth=2, color='r')
ax3.grid(True)

fig.show()
plt.show()
