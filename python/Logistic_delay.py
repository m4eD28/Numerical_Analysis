import matplotlib.pyplot as plt
import numpy as np

T_delay = 1
b = 20
T = b*T_delay
N_delay = 200
N_T = b*N_delay
h = T_delay/N_delay
r = 1.5
K = 10

tt = np.linspace(-T_delay, T, N_delay + N_T + 1)
x0 = 1
X = [0]*(N_delay + N_T + 1)

for i in range(N_delay+1):
    X[i] = x0

for i in range(N_delay, N_delay + N_T):
    X[i+1] = X[i] + h*r*(1 - X[i-N_delay]/K)*X[i]

plt.plot(tt, X, color='b')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.show()
