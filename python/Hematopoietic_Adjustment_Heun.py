import matplotlib.pyplot as plt
import math
import numpy as np

T_delay = 2
b = 20
T = b*T_delay
N_delay = 100
N_T = b*N_delay
h = T_delay/N_delay
M = 100

#
g = 1
a = 1
lam = 2
c0 = 0.6
m = 7


def dcdt(Ci, Ci_m):
    numerator = lam * (a**m) * Ci_m
    denominator = a**m + Ci_m**m
    return numerator / denominator - g * Ci


tt = np.linspace(-T_delay, T, N_delay + N_T + 1)
C = [0]*(N_delay + N_T + 1)

for i in range(N_delay+1):
    C[i] = c0

for i in range(N_delay, N_delay + N_T):
    C_mid = C[i] + h * (dcdt(C[i], C[i-M]))
    C[i+1] = C[i] + h / 2 * (dcdt(C[i], C[i-M]) + (dcdt(C_mid, C[i-M+1])))

plt.plot(tt, C, color='b')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.show()
print('C21_M = %.8e' % C[-1])
