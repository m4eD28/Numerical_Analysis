#!/usr/env/bin python3
from math import exp
from math import log10
import matplotlib.pyplot as plt
from typing import List
# (造血モデル)のパラメータ
g = 1
a = 1
lam = 2
c0 = 0.6
m = 7
# heun法のパラメーター
td = 2
M = 100
h = td/M


def dcdt(Ci: float, Ci_m: float):
    numerator: float = lam * (a**m) * Ci_m
    denominator: float = a**m + Ci_m**m
    return numerator/denominator - g*Ci


def main()->None:
    C_list: List[float] = [c0]
    time_list: List[float] = [-td]
    for i in range(1, M+1):
        t = time_list[-1] + h
        time_list.append(t)
        C_list.append(c0)
    for i in range(M, 21*M):
        t = time_list[-1] + h
        c_mid = C_list[-1] + h * dcdt(C_list[i], C_list[i-M])
        nextC = C_list[-1] + h/2.0*(dcdt(C_list[i], C_list[i-M]) + dcdt(c_mid, C_list[i-M+1]))
        C_list.append(nextC)
        time_list.append(t)
    plt.plot(time_list[M:21*M], C_list[M:21*M], color='k')
    plt.title('g = %d, a = %d, lambda = %d, c0 = %.2e, m = %.2e' % (g, a, lam, c0, m))
    plt.show()
    print('C_21M = %.8e' % C_list[-1])


if __name__ == '__main__':
    main()
