#!/usr/env/bin python3
from math import exp
from math import log10
import matplotlib.pyplot as plt
from typing import List
td = 1
M = 10
h = td/M

c0 = 80
K = 100
r = 1
p = 20



def dcdt(Ci: float, Ci_m: float):
    return r*Ci*(1 - Ci_m/K) - p


def main()->None:
    C_list: List[float] = [c0]
    time_list: List[float] = [-td]
    for i in range(1, M+1):
        t = time_list[-1] + h
        time_list.append(t)
        C_list.append(c0)
    for i in range(M, 10*M):
        t = time_list[-1] + h
        c_mid = C_list[-1] + h * dcdt(C_list[i], C_list[i-M])
        nextC = C_list[-1] + h/2.0*(dcdt(C_list[i], C_list[i-M]) + dcdt(c_mid, C_list[i-M+1]))
        C_list.append(nextC)
        time_list.append(t)
    plt.plot(time_list[M:10*M], C_list[M:10*M], color='k')
    # plt.title('g = %d, a = %d, lambda = %d, c0 = %.2e, m = %.2e' % (g, a, lam, c0, m))
    plt.show()
    print('C_10M = %.8e' % C_list[-1])


if __name__ == '__main__':
    main()
