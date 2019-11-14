#!/usr/env/bin python3
from math import exp
from math import log10
import matplotlib.pyplot as plt
from typing import List
td = 1
nd = 100
h = td/nd
nt = 20*nd

c0 = 0.1

a = 1


def dcdt(Ci: float, Ci_m: float):
    return a*Ci*(1-Ci_m)*(1+Ci_m)


def main()->None:
    C_list: List[float] = [c0]
    time_list: List[float] = [-td]
    for i in range(1, nd+1):
        t = time_list[-1] + h
        time_list.append(t)
        C_list.append(c0)
    for i in range(nd, nd+20*nd):
        t = time_list[-1] + h
        c_mid = C_list[-1] + h * dcdt(C_list[i], C_list[i-nd])
        nextC = C_list[-1] + h/2.0*(dcdt(C_list[i], C_list[i-nd]) + dcdt(c_mid, C_list[i-nd+1]))
        C_list.append(nextC)
        time_list.append(t)
    plt.plot(time_list[0:-1], C_list[0:-1], color='k')
    # plt.title('g = %d, a = %d, lambda = %d, c0 = %.2e, m = %.2e' % (g, a, lam, c0, m))
    plt.show()
    print('C_Nt= %.8e' % C_list[nt])


if __name__ == '__main__':
    main()
