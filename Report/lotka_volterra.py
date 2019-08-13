import matplotlib.pyplot as plt
import math
import numpy as np
import random
import copy
import sys


def euler_lotka_volterra(u0, v0, a1, b1, c1, a2, b2, c2, T, N, color='k'):
    tau = T/N
    uu = np.array([], dtype='float128')
    uu = np.append(uu, u0)
    vv = np.array([], dtype='float128')
    vv = np.append(vv, v0)
    # a1 = random.uniform(-1.0, 1.0)
    # a2 = random.uniform(-1.0, 1.0)
    # b1 = random.uniform(-1.0, 1.0)
    # b2 = random.uniform(-1.0, 1.0)
    # c1 = random.uniform(-1.0, 1.0)
    # c2 = random.uniform(-1.0, 1.0)
    u, v, t = u0, v0, 0
    print(uu[0])
    print(vv[0])
    for i in range(N):
        u_new = u + tau*u*(a1 + b1*u + c1*v)
        v_new = v + tau*v*(a2 + b2*u + c2*v)
        uu = np.append(uu, u_new)
        vv = np.append(vv, v_new)
        u = u_new
        v = v_new
        sys.stdout.write("\r%f" % float(i/N*100))
        sys.stdout.flush()

    plt.figure(figsize=(14, 10), dpi=300)
    tt = np.linspace(0, T, N+1)
    plt.subplot2grid((2, 4), (0, 0), colspan=2, rowspan=2)
    plt.plot(uu, vv, color=color)
    plt.xlabel('u(t)')
    plt.ylabel('v(t)')
    # plt.title('u0 = %.2f, v0 = %.2f, a1 = %.2f, b1 = %.2f c1 = %.2f, a2 = %.2f, b2 = %.2f, c2 = %.2f, T = %.2e, N = %.2e' % (u0, v0, a1, b1, c1, a2, b2, c2, T, N))
    plt.subplot2grid((2, 4), (0, 2), colspan=2, rowspan=1)
    plt.plot(tt, uu, color=color)
    plt.xlabel('t')
    plt.ylabel('u(t)')
    plt.subplot2grid((2, 4), (1, 2), colspan=2, rowspan=1)
    plt.plot(tt, vv, color=color)
    plt.xlabel('t')
    plt.ylabel('v(t)')
    plt.tight_layout()
    str = ('u%.1f v%.1f a1%.1f b1%.1f c1%.1f a2%.1f b2%.1f c2%.1f t%.2e n%.2e' % (u0, v0, a1, b1, c1, a2, b2, c2, T, N))
    str = str.replace('.', ',')
    str = str.replace(' ', '')
    # str = str.replace('-', '')
    plt.savefig(str+'.png')
    return uu, vv

# def error(T, N):
#     uu = euler_logistic(1, 1, 10, T, N)
#     uu_analystic = logistic_analysitc(1, 1, 10, T)
#     error = abs(uu[N] - uu_analystic)
#     h = T/N
#     # print('1 = %.10e' % uu_analystic)
#     # print('2 = %.10e' % uu[N])
#     print('N = %d, h = %.e, |u(N) - u(T)| = %.2e' % (N, h, error))


if __name__ == '__main__':
    # x1, y1 = euler_lotka_volterra(3, 3, 1, 0, -1, -1, 0, 0, 100, 1000, )
    # # plt.show()
    # x1 = copy.deepcopy(x1)
    # y1 = copy.deepcopy(y1)
    # x2, y2 = euler_lotka_volterra(3, 3, 1, 0, -1, -1, 0, 0, 100, 10000, )
    # print('(xの誤差) = %f, (yの誤差) = %f' % (abs(x1[-1]-x2[-1]), abs(y1[-1]-y2[-1])))

    # N = 1000
    # euler_lotka_volterra(3, 3, 22, -0.2, -0.4, 17, -0.5, -0.3, 10, 1000, )
    # euler_lotka_volterra(3, 3, 14, -1.1, -0.9, 20, -0.5, -0.4, 10, 1000, )
    # euler_lotka_volterra(3.0, 3.0, 1.9, 0, -1.901, 1.9, 0, 1.9, 10, 1000, )
    # euler_lotka_volterra(3.0, 3.0, 25, -0.1, -0.4, 20, -0.2, -0.1, 10, 1000, )
    # euler_lotka_volterra(3.0, 3.0, 23.2, -0.6, -0.401, 18.8, -0.5, -0.3, 10, 1000, )
    # euler_lotka_volterra(3.0, 3.0, 1, 0, -1, -1, 1, 0, 10, 1000, )

    # N = 2000
    # euler_lotka_volterra(3, 3, 22, -0.2, -0.4, 17, -0.5, -0.3, 10, 2000, )
    # euler_lotka_volterra(3, 3, 14, -1.1, -0.9, 20, -0.5, -0.4, 10, 2000, )
    # euler_lotka_volterra(3.0, 3.0, 1.9, 0, -1.901, 1.9, 0, 1.9, 10, 2000, )
    # euler_lotka_volterra(3.0, 3.0, 25, -0.1, -0.4, 20, -0.2, -0.1, 10, 2000, )
    # euler_lotka_volterra(3.0, 3.0, 23.2, -0.6, -0.401, 18.8, -0.5, -0.3, 10, 2000, )
    # euler_lotka_volterra(3.0, 3.0, 1, 0, -1, -1, 1, 0, 10, 2000, )
