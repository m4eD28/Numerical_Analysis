import matplotlib.pyplot as plt
import math
import numpy as np


def euler_lotka_volterra(u0, v0, a1, b1, c1, a2, b2, c2, T, N, color='k'):
    tau = T/N
    uu = np.array([], dtype='float128')
    uu = np.append(uu, u0)
    vv = np.array([], dtype='float128')
    vv = np.append(vv, v0)
    u, v, t = u0, v0, 0
    for i in range(N):
        u_new = u + tau*u*(a1 + b1*u + c1*v)
        v_new = v + tau*v*(a2 + b2*u + c2*v)
        uu = np.append(uu, u_new)
        vv = np.append(vv, v_new)
        u = u_new
        v = v_new

    # plt.plot(uu, vv, color=color)
    # plt.title('x0 = %.2e, y0 = %.2e, mu = %.2e, omega = %.2e' % (x0, y0, mu, omega))

    plt.figure(figsize=(14, 10), dpi=300)
    tt = np.linspace(0, T, N+1)
    plt.subplot2grid((2, 4), (0, 0), colspan=2, rowspan=2)
    plt.plot(uu, vv, color=color)
    plt.xlabel('u(t)')
    plt.ylabel('v(t)')
    plt.title('u0 = %.2f, v0 = %.2f, a1 = %.2f, b1 = %.2f c1 = %.2f, a2 = %.2f, b2 = %.2f, c2 = %.2f' % (u0, v0, a1, b1, c1, a2, b2, c2))
    plt.subplot2grid((2, 4), (0, 2), colspan=2, rowspan=1)
    plt.plot(tt, uu, color=color)
    plt.xlabel('t')
    plt.ylabel('u(t)')
    plt.subplot2grid((2, 4), (1, 2), colspan=2, rowspan=1)
    plt.plot(tt, vv, color=color)
    plt.xlabel('t')
    plt.ylabel('v(t)')
    plt.tight_layout()
    str = ('u%.1f v%.1f a1%.1f b1%.1f c1%.1f a2%.1f b2%.1f c2%.1f' % (u0, v0, a1, b1, c1, a2, b2, c2))
    str = str.replace('.', '_')
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
    # (1) -------------------------------------
    # uu = euler_logistic(1, 1, 10, 8, 400, )
    # plt.show()
    # print('u(n) = %.10e' % (uu[400]))

    # (2) ---------------------------------------
    # plt.subplots_adjust(wspace=1.6, hspace=0.6)
    # plt.figure(figsize=(20, 15), dpi=50)
    # plt.subplot(4, 3, 1)

    euler_lotka_volterra(10, 10, 2.25, -0.2, -0.5, 2.35, -0.15, -0.45, 1000, 1000, )

    # plt.tight_layout()
    # plt.show()

    # (3) --------------------------------
    # error(8, 100)
    # error(8, 200)
    # error(8, 400)
    # error(8, 800)
    # error(8, 1600)
