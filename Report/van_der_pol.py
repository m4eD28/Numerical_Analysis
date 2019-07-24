import matplotlib.pyplot as plt
import math
import numpy as np


def euler_van_der_pl(x0, y0, mu, omega, T, N, color='k'):
    tau = T/N
    xx = np.array([], dtype='float128')
    # xx = []
    xx = np.append(xx, x0)
    yy = np.array([], dtype='float128')
    # yy = []
    yy = np.append(yy, y0)
    x, y, t= x0, y0, 0
    for i in range(N):
        x_new = x + tau*y
        y_new = y + tau*(mu * (1 - x**2) * y - omega**2 * x)
        xx = np.append(xx, x_new)
        yy = np.append(yy, y_new)
        x = x_new
        y = y_new

    # plt.plot(xx, yy, color=color)
    # plt.title('x0 = %.2e, y0 = %.2e, mu = %.2e, omega = %.2e' % (x0, y0, mu, omega))

    plt.figure(figsize=(14, 10), dpi=300)
    tt = np.linspace(0, T, N+1)
    plt.subplot2grid((2, 4), (0, 0), colspan=2, rowspan=2)
    plt.plot(xx, yy, color=color)
    plt.xlabel('x(t)')
    plt.ylabel('y(t)')
    plt.title('x0 = %.2f, y0 = %.2f, mu = %.2f, omega = %.2f' % (x0, y0, mu, omega))
    plt.subplot2grid((2, 4), (0, 2), colspan=2, rowspan=1)
    plt.plot(tt, xx, color=color)
    plt.xlabel('t')
    plt.ylabel('x(t)')
    plt.subplot2grid((2, 4), (1, 2), colspan=2, rowspan=1)
    plt.plot(tt, yy, color=color)
    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.tight_layout()
    str = ('x%.1fy%.1fmu%.1fomega%.1f' % (x0, y0, mu, omega))
    str = str.replace('.', '_')
    # str = str.replace('-', '')
    plt.savefig(str+'.png')
    return xx, yy

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

    euler_van_der_pl(1, 1, 1, 1, 100, 1000, )

    euler_van_der_pl(1, 1, -1, 1, 100, 1000, )
    euler_van_der_pl(1, 1, 0.25, 1, 100, 1000, )
    euler_van_der_pl(1, 1, 1.25, 1, 100, 1000, )
    euler_van_der_pl(1, 1, 2, 1, 100, 1000, )

    euler_van_der_pl(1, 1, 1, -1, 100, 1000, )
    euler_van_der_pl(1, 1, 1, 0.25, 100, 1000, )
    euler_van_der_pl(1, 1, 1, 1.25, 100, 1000, )
    euler_van_der_pl(1, 1, 1, 2, 100, 1000, )

    euler_van_der_pl(3, 3, 1, 1, 100, 1000, )
    euler_van_der_pl(-3, -3, 1, 1, 100, 1000, )

    # plt.tight_layout()
    # plt.show()

    # (3) --------------------------------
    # error(8, 100)
    # error(8, 200)
    # error(8, 400)
    # error(8, 800)
    # error(8, 1600)
