import matplotlib.pyplot as plt
import numpy as np
import sys
import copy


def euler_van_der_pl(x0, y0, mu, omega, T, N, color='k'):
    tau = T/N
    xx = np.array([], dtype='float128')
    xx = np.append(xx, x0)
    yy = np.array([], dtype='float128')
    yy = np.append(yy, y0)
    x, y = x0, y0
    for i in range(N):
        x_new = x + tau*y
        y_new = y + tau*(mu * (1 - x**2) * y - omega**2 * x)
        xx = np.append(xx, x_new)
        yy = np.append(yy, y_new)
        x = x_new
        y = y_new
        sys.stdout.write("\r%f" % float(i/N*100))
        sys.stdout.flush()

    plt.figure(figsize=(14, 10), dpi=300)
    tt = np.linspace(0, T, N+1)
    plt.subplot2grid((2, 4), (0, 0), colspan=2, rowspan=2)
    plt.plot(xx, yy, color=color)
    plt.xlabel('x(t)')
    plt.ylabel('y(t)')
    # plt.title('x0 = %.2f, y0 = %.2f, mu = %.2f, omega = %.2f, T = %.2e, N = %.2e' % (x0, y0, mu, omega, T, N))
    plt.subplot2grid((2, 4), (0, 2), colspan=2, rowspan=1)
    plt.plot(tt, xx, color=color)
    plt.xlabel('t')
    plt.ylabel('x(t)')
    plt.subplot2grid((2, 4), (1, 2), colspan=2, rowspan=1)
    plt.plot(tt, yy, color=color)
    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.tight_layout()
    str = ('x%.1fy%.1fmu%.1fomega%.1ft%.2en%.2e' % (x0, y0, mu, omega, T, N))
    str = str.replace('.', ',')
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
    # T = 100
    # x1, y1 = euler_van_der_pl(1, 1, 1, 1, 100, 1000, )
    # plt.show()
    # x1 = copy.deepcopy(x1)
    # y1 = copy.deepcopy(y1)
    # x2, y2 = euler_van_der_pl(1, 1, 1, 1, 100, 10000, )
    # print('(xの誤差) = %f, (yの誤差) = %f' % (abs(x1[-1]-x2[-1]), abs(y1[-1]-y2[-1])))
    # plt.show()
    # # plt.subplot2grid((2, 4), (0, 0), colspan=2, rowspan=2)
    # # plt.plot(x1, y1)
    # # plt.subplot2grid((2, 4), (0, 2), colspan=2, rowspan=2)
    # # plt.plot(x2, y2)
    # # plt.tight_layout()
    # # plt.show()
    # # plt.savefig('van_der_pol_error.png')

    # N = 1000
    # euler_van_der_pl(1, 1, 1, 1, 100, 1000, )

    # euler_van_der_pl(1, 1, -1, 1, 100, 1000, )
    euler_van_der_pl(1, 1, 0, 1, 100, 1000, )
    # # euler_van_der_pl(1, 1, 0.25, 1, 100, 1000, )
    # # euler_van_der_pl(1, 1, 1.25, 1, 100, 1000, )
    # euler_van_der_pl(1, 1, 2, 1, 100, 1000, )

    # euler_van_der_pl(1, 1, 1, -1, 100, 1000, )
    # euler_van_der_pl(1, 1, 1, 0, 100, 1000, )
    # # euler_van_der_pl(1, 1, 1, 0.25, 100, 1000, )
    # # euler_van_der_pl(1, 1, 1, 1.25, 100, 1000, )
    # euler_van_der_pl(1, 1, 1, 2, 100, 1000, )

    # euler_van_der_pl(3, 3, 1, 1, 100, 1000, )
    # euler_van_der_pl(-3, -3, 1, 1, 100, 1000, )

    # N = 2000
    # euler_van_der_pl(1, 1, 1, 1, 100, 2000, )

    # euler_van_der_pl(1, 1, -1, 1, 100, 2000, )
    # euler_van_der_pl(1, 1, 0, 1, 100, 2000, )
    # # euler_van_der_pl(1, 1, 0.25, 1, 100, 2000, )
    # # euler_van_der_pl(1, 1, 1.25, 1, 100, 2000, )
    # euler_van_der_pl(1, 1, 2, 1, 100, 2000, )

    # euler_van_der_pl(1, 1, 1, -1, 100, 2000, )
    # euler_van_der_pl(1, 1, 1, 0, 100, 2000, )
    # # euler_van_der_pl(1, 1, 1, 0.25, 100, 2000, )
    # # euler_van_der_pl(1, 1, 1, 1.25, 100, 2000, )
    # euler_van_der_pl(1, 1, 1, 2, 100, 2000, )

    # euler_van_der_pl(3, 3, 1, 1, 100, 2000, )
    # euler_van_der_pl(-3, -3, 1, 1, 100, 2000, )

    # N = 10000
    # euler_van_der_pl(1, 1, 1, 1, 100, 10000, )

    # euler_van_der_pl(1, 1, -1, 1, 100, 10000, )
    # euler_van_der_pl(1, 1, 0, 1, 100, 10000, )
    # # euler_van_der_pl(1, 1, 0.25, 1, 100, 10000, )
    # # euler_van_der_pl(1, 1, 1.25, 1, 100, 10000, )
    # euler_van_der_pl(1, 1, 2, 1, 100, 10000, )

    # euler_van_der_pl(1, 1, 1, -1, 100, 10000, )
    # euler_van_der_pl(1, 1, 1, 0, 100, 10000, )
    # # euler_van_der_pl(1, 1, 1, 0.25, 100, 10000, )
    # # euler_van_der_pl(1, 1, 1, 1.25, 100, 10000, )
    # euler_van_der_pl(1, 1, 1, 2, 100, 10000, )

    # euler_van_der_pl(3, 3, 1, 1, 100, 10000, )
    # euler_van_der_pl(-3, -3, 1, 1, 100, 10000, )
