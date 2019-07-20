import matplotlib.pyplot as plt
import math
import numpy as np

def logistic_analysitc(u0, r, K, T):
    return K / (1 + (K/u0 - 1) * math.exp(-1*r*T))

def euler_logistic(u0, r, K, T, N, color='k'):
    tau = T/N
    tt = np.linspace(0, T, N+1)
    uu = []
    uu = np.append(uu, u0)
    u, t= u0, 0
    for i in range(N):
        u_new = u + tau*r*u*(1 - u/K)
        t = t + tau
        uu = np.append(uu, u_new)
        u = u_new

    plt.plot(tt, uu, color=color)
    plt.xlabel('t')
    plt.ylabel('u(t)')
    plt.title('u(t), u0 = %d, r = %d, K = %d' % (u0, r, K))
    return uu

def error(T, N):
    uu = euler_logistic(1, 1, 10, T, N)
    uu_analystic = logistic_analysitc(1, 1, 10, T)
    error = abs(uu[N] - uu_analystic)
    h = T/N
    # print('1 = %.10e' % uu_analystic)
    # print('2 = %.10e' % uu[N])
    print('N = %d, h = %.e, |u(N) - u(T)| = %.2e' % (N, h, error))

if __name__ == '__main__':
    # (1) -------------------------------------
    uu = euler_logistic(1, 1, 10, 8, 400, )
    plt.show()
    print('u(n) = %.10e' % (uu[400]))

    # (2) ---------------------------------------
    # plt.subplots_adjust(wspace=1.6, hspace=0.6)
    plt.figure(figsize=(20, 15), dpi=50)
    plt.subplot(4, 3, 1)
    euler_logistic(5, -1, 7, 8, 400, )
    plt.subplot(4, 3, 2)
    euler_logistic(5, 0.5, 7, 8, 400, )
    plt.subplot(4, 3, 3)
    euler_logistic(5, 3, 7, 8, 400, )

    plt.subplot(4, 3, 4)
    euler_logistic(5, -1, 20, 8, 400, )
    plt.subplot(4, 3, 5)
    euler_logistic(5, 0.5, 20, 8, 400, )
    plt.subplot(4, 3, 6)
    euler_logistic(5, 3, 20, 8, 400, )

    plt.subplot(4, 3, 7)
    euler_logistic(30, -1, 7, 8, 400, )
    plt.subplot(4, 3, 8)
    euler_logistic(30, 0.5, 7, 8, 400, )
    plt.subplot(4, 3, 9)
    euler_logistic(30, 3, 7, 8, 400, )

    plt.subplot(4, 3, 10)
    euler_logistic(30, -1, 20, 8, 400, )
    plt.subplot(4, 3, 11)
    euler_logistic(30, 0.5, 20, 8, 400, )
    plt.subplot(4, 3, 12)
    euler_logistic(30, 3, 20, 8, 400, )
    # plt.xlabel('t')
    # plt.ylabel('u(t)')
    plt.tight_layout()
    plt.show()

    # (3) --------------------------------
    error(8, 100)
    error(8, 200)
    error(8, 400)
    error(8, 800)
    error(8, 1600)

