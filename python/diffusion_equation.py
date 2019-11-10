'''One dimention diffusion equation'''
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as linalg
import math
# from typing import List


def analysis(x):
    return x * (1 - x)


def main() -> None:
    T: float = 0.25
    alpha: float = 0
    beta: float = 0
    N_t: int = 5000
    N_x: int = 50
    dt: float = T / N_t
    h: float = 1 / N_x
    u = [0] * (N_x + 1)
    u_pre = [0] * (N_x + 1)
    u0 = [0] * (N_x + 1)
    f_k = [0] * (N_x + 1)
    x = np.linspace(0, 1, N_x + 1)

    t = 0
    for i in range(N_x):
        u0[i] = 200 * pow(math.e, -100 * pow(i * h - 0.5, 2))
        # u0[i] = math.pow(i*h, 2) / 5

    plt.plot(x, u0, color='b')
    plt.xlabel('x')
    plt.ylabel('y(x,t)')
    plt.title('t = {0}'.format(0))
    plt.show()

    u_pre = u0
    for k in range(N_t):
        t = t + dt
        u[0] = alpha
        u[N_x] = beta
        for i in range(N_x):
            f_k[i] = 0
        for i in range(1, N_x):
            u[i] = (dt / pow(h, 2)) * u_pre[i-1] + (1 - 2 * (dt / pow(h, 2))) * u_pre[i] + (dt / pow(h, 2)) * u_pre[i+1] + f_k[i] * dt

        if (k + 1) % 100 == 0:
            plt.plot(x, u, color='b')
            plt.xlabel('x')
            plt.ylabel('u(x)')
            plt.title('t = {0}'.format(t))
        u_pre = u

    plt.show()

if __name__ == '__main__':
    main()
