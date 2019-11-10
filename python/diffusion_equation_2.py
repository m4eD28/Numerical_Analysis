"""One dimention diffusion quation"""
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    """
    main
    """
    time_end: float = 0.13
    alpha: float = 0
    beta: float = 1
    n_t: int = 40000
    n_x: int = 100
    step: float = 1 / n_x
    tau: float = time_end / n_t
    x_line = np.linspace(0, 1, n_x + 1)
    u_pre: np.array = np.full(n_x + 1,
                              (200 * np.exp(-100 * np.power((x_line-0.5), 2))),
                              dtype='float')
    # u_pre: np.array = np.full(n_x + 1,
    #                           np.power(x_line, 2) / 5,
    #                           dtype='float')
    u_array: np.array = np.zeros_like(u_pre, dtype='float')
    f_k: np.array = np.full_like(u_pre, 0, dtype='float')
    # time: float = np.linspace(0, time_end, n_x)
    time: float = 0

    for k in range(n_t):
        time = time + tau
        u_array[0] = alpha
        u_array[-1] = beta
        for i in range(n_x):
            f_k[i] = np.power(u_pre[i], 2) / 5
        for i in range(1, len(u_array)-1):
            u_array[i] = u_pre[i] +\
                         tau / step**2 *\
                         (u_pre[i+1] - 2 * u_pre[i] + u_pre[i-1]) +\
                         tau * f_k[i]
        if (k + 1) % 100 == 0:
            plt.plot(x_line, u_array, color='b')
            plt.xlabel('x')
            plt.ylabel('u(x)')
            plt.title('t = {0}'.format(time))
                # plt.show()
        u_pre = u_array
    plt.show()

if __name__ == '__main__':
    main()
