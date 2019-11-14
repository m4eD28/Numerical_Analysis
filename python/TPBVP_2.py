# Two point boundary value probmel
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as linalg
import math
# from typing import List


def analysis(x):
    # return x * (2 - x)
    # return -2 / (1 - np.exp(-4)) * np.exp(-4 * x) - x + 2 / (1 - np.exp(-4))
    return - math.e * (np.exp(-x) - 1) - x


def main() -> None:
    N: int = 10
    h: float = 1 / N
    alpha: float = 0
    beta: float = 0
    p = np.full(N, 1, dtype='float64')
    q = np.full(N, 0, dtype='float64')
    r = np.full(N, -1, dtype='float64')

    df_eq = np.zeros((N, N), dtype='float64')
    for i in range(len(df_eq)):
        for j in range(len(df_eq)):
            if i == j:
                df_eq[i][j] = h**2 * q[i] - 2
            elif i + 1 == j:
                df_eq[i][j] = 1 + (h * p[i]) / 2
            elif i - 1 == j:
                df_eq[i][j] = 1 - (h * p[i]) / 2
    df_eq[-1][-1] = 1
    df_eq[-1][-2] = -1

    ans_vector = np.zeros(N, dtype='float64')
    for i in range(len(ans_vector)):
        if i == 0:
            ans_vector[i] = h**2 * r[i] - alpha * (1 - (h * p[i]) / 2)
        elif i == len(ans_vector)-1:
            ans_vector[i] = h * beta
        else:
            ans_vector[i] = h**2 * r[i]
    print(df_eq)
    print(ans_vector)

    lu_facotr = linalg.lu_factor(df_eq)
    y = linalg.lu_solve(lu_facotr, ans_vector)
    y = np.insert(y, 0, alpha)

    x = np.linspace(0, 1, N+1)
    xx = np.linspace(0, 1)
    # plt.plot(xx, analysis(xx), color='r')
    plt.plot(x, y, color='b')
    title: str = 'TPBVP N = %d' % N
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

    error1: float = np.amax(abs(y - analysis(x)))
    error2: float = abs(y[-1] - analysis(1))
    print('error1 = %.3e' % error1)
    print('error2 = %.3e' % error2)
    print('Y_N ~ %.8e' % y[-1])


if __name__ == '__main__':
    main()
