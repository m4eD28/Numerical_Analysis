'''TPBVP'''
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as linalg


# 解析解
def analysis(x) -> np.array:
    return -2 / (1 - np.exp(-4)) * np.exp(-4 * x) - x + 2 / (1 - np.exp(-4))


# main
def main() -> None:

    N: int = 24
    h: float = 1 / N

    alpha: float = 0
    beta: float = 1
    p: np.array = np.full(N, 4, dtype='float64')
    q: np.array = np.full(N, 0, dtype='float64')
    r: np.array = np.full(N, -4, dtype='float64')

    df_eq: np.array = np.zeros((N-1, N-1), dtype='float64')
    for i in range(len(df_eq)):
        for j in range(len(df_eq)):
            if i == j:
                df_eq[i][j] = h**2 * q[i] - 2
            elif i + 1 == j:
                df_eq[i][j] = 1 + (h * p[i]) / 2
            elif i - 1 == j:
                df_eq[i][j] = 1 - (h * p[i]) / 2

    ans_vector: np.array = np.zeros(N-1, dtype='float64')
    for i in range(len(ans_vector)):
        if i == 0:
            ans_vector[i] = h**2 * r[i] - alpha * (1 - (h * p[i]) / 2)
        elif i == len(ans_vector)-1:
            print('hello')
            ans_vector[i] = h**2 * r[i] - beta * (1 + (h * p[i]) / 2)
        else:
            ans_vector[i] = h**2 * r[i]

    print(df_eq)
    print(ans_vector)

    lu_facotr: np.array = linalg.lu_factor(df_eq)
    y: np.array = linalg.lu_solve(lu_facotr, ans_vector)
    y = np.insert(y, 0, alpha)
    y = np.append(y, beta)

    x = np.linspace(0, 1, N+1)
    xx = np.linspace(0, 1)
    plt.plot(xx, analysis(xx), color='r', label='Y')
    plt.plot(x, y, color='b', label='y(x)')
    title: str = 'TPBVP N = %d' % N
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(loc='upper left')
    plt.show()

    error: float = np.amax(abs(y - analysis(x)))
    print('error = %.3e' % error)
    print('h = %.3e' % h)
    print('Y_1 ~ %.5e' % y[1])


if __name__ == '__main__':
    main()
