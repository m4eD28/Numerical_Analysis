import math
from matplotlib import pyplot as plt

def f(x):
    return (x-1)**2 * (x+1)

def df(x):
    return 2*(x-1)*(x+1) + (x-1)**2

if __name__ == '__main__':
    print('Newton')
    x0 = 7; xs = 1; eps = 1e-7; Max = 60
    x = x0; N = 0; x_N = 0;
    Error = [];

    print('step \t x_k \t error \t order')
    for k in range(Max+1):
        if(abs(f(x))<eps):
            N = k+1; x_N = x;
            print('N={}, x_N = {:.12e}'.format(N,x_N))
            break
        x = x0 - f(x0)/df(x0)
        print('{} \t {:.12e} \t {:.2e} \t {:.2e}'.format(k+1,x,x-2,xs))
        x0 = x;
        Error.append(abs(x-xs));

    plt.plot(range(1, N), Error, marker='o')
    plt.yscale('log')
    plt.show()
