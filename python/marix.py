import numpy as np
from scipy.sparse import dia_matrix identity

# construct matrix A
data = np.array([np.ones(n), -2.0*np.ones(n), np.ones(n)])
offsets = np.array([-1, 0, 1])
B = dia_matrix((data, offsets), shape=(n, n))
E = identity(n)
A = E + r * B
