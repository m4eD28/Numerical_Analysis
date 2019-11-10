from typing import List
import numpy as np


def Forward_erase(A: np.array, B: np.array) -> np.array:
    alpha: float
    for k in range(len(A)-1):
        for i in range(k+1, len(A)):
            alpha = A[i][k] / A[k][k]
            for j in range(k+1, len(A)):
                A[i][j] = A[i][j] - alpha * A[k][j]
            B[i] = B[i] - alpha * B[k]
    return A


def Backward_sub(A: np.array, b: np.array) -> np.array:
    for k in range(len(a)-1, 0):
        for j in range(k+1, )
