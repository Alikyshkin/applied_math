import numpy as np
from scipy.sparse import lil_matrix
from lu_decomposition import *


# A - матрица хранящаяся в разреженном виде
def get_inverse_sparse(A):
    n = A.shape[0]

    L, U = get_LU_sparse(A)

    Y = []

    for i in range(n):
        y = solve_triangular_system_sparse(
            L, np.eye(100, 1, k=-i), lower=True)
        Y.append(y)

    AI = lil_matrix(A.shape)
    for i in range(n):
        AI[:, i] = solve_triangular_system_sparse(
            U, Y[i], lower=False)

    return AI
