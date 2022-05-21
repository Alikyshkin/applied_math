from inverse import get_inverse_sparse
from scipy.sparse import dok_matrix


# Итерационный метод Зейделя
def zeidel(A, B, EPS):
    U = dok_matrix(A.shape)
    L = dok_matrix(A.shape)
    for i in range(A.shape[0]):
        for j in range(A.shape[0]):
            if i < j:
                U[i, j] = A[i, j]
            else:
                L[i, j] = A[i, j]
    Linv = get_inverse_sparse(L)
    T = -Linv @ U
    C = Linv @ B
    x = dok_matrix(B.shape)
    prev = dok_matrix(B.shape)
    for i in range(B.shape[0]):
        x[i, 0] = 1
    count = 0
    while isNotEqual(x, prev, EPS):
        prev = x
        x = T @ x + C
        count += 1
    return x, count


# Сравнение
def isNotEqual(A, B, EPS):
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if abs(A[i, j] - B[i, j]) > EPS:
                return True
    return False
