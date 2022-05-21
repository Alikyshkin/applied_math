from scipy.sparse import dok_matrix


# Найти LU-разложение матрицы A, хранящейся в разреженном виде
def get_LU_sparse(A):
    EPS = 1e-5

    L = dok_matrix(A.shape)
    U = dok_matrix(A.shape)

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if i <= j:
                t = A[i, j] - (L[i, :i] @ U[:i, j]).sum()

                if abs(t) > EPS:
                    U[i, j] = t

                if i == j:
                    L[i, i] = 1.0

            else:
                if U[j, j] == 0.0:
                    raise ValueError('LU decomposition does not exist!')

                t = (A[i, j] - (L[i, :j] @ U[:j, j]).sum()) / U[j, j]
                if abs(t) > EPS:
                    L[i, j] = t

    return L.tocsr(), U.tocsr()


def get_LU_sparse_count(A, count):
    EPS = 1e-100

    L = dok_matrix(A.shape)
    U = dok_matrix(A.shape)

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if i <= j:
                t = A[i, j] - (L[i, :i] @ U[:i, j]).sum()
                count = count + 1

                if abs(t) > EPS:
                    U[i, j] = t

                if i == j:
                    L[i, i] = 1.0

            else:
                if U[j, j] == 0.0:
                    raise ValueError('LU decomposition does not exist!')

                t = (A[i, j] - (L[i, :j] @ U[:j, j]).sum()) / U[j, j]
                count = count + 1
                if abs(t) > EPS:
                    L[i, j] = t

    return L.tocsr(), U.tocsr(), count


# Найти тривиальное решение уравнения Ax=B,
# где A - нижнетреугольная матрица, хранящаяся в разреженном виде
# и B - вектор в правой части
# или где A - верхнетреугольная матрица, хранящаяся в разреженном виде
# и B - вектор в правой части
def solve_triangular_system_sparse(A, b, lower):
    n = A.shape[0]

    x = dok_matrix((n, 1))

    if lower:
        for i in range(n):
            x[i] = (b[i].sum() - (A[i, :i] @ x[:i]).sum()) / A[i, i]
    else:
        for i in range(n - 1, -1, -1):
            x[i] = (b[i].sum() - (A[i, i + 1:] @ x[i + 1:]).sum()) / A[i, i]
    return x


def solve_triangular_system_sparse_count(A, b, count, lower):
    n = A.shape[0]

    x = dok_matrix((n, 1))

    if lower:
        for i in range(n):
            x[i] = (b[i].sum() - (A[i, :i] @ x[:i]).sum()) / A[i, i]
            count = count + 1
    else:
        for i in range(n - 1, -1, -1):
            x[i] = (b[i].sum() - (A[i, i + 1:] @ x[i + 1:]).sum()) / A[i, i]
            count = count + 1
    return x, count


# Найти решение системы Ax=B,
# где A - невырожденная матрица, хранящаяся в разреженном виде
# и B - вектор в правой части
def solve_system_sparse(A, b):
    n = A.shape[0]

    L, U = get_LU_sparse(A)

    y = solve_triangular_system_sparse(L, b, lower=True)
    x = solve_triangular_system_sparse(U, y, lower=False)

    return x


def solve_system_sparse_count(A, b):
    n = A.shape[0]
    count = 0
    L, U, count = get_LU_sparse_count(A, count)

    y, count = solve_triangular_system_sparse_count(L, b, count, lower=True)
    x, count = solve_triangular_system_sparse_count(U, y, count, lower=False)

    return x, count
