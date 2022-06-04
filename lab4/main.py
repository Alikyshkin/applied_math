import time

from scipy.sparse import csr_matrix, dok_matrix
from numpy import identity, flip, round
from math import sqrt
from random import uniform


# Создать матрицу с диагональным преобладанием
def create_diagonal_matrix(k, n):
    A = dok_matrix((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i, j] = uniform(0.0, 4.0) - 4.0
    for i in range(n):
        t1 = -sum(A[i, k] for k in range(i))
        t2 = -sum(A[i, k] for k in range(i + 1, n))
        t = t1 + t2
        A[i, i] = t + pow(10.0, -k)
    return A


# Создание матрицы Гилберта
def create_gilbert_matrix(n):
    M = dok_matrix((n, n))
    for i in range(n):
        for j in range(n):
            M[i, j] = 1 / (i + j + 1)
    return M


# С диагональным ли преобладанием матрица
def isDDM(m, n):
    # for each row
    sum = 0
    sum_diags = 0
    for i in range(0, n):
        for j in range(0, n):
            sum = sum + abs(m[i, j])

        sum = sum - abs(m[i, i])
        sum_diags += abs(m[i, i])

    if (sum_diags < sum):
        return False

    return True


# Метод Якоби (только для csr матриц)
def jacobi_for_csr_matrix(a, tol=1.0e-9):
    def maxElem1(a):
        n = a.shape[0]
        aMax = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                if abs(a[i, j]) >= aMax:
                    aMax = abs(a[i, j])
                    k = i
                    l = j
        return aMax, k, l

    def rotate(a, p, k, l):
        n = a.shape[0]

        aDiff = a[l, l] - a[k, k]
        if abs(a[k, l]) < abs(aDiff) * 1.0e-36:
            t = a[k, l] / aDiff
        else:
            phi = aDiff / (2.0 * a[k, l])
            t = 1.0 / (abs(phi) + sqrt(phi ** 2 + 1.0))
            if phi < 0.0: t = -t
        c = 1.0 / sqrt(t ** 2 + 1.0)
        s = t * c
        tau = s / (1.0 + c)
        temp = a[k, l]
        a[k, l] = 0.0
        a[k, k] = a[k, k] - t * temp
        a[l, l] = a[l, l] + t * temp
        for i in range(k):  # если i < k
            temp = a[i, k]
            a[i, k] = temp - s * (a[i, l] + tau * temp)
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])
        for i in range(k + 1, l):  # если k < i < l
            temp = a[k, i]
            a[k, i] = temp - s * (a[i, l] + tau * a[k, i])
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])
        for i in range(l + 1, n):  # если i > l
            temp = a[k, i]
            a[k, i] = temp - s * (a[l, i] + tau * temp)
            a[l, i] = a[l, i] + s * (temp - tau * a[l, i])
        for i in range(n):  # обновляем матрицу преобразований
            temp = p[i, k]
            p[i, k] = temp - s * (p[i, l] + tau * p[i, k])
            p[i, l] = p[i, l] + s * (temp - tau * p[i, l])

    n = a.shape[0]
    maxRot = 5 * (n ** 2)  # максимальное число вращений
    p = identity(n) * 1.0  # матрица преобразования
    for i in range(maxRot):  # цикл вращений метода Якоби
        aMax, k, l = maxElem1(a)
        if aMax < tol: return csr_matrix.diagonal(flip(a)), round(p), i
        rotate(a, p, k, l)
    print('Jacobi method did not converge')


# Метод Якоби (для матриц Гильберта и с диагональным преобладанием)
def jacobi_for_dok_matrix(a, tol=1.0e-9):
    def maxElem1(a):
        n = a.shape[0]

        aMax = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                if abs(a[i, j]) >= aMax:
                    aMax = abs(a[i, j])
                    k = i
                    l = j
        return aMax, k, l

    def rotate(a, p, k, l):
        n = a.shape[0]

        aDiff = a[l, l] - a[k, k]
        if abs(a[k, l]) < abs(aDiff) * 1.0e-36:
            t = a[k, l] / aDiff
        else:
            phi = aDiff / (2.0 * a[k, l])
            t = 1.0 / (abs(phi) + sqrt(phi ** 2 + 1.0))
            if phi < 0.0: t = -t
        c = 1.0 / sqrt(t ** 2 + 1.0)
        s = t * c
        tau = s / (1.0 + c)
        temp = a[k, l]
        a[k, l] = 0.0
        a[k, k] = a[k, k] - t * temp
        a[l, l] = a[l, l] + t * temp
        for i in range(k):  # если i < k
            temp = a[i, k]
            a[i, k] = temp - s * (a[i, l] + tau * temp)
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])
        for i in range(k + 1, l):  # если k < i < l
            temp = a[k, i]
            a[k, i] = temp - s * (a[i, l] + tau * a[k, i])
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])
        for i in range(l + 1, n):  # если i > l
            temp = a[k, i]
            a[k, i] = temp - s * (a[l, i] + tau * temp)
            a[l, i] = a[l, i] + s * (temp - tau * a[l, i])
        for i in range(n):  # обновляем матрицу преобразований
            temp = p[i, k]
            p[i, k] = temp - s * (p[i, l] + tau * p[i, k])
            p[i, l] = p[i, l] + s * (temp - tau * p[i, l])

    n = a.shape[0]
    maxRot = 5 * (n ** 2)  # максимальное число вращений
    p = identity(n) * 1.0  # матрица преобразования
    for i in range(maxRot):  # цикл вращений метода Якоби
        aMax, k, l = maxElem1(a)
        if aMax < tol: return dok_matrix.diagonal(flip(a)), round(p), i
        rotate(a, p, k, l)
    print('Jacobi method did not converge')


# Задание 1
print("Задание 1")
# Создаем квадратную матрицу в разреженно-столбцовом формате
A = csr_matrix(([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 3, 6, 9]), shape=(3, 3))
print("Matrix A")
print(A, "\n")

start_time = time.time()
lam, x, i = jacobi_for_csr_matrix(A)
print("n = ", 3, "; вращения: ", str(i + 1))
print("--- %s seconds ---" % (time.time() - start_time))
print("eigenvalues in vector:")
print(lam, "\n")
print("eigenvectors as columns of matrix:")
print(x, "\n")

# Задание 2
print("Задание 2")
# Создаем матрицу с диагональным преобладанием
n = 10  # размерность матрицы
A = create_diagonal_matrix(1, n)
# print("Matrix A")
# print(A, "\n")

# проверяем что сгенерировалась правильная матрица
if isDDM(A, n):
    print("YES, it is a DDM")
else:
    print("NO, it is not a DDM")

for count in range(0, 20):
    A = create_diagonal_matrix(count, n)
    start_time = time.time()
    lam, x, i = jacobi_for_dok_matrix(A)
    print("k = ", count, "вращения: ", str(i + 1))
    print("--- %s seconds ---" % (time.time() - start_time))

# Задание 3
print("Задание 3")
# Создаем матрицу Гильберта
for count in range(2, 26):
    A = create_gilbert_matrix(count)
    start_time = time.time()
    lam, x, i = jacobi_for_dok_matrix(A)
    print("n = ", count, "вращения: ", str(i + 1))
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    print('')
