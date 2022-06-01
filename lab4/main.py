import numpy as numpy
from scipy.sparse import csr_matrix, dok_matrix, lil_matrix

import time
from jacobi1_method import *

# '''
# Задание 1
print("Задание 1")
# Создаем квадратную матрицу в разреженно-столбцовом формате
# A = np.array([[1.0,8.0,7.0],[3.0,7.0,5.0],[4.0,6.0,5.0]])
A = csr_matrix(([1.0, 8.0, 7.0, 3.0, 7.0, 5.0, 4.0, 6.0, 5.0], [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 3, 6, 9]), shape=(3, 3))
# A = np.array([[1.0,1.0,1.0],[1.0,1.0,1.0],[1.0,1.0,1.0]])
# A = csr_matrix(([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 3, 6, 9]), shape=(3, 3))
# A = create_diagonal_preobl_matrix(1, 3)  # k = 1, n = 3
# A = create_gilbert_matrix(3)
print("Matrix A\n", A)

print("\nОтвет который должен получиться: ")
# print("Полином: -Л^3 + 13Л^2 + 35Л - 25")
print("Cобств числа: 15.1951, -2.78571, 0.59061")
print("Cобств вектора: \n[1.05454, 0.996158, 1] \n[-2.1846, 0.158783, 1] \n[0.227625, -0.886648, 1]")
# print("Cобств числа: 3, 0, 0")
# print("Cобств вектора: \n[1, 1, 1] \n[-1, 0, 1] \n[-1, 1, 0]")
print("\nОтвет программы: ")

print("JACOBI 1")
start_time = time.time()
lam, x, i = jacobi1(A)
print("вращения: ", str(i + 1))
print("--- %s seconds ---" % (time.time() - start_time))
print("eigenvalues in vector:")
print(lam, "\n")
print("eigenvectors as columns of matrix:")
print(x, "\n")

if __name__ == '__main__':
    print('')
