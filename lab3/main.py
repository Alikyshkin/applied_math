from scipy.sparse import csr_matrix, dok_matrix, lil_matrix
import time
from creators import *
from iter_method import *
from lu_decomposition import *

'''
# Задание 1
# Создаем матрицу в разреженно-столбцовом формате
A = csr_matrix(([1, 8, 7, 3, 7, 5, 4, 6, 5], [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 3, 6, 9]), shape=(3, 3))
print("Matrix A")
print(A, "\n")

b = get_inverse_sparse(A)
print("Inversed matrix A")
print(b, "\n")

# Lu = get_LU_sparse(A)
array = get_LU_sparse(A)

# Нижняя
L = array[0]
print("L sparse of matrix A")
print(L, "\n")

# Верхняя
U = array[1]
print("U sparse of matrix A")
print(U, "\n")
'''

'''
# Задание 2
print('_______________________________________________________________________________________________________________')
print('Метод Зейделя: ')

A = create_diagonal_preobl_matrix(2, 4)
B = csr_matrix(createArray(4, 1))
EPS = 1e-5

print('Мaтрица А: ')
print(A)
print('Мaтрица В: ')
print(B)


C, count = zeidel(A, B, EPS)
print('Мaтрица C: ')
print(C)
print('count: ')
print(count)
'''

# '''
# Задание 3
print('_______________________________________________________________________________________________________________')
print('Матрица:  create_diagonal_preobl_matrix(1, 20)')
print('LU: ')
# A = csr_matrix(createArray(20, 20))
# B = csr_matrix(createArray(20, 1))
A = create_diagonal_preobl_matrix(1, 20)
B = create_diagonal_preobl_matrix(1, 20)
start_time = time.time()
C = (solve_system_sparse(A, B))
print("--- %s seconds ---" % (time.time() - start_time))

print('Метод Зейделя: ')
EPS = 1e-5
# A = create_diagonal_preobl_matrix(13, 20)
# B = create_diagonal_preobl_matrix(13, 20)
start_time1 = time.time()
C, count = zeidel(A, B, EPS)
print("--- %s seconds ---" % (time.time() - start_time1))

# '''

'''
# Задание 4
print('_______________________________________________________________________________________________________________')
print('LU: ')
A = create_gilbert_matrix(10)
B = create_gilbert_matrix(10)
start_time = time.time()
C = (solve_system_sparse(A, B))

print("--- %s seconds ---" % (time.time() - start_time))

print('Метод Зейделя: ')
EPS = 1e-5
A = create_gilbert_matrix(10)
B = create_gilbert_matrix(10)

start_time = time.time()
C, count = zeidel(A, B, EPS)
print("--- %s seconds ---" % (time.time() - start_time))
print('count: ')
print(count)

'''
'''

# Задание 5
print('_______________________________________________________________________________________________________________')
print('LU: ')
A = csr_matrix(createArray(100, 100))
B = csr_matrix(createArray(100, 100))
start_time = time.time()
C, count = (solve_system_sparse_count(A, B))
print('count: ')
print(count)
print("--- %s seconds ---" % (time.time() - start_time))

print('Метод Зейделя: ')
EPS = 1e-5
A = csr_matrix(createArray(10, 10))
B = csr_matrix(createArray(10, 10))

start_time = time.time()
C, count = zeidel(A, B, EPS)
print("--- %s seconds ---" % (time.time() - start_time))
print('count: ')
print(count)
'''

# A = csr_matrix(([10, -1, 2, -1, 11, -1, 3, 2, -1, 10, -1, 3, -1, 8], [0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3],[0, 3, 7, 11, 14]),shape=(4, 4))
# B = csr_matrix(([6, 25, -11, 15], [0, 0, 0, 0], [0, 1, 2, 3, 4]), shape=(4, 1))
# A = csr_matrix(([10, 1, -1, 1, 10, -1, -1, 1, 10], [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 3, 6, 9]), shape=(3, 3))
# B = csr_matrix(([11, 10, 10], [0, 0, 0], [0, 1, 2, 3]), shape=(3, 1))

# A = csr_matrix(([2, 3, 5, 7], [0, 1, 0, 1], [0, 2, 4]), shape=(2, 2))
# B = csr_matrix(([11, 13], [0, 0], [0, 1, 2]), shape=(2, 1))

if __name__ == '__main__':
    print('')
