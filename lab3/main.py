from scipy.sparse import csr_matrix, dok_matrix, lil_matrix
import time
from creators import *
from iter_method import*
from lu_decomposition import *

A = csr_matrix(([1, 8, 7, 3, 7, 5, 4, 6, 5], [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 3, 6, 9]), shape=(3, 3))
print(A, "\n")


b = get_inverse_sparse(A)
Lu = get_LU_sparse(A)
array = get_LU_sparse(A)
L = array[0]
U = array[1]
print(L, "\n")
print(U, "\n")
print(b, "\n")


# A = csr_matrix(([10, -1, 2, -1, 11, -1, 3, 2, -1, 10, -1, 3, -1, 8], [0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3],[0, 3, 7, 11, 14]),shape=(4, 4))
# B = csr_matrix(([6, 25, -11, 15], [0, 0, 0, 0], [0, 1, 2, 3, 4]), shape=(4, 1))
A = csr_matrix(([10, 1, -1, 1, 10, -1, -1, 1, 10], [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 3, 6, 9]), shape=(3, 3))
B = csr_matrix(([11, 10, 10], [0, 0, 0], [0, 1, 2, 3]), shape=(3, 1))


# A = csr_matrix(([2, 3, 5, 7], [0, 1, 0, 1], [0, 2, 4]), shape=(2, 2))
# B = csr_matrix(([11, 13], [0, 0], [0, 1, 2]), shape=(2, 1))





# A = csr_matrix(createArray(20, 20))
# B = csr_matrix(createArray(20, 1))
# A = create_gilbert_matrix(4)
A = create_diagonal_preobl_matrix(2, 4)
# print(A)
# A = csr_matrix(([-3.999999, -1, -3, -3, -4.999999, -2, -1, -0.999999],[0,1,2,0,1,2,1,2],[0,3,6,8]),shape=(3,3))
B = csr_matrix(createArray(4, 1))
# B = csr_matrix(([4, 2, 1], [0, 0, 0], [0, 1, 2, 3]), shape=(3, 1))
EPS = 1e-5
start_time = time.time()
print(A, "\n", B)
C, d = zeidel(A, B, EPS)
print(C, "\n", d)
print("--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()
print(solve_system_sparse(A, B))
print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
    print('')

