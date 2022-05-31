import numpy as numpy
from scipy.sparse import csr_matrix


def jacobi3(A, tol=1.0e-9):
    a = A.copy()

    def maxElem(a):
        '''
            Find largest off-diag. element a[k,l]
        '''
        n = len(a)
        # n = a.shape[0]
        aMax = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                if abs(a[i, j]) >= aMax:
                    aMax = abs(a[i, j])
                    k = i
                    l = j
        return aMax, k, l

    def rotate(a, p, k, l):
        '''
            Rotate to make a[k,l] = 0
        '''
        n = len(a)
        # n = a.shape[0]
        aDiff = a[l, l] - a[k, k]

        if abs(a[k, l]) < abs(aDiff) * 1.0e-36:
            t = a[k, l] / aDiff
        else:
            phi = aDiff / (2.0 * a[k, l])
            t = 1.0 / (abs(phi) + numpy.sqrt(phi ** 2 + 1.0))
            if phi < 0.0:
                t = -t

        c = 1.0 / numpy.sqrt(t ** 2 + 1.0)
        s = t * c
        tau = s / (1.0 + c)
        temp = a[k, l]
        a[k, l] = 0.0
        a[k, k] = a[k, k] - t * temp
        a[l, l] = a[l, l] + t * temp

        for i in range(k):  # Case of i < k
            temp = a[i, k]
            a[i, k] = temp - s * (a[i, l] + tau * temp)
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])

        for i in range(k + 1, l):  # Case of k < i < l
            temp = a[k, i]
            a[k, i] = temp - s * (a[i, l] + tau * a[k, i])
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])

        for i in range(l + 1, n):  # Case of i > l
            temp = a[k, i]
            a[k, i] = temp - s * (a[l, i] + tau * temp)
            a[l, i] = a[l, i] + s * (temp - tau * a[l, i])

        for i in range(n):  # Update transformation matrix
            temp = p[i, k]
            p[i, k] = temp - s * (p[i, l] + tau * p[i, k])
            p[i, l] = p[i, l] + s * (temp - tau * p[i, l])

    n = len(a)
    # n = a.shape[0]
    maxRot = 5 * (n ** 2)  # Set rotation number limit
    p = numpy.identity(n) * 1.0  # Initialize transformation matrix

    for i in range(maxRot):  # Jacobi rotation loop
        aMax, k, l = maxElem(a)

        if aMax < tol:
            print("Number of rotations: " + str(i + 1))
            return numpy.diagonal(numpy.flip(a)), p
            # return csr_matrix.diagonal(a), p

        rotate(a, p, k, l)

    print('Jacobi method did not converge')