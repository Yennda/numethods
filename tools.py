import math as m
import numpy as np
import copy


def integrate(f, a, b, n):
    if n % 2 != 0:
        raise Exception('n must be even')
    h = (b - a) / n
    return h * sum([f(a + i * h) for i in range(1, n)]) + h * (f(a) + f(b)) / 2


def integrate2(f, a, b, n):
    if n % 2 != 0:
        raise Exception('n must be even')
    h = (b - a) / n

    s = (f(a) + f(b)) / 2
    for i in range(1, n):
        s += f(a + i * h)

    return h * s


def is_strong_reg(H):
    for i in range(1, len(H)):
        nH = np.matrix(H)
    return not m.isclose(np.linalg.det(nH[:i, :i]), 0, rel_tol=0.001)


def one(n):
    o = []
    for i in range(n):
        r = n * [0]
        r[i] = 1
        o.append(r)
    return o


def zero(n):
    o = []
    for i in range(n):
        r = n * [0]
        o.append(r)
    return o


def dot(a, b):
    if len(a[0]) != len(b):
        raise Exception('not suitable matrix dimensions')
    C = []
    for i in range(len(a)):
        r = []
        for j in range(len(b)):
            r.append(sum([a[i][k] * b[k][j] for k in range(len(a[0]))]))
        C.append(r)

    return C


def lu(A):
    n = len(A)
    L = one(n)
    U = zero(n)
    U[0][:] = A[0][:]
    for i in range(n):
        L[i][0] = A[i][0] / U[0][0]

    def matrixl(i, j, A):
        return (A[i][j] - sum([L[i][k] * U[k][j] for k in range(j)])) / U[j][j]

    def matrixu(i, j, A):
        return A[i][j] - sum([L[i][k] * U[k][j] for k in range(i)])

    for j in range(1, n):
        U[j][j] = matrixu(j, j, A)
        for i in range(j + 1, n):
            L[i][j] = matrixl(i, j, A)
        for i in range(j, n):
            U[j][i] = matrixu(j, i, A)
    return L, U


def minus(a, b):
    if len(a) != len(b):
        raise Exception('different lengths of arrays')
    return [a[i] - b[i] for i in range(len(a))]


def eigenvalues(A, n):
    L, U = lu(A)

    for i in range(n):
        if i == n // 10 * 9:
            B = A
        L, U = lu(A)
        A = dot(U, L)
    return [A[i][i] for i in range(len(A))], [B[i][i] for i in range(len(B))]


def print_mat(self, M):
    for r in M:
        print(['{:02.2f}'.format(i) for i in r])

# def list_to_csv(data, path):
#     file=open(path, w)
#     for r in data:
#         file.write()
