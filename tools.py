import math as m
import numpy as np


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


def ludec(A):
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

