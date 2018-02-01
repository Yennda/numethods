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


def is_small(A):
    for a in A:
        if m.fabs(a) > 1e-10:
            return False
    else:
        return True


def mmax(A):
    B = copy.copy(A)
    lA = []
    n = len(B)

    for i in range(n):
        lA.extend(B[i])
        lA.remove(B[i][i])
    index = lA.index(max(lA))
    i = index // (n - 1)
    j = index % (n - 1)
    if j >= i:
        j += 1
    return i, j


def S_matrix(n, p, q, a):
    S = one(n)
    S[p][p] = m.cos(a)
    S[q][q] = m.cos(a)
    S[p][q] = m.sin(a)
    S[q][p] = -m.sin(a)
    return S


def trans(A):
    AT = []
    for i in range(len(A)):
        r = [A[j][i] for j in range(len(A))]
        AT.append(r)
    return AT


def eigenLU(A, n):
    L, U = lu(A)
    for i in range(n):
        lA = A
        L, U = lu(lA)
        A = dot(U, L)
        if is_small([A[i][i] - lA[i][i] for i in range(len(A))]):
            return [A[i][i] for i in range(len(A))]
    return [A[i][i] for i in range(len(A))]


def eigenJacobi(A, n):
    d = len(A)
    for i in range(n):
        lA = A
        p, q = mmax(A)
        if (A[q][q] - A[p][p]) == 0:
            a = m.pi / 2
        else:
            a = 0.5 * m.atan(2 * A[p][q] / (A[q][q] - A[p][p]))
        S = S_matrix(d, p, q, a)
        A = dot(trans(S), dot(lA, S))

        if is_small([A[i][i] - lA[i][i] for i in range(len(A))]):
            # print(i)
            return [A[i][i] for i in range(len(A))]
    return [A[i][i] for i in range(len(A))]


def print_mat(M):
    for r in M:
        print(['{:02.2f}'.format(i) for i in r])
