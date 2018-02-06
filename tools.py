import math as m
import numpy as np
import copy


def integrate(f, a, b, n):
    h = (b - a)
    Io = (f(a) + f(b)) / 2 * h

    h /= 2
    In = f(a + h)
    I = Io / 2 + In * h
    k = 1

    while k <= n and not (m.fabs(I - Io) < 1e-6 and k > 4):
        Io = I
        h /= 2
        In = sum([f(a + h + i * 2 * h) for i in range(0, 2 ** k)])
        I = Io / 2 + In * h
        k += 1
    return I


def integrate_old(f, a, b, n):
    if n % 2 != 0:
        raise Exception('n must be even')
    h = (b - a) / n
    return h * sum([f(a + i * h) for i in range(1, n)]) + h * (f(a) + f(b)) / 2


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
        if m.fabs(a) > 1e-19:
            return False
    else:
        return True


def matrix_max(A):
    B = copy.copy(A)
    lA = []
    n = len(B)

    for i in range(n):
        lA.extend(B[i])
        lA.remove(B[i][i])
    lA = [m.fabs(l) for l in lA]
    index = lA.index(max(lA))
    i = index // (n - 1)
    j = index % (n - 1)
    if j >= i:
        j += 1
    return i, j


def S_matrix(n, p, q, a):
    S = one(n)
    St = one(n)

    S[p][p] = m.cos(a)
    S[q][q] = m.cos(a)
    S[p][q] = m.sin(a)
    S[q][p] = -m.sin(a)

    St[p][p] = m.cos(a)
    St[q][q] = m.cos(a)
    St[p][q] = -m.sin(a)
    St[q][p] = m.sin(a)
    return S, St


def trans(A):
    AT = []
    for i in range(len(A)):
        r = [A[j][i] for j in range(len(A))]
        AT.append(r)
    return AT


def eigenLU(A, n):
    for i in range(n):
        lA = A[-1][-1]
        L, U = lu(A)
        A = dot(U, L)

        # if m.fabs(A[-1][-1] - lA) < 1e-19:
        #     print('sooner {}'.format(i))
        #     return [A[i][i] for i in range(len(A))]
    return [A[i][i] for i in range(len(A))]


def eigenJacobi(A, n):
    d = len(A)
    for i in range(n):
        if i == n // 10 * 9:
            B = A
        lA = A
        p, q = matrix_max(A)
        if (A[q][q] - A[p][p]) == 0:
            a = m.pi / 2
        else:
            a = 0.5 * m.atan(2 * A[p][q] / (A[q][q] - A[p][p]))
        S = S_matrix(d, p, q, a)
        A = dot(S[1], dot(lA, S[0]))

        # if is_small([A[i][i] - lA[i][i] for i in range(len(A))]):
        #     print('sooner {}'.format(i))
        #     return [A[i][i] for i in range(len(A))]

    return [A[i][i] for i in range(len(A))]


def print_mat(M):
    for r in M:
        print(['{:02.4f}'.format(i) for i in r])
