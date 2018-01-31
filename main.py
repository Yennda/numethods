import math as m
import tools as t
import time
import numpy as np

a = -3
b = 5
N = 1000
GAMMA = 16.65

"""
-rewrite b-a
- mayebe I could write rewrite the whole f-n without calling subf-ns
"""


def phi(k, x):
    if k < 1:
        raise Exception('k = 1, 2, ...')
    return (2 / (b - a)) ** 0.5 * m.sin(k * m.pi * (x - a) / (b - a))


def morse(x):
    return m.exp(-2 * x) - 2 * m.exp(-x)


def potential_matrix_fn(k, l, x):
    return phi(k, x) * phi(l, x) * morse(x)


def potential_matrix(k, l):
    return t.integrate(lambda x: potential_matrix_fn(k, l, x), a, b, N)


def kinetic_matrix(k, l):
    if k != l:
        return 0
    return (1 / GAMMA * (k * m.pi / (b - a))) ** 2


def hamilton_matrix(r):
    matrix = []
    for i in range(1, r):
        row = []
        for j in range(1, r):
            row.append(potential_matrix(i, j) + kinetic_matrix(i, j))
        matrix.append(row)
    return matrix


def print_mat(M):
    for r in M:
        print(['{:02.2f}'.format(i) for i in r])


def compute(n):
    lt = time.time()

    H = hamilton_matrix(n)
    print('\nT {:.3f} s'.format(time.time() - lt))

    print('Matrix of hamiltonian H:')
    if n > 10:
        print('too large')
    else:
        print_mat(H)

    print('H is strongly regular: {}'.format(t.is_strong_reg(H)))

    res = t.eigenvalues(H, 500)
    print(
        'Eigenvalues:\n {}\nEigenvalues at 0.9n:\n {}'.format(res[0], res[1]))

    print('\nT {:.3f} s'.format(time.time() - lt))


dt = time.time()
# compute(50)

print('\nTotal time: {:.3f} s'.format(time.time() - dt))
