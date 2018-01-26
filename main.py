import math as m
import tools as t
import time
import scipy

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


def potential_matrix_fn(k, l, x):
    return phi(k, x) * phi(l, x) * (m.exp(-2 * x) - 2 * m.exp(-x))


def potential_matrix(k, l):
    return t.integrate(lambda x: potential_matrix_fn(k, l, x), a, b, N)


def matrix_hamilton(r, c):
    matrix = []
    for i in range(1, r):
        row = []
        for j in range(1, c):
            row.append(kinetic_matrix(i, j) + potential_matrix(i, j))
        matrix.append(row)
    return matrix


def print_mat(M):
    for r in M:
        print(r)


def kinetic_matrix(k, l):
    if k != l:
        return 0
    return GAMMA ** -2 * (k * m.pi / (b - a)) ** 2


dt = time.time()
print('Tschus Welt!')

H = matrix_hamilton(10, 10)
print_mat(H)

lt = time.time()
print_mat(t.ludec(H)[0])
print_mat(t.ludec(H)[1])

print('\nElapsed time: {}'.format(time.time() - lt))

print('\nElapsed time: {}'.format(time.time() - dt))
