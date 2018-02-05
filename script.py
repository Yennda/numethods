import math as m
import tools as t
import time
import scipy.integrate as it

# a = -3
# b = 5
# N = 1000
GAMMA = 16.65

"""
-rewrite b-a
- mayebe I could write rewrite the whole f-n without calling subf-ns
"""


class DVR():

    def __init__(self, a, b, N):
        self.a = a
        self.b = b
        self.N = N

    def phi(self, k, x):
        if k < 1:
            raise Exception('k = 1, 2, ...')
        return (2 / (self.b - self.a)) ** 0.5 * m.sin(
            k * m.pi * (x - self.a) / (self.b - self.a))

    def morse(self, x):
        return m.exp(-2 * x) - 2 * m.exp(-x)

    def potential_matrix_fn(self, k, l, x):
        return self.phi(k, x) * self.phi(l, x) * self.morse(x)

    def potential_matrix(self, k, l):
        return t.integrate_old(lambda x: self.potential_matrix_fn(k, l, x),
                               self.a, self.b, self.N)

        # return it.quad(lambda x: self.potential_matrix_fn(k, l, x), self.a,  # self.b)[0]

    def kinetic_matrix(self, k, l):
        if k != l:
            return 0
        return (1 / GAMMA * (k * m.pi / (self.b - self.a))) ** 2

    def hamilton_matrix(self, r):
        matrix = []
        for i in range(1, r):
            row = []
            for j in range(1, r):
                row.append(
                    self.potential_matrix(i, j) + self.kinetic_matrix(i, j))
            matrix.append(row)
        return matrix
