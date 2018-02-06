import math as m
import tools as t

GAMMA = 16.65


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

    def potential_matrix_fn_old(self, k, l, x):
        return self.phi(k, x) * self.phi(l, x) * self.morse(x)

    def potential_matrix_fn(self, k, l, x):
        return self.phi(k, x) * self.phi(l, x) * self.morse(x)

    def potential_matrix(self, k, l):
        return t.integrate_trapez(lambda x: self.potential_matrix_fn(k, l, x),
                                  self.a, self.b, self.N)

    def kinetic_matrix(self, k, l):
        if k != l:
            return 0
        return (1 / GAMMA * (k * m.pi / (self.b - self.a))) ** 2

    def hamilton_matrix(self, r):
        matrix = t.zero(r - 1)
        for i in range(r - 1):
            for j in range(i + 1):
                element = self.potential_matrix(i + 1,
                                                j + 1) + self.kinetic_matrix(
                    i + 1, j + 1)
                matrix[i][j] = self.potential_matrix(i + 1,
                                                     j + 1) + self.kinetic_matrix(
                    i + 1, j + 1)
                matrix[j][i] = matrix[i][j]
        return matrix
