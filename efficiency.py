import timeit

print(timeit.timeit('t.integrate(lambda x: x ** 2, 0, 1, 100)', setup='import tools as t', number=100))
print(timeit.timeit('t.integrate2(lambda x: x ** 2, 0, 1, 100)', setup='import tools as t', number=100))
print('\n')
print(timeit.timeit('t.integrate(lambda x: x ** 2, 0, 1, 100)', setup='import tools as t', number=10000))
print(timeit.timeit('t.integrate2(lambda x: x ** 2, 0, 1, 100)', setup='import tools as t', number=10000))

print('LU x Jacobi')
A = [[2, 0, 0], [0, 3, 4], [0, 4, 9]]
B = [[4, -5, 7], [1, -4, 9], [-4, 0, 5]]
# print(timeit.timeit('t.eigenLU([[4, -1, 3], [8, -2, 6], [-4, 1, -3]], 500)', setup='import tools as t', number=10000))
print(
    timeit.timeit('t.eigenJacobi([[4, -1, 3], [8, -2, 6], [-4, 1, -3]], 500)', setup='import tools as t', number=10000))
