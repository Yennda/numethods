import timeit

print(timeit.timeit('t.integrate(lambda x: x ** 2, 0, 1, 100)', setup='import tools as t', number=100))
print(timeit.timeit('t.integrate2(lambda x: x ** 2, 0, 1, 100)', setup='import tools as t', number=100))
print('\n')
print(timeit.timeit('t.integrate(lambda x: x ** 2, 0, 1, 100)', setup='import tools as t', number=10000))
print(timeit.timeit('t.integrate2(lambda x: x ** 2, 0, 1, 100)', setup='import tools as t', number=10000))
