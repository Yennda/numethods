import math as m
import tools as t
import time
from script import DVR


def compute(n):
    lt = time.time()

    H = dvr.hamilton_matrix(n)
    print('\nT {:.3f} s'.format(time.time() - lt))

    print('Matrix of hamiltonian H:')
    if n > 10:
        print('too large')
    else:
        t.print_mat(H)

    res = t.eigenLU(H, 500)
    print('Eigenvalues:\n {}'.format(res))

    print('\nT {:.3f} s'.format(time.time() - lt))


def compute_simple(n):
    H = dvr.hamilton_matrix(n)
    return sorted(t.eigenLU(H, 1500))[0:3]


def err_integral(dim):
    lis = []
    times = []
    for i in range(0, 1000, 100):
        dvr.N += i

        rt = time.time()
        lis.append(compute_simple(dim))
        times.append(time.time() - rt)

        print('Loop {}: {}'.format(i, time.time() - rt))

    print(
        'value_ig[0].extend({})\nvalue_ig[1].extend({})\nvalue_ig[2].extend({})\n'.format(
            [lis[i][0] for i in range(len(lis))],
            [lis[i][1] for i in range(len(lis))],
            [lis[i][2] for i in range(len(lis))]))

    print('times_ig.extend({})'.format(times))


def err_dimension():
    lis = []
    times = []
    for i in range(30, 50):
        rt = time.time()
        lis.append(compute_simple(i))
        times.append(time.time() - rt)

        print('Loop {}: {}'.format(i, time.time() - rt))

    print(
        'value[0].extend({})\nvalue[1].extend({})\nvalue[2].extend({})\n'.format(
            [lis[i][0] for i in range(len(lis))],
            [lis[i][1] for i in range(len(lis))],
            [lis[i][2] for i in range(len(lis))]))

    print('times.extend({})'.format(times))


def err_interval(dim):
    lis = []
    times = []
    for i in range(0, 10):
        dvr.a -= i / 10
        dvr.b += i / 10

        rt = time.time()
        lis.append(compute_simple(dim))
        times.append(time.time() - rt)

        print('Loop {}: {}'.format(i, time.time() - rt))

    print(
        'value_iv[0].extend({})\nvalue_iv[1].extend({})\nvalue_iv[2].extend({})\n'.format(
            [lis[i][0] for i in range(len(lis))],
            [lis[i][1] for i in range(len(lis))],
            [lis[i][2] for i in range(len(lis))]))

    print('times_iv.extend({})'.format(times))


dt = time.time()
rt = time.time()

dvr = DVR(a=-3, b=5, N=100)
dim = 70

# err_dimension()
# err_interval(dim)
err_integral(dim)

print('\nTotal time: {:.3f} s'.format(time.time() - dt))
