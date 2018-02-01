import math as m
import tools as t
import time
from script import DVR

dvr = DVR(a=-3, b=5, N=1000)


def compute(n):
    lt = time.time()

    H = dvr.hamilton_matrix(n)
    print('\nT {:.3f} s'.format(time.time() - lt))

    print('Matrix of hamiltonian H:')
    if n > 10:
        print('too large')
    else:
        t.print_mat(H)

    # print('H is strongly regular: {}'.format(t.is_strong_reg(H)))

    res = t.eigenLU(H, 500)
    print('Eigenvalues:\n {}'.format(res))

    print('\nT {:.3f} s'.format(time.time() - lt))


def compute_simple(n):
    H = dvr.hamilton_matrix(n)
    return t.eigenLU(H, 500)[0:3]


def err_integral(n, dN):
    eigen = t.eigenLU(dvr.hamilton_matrix(n), 500)[0:3]
    dvr.N += dN

    diff = t.minus(eigen, t.eigenLU(dvr.hamilton_matrix(n), 500)[0:3])
    return [d / dN for d in diff]


def err_dimension(n, dn):
    eigen = t.eigenLU(dvr.hamilton_matrix(n), 500)[0:3]

    diff = t.minus(eigen, t.eigenLU(dvr.hamilton_matrix(n + dn), 500)[0:3])
    return [d / dn for d in diff]


def err_interval(n, a=0, b=0):
    eigen = t.eigenLU(dvr.hamilton_matrix(n), 500)[0:3]
    dvr.a += a
    dvr.b += b

    diff = t.minus(eigen, t.eigenLU(dvr.hamilton_matrix(n), 500)[0:3])
    return [d / (b - a) for d in diff]


dt = time.time()
rt = time.time()

# l = []
# for i in range(1, 10):
#     b = i / 10
#     print(b)
#     print('T {}'.format(time.time() - dt))
#     l.append(err_interval(20, b=b))
#
# print(l)

dvr.N = 500

# H = dvr.hamilton_matrix(10)
# t.print_mat(H)
#
# print(t.eigenLU(H, 500)[0:3])
#
# print(t.eigenJacobi(H, 500))

A = [[-1, 1, 1], [1, 1, -1], [1, -1, 1]]

print(t.eigenLU(A, 500))

print(t.eigenJacobi(A, 500))
# lis = []
# times = []
# for i in [100, 110, 130]:
#     lis.append(compute_simple(i))
#     print('{}: {:.3f} s'.format(i, time.time() - rt))
#     times.append(time.time() - rt)
#     rt = time.time()
#
# print([lis[i][0] for i in range(len(lis))])
# print([lis[i][1] for i in range(len(lis))])
# print([lis[i][2] for i in range(len(lis))])
# print(times)

print('\nTotal time: {:.3f} s'.format(time.time() - dt))
