import math as m
import tools as t
import time
from script import DVR
import numpy as np


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


dvr = DVR(a=-3, b=5, N=500)

dt = time.time()
rt = time.time()

# print(sorted(t.eigenLU(H, 500), reverse=True))
# print(sorted(t.eigenJacobi(H, 500), reverse=True))
# print(time.time() - rt)
#
# # rt = time.time()
# # print(sorted(t.eigenJacobi(H, 500), reverse=True))
# # print(time.time() - rt)
#
# rt = time.time()
# print(np.linalg.eig(H)[0])

# print(time.time() - rt)

lis = []
times = []
for i in [60, 70, 80]:
    # lis.append(compute_simple(i))

    H = dvr.hamilton_matrix(i)
    # t.print_mat(H)
    ev = sorted(np.linalg.eig(H)[0], reverse=True)
    lis.append(ev[:3])

    times.append(time.time() - rt)
    print('{}: {}'.format(i, time.time() - rt))
    print(ev)


print('value={}'.format([[lis[i][0] for i in range(len(lis))],
                            [lis[i][1] for i in range(len(lis))],
                            [lis[i][2] for i in range(len(lis))]]))

print('times={}'.format(times))

print('\nTotal time: {:.3f} s'.format(time.time() - dt))
