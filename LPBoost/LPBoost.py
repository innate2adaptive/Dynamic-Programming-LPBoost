__author__ = 'yuxinsun'

import numpy as np
from cvxopt import *

# Linear Programming Boosting
# Required Library:
# numpy
# cvxopt

# Convex optimisation with cvxopt
# Input:
# x: hypothesis space - m*n matrix
# y: desired labels - m*1 vector
# D: D in LPBoost, usually D = 1/(m*nu) - scalar
# Output:
# Optimisation results:
# u: primal classification cost - m*1 vector
# beta: primal variable - scalar
# c4.multiplier: dual weights - m*1 vector

def LPcvx(z, y, D):

    m = y.size

    # Initialise
    u = modeling.variable(m, 'u')
    u.value = matrix(np.ones(m)/m)
    beta = modeling.variable(1, 'beta')
    beta.value = 0

    # Constraints
    c1 = (modeling.sum(u) == 1)
    c2 = (u <= D)
    c3 = (u >= 0)
    c4 = (modeling.dot(matrix(z), u) <= beta)

    # Solve problems
    lp = modeling.op(beta, [c1, c2, c3, c4])
    solvers.options['show_progress']=False
    sol = lp.solve()

    return u.value, beta.value, c4.multiplier.value

def LPBoostAg(x, y, D, maxIter):

    m = y.size
    u = np.ones(m)/m

    temp = np.multiply(u, y)
    h = np.dot(temp, x)


    ind = np.argmax(h)
    crit = np.max(h)
    beta = 0
    F = []
    idx = []

    counter = 1
    while crit >= beta+10**(-6) and counter <= maxIter:

        F.append(x[:, ind]*y)
        opt = np.asarray(F)
        opt = np.transpose(opt)

        idx.append(ind)

        [u, beta, a] = LPcvx(opt, y, D)
        u = np.squeeze(np.asarray(u))

        temp = np.multiply(u, y)
        h = np.dot(temp, x)

        ind = np.argmax(h)
        crit = np.max(h)

        beta = np.asarray(beta)
        beta = beta[0]

        print('Iteration: %d, beta: %f, criterion: %f' % (counter, beta, crit))

        counter += 1

    return u, beta, a, idx


