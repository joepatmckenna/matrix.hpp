import numpy as np


# line fcn
# f(x) = c0 + c1*x1 + ... + cn*xn
def f(x):
    return c[0] + np.dot(c[1:], x)


# n: dimension
for n in range(2, 100):
    # line coeffs
    c = np.random.random(n + 1)
    np.savetxt('data/c_%i.txt' % (n, ), c)

    # num data pts
    m = 10 * n

    for dataset in range(10):
        # independent vectors x
        x = np.array([2 * (np.random.random(n) - .5) for i in range(m)])

        # dependent values y, with noise
        y = np.array([f(xi) + .2 * (np.random.random() - .5) for xi in x])

        # zip x and y into array z
        z = np.empty([m, n + 1])
        for i in range(m):
            z[i, :n] = x[i]
            z[i, -1] = y[i]
        np.savetxt('data/xy_%i_%i.txt' % (n, dataset), z)
