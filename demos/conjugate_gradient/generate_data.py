import numpy as np

# dimension
n = 100

# line coeffs
c = np.random.random(n + 1)

# line fcn
# f(x) = c0 + c1*x1 + ... + cn*xn
def f(x):
    return c[0] + np.dot(c[1:], x)


# num data pts
m = 500

# independent vectors x
x = np.array([np.random.random(n) for i in range(m)])

# dependent values y, with noise
y = np.array([f(xi) + np.random.random() for xi in x])

# zip x and y into array z
z = np.empty([m, n + 1])
for i in range(m):
    z[i, :n] = x[i]
    z[i, -1] = y[i]

# write to file
np.savetxt('conj_grad_data.txt', z)
