import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# line fcn
# f(x) = c0 + c1*x
def f(x):
    return c[0] + c[1] * x


z = np.loadtxt('data/xy_1_0.txt')
x = z[:, 0]
y = z[:, 1]

times = np.loadtxt('data/times.txt')

fig, ax = plt.subplots(1, 2, figsize=(6, 2.5))

x_grid = np.linspace(-1, 1)
# c = np.loadtxt('data/c_1.txt')
# ax[0].plot(x_grid, f(x_grid), 'r-')
c = np.loadtxt('data/c_1_min.txt')
ax[0].plot(x_grid, f(x_grid), 'k-')
ax[0].plot(x, y, 'ko')

ax[0].set_xlabel('$x$')
ax[0].set_ylabel('$y$')

ax[1].loglog(times[1:], 'k-')

ax[1].set_xlabel('$n$')
ax[1].set_ylabel('time (ms)')

plt.tight_layout()
plt.savefig('fig.png')
