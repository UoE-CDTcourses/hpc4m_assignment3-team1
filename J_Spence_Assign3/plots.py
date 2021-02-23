import numpy as np
import matplotlib.pyplot as plt

M = 2306
fs = 14
x = np.linspace(-1, 1, M)
X, Y = np.meshgrid(x, x)
# plot strong error 
times = np.loadtxt('runtimes-hor.txt')
timesvert = np.loadtxt('runtimes-vert.txt')
timessq = np.loadtxt('runtimes-sq.txt')
fig, ax = plt.subplots(1, 1, figsize = (6,6))
ax.plot(times[:,0], times[times[:,0] == 1,1]/times[:,1], 'r-o', label = 'horizontal')
ax.plot(timesvert[:,0], timesvert[timesvert[:,0] == 1,1]/timesvert[:,1], 'b-o', label = 'vertical')
ax.plot(timessq[:,0], timessq[timessq[:,0] == 1,1]/timessq[:,1], 'g-o', label = 'square')
ax.plot(np.linspace(1, 32), np.linspace(1, 32), 'k--')
ax.set_xlabel('Number of processors', fontsize = fs)
ax.set_ylabel('Speed-up', fontsize = fs)
ax.legend()
fig.savefig('strong_scaling.png', dpi = 512)

fig, ax = plt.subplots(1, 1, figsize = (6,6))
ax.plot(times[:,0], times[times[:,0] == 1,1]/(times[:,0]*times[:,1]), 'r-o', label = 'Horizontal')
ax.plot(timesvert[:,0], timesvert[timesvert[:,0] == 1,1]/(timesvert[:,0]*timesvert[:,1]), 'b-o', label = 'vertical')
ax.plot(timessq[:,0], timessq[timessq[:,0] == 1,1]/(timessq[:,0]*timessq[:,1]), 'g-o', label = 'square')
ax.plot(np.linspace(1, 32,2), np.ones(2), 'k--')
ax.set_xlabel('Number of processors', fontsize = fs)
ax.set_ylabel('Parallel Efficiency', fontsize = fs)
ax.legend()
fig.savefig('strong_efficiency.png', dpi = 512)




# Plot initial condition 
U = np.loadtxt('initial-cond.txt')
fig, ax = plt.subplots(subplot_kw = {"projection":"3d"})
ax.plot_surface(X, Y, U, rstride = 10, cstride = 10, edgecolors = 'k', linewidth = 0.5)
ax.set_xlabel(r'$y$', fontsize = fs)
ax.set_ylabel(r'$x$', fontsize = fs)
ax.set_zlim((-0.5, 1))
fig.savefig('U_init.png', dpi = 512)

for i in range(3):
    # load data
    U = np.loadtxt('data' + str(i+1) + '.txt')
    fig, ax = plt.subplots(subplot_kw = {"projection":"3d"})
    ax.plot_surface(X, Y, U, rstride = 10, cstride = 10, edgecolors = 'k', linewidth = 0.5)
    ax.set_title(r'$t = {}$'.format(str((i+1)/3)[:7]))
    ax.set_xlabel(r'$x$', fontsize = fs)
    ax.set_ylabel(r'$y$', fontsize = fs)
    ax.set_zlim((-0.5, 1))
    fig.savefig('U_' + str(i) + '.png', dpi = 512)


