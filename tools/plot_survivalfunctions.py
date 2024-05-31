import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats
from mpl_toolkits.axisartist.axislines import SubplotZero


fig = plt.figure(figsize=(5, 3))
ax = SubplotZero(fig, 111)
fig.add_subplot(ax)

helmholtzlightblue = (20/255, 200/255, 255/255)

x = np.linspace(0, 4.4, 1000)
ax.plot(x, 1-scipy.stats.lognorm.cdf(x, s=0.5, scale=1),
        color=helmholtzlightblue)
ax.set_xlabel(r'$\tau $ ', fontsize=50)
ax.set_ylabel(r'$\gamma_I^{\,R}(\tau)$', fontsize=50)
ax.set_xlim(xmin=0, xmax=4.5)
ax.set_ylim(ymin=0, ymax=1.1)
ax.set_xticks([0, 1, 2, 3, 4])
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

for direction in ["xzero", "yzero"]:
    # adds arrows at the ends of each axis
    ax.axis[direction].set_axisline_style("-|>")

    # adds X and Y-axis from the origin
    ax.axis[direction].set_visible(True)

for direction in ["left", "right", "bottom", "top"]:
    # hides borders
    ax.axis[direction].set_visible(False)


fig.savefig('plots/survivalfunctions/lognorm.png',
            bbox_inches='tight', dpi=500)
plt.show()

fig = plt.figure(figsize=(5, 3))
ax = SubplotZero(fig, 111)
fig.add_subplot(ax)

x = np.linspace(0, 4.4, 1000)
ax.plot(x, 1-scipy.stats.expon.cdf(x, scale=1),
        color=helmholtzlightblue)
ax.set_xlabel(r'$\tau $ ', fontsize=50)
ax.set_ylabel(r'$\gamma_I^{\,R}(\tau)$', fontsize=50)
ax.set_xlim(xmin=0, xmax=4.5)
ax.set_ylim(ymin=0, ymax=1.1)
ax.set_xticks([0, 1, 2, 3, 4])
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

for direction in ["xzero", "yzero"]:
    # adds arrows at the ends of each axis
    ax.axis[direction].set_axisline_style("-|>")

    # adds X and Y-axis from the origin
    ax.axis[direction].set_visible(True)

for direction in ["left", "right", "bottom", "top"]:
    # hides borders
    ax.axis[direction].set_visible(False)


fig.savefig('plots/survivalfunctions/exponential.png',
            bbox_inches='tight', dpi=500)
plt.show()

fig = plt.figure(figsize=(5, 3))
ax = SubplotZero(fig, 111)
fig.add_subplot(ax)

x = np.linspace(0, 4.4, 1000)
ax.plot(x, 1-scipy.stats.gamma.cdf(x, a=2, scale=0.8),
        color=helmholtzlightblue)
ax.set_xlabel(r'$\tau $ ', fontsize=50)
ax.set_ylabel(r'$\gamma_I^{\,R}(\tau)$', fontsize=50)
ax.set_xlim(xmin=0, xmax=4.5)
ax.set_ylim(ymin=0, ymax=1.1)
ax.set_xticks([0, 1, 2, 3, 4])
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

for direction in ["xzero", "yzero"]:
    # adds arrows at the ends of each axis
    ax.axis[direction].set_axisline_style("-|>")

    # adds X and Y-axis from the origin
    ax.axis[direction].set_visible(True)

for direction in ["left", "right", "bottom", "top"]:
    # hides borders
    ax.axis[direction].set_visible(False)


fig.savefig('plots/survivalfunctions/gamma.png',
            bbox_inches='tight', dpi=500)
plt.show()
