"""
script to plot the density of the Erlang distribution and the survival function for different n or variances
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.stats as stats


def plot_erlang_survival():
    """ Plots the survival functions of the Erlang distribution with different numbers n.
        The plot is saved in a folder "Plots/".
        The expected value of the distributions is set to 10.."""

    gamma = 1.0/10.0
    t = np.arange(0, 30, 1 / 30)

    plt.figure('Survivalfunction')

    for n in [1, 5, 10, 20, 150, 300]:
        # a is the shape parameter, scale the inverse rate parameter
        y = 1-stats.gamma.cdf(t, a=n, scale=1/(n*gamma))
        plt.plot(t, y, linewidth=1.2, label="n="+str(n))

    plt.legend(fontsize=14)
    plt.xlim(left=0, right=30)
    plt.ylim(bottom=0)
    plt.xlabel('Time $\\tau$ (in days)', fontsize=13)
    plt.ylabel('$1-F_{n\\, /\\, T,n}(\\tau)$', fontsize=14)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if not os.path.isdir('plots'):
        os.makedirs('plots')
    plt.savefig('plots/survival.png', bbox_inches='tight', dpi=500)


def plot_erlang_density():
    """ Plots the density of the Erlang distribution with different numbers n or variances.
        The plot is saved in a "Plots/" folder.
        The expected value of the distributions is set to 10."""

    gamma = 1.0/10.0
    t = np.arange(0, 30, 1 / 30)

    plt.figure('density')

    for n in [1, 5, 10, 20, 150, 300]:
        y = stats.gamma.pdf(t, a=n, scale=1/(n*gamma))
        plt.plot(t, y, linewidth=1.2, label="n="+str(n))

    plt.legend(fontsize=14)
    plt.xlim(left=0, right=30)
    plt.ylim(bottom=0)
    plt.xlabel('Time $\\tau$ (in days)', fontsize=13)
    plt.ylabel('$f_{n\\, /\\, T,n}(\\tau)$', fontsize=14)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if not os.path.isdir('plots'):
        os.makedirs('plots')
    plt.savefig('plots/density.png', bbox_inches='tight', dpi=500)


if __name__ == '__main__':
    plot_erlang_survival()
    plot_erlang_density()
