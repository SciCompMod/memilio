"""
script to plot the density of the Erlang distribution and the survival function for different n or variances
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.stats as stats


def plot_erlang_survival(ns):
    """ Plots the survival functions of the Erlang distribution with different numbers n.
        The plot is saved in a folder "Plots/".
        The expected value of the distributions is set to 10.."""

    gamma = 1.0/10.0
    t = np.arange(0, 30, 1 / 30)

    plt.figure('Survivalfunction')

    for n in ns:
        # a is the shape parameter, scale the inverse rate parameter
        y = 1-stats.gamma.cdf(t, a=n, scale=1/(n*gamma))
        plt.plot(t, y, linewidth=1.2, label="n="+str(n))

    plt.legend(fontsize=14)
    plt.xlim(left=0, right=30)
    plt.ylim(bottom=0)
    plt.xlabel('Time $\\tau$ [days]', fontsize=13)
    plt.ylabel('$1-F_{n\\, /\\, T,n}(\\tau)$', fontsize=14)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/survival.png', bbox_inches='tight', dpi=500)
    plt.close()


def plot_erlang_density(ns):
    """ Plots the density of the Erlang distribution with different numbers n or variances.
        The plot is saved in a "Plots/" folder.
        The expected value of the distributions is set to 10."""

    gamma = 1.0/10.0
    t = np.arange(0, 30, 1 / 30)

    plt.figure('density')

    for n in ns:
        y = stats.gamma.pdf(t, a=n, scale=1/(n*gamma))
        plt.plot(t, y, linewidth=1.2, label="n="+str(n))

    plt.legend(fontsize=14)
    plt.xlim(left=0, right=30)
    plt.ylim(bottom=0)
    plt.xlabel('Time $\\tau$ [days]', fontsize=13)
    plt.ylabel('$f_{n\\, /\\, T,n}(\\tau)$', fontsize=14)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/density.png', bbox_inches='tight', dpi=500)
    plt.close()


def plot_contact_matrix():
    """ """
    contact_pattern = np.matrix([[3.9547, 1.1002, 2.9472,  2.05, 0.3733, 0.0445],
                                 [0.3327, 3.5892, 1.236, 1.9208, 0.2681, 0.0161],
                                 [0.246, 0.7124, 5.6518, 3.2939, 0.2043, 0.0109],
                                 [0.1742, 0.8897, 3.3124, 4.5406, 0.4262, 0.0214],
                                 [0.0458, 0.1939, 0.5782, 1.3825, 1.473, 0.0704],
                                 [0.1083, 0.1448, 0.4728, 0.9767, 0.6266, 0.1724]])
    plt.imshow(contact_pattern, vmin=0, cmap='plasma')
    plt.colorbar()
    plt.xticks(ticks=range(6), labels=[
               '0–4', '5–14', '15–34', '35–59', '60–79', '80+'])
    plt.gca().xaxis.set_ticks_position('top')
    plt.yticks(ticks=range(6), labels=[
               '0–4', '5–14', '15–34', '35–59', '60–79', '80+'])

    for i in range(6):
        for j in range(6):
            plt.text(j, i, round(contact_pattern[i, j], 2),
                     ha="center", va="center", color="w")

    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/contact_pattern.png', bbox_inches='tight', dpi=500)


def deviations_final_sizes():
    # Define final sizes for different reproduction numbers.
    # Order: ODE, LCT3, LCT10, LCT50
    order = [1, 3, 10, 50]
    final_size2 = [66187889.688620, 66177703.545307,
                   66173558.070585, 66172050.325873]
    final_size4 = [81489437.728607, 81487771.258965,
                   81487273.028473, 81487137.117737]
    final_size10 = [83151138.095973, 83151130.429073,
                    83151128.860178, 83151128.506111]
    print("Final size deviations for R=2:")
    for i in range(1, len(order)):
        print(str(i)+": ", end=" ")
        result = 100*(final_size2[i]-final_size2[0])/final_size2[0]
        print(f"{result:.5f}")

    print("Final size deviations for R=4:")
    for i in range(1, len(order)):
        print(str(i)+": ", end=" ")
        result = 100*(final_size4[i]-final_size4[0])/final_size4[0]
        print(f"{result:.5f}")

    print("Final size deviations for R=10:")
    for i in range(1, len(order)):
        print(str(i)+": ", end=" ")
        result = 100*(final_size10[i]-final_size10[0])/final_size10[0]
        print(f"{result:.8f}")


if __name__ == '__main__':
    deviations_final_sizes()
    # ns = list([1, 3, 5, 10, 20, 50, 150, 300])
    # plot_erlang_survival(ns)
    # plot_erlang_density(ns)
    # plot_contact_matrix()
