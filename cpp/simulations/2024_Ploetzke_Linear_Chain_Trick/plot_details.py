#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Lena Ploetzke
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
"""@plot_details.py
Functions to create plots that do not use simulation results.
There is one function to visualize the age resolved contact pattern.
Moreover, there is functionality to plot the density and the survival function of
the Erlang-distribution with different parameters.
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.stats as stats

plotfolder = 'Plots'


def plot_erlang_survival(ns, file_name=''):
    """ Plots the survival functions of Erlang distributions with rate parameter n/mean
     and shape parameter n. Different choices for n are plotted.
    The mean value of the distribution is set to 10 in this function.

    @param[in] ns: Vector with values for n. For each entry, one Erlang distribution is added to the plot.
    @param[in] file_name: The name of the file where the plot will be saved.
            If an empty string is provided, the plot will not be saved.
    """
    mean = 10.0
    t = np.arange(0, 30, 1 / 30)

    plt.figure(file_name)

    for n in ns:
        # a is the shape parameter, scale the inverse rate parameter.
        y = 1-stats.gamma.cdf(t, a=n, scale=mean/n)
        plt.plot(t, y, linewidth=1.2, label="n="+str(n))

    plt.legend(fontsize=14)
    plt.xlim(left=0, right=30)
    plt.ylim(bottom=0)
    plt.xlabel('Time $\\tau$ [days]', fontsize=13)
    plt.ylabel('$1-F_{n\\, /\\, T,n}(\\tau)$', fontsize=14)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if file_name:
        plt.savefig(plotfolder+'/'+file_name+'.png',
                    bbox_inches='tight', dpi=500)
    plt.close()


def plot_erlang_density(ns, file_name=''):
    """ Plots the density of Erlang distributions with rate parameter n/mean
     and shape parameter n. Different choices for n are plotted.
    The mean value of the distribution is set to 10 in this function.

    @param[in] ns: Vector with values for n. For each entry, one Erlang distribution is added to the plot.
    @param[in] file_name: The name of the file where the plot will be saved.
            If an empty string is provided, the plot will not be saved.
    """
    mean = 10.0
    t = np.arange(0, 30, 1 / 30)

    plt.figure(file_name)

    for n in ns:
        # a is the shape parameter, scale the inverse rate parameter.
        y = stats.gamma.pdf(t, a=n, scale=mean/n)
        plt.plot(t, y, linewidth=1.2, label="n="+str(n))

    plt.legend(fontsize=14)
    plt.xlim(left=0, right=30)
    plt.ylim(bottom=0)
    plt.xlabel('Time $\\tau$ [days]', fontsize=13)
    plt.ylabel('$f_{n\\, /\\, T,n}(\\tau)$', fontsize=14)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if file_name:
        plt.savefig(plotfolder+'/'+file_name+'.png',
                    bbox_inches='tight', dpi=500)
    plt.close()


def plot_contact_matrix(file_name=''):
    """ Visualizes the age resolved contact pattern for Germany. For each RKI age group represented in the rows,
    the average number of daily contacts with the respective age group is provided in the columns.

    @param[in] file_name: The name of the file where the plot will be saved.
            If an empty string is provided, the plot will not be saved.
    """
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

    if file_name:
        plt.savefig(plotfolder+'/'+file_name+'.png',
                    bbox_inches='tight', dpi=500)
    plt.close()


def main():
    if not os.path.isdir(plotfolder):
        os.makedirs(plotfolder)
    ns = list([1, 3, 5, 10, 20, 50, 150, 300])
    plot_erlang_survival(ns, "survival")
    plot_erlang_density(ns, "density")
    plot_contact_matrix("contact_pattern")


if __name__ == "__main__":
    main()
