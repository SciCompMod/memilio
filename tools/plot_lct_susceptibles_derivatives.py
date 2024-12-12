#############################################################################
# Copyright (C) 2020-2024 MEmilio
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

import h5py
import os
import math
import matplotlib.pyplot as plt
import numpy as np

import memilio.epidata.getDataIntoPandasDataFrame as gd


"""@plot_lct_susceptibles_derivatives.py
Function to create a plot with the number of susceptibles as well as the first 
and the second derivative.
"""

# Define color and style to be used while plotting for different models to make plots consistent.
color_dict = {'ODE': '#1f77b4',
              'LCT3': '#2ca02c',
              'LCT10': '#ff7f0e',
              'LCT20': '#9467bd',
              'LCT50':  '#e377c2',
              }
fontsize_labels = 14
fontsize_legends = 14
plotfolder = 'Plots_fictional'


def derivatives_susceptibles(files, legendplot,  filename_plot="derivatives_susceptibles"):
    fig, axs = plt.subplots(
        3, 1, sharex='all', num=filename_plot, tight_layout=True)

    # Set the ylabels manually to ensure that the important information can be seen.
    ylabelsusbottom = 1e14
    ylabelsustop = 0
    ylabelderivbottom = 10000
    ylabelderivtop = 0
    ylabel2ndderivbottom = 10000
    ylabel2ndderivtop = 0

    # Add simulation results to plot.
    for file in range(len(files)):
        # Load data.
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result.
        total = data['Total'][:, :]
        if (total.shape[1] != 8):
            raise gd.DataError("Expected a different number of compartments.")

        incidence = (total[:-1, 0]-total[1:, 0])/(dates[1:]-dates[:-1])
        incidencederivative = (incidence[1:] -
                               incidence[:-1])/(dates[2:]-dates[1:-1])
        index = np.where(dates > 2.05)[0][0]
        if (np.min(total[index:, 0]) < ylabelsusbottom):
            ylabelsusbottom = np.min(total[index:, 0])
        if (np.max(total[index:, 0]) > ylabelsustop):
            ylabelsustop = np.max(total[index:, 0])
        if (np.min(incidence[index:]) < ylabelderivbottom):
            ylabelderivbottom = np.min(incidence[index:])
        if (np.max(incidence[index:]) > ylabelderivtop):
            ylabelderivtop = np.max(incidence[index:])
        if (np.min(incidencederivative[index:]) < ylabel2ndderivbottom):
            ylabel2ndderivbottom = np.min(incidencederivative[index:])
        if (np.max(incidencederivative[index:]) > ylabel2ndderivtop):
            ylabel2ndderivtop = np.max(incidencederivative[index:])

        # Plot result.
        if legendplot[file] in color_dict:
            axs[0].plot(dates, total[:, 0], label=legendplot[file], linewidth=1.2,
                        linestyle="solid",
                        color=color_dict[legendplot[file]])
            axs[1].plot(dates[1:], incidence, label=legendplot[file], linewidth=1.2,
                        linestyle="solid",
                        color=color_dict[legendplot[file]])
            axs[2].plot(dates[2:], incidencederivative, label=legendplot[file], linewidth=1.2,
                        linestyle="solid",
                        color=color_dict[legendplot[file]])
        else:
            axs[0].plot(dates, total[:, i],
                        label=legendplot[file], linewidth=1.2)
            axs[1].plot(dates[1:], incidence,
                        label=legendplot[file], linewidth=1.2)
            axs[2].plot(dates[2:], incidencederivative,
                        label=legendplot[file], linewidth=1.2)

        h5file.close()

    # Define some characteristics of the plot.
    axs[0].set_title("Susceptibles", fontsize=8, pad=3)
    axs[0].set_ylabel("S(t)", fontsize=8)
    axs[0].set_ylim(bottom=ylabelsusbottom-200, top=ylabelsustop+200)

    axs[1].set_title("New transmissions", fontsize=8, pad=3)
    axs[1].set_ylabel("-S'(t)", fontsize=8)
    axs[1].set_ylim(bottom=ylabelderivbottom*0.96, top=ylabelderivtop*1.04)

    axs[2].set_title("Derivative of new transmissions", fontsize=8, pad=3)
    axs[2].set_ylabel("-S''(t)", fontsize=8)
    axs[2].set_ylim(bottom=ylabel2ndderivbottom-10, top=ylabel2ndderivtop+10)

    for i in range(3):
        axs[i].set_xlim(left=2, right=dates[-1])
        axs[i].grid(True, linestyle='--')
        axs[i].tick_params(axis='y', labelsize=7)
        axs[i].tick_params(axis='x', labelsize=7)

    fig.supxlabel('Simulation time [days]', fontsize=9)

    lines, labels = axs[0].get_legend_handles_labels()
    lgd = fig.legend(lines, labels, ncol=len(legendplot),  loc='outside lower center',
                     fontsize=8, bbox_to_anchor=(0.5, - 0.065), bbox_transform=fig.transFigure)
    fig.set_size_inches(6, 5.5)
    fig.tight_layout(pad=0, w_pad=0.3, h_pad=0.4)
    fig.subplots_adjust(bottom=0.09)

    # Save result.
    fig.savefig(plotfolder+'/'+filename_plot+'.png',
                bbox_extra_artists=(lgd,),  bbox_inches='tight', dpi=500)
    plt.close()


if __name__ == '__main__':
    if not os.path.isdir(plotfolder):
        os.makedirs(plotfolder)
    # Susceptibles with derivatives for rise R0 short and drop.

    # Simulation results should be stored in this folder.
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "data", "simulation_lct_noage")
    folder = os.path.join(data_dir, "riseR0short")
    derivatives_susceptibles([os.path.join(folder, "fictional_lct_2.0_1"), os.path.join(folder, "fictional_lct_2.0_3"),
                              os.path.join(folder, "fictional_lct_2.0_10"), os.path.join(folder, "fictional_lct_2.0_50")],
                             legendplot=list(["ODE", "LCT3", "LCT10", "LCT50"]), filename_plot="susceptibles_derivatives_rise")
    folder = os.path.join(data_dir, "dropR0")
    derivatives_susceptibles([os.path.join(folder, "fictional_lct_0.5_1"), os.path.join(folder, "fictional_lct_0.5_3"),
                              os.path.join(folder, "fictional_lct_0.5_10"), os.path.join(folder, "fictional_lct_0.5_50")],
                             legendplot=list(["ODE", "LCT3", "LCT10", "LCT50"]), filename_plot="susceptibles_derivatives_drop")
