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

"""@plot_results_lct_secir_fictional.py
Functions to plot and compare results of simulations with different kind of models,
eg LCT, IDE or ODE SECIR models without division in agegroups.
To compare with real data, use the script plot_results_lct_secir_real.py.

The data to be plotted should be stored in a '../data/simulation_lct' folder as .h5 files.
Data could be generated eg by executing the file ./cpp/examples/lct_secir_fictional_scenario.cpp.
"""

import h5py
import os
import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import Normalize

import memilio.epidata.getDataIntoPandasDataFrame as gd

# Define compartments.
secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
              5: 'ICU', 6: 'Recovered', 7: 'Dead'}

# Define color and style to be used while plotting for different models to make plots consistent.
color_dict = {'ODE': '#1f77b4',
              'LCT3': '#2ca02c',
              'LCT10': '#ff7f0e',
              'LCT20': '#9467bd',
              'LCT50':  '#e377c2',
              }
linestyle_dict = {'ODE': 'solid',
                  'LCT3': 'solid',
                  'LCT10': 'solid',
                  'LCT20': 'solid',
                  'LCT50': 'solid'
                  }


def compare_single_compartment(files, legendplot,  compartment_idx=1, filename_plot="compare_single_compartment"):
    """ Plots the result of a simulation in a single Plot for specified compartments.
        The result should consist of accumulated numbers for subcompartments.

    @param[in] files: paths of the files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] compartment_idx: index of the compartment one wants to compare.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """
    # define plot
    plt.figure(filename_plot)

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
            raise gd.DataError(
                "Expected a different number of compartments.")

        # Plot result.
        if legendplot[file] in linestyle_dict:
            plt.plot(dates, total[:, compartment_idx], linewidth=1.2,
                     linestyle=linestyle_dict[legendplot[file]], color=color_dict[legendplot[file]])
        else:
            plt.plot(dates, total[:, compartment_idx], linewidth=1.2)

        h5file.close()

    plt.xlabel('Simulation Time [days]', fontsize=16)
    plt.yticks(fontsize=9)
    plt.ylabel('Number of ' +
               secir_dict[compartment_idx] + " Individuals", fontsize=14)
    # plt.ylim(bottom=0)
    plt.xlim(left=0, right=dates[-1])
    plt.legend(legendplot, fontsize=14, framealpha=0.5)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # Save result.
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


def compare_compartments_horizontal(files, legendplot,  filename_plot="compare_compartments", compartment_indices=range(8)):
    """ Creates a Plot with one subplot per compartment and one line per result one wants to compare.

    @param[in] files: paths of the files (without file extension .h5) with the simulation results
        that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models).
        Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """

    fig, axs = plt.subplots(2,
                            math.ceil(len(compartment_indices)/2), sharex='all', num=filename_plot, tight_layout=True)

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
        # Plot result.
        left_right = 0
        up_down = 0
        for i in compartment_indices:
            if legendplot[file] in linestyle_dict:
                axs[up_down, left_right].plot(dates,
                                              total[:, i], label=legendplot[file], linewidth=1.2,
                                              linestyle=linestyle_dict[legendplot[file]],
                                              color=color_dict[legendplot[file]])
            else:
                axs[up_down, left_right].plot(dates,
                                              total[:, i], label=legendplot[file], linewidth=1.2)
            axs[up_down, left_right].set_title(
                secir_dict[i], fontsize=8, pad=3)
            if (left_right < math.ceil(len(compartment_indices)/2)-1):
                left_right += 1
            else:
                left_right = 0
                up_down += 1
        h5file.close()

    # Define some characteristics of the plot.
    for i in range(math.ceil(len(compartment_indices)/2)*2):
        axs[i % 2, int(i/2)].set_xlim(left=0, right=dates[-1])
        if legendplot[0] == "ODE":
            axs[i % 2, int(i/2)].set_ylim(bottom=0)
        axs[i % 2, int(i/2)].grid(True, linestyle='--')
        axs[i % 2, int(i/2)].tick_params(axis='y', labelsize=7)
        axs[i % 2, int(i/2)].tick_params(axis='x', labelsize=7)

    fig.supxlabel('Simulation Time [days]', fontsize=9)

    lines, labels = axs[0, 0].get_legend_handles_labels()
    lgd = fig.legend(lines, labels, ncol=len(legendplot),  loc='outside lower center',
                     fontsize=8, bbox_to_anchor=(0.5, - 0.065), bbox_transform=fig.transFigure)
    fig.set_size_inches(7.5/3*math.ceil(len(compartment_indices)/2), 10.5/2.5)
    fig.tight_layout(pad=0, w_pad=0.3, h_pad=0.4)
    fig.subplots_adjust(bottom=0.09)

    # Save result.
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    fig.savefig('Plots/'+filename_plot+'.png',
                bbox_extra_artists=(lgd,),  bbox_inches='tight', dpi=500)


def compare_all_compartments(files, legendplot,  filename_plot="compare_compartments", compartment_indices=range(8)):
    """ Creates a 4x2 Plot with one subplot per compartment and one line per result one wants to compare.

    @param[in] files: paths of the files (without file extension .h5) with the simulation results
        that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models).
        Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """

    fig, axs = plt.subplots(
        math.ceil(len(compartment_indices)/2), 2, sharex='all', num=filename_plot, tight_layout=True)

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
        # Plot result.
        left_right = 0
        index = 0
        if legendplot[file] in linestyle_dict:
            for i in compartment_indices:
                axs[int(index/2), left_right % 2].plot(dates,
                                                       total[:, i], label=legendplot[file], linewidth=1.2,
                                                       linestyle=linestyle_dict[legendplot[file]],
                                                       color=color_dict[legendplot[file]])
                axs[int(index/2), left_right %
                    2].set_title(secir_dict[i], fontsize=8)
                left_right += 1
                index += 1
        else:
            for i in compartment_indices:
                axs[int(index/2), left_right % 2].plot(dates,
                                                       total[:, i], label=legendplot[file], linewidth=1.2)
                axs[int(index/2), left_right %
                    2].set_title(secir_dict[i], fontsize=8)
                left_right += 1
                index += 1
        h5file.close()

    # Define some characteristics of the plot.
    for i in range(math.ceil(len(compartment_indices)/2)*2):
        axs[int(i/2), i % 2].set_xlim(left=0, right=dates[-1])
        axs[int(i/2), i % 2].set_ylim(bottom=0)
        axs[int(i/2), i % 2].grid(True, linestyle='--')
        axs[int(i/2), i % 2].tick_params(axis='y', labelsize=7)
        axs[int(i/2), i % 2].tick_params(axis='x', labelsize=7)
        # axs[int(i/2), i % 2].xaxis.set_ticks(np.arange(0, dates[-1]+1, 5))

    fig.supxlabel('Simulation Time [days]', fontsize=9)

    lines, labels = axs[0, 0].get_legend_handles_labels()
    lgd = fig.legend(lines, labels, ncol=len(legendplot),  loc='outside lower center',
                     fontsize=10, bbox_to_anchor=(0.5, - 0.06), bbox_transform=fig.transFigure)

    plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.1)
    plt.subplots_adjust(bottom=0.09)

    # Save result.
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    fig.savefig('Plots/'+filename_plot+'.png',
                bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=500)


def plot_subcompartments3D(file, subcompartments, compartment_idx, first_time_idx, filename_plot="compare_new_infections"):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    h5file = h5py.File(str(file) + '.h5', 'r')

    if not ((compartment_idx > 0) and (compartment_idx < 7)):
        raise Exception("Use valid compartment index.")
    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]
    dates = data['Time'][:]
    # As there should be only one Group, total is the simulation result.
    total = data['Total'][:, :]
    first_idx_comp = 1+(compartment_idx-1)*subcompartments

    # Make data
    x = [i-0.5 for i in range(1, subcompartments+1)]
    x, y = np.meshgrid(x, dates[first_time_idx:]-0.5)
    flat_x, flat_y = x.ravel(), y.ravel()
    comp_total = 0
    z = np.zeros(x.shape)
    for i in range(x.shape[0]):
        comp_total = np.sum(
            total[first_time_idx+i, first_idx_comp:first_idx_comp+subcompartments])
        for j in range(x.shape[1]):
            z[i, j] = total[first_time_idx+i, first_idx_comp+j]
    flat_z = z.ravel()
    bottom = np.zeros_like(flat_x)

    cmap = cm.get_cmap('plasma')
    norm = Normalize(vmin=min(flat_z), vmax=max(flat_z))
    colors = cmap(norm(flat_z))

    ax.bar3d(flat_x, flat_y, 0, 1, 1, flat_z, shade=True, color=colors)

    ax.set_xlabel('Index of Subcompartment', labelpad=-5, fontsize=9)
    a = list(range(0, subcompartments+1, 5))
    a[0] += 1
    ax.set_xticks(a)
    ax.set_xlim(left=0.5, right=subcompartments+0.5)
    ax.tick_params(axis='x', pad=-5, labelsize=9)
    ax.set_xticks(range(1, subcompartments), minor=True)

    ax.set_ylabel('Simulation Time [days]', labelpad=-5, fontsize=9)
    ax.set_ylim(bottom=first_time_idx-0.5, top=dates[-1]+0.5)
    ax.set_yticks(range(int(dates[first_time_idx]),
                  int(dates[-1])+1, 2))
    ax.set_yticklabels(range(int(dates[first_time_idx]),
                             int(dates[-1])+1, 2), ha='left', va='center')
    ax.set_yticks(range(int(dates[first_time_idx]),
                        int(dates[-1])+1), minor=True)
    ax.tick_params(axis='y', pad=-4, labelsize=9)

    ax.set_zlim(bottom=0)
    ax.tick_params(axis='z', pad=1, labelsize=9)

    sc = cm.ScalarMappable(cmap=cmap, norm=norm)
    sc.set_array([])
    cbar = plt.colorbar(sc, ax=ax, shrink=0.75)

    h5file.close()
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


def plot_maxpeak_incidence(func_get_file_name, R0s, subcomps, filename_plot="maxpeak2d"):
    fig = plt.figure()

    for R0 in R0s:
        y = np.zeros(len(subcomps))
        for i in range(len(subcomps)):
            y[i] = get_maxpeak_incidence(
                func_get_file_name(R0=R0, subcompartment=subcomps[i]))
        plt.plot(subcomps, y, 'o-', linewidth=1.2,
                 label=str(R0))

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10,
               title="$\\mathcal{R}_{\\text{eff}}\\approx$", title_fontsize=14)
    plt.ylim(bottom=0)
    plt.xlabel('Number of Subcompartments', fontsize=13)
    plt.ylabel('Maximum daily new Transmissions', fontsize=13)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


def get_maxpeak_incidence(file):
    # Load data.
    h5file = h5py.File(str(file) + '.h5', 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]
    dates = data['Time'][:]
    # As there should be only one Group, total is the simulation result.
    total = data['Total'][:, :]
    incidence = (total[:-1, 0]-total[1:, 0])/(dates[1:]-dates[:-1])
    max = np.max(incidence)
    h5file.close()
    return max


def plot_day_peak_incidence(func_get_file_name, R0s, subcomps, filename_plot="maxpeak2d"):
    fig = plt.figure()

    for R0 in R0s:
        y = np.zeros(len(subcomps))
        for i in range(len(subcomps)):
            y[i] = get_day_peak_incidence(
                func_get_file_name(R0=R0, subcompartment=subcomps[i]))
        plt.plot(subcomps, y, 'x-', linewidth=1.2,
                 label=str(R0))

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10,
               title="$\\mathcal{R}_{\\text{eff}}\\approx$", title_fontsize=14)
    plt.ylim(bottom=0)
    plt.xlabel('Number of Subcompartments', fontsize=13)
    plt.ylabel('Simulation Day of Peak [days]', fontsize=13)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


def get_day_peak_incidence(file):
    # Load data.
    h5file = h5py.File(str(file) + '.h5', 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]
    dates = data['Time'][:]
    # As there should be only one Group, total is the simulation result.
    total = data['Total'][:, :]
    incidence = (total[:-1, 0]-total[1:, 0])/(dates[1:]-dates[:-1])
    argmax = np.argmax(incidence)
    h5file.close()
    return dates[argmax]


def plot_new_infections(files, ylim, legendplot, filename_plot="compare_new_infections", tmax=0):
    """ Single plot to compare the incidence of different results. Incidence means the number of people leaving the
        susceptible class per day.

    @param[in] files: paths of the files (without file extension .h5) with the simulation results
         that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models).
        Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] ylim: upper limit for the y-axis of the plot.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """
    plt.figure(filename_plot)

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
        incidence = (total[:-1, 0]-total[1:, 0])/(dates[1:]-dates[:-1])
        # Plot result.
        if legendplot[file] in linestyle_dict:
            plt.plot(dates[1:], incidence, linewidth=1.2,
                     linestyle=linestyle_dict[legendplot[file]], color=color_dict[legendplot[file]])
        else:
            plt.plot(dates[1:], incidence, linewidth=1.2)

        h5file.close()

    plt.xlabel('Simulation Time [days]', fontsize=16)
    plt.yticks(fontsize=9)
    plt.ylabel('Daily new Transmissions', fontsize=14)
    plt.ylim(bottom=0, top=ylim)
    if tmax > 0:
        plt.xlim(left=0, right=tmax)
    else:
        plt.xlim(left=0, right=dates[-1])
    plt.legend(legendplot, fontsize=14, framealpha=0.5)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # Save result.
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


if __name__ == '__main__':
    # simulation results should be stored in this folder.
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "data", "simulation_lct_noage")

    cases = [0, 1, 2, 3, 4]
    for case in cases:
        if case == 0:
            # rise R0 short
            folder = os.path.join(data_dir, "riseR0short")
            plot_new_infections([os.path.join(folder, "fictional_lct_2.0_1"), os.path.join(folder, "fictional_lct_2.0_3"),
                                os.path.join(folder, "fictional_lct_2.0_10"),
                                os.path.join(folder, "fictional_lct_2.0_50")],
                                20000, legendplot=list(["ODE", "LCT3", "LCT10", "LCT50"]),
                                filename_plot="new_infections_rise2.0")
            compare_single_compartment([os.path.join(folder, "fictional_lct_2.0_1"), os.path.join(folder, "fictional_lct_2.0_3"),
                                        os.path.join(
                                            folder, "fictional_lct_2.0_10"),
                                        os.path.join(folder, "fictional_lct_2.0_50")],
                                       legendplot=list(
                ["ODE", "LCT3", "LCT10", "LCT50"]),
                compartment_idx=2, filename_plot="carrier_compartment_rise2.0")
            compare_single_compartment([os.path.join(folder, "fictional_lct_2.0_1"), os.path.join(folder, "fictional_lct_2.0_3"),
                                        os.path.join(
                                            folder, "fictional_lct_2.0_10"),
                                        os.path.join(folder, "fictional_lct_2.0_50")],
                                       legendplot=list(
                ["ODE", "LCT3", "LCT10", "LCT50"]),
                compartment_idx=3, filename_plot="infected_compartment_rise2.0")

        elif case == 1:
            # drop R0 short
            folder = os.path.join(data_dir, "dropR0")
            plot_new_infections([os.path.join(folder, "fictional_lct_0.5_1"), os.path.join(folder, "fictional_lct_0.5_3"),
                                os.path.join(folder, "fictional_lct_0.5_10"),
                                os.path.join(folder, "fictional_lct_0.5_50")],
                                4100, legendplot=list(["ODE", "LCT3", "LCT10", "LCT50"]),
                                filename_plot="new_infections_drop0.5")
            compare_single_compartment([os.path.join(folder, "fictional_lct_0.5_1"),
                                        os.path.join(
                                            folder, "fictional_lct_0.5_3"),
                                        os.path.join(
                                            folder, "fictional_lct_0.5_10"),
                                        os.path.join(folder, "fictional_lct_0.5_50")],
                                       legendplot=list(
                ["ODE", "LCT3", "LCT10", "LCT50"]),
                compartment_idx=2, filename_plot="carrier_compartment_drop0.5")
            compare_single_compartment([os.path.join(folder, "fictional_lct_0.5_1"),
                                        os.path.join(folder,
                                                     "fictional_lct_0.5_3"),
                                        os.path.join(folder,
                                                     "fictional_lct_0.5_10"),
                                        os.path.join(folder,
                                                     "fictional_lct_0.5_50")
                                        ],
                                       legendplot=list(
                ["ODE", "LCT3", "LCT10", "LCT50"]),
                compartment_idx=3, filename_plot="infected_compartment_drop0.5")
        elif case == 2:
            # rise r0 2.0 long
            folder = os.path.join(data_dir, "riseR0long")
            plot_new_infections([os.path.join(folder, "fictional_lct_2.0_1"),
                                os.path.join(folder,
                                             "fictional_lct_2.0_3"),
                                os.path.join(folder,
                                             "fictional_lct_2.0_10"),
                                os.path.join(folder,
                                             "fictional_lct_2.0_50")
                                 ],
                                2.1e6, legendplot=list(["ODE", "LCT3", "LCT10", "LCT50"]),
                                filename_plot="new_infections_rise2.0_long", tmax=150)
            compare_compartments_horizontal([os.path.join(folder, "fictional_lct_2.0_1"),
                                            os.path.join(folder,
                                                         "fictional_lct_2.0_3"),
                                            os.path.join(folder,
                                                         "fictional_lct_2.0_10"),
                                            os.path.join(folder,
                                                         "fictional_lct_2.0_50")
                                             ],
                                            legendplot=list(
                ["ODE", "LCT3", "LCT10", "LCT50"]),
                filename_plot="compartments_rise2.0_long", compartment_indices=[0, 1, 2, 3, 4, 5, 6, 7])
        elif case == 3:
            # rise r0 4.0 long
            folder = os.path.join(data_dir, "riseR0long")
            plot_new_infections([os.path.join(folder, "fictional_lct_4.0_1"),
                                os.path.join(folder,
                                             "fictional_lct_4.0_3"),
                                os.path.join(folder,
                                             "fictional_lct_4.0_10"),
                                os.path.join(folder,
                                             "fictional_lct_4.0_50")
                                 ],
                                6e6, legendplot=list(["ODE", "LCT3", "LCT10", "LCT50"]),
                                filename_plot="new_infections_rise4.0_long", tmax=70)
        elif case == 4:
            data_dir = os.path.join(os.path.dirname(
                __file__), "..", "data", "simulation_lct_agevsnoage")
            # youngvsold
            compare_compartments_horizontal([os.path.join(data_dir, "fictional_lct_ageres_10_agegroupinit_2"),
                                            os.path.join(
                                                data_dir, "fictional_lct_ageres_10_agegroupinit_5"),
                                            os.path.join(data_dir, "fictional_lct_notageres_10")],
                                            legendplot=list(
                ["University Scenario", "Retirement Home Scenario", "Non age-resolved Scenario"]),
                filename_plot="compartments_agevsnoage", compartment_indices=[0, 1, 2, 3, 4, 5, 6, 7])
