#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
"""@plot_results_lct_secir.py
functions to plot results of a simulation with a LCT SECIR model with subcompartments.
There is also a method to compare different results of different models.

The data to be plotted should be stored in a '../data/simulation_lct' folder as .h5 files.
Data could be generated eg by executing the file ./cpp/examples/lct_secir_compare.cpp.
"""

import h5py
import os
import pandas as pd
import matplotlib.pyplot as plt

import memilio.epidata.getDataIntoPandasDataFrame as gd

# Define compartments
secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
              5: 'ICU', 6: 'Recovered', 7: 'Death'}


def plot_lct_subcompartments(file, vec_subcompartments=[1, 20, 20, 20, 20, 20, 1, 1], filename_plot="LCT_subcompartments"):
    """ Plots the result of a simulation with an LCT SECIR model in an 4x2 Plot.
        In each subplot, one compartment with one line per subcompartment is plotted.

    @param[in] file: path of the file (without file extension .h5) with the simulation result of an lct model with subcompartments.
        The number of subcompartments should be fitting to the ones specified in get_subcompartments().
    @param[in] vec_subcompartments: vector with the number of subcompartments per compartment used for simulation.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """
    secir_dict_short = {0: 'S', 1: 'E', 2: 'C', 3: 'I', 4: 'H',
                        5: 'U', 6: 'R', 7: 'D'}
    # load data
    h5file = h5py.File(str(file) + '.h5', 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]
    dates = data['Time'][:]
    # As there should be only one Group, total is the simulation result
    total = data['Total'][:, :]
    if (total.shape[1] != sum(vec_subcompartments)):
        print(total.shape[1])
        print(sum(vec_subcompartments))
        raise gd.DataError("Expected a different number of subcompartments.")
    h5file.close()

    # Plot result.
    plt.style.use('seaborn-v0_8-colorblind')
    fig, axs = plt.subplots(
        4, 2, sharex='all', num=filename_plot)
    for i in range(len(vec_subcompartments)):
        idx_start = sum(vec_subcompartments[0:i])
        axs[int(i/2), i % 2].plot(dates,
                                  total[:, idx_start:idx_start+vec_subcompartments[i]], linewidth=1.5)
        axs[int(i/2), i % 2].grid(True, linestyle='--')
        legendplot = []
        if (vec_subcompartments[i] == 1):
            legendplot.append(secir_dict_short[i])
        elif (vec_subcompartments[i] < 5):
            for j in range(vec_subcompartments[i]):
                legendplot.append(secir_dict_short[i]+str(j+1))

        axs[int(i/2), i % 2].legend(legendplot, fontsize=10)
        axs[int(i/2), i % 2].set_xlim(left=0)
    fig.supxlabel('Time')
    fig.supylabel('Number of persons')
    plt.tight_layout()

    # save plot
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    fig.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


def plot_single_result(file, compartment_idx=range(8), filename_plot="LCT_compartments"):
    """ Plots the result of a simulation with an LCT SECIR model in a single Plot for specified comparments.
        The result should consist of accumulated numbers for subcompartments.

    @param[in] file: path of the file (without file extension .h5) with the simulation result of an lct
        model with accumulated numbers for subcompartments.
    @param[in] compartment_idx: indexes of the compartments that should be plot, eg in form of a list.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """

    # load data
    h5file = h5py.File(str(file) + '.h5', 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]
    dates = data['Time'][:]
    # As there should be only one Group, total is the simulation result
    total = data['Total'][:, compartment_idx]
    h5file.close()

    # plot result
    plt.figure(filename_plot)
    plt.plot(dates, total, linewidth=1.5)

    legendplot = []
    for i in compartment_idx:
        legendplot.append(secir_dict[i])
    plt.legend(legendplot, fontsize=14)

    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Number of persons', fontsize=10)
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    plt.grid(True, linestyle='--')
    plt.style.use('seaborn-v0_8-colorblind')
    plt.tight_layout()

    # save result
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


def compare_single_compartment(files, legendplot,  compartment_idx=1, filename_plot="compare_single_compartment"):
    """ Plots the result of a simulation with an LCT SECIR model in a single Plot for specified comparments.
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
    plt.style.use('seaborn-v0_8-colorblind')

    # add results to plot
    for file in range(len(files)):
        # load data
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result
        total = data['Total'][:, :]
        if (total.shape[1] != 8):
            raise gd.DataError(
                "Expected a different number of compartments.")

        # plot result
        plt.plot(dates, total[:, compartment_idx], linewidth=1.5)
        h5file.close()

    plt.xlabel('Zeit', fontsize=14)
    plt.ylabel('Anzahl der Individuen in ' +
               secir_dict[compartment_idx], fontsize=14)
    plt.xlim(left=0)
    plt.legend(legendplot, fontsize=10)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # save result
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


def compare_all_compartments(files, legendplot, filename_plot="compare_compartments"):
    """ Creates a 4x2 Plot with one subplot per compartment and one line per result one wants to compare.

    @param[in] files: paths of the files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """

    plt.style.use('seaborn-v0_8-colorblind')
    fig, axs = plt.subplots(
        4, 2, sharex='all', num=filename_plot, tight_layout=True)

    # add results to plot
    for file in range(len(files)):
        # load data
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result
        total = data['Total'][:, :]
        if (total.shape[1] != 8):
            raise gd.DataError(
                "Expected a different number of compartments.")
        # plot data
        for i in range(8):
            axs[int(i/2), i % 2].plot(dates,
                                      total[:, i], label=legendplot[file], linewidth=1.5)

        h5file.close()

    # define some characteristics of the plot
    for i in range(8):
        axs[int(i/2), i % 2].set_title(secir_dict[i], fontsize=10)
        axs[int(i/2), i % 2].set_xlim(left=0)
        axs[int(i/2), i % 2].grid(True, linestyle='--')
        axs[int(i/2), i % 2].tick_params(axis='y', labelsize=8)
        axs[int(i/2), i % 2].tick_params(axis='x', labelsize=8)

    fig.supxlabel('Zeit', fontsize=12)
    fig.supylabel('Anzahl an Personen', fontsize=12)

    lines, labels = axs[0, 0].get_legend_handles_labels()
    lgd = fig.legend(lines, labels, ncol=len(legendplot),  loc='outside lower center',
                     fontsize=10, bbox_to_anchor=(0.5, - 0.06), bbox_transform=fig.transFigure)

    plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.1)
    plt.subplots_adjust(bottom=0.09)

    # save result
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    fig.savefig('Plots/'+filename_plot+'.png',
                bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=500)


def plot_new_infections(files, legendplot, filename_plot="compare_new_infections"):
    """ Single plot to compare the incidence of different results. Incidence means the number of people leaving the susceptible class per day.

    @param[in] files: paths of the files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
     @param[in] filename_plot: name to use as the file name for the saved plot.
    """
    # define plot
    plt.figure(filename_plot)
    plt.style.use('seaborn-v0_8-colorblind')

    # add results to plot
    for file in range(len(files)):
        # load data
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result
        total = data['Total'][:, :]
        if (total.shape[1] != 8):
            raise gd.DataError(
                "Expected a different number of compartments.")
        incidence = (total[:-1, 0]-total[1:, 0])/(dates[1:]-dates[:-1])

        # plot result
        plt.plot(dates[1:], incidence, linewidth=1.5)

        h5file.close()

    plt.xlabel('Zeit', fontsize=14)
    plt.ylabel('Neuansteckungen pro Tag', fontsize=14)
    plt.ylim(bottom=0)
    plt.xlim(left=0, right=dates[-1])
    plt.legend(legendplot, fontsize=10)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # save result
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


if __name__ == '__main__':
    # simulation results should be stored in folder "../data/simulation_lct"
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "data", "simulation_lct")

    # plot results
    """plot_lct_subcompartments(os.path.join(data_dir, "result_lct_subcompartments_fictional_3"), vec_subcompartments=[
                             1, 3, 3, 3, 3, 3, 1, 1], filename_plot="LCT_Subcompartments_3")
    compartment_index = [1, 2, 3]
    plot_single_result(os.path.join(data_dir, "result_lct_fictional_3"),
                       compartment_idx=compartment_index, filename_plot="LCT_Compartments_3")"""

    # compare lct and ode model
    # plot_new_infections([os.path.join(data_dir, "init", "lct_init_transitions"), os.path.join(data_dir, "init", "lct_init_mean"), os.path.join(data_dir, "init", "lct_init_first")],
    #                  legendplot=list(["Übergänge", "Mittelwerte", "Erstes Subkompartiment"]), filename_plot="compare_new_infections_initialization")

    plot_new_infections([os.path.join(data_dir, "dropR0", "fictional_drop_ode"),  os.path.join(data_dir, "dropR0", "fictional_drop_lct3"), os.path.join(data_dir, "dropR0", "fictional_drop_lct10"), os.path.join(data_dir, "dropR0", "fictional_drop_ide10")],
                        legendplot=list(["ODE", "LCT3", "LCT10", "IDE10"]), filename_plot="compare_drop")
    compare_all_compartments([os.path.join(data_dir, "dropR0", "fictional_drop_ode"), os.path.join(data_dir, "dropR0", "fictional_drop_lctvar"), os.path.join(data_dir, "dropR0", "fictional_drop_idevar")],
                             legendplot=list(["ODE", "LCTvar", "IDEvar"]), filename_plot="compare_drop_compartments")
    """plot_new_infections([os.path.join(data_dir, "result_ide_fictional"), os.path.join(data_dir, "result_lct_fictional_var")],
                        legendplot=list(["IDE", "LCTvar"]), filename_plot="IDEvsLCT_var_new_infections")
    compare_all_compartments([os.path.join(data_dir, "result_ide_fictional"), os.path.join(data_dir, "result_lct_fictional_var")],
                             legendplot=list(["IDE", "LCTvar"]), filename_plot="IDEvsLCT_var")
    compare_single_compartment([os.path.join(data_dir, "result_ide_fictional"), os.path.join(data_dir, "result_lct_fictional_var")],
                               legendplot=list(["IDE", "LCTvar"]), compartment_idx=2, filename_plot="carrier_IDEvsLCT_var")"""
