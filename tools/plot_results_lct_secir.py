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

The data to be plotted should be stored in the current directory as .h5 files. 
Data could be generated eg by executing the file ./cpp/examples/lct_secir_compare.cpp.
"""

import h5py
import os
import pandas as pd
import matplotlib.pyplot as plt

import memilio.epidata.getDataIntoPandasDataFrame as gd


def get_subcompartments():
    """ Defines the used subcompartments for simulation of the LCT model.

    @return Vector with the number of subcompartments per compartment and a dictionary giving names 
        for the subcompartments.
    """
    vec_subcompartments = [1, 2, 3, 1, 1, 5, 1, 1]
    lct_secir_dict = {0: 'S', 1: 'E1', 2: 'E2', 3: 'C1', 4: 'C2',
                      5: 'C3', 6: 'I', 7: 'H', 8: 'U1', 9: 'U2', 10: 'U3', 11: 'U4', 12: 'U5', 13: 'R', 14: 'D'}
    return (vec_subcompartments, lct_secir_dict)


def plot_lct_subcompartments(file, save=True):
    """ Plots the result of a simulation with an LCT SECIR model in an 4x2 Plot. 
        In each subplot, one compartment with one line per subcompartment is plotted.

    @param[in] file: name of the file (without file extension .h5) with the simulation result of an lct model with subcompartments.
        The number of subcompartments should be fitting to the ones specified in get_subcompartments().
    @param[in] save: if save is True, the plot is saved in a folder named Plots.
    """
    (vec_subcompartments, lct_secir_dict) = get_subcompartments()

    # load data
    input_file = os.path.join(os.getcwd(), str(file))
    h5file = h5py.File(input_file + '.h5', 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]
    dates = data['Time'][:]
    # As there should be only one Group, total is the simulation result
    total = data['Total'][:, :]
    if (total.shape[1] != sum(vec_subcompartments)):
        raise gd.DataError("Expected a different number of subcompartments.")

    # Plot result.
    fig, axs = plt.subplots(
        4, 2, sharex='all', num='Subcompartments LCT SECIR model fictional scenario')
    for i in range(len(vec_subcompartments)):
        idx_start = sum(vec_subcompartments[0:i])
        axs[int(i/2), i % 2].plot(dates,
                                  total[:, idx_start:idx_start+vec_subcompartments[i]])
        axs[int(i/2), i % 2].grid(True, linestyle='--')
        legendplot = []
        for j in range(idx_start, idx_start+vec_subcompartments[i]):
            legendplot.append(lct_secir_dict[j])
        axs[int(i/2), i % 2].legend(legendplot, fontsize=10)
        axs[int(i/2), i % 2].set_ylim(bottom=0)
        axs[int(i/2), i % 2].set_xlim(left=0)
    fig.supxlabel('Time')
    fig.supylabel('Number of persons')

    # save plot
    if save:
        if not os.path.isdir('Plots'):
            os.makedirs('Plots')
        fig.savefig('Plots/lct_secir_fictional_subcompartments.png',
                    bbox_inches='tight', dpi=500)
    plt.show()
    h5file.close()


def plot_lct_result(file, compartment_idx=range(8), save=True):
    """ Plots the result of a simulation with an LCT SECIR model in a single Plot for specified comparments. 
        The result should consist of accumulated numbers for subcompartments.

    @param[in] file: name of the file (without file extension .h5) with the simulation result of an lct 
        model with accumulated numbers for subcompartments.
    @param[in] compartment_idx: indexes of the compartments that should be plot, eg in form of a list.
    @param[in] save: if save is True, the plot is saved in a folder named Plots.
    """
    # Define compartments
    secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                  5: 'ICU', 6: 'Recovered', 7: 'Death'}

    # load data
    input_file = os.path.join(os.getcwd(), str(file))
    h5file = h5py.File(input_file + '.h5', 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]
    dates = data['Time'][:]
    # As there should be only one Group, total is the simulation result
    total = data['Total'][:, compartment_idx]

    # plot result
    subcompartments = get_subcompartments()[0]
    heading_sub = ', '.join(str(sub) for sub in subcompartments)
    plt.figure(
        'SECIR LCT model, subcompartments:['+heading_sub+'], fictional scenario')
    plt.plot(dates, total, linewidth=1.0)

    legendplot = []
    for i in compartment_idx:
        legendplot.append(secir_dict[i])
    plt.legend(legendplot, fontsize=14)

    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Number of persons', fontsize=10)
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    plt.grid(True, linestyle='--')

    # save result
    if save:
        if not os.path.isdir('Plots'):
            os.makedirs('Plots')
        plt.savefig('Plots/lct_secir_fictional.png',
                    bbox_inches='tight', dpi=500)
    plt.show()
    h5file.close()


def compare_results(files, legendplot, save=True):
    """ Creates a 4x2 Plot with one subplot per compartment and one line per result one wants to compare.

    @param[in] files: names of the files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] save: if save is True, the plot is saved in a folder named Plots.
    """
    # Define compartments
    secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                  5: 'ICU', 6: 'Recovered', 7: 'Death'}

    fig, axs = plt.subplots(4, 2, sharex='all', num='Compare files')

    # add results to plot
    for file in range(len(files)):
        # load data
        input_file = os.path.join(os.getcwd(), str(files[file]))
        h5file = h5py.File(input_file + '.h5', 'r')

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
                "Expected a different number of subcompartments.")
        # plot data
        for i in range(8):
            axs[int(i/2), i % 2].plot(dates,
                                      total[:, i], label=legendplot[file])

        h5file.close()

    # define some characteristics of the plot
    for i in range(8):
        axs[int(i/2), i % 2].set_title(secir_dict[i], fontsize=8)
        axs[int(i/2), i % 2].set_ylim(bottom=0)
        axs[int(i/2), i % 2].set_xlim(left=0)
        axs[int(i/2), i % 2].grid(True, linestyle='--')
        axs[int(i/2), i % 2].legend(fontsize=8)

    fig.supxlabel('Time')
    fig.supylabel('Number of persons')

    # save result
    if save:
        if not os.path.isdir('Plots'):
            os.makedirs('Plots')
        fig.savefig('Plots/lct_secir_fictional_subcompartments.png',
                    bbox_inches='tight', dpi=500)
    plt.show()


if __name__ == '__main__':
    # plot lct result
    arr = list(range(1, 5))
    arr.append(7)
    plot_lct_result("result_lct", arr)
    # compare lct and ode model
    compare_results(["result_lct", "result_ode"],
                    legendplot=list(["LCT", "ODE"]), save=True)
