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

# import memilio.epidata.getDataIntoPandasDataFrame as gd

# Define compartments
secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
              5: 'ICU', 6: 'Recovered', 7: 'Dead'}


def compare_results(timestep, setting, legendplot, save=True):
    """ Creates a 4x2 Plot with one subplot per compartment and one line per result one wants to compare.
    @param[in] files: paths of the files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] save: if save is True, the plot is saved in a folder named Plots.
    """
    files = [os.path.join(data_dir, f"result_ode_dt={timestep}_setting{setting}"), os.path.join(
        data_dir, f"result_ide_dt={timestep}_setting{setting}")]

    fig, axs = plt.subplots(4, 2, sharex='all', num='Compare files')
    # helmholtzdarkblue, helmholtzclaim
    colors = [(0, 40/255, 100/255), (20/255, 200/255, 255/255)]
    linestyles = ['-', '--']
    # add results to plot
    for file in range(len(files)):
        # load data
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]

        if len(data['Total'][0]) == 8:
            # As there should be only one Group, total is the simulation result
            total = data['Total'][:, :]
        elif len(data['Total'][0]) == 10:
            # in ODE there are two compartments we don't use, throw these out
            total = data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]]

        dates = data['Time'][:]

        if (total.shape[1] != 8):
            raise gd.DataError(
                "Expected a different number of compartments.")
        # plot data
        for i in range(8):
            axs[int(i/2), i % 2].plot(dates,
                                      total[:, i], label=legendplot[file], color=colors[file], linestyle=linestyles[file])

        h5file.close()

    # define some characteristics of the plot
    for i in range(8):
        axs[int(i/2), i % 2].set_title(secir_dict[i], fontsize=8)
        # axs[int(i/2), i % 2].set_ylim(bottom=0)
        axs[int(i/2), i % 2].set_xlim(left=0)
        axs[int(i/2), i % 2].grid(True, linestyle='--')
        axs[int(i/2), i % 2].legend(fontsize=8)

    fig.supxlabel(' Time')
    fig.supylabel('Number of persons')
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=None, hspace=0.6)
    # plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.3)

    # save result
    if save:
        if not os.path.isdir('plots'):
            os.makedirs('plots')
        plt.savefig(f'plots/ide_ode_compare_dt={timestep}_setting{setting}.png',
                    bbox_inches='tight', dpi=500)
    # plt.show()


if __name__ == '__main__':
    # Path to simulation results
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "results")

    timestep = '1e-2'
    setting = 4
    # Plot comparison of ODE and IDE models
    compare_results(timestep, setting,
                    legendplot=list(["ODE", "IDE"]), save=True)
