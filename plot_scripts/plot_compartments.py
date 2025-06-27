#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Anna Wendler
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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from memilio.epidata import getDataIntoPandasDataFrame as gd


def plot_susceptibles(files, fileending, save_dir=""):
    """
    Plots simulation results of Susceptibles.

    @param[in] files Expects list of two files with ODE and IDE simulation results for compartments, respectively, in
        this order.
    @param[in] fileending Determines file ending of saved plot. Default is an empty string leading to no further
        specification.
    @param[in] save_dir Directory where plot will be stored. Default is an empty string leading to the plot not being
        saved.
    """
    # Define compartments
    secir_dict = {0: 'Susceptible',  1: 'Infected', 2: 'Recovered'}

    # Define plot.
    num_plots = 1
    fig, axs = plt.subplots(1, num_plots, sharex='all', num='Compare files')

    colors = ["C0", "limegreen"]
    linestyles = ['-', '--']
    linewidth = 1
    labels = ["Groundtruth", "IDE result"]

    # Add results to plot.
    for file in range(len(files)):
        # Load data.
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]

        # As there should be only one Group, total is the simulation result.
        total = data['Total'][:, :]

        dates = data['Time'][:]

        # Plot data.
        for i in range(num_plots):
            axs.plot(dates,
                     total[:, i], label=labels[file],  linestyle=linestyles[file], color=colors[file], linewidth=linewidth)

        h5file.close()

    # Define some characteristics of the plot
    for i in range(num_plots):
        axs.set_title(secir_dict[i], fontsize=8)
        axs.set_xlim(left=dates[0], right=dates[-1])
        axs.grid(True, linestyle='--', alpha=0.5)
        axs.ticklabel_format(axis='y',
                             style='sci', scilimits=(0, 0))

    fig.legend(labels, bbox_to_anchor=(0.1, -0.73, 0.8, 0.8),
               fancybox=False, shadow=False, ncol=1)

    fig.supxlabel('Simulation time [days]')
    fig.supylabel('Number of individuals')
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=None, hspace=0.6)

    plt.tight_layout()

    # Save result.
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f"compare_compartments_{fileending}.png",
                    bbox_inches='tight', dpi=500)

    plt.clf()


if __name__ == '__main__':

    dir_name = "exponential_paper_example"

    # Path where simulation results (generated with ide_changepoints.cpp) are stored.
    result_dir = os.path.join(os.path.dirname(
        __file__),  f"../simulation_results/{dir_name}/")
    # Path where plots will be stored.
    plot_dir = os.path.join(os.path.dirname(
        __file__),  f"../plots/{dir_name}/")

    gregory_orders = ["1", "2", "3"]
    ide_exponent = "0"
    groundtruth_exponent = "5"

    # for gregory_order in gregory_orders:
    #     plot_susceptibles([os.path.join(result_dir, f"result_ide_dt=1e-4_gregoryorder=3_finitedifforder=1"),
    #                        os.path.join(result_dir, f"result_ide_dt=1e-{dt_exp}_gregoryorder={gregory_order}_finitedifforder=1")],
    #                       fileending=f"dt=1e-{dt_exp}_gregory={gregory_order}", save_dir=plot_dir)

    # for gregory_order in gregory_orders:
    #     plot_susceptibles([os.path.join(result_dir, f"result_ide_dt=1e-{groundtruth_exponent}_gregoryorder={3}"),
    #                        os.path.join(result_dir, f"result_ide_dt=1e-{ide_exponent}_gregoryorder={gregory_order}")],
    #                       fileending=f"dt=1e-{ide_exponent}_gregory={gregory_order}", save_dir=plot_dir)
    # for ODE

    groundtruth_exponent = "6"
    for gregory_order in gregory_orders:
        plot_susceptibles([os.path.join(result_dir, f"result_ode_dt=1e-{groundtruth_exponent}"),
                           os.path.join(result_dir, f"result_ide_dt=1e-{ide_exponent}_gregoryorder={gregory_order}")],
                          fileending=f"dt=1e-{ide_exponent}_gregory={gregory_order}", save_dir=plot_dir)
