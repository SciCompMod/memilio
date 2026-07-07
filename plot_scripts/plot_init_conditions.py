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
    num_plots = 3
    fig, axs = plt.subplots(1, num_plots, sharex='all', num='Compare files')

    colors = ["C0", "limegreen", "Orange"]
    linestyles = ['-', '--', ':']
    linewidth = 1
    labels = ["Groundtruth", "Detailed", "Simple"]

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
            axs[i].plot(dates,
                        total[:, i], label=labels[file],  linestyle=linestyles[file], color=colors[file], linewidth=linewidth)

        h5file.close()

    # Define some characteristics of the plot
    for i in range(num_plots):
        axs[i].set_title(secir_dict[i], fontsize=8)
        axs[i].set_xlim(left=0, right=dates[-1])
        axs[i].grid(True, linestyle='--', alpha=0.5)
        axs[i].ticklabel_format(axis='y',
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

    plt.close()


def plot_flow_S_to_I(files, fileending, save_dir=""):
    # Define compartments
    flow_dict = {0: 'SusceptibleToExposed'}

    # Define plot.
    num_plots = 1
    fig, axs = plt.subplots(1, num_plots, sharex='all', num='Compare files')

    colors = ["C0", "limegreen", "Orange"]
    linestyles = ['-', '--', ':']
    linewidth = 1
    labels = ["Groundtruth", "Detailed", "Simple"]

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
        # if file != 1:
        if file == 0:
            axs.plot(dates,
                     total[:, 0], label=labels[file],  linestyle=linestyles[file], color=colors[file], linewidth=linewidth)

        h5file.close()

    # Define some characteristics of the plot

    axs.set_title(flow_dict[0], fontsize=8)
    # axs.set_xlim(left=0, right=dates[-1])
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
        plt.savefig(save_dir + f"compare_flow_{fileending}.png",
                    bbox_inches='tight', dpi=500)

    plt.close()


def plot_infectionage_distribution(files, fileending, save_dir=""):
    # Define compartments
    flow_dict = {0: 'Infected per infection age'}

    # Define plot.
    num_plots = 1
    fig, axs = plt.subplots(1, num_plots, sharex='all', num='Compare files')

    colors = ["C0", "limegreen", "Orange"]
    linestyles = ['-', '--', ':']
    linewidth = 1
    labels = ["Groundtruth", "Detailed", "Simple"]

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
        # if file != 1:
        s = 1
        if file == 2:
            s = 10

        axs.scatter(dates,
                    total[:], label=labels[file],  linestyle=linestyles[file], color=colors[file], linewidth=linewidth, s=s)

        h5file.close()

    # Define some characteristics of the plot

    axs.set_title(flow_dict[0], fontsize=8)
    # axs.set_xlim(left=0, right=dates[-1])
    axs.grid(True, linestyle='--', alpha=0.5)
    axs.ticklabel_format(axis='y',
                         style='sci', scilimits=(0, 0))

    fig.legend(labels, bbox_to_anchor=(0.1, -0.73, 0.8, 0.8),
               fancybox=False, shadow=False, ncol=1)

    fig.supxlabel('Infection age [days]')
    fig.supylabel('Number of individuals')
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=None, hspace=0.6)

    plt.tight_layout()

    # Save result.
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f"infected_per_infection_age_{fileending}.png",
                    bbox_inches='tight', dpi=500)

    plt.close()


def subfolders_scandir(path):
    # path = os.path.dirname(path)
    print(path)
    with os.scandir(path) as it:
        return [entry.name for entry in it if entry.is_dir()]


if __name__ == '__main__':

    root_dir = os.path.join(os.path.dirname(
        __file__), "../simulation_results")

    main_dir = "2026-07-06/diff_groundtruth_init"

    relevant_dir = os.path.join(root_dir, main_dir)
    sub_dirs = subfolders_scandir(relevant_dir)
    # sub_dirs = [
    #     "nonconst_contacts_t0=0_tinit=20_tmax=30_scalingtime=25_damping=0.5"]

    gregory_order = 3

    for sub_dir in sub_dirs:
        print(main_dir + "/" + sub_dir)

        # Path where simulation results are stored.
        result_dir = os.path.join(os.path.dirname(
            __file__),  f"../simulation_results/{main_dir}/{sub_dir}")
        # Path where plots will be stored.
        plot_dir = os.path.join(os.path.dirname(
            __file__),  f"../plots/{main_dir}/{sub_dir}/")

        files = os.listdir(result_dir)

        ide_exponent = 2
        if f'detailed_dt=1e-{ide_exponent}_gregoryorder=3.h5' in files:
            print(ide_exponent)

            plot_susceptibles([os.path.join(result_dir, f"groundtruth_dt=1e-{ide_exponent}_gregoryorder={gregory_order}"),
                               os.path.join(
                result_dir, f"detailed_dt=1e-{ide_exponent}_gregoryorder={gregory_order}"),
                os.path.join(result_dir, f"simple_dt=1e-{ide_exponent}_gregoryorder={gregory_order}")],
                fileending=f"dt=1e-{ide_exponent}", save_dir=plot_dir)

            plot_flow_S_to_I([os.path.join(result_dir, f"groundtruth_flows_dt=1e-{ide_exponent}_gregoryorder={gregory_order}"),
                              os.path.join(
                result_dir, f"detailed_flows_dt=1e-{ide_exponent}_gregoryorder={gregory_order}"),
                os.path.join(result_dir, f"simple_flows_dt=1e-{ide_exponent}_gregoryorder={gregory_order}")],
                fileending=f"dt=1e-{ide_exponent}", save_dir=plot_dir)

            plot_infectionage_distribution([os.path.join(result_dir, f"groundtruth_infectionagedistribution_dt=1e-{ide_exponent}_gregoryorder={gregory_order}"),
                                            os.path.join(
                result_dir, f"detailed_infectionagedistribution_dt=1e-{ide_exponent}_gregoryorder={gregory_order}"),
                os.path.join(result_dir, f"simple_infectionagedistribution_dt=1e-{ide_exponent}_gregoryorder={gregory_order}")],
                fileending=f"dt=1e-{ide_exponent}", save_dir=plot_dir)
