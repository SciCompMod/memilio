import h5py
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from memilio.epidata import getDataIntoPandasDataFrame as gd

from plotting_settings import *

model_colors = [colors['Light green'], colors['Rose'],
                colors['Light blue'], colors['Teal']]


def plot_model_comparison_all_compartments(files, save_dir=""):
    """
    Plots the result of the simulation.

    @param[in] files Expects list of two files with ODE and IDE simulation results for compartments, respectively, in
        this order.
    @param[in] fileending Determines file ending of saved plot. Default is an empty string leading to no further
        specification.
    @param[in] save_dir Directory where plot will be stored. Default is an empty string leading to the plot not being
        saved.
    """
    # Define compartments
    secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                  5: 'ICU', 6: 'Recovered', 7: 'Dead'}

    # Define plot.
    fig, axs = plt.subplots(4, 2, sharex='all')
    num_plots = 8
    linewidth = 1
    labels = ['ODE', 'LCT', 'IDE', 'ABM']

    linestyles = ['-', '--', ':', '-']
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
        # ODE
        if file == 0:
            # Define indices of compartmens in ODE results separately because results include two compartments for
            # confirmed cases that we do not consider.
            ode_secir_indices = [0, 1, 2, 4, 6, 7, 8, 9]
            for counter, ode_index in enumerate(ode_secir_indices):
                axs[int(counter/2), counter % 2].plot(dates,
                                                      total[:, ode_index], label=labels[file], color=model_colors[file], linestyle=linestyles[file], linewidth=linewidth)

        # LCT and IDE
        elif file == 1 or file == 2:
            for i in range(num_plots):
                axs[int(i/2), i % 2].plot(dates,
                                          total[:, i], label=labels[file], color=model_colors[file], linestyle=linestyles[file], linewidth=linewidth)

        # ABM
        # TODO
        elif file == 3:
            pass

        h5file.close()

    # Define some characteristics of the plot
    for i in range(num_plots):
        axs[int(i/2), i % 2].set_title(secir_dict[i], fontsize=8)
        axs[int(i/2), i % 2].set_xlim(left=dates[0], right=dates[-1])
        axs[int(i/2), i % 2].grid(True, linestyle='--', alpha=0.5)
        axs[int(i/2), i % 2].ticklabel_format(axis='y',
                                              style='sci', scilimits=(0, 0))

    fig.legend(labels, bbox_to_anchor=(0.1, -0.73, 0.8, 0.8),
               fancybox=False, shadow=False, ncol=1)

    fig.supxlabel('Simulation time [days]')
    fig.supylabel('Number of individuals')
    # plt.subplots_adjust(left=None, bottom=None, right=None,
    #                     top=None, wspace=None, hspace=0.6)

    # plt.tight_layout()

    # Save result.
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f"all_compartments.png",
                    bbox_inches='tight', dpi=500)

    plt.clf()


def plot_model_comparison_one_compartment(files, compartment_index, exponential_scenario=False, save_dir=""):
    """
    Plots the result of the simulation.

    @param[in] files Expects list of two files with ODE and IDE simulation results for compartments, respectively, in
        this order.
    @param[in] fileending Determines file ending of saved plot. Default is an empty string leading to no further
        specification.
    @param[in] save_dir Directory where plot will be stored. Default is an empty string leading to the plot not being
        saved.
    """
    # Define compartments
    secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                  5: 'ICU', 6: 'Recovered', 7: 'Dead'}

    # Define plot.
    fig, axs = plt.subplots(1, 1, sharex='all', figsize=(5, 4))
    num_plots = 1
    linewidth = 3
    labels = ['ODE', 'LCT', 'IDE', 'ABM']

    # colors = ["C0", "limegreen"]
    linestyles = ['-', '--', ':', '-']
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
        # ODE
        if file == 0:
            # Define indices of compartmens in ODE results separately because results include two compartments for
            # confirmed cases that we do not consider.
            ode_secir_indices = [0, 1, 2, 4, 6, 7, 8, 9]
            compartment_index_ode = ode_secir_indices[compartment_index]

            axs.plot(dates,
                     total[:, compartment_index_ode], label=labels[file],  color=model_colors[file], linestyle=linestyles[file], linewidth=linewidth)

        # LCT and IDE
        elif file == 1 or file == 2:
            axs.plot(dates,
                     total[:, compartment_index], label=labels[file],  color=model_colors[file], linestyle=linestyles[file], linewidth=linewidth)

        # ABM
        # TODO
        elif file == 3:
            axs.plot(0, 0, label=labels[file],  color=model_colors[file],
                     linestyle=linestyles[file], linewidth=linewidth)

        h5file.close()

    # Define some characteristics of the plot
    if exponential_scenario:
        axs.set_title("Scenario 1")
    else:
        axs.set_title("Scenario 2")

    axs.set_xlim(left=dates[0], right=dates[-1])
    # axs.grid(True, linestyle='--', alpha=0.5)
    axs.ticklabel_format(axis='y',
                         style='sci', scilimits=(0, 0))

    axs.set_xlabel('Simulation time [days]', labelpad=10)
    if compartment_index == 3:
        axs.set_ylabel('Mildly symptomatic individuals', labelpad=10)
    if compartment_index == 5:
        axs.set_ylabel('Individuals on ICU', labelpad=10)

    fig.legend(labels, bbox_to_anchor=(0.55, 0.85),
               fancybox=False, shadow=False, ncol=2)

    # plt.subplots_adjust(left=0.5, bottom=None, right=0.7,
    #                     top=None, wspace=None, hspace=None)

    set_fontsize()

    plt.tight_layout()

    # Save result.
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f"model_comparison_{secir_dict[compartment_index]}.png",
                    bbox_inches='tight', dpi=dpi)

    plt.clf()


if __name__ == '__main__':

    exponential_scenario = False

    if exponential_scenario:
        # Path where simulation results are stored.
        result_dir = os.path.join(os.path.dirname(
            __file__), "../../..", "simulation_results/compare_abm_ide_lct_ode/exponential/")
        # Path where plots will be stored.
        plot_dir = os.path.join(os.path.dirname(
            __file__), "../../..", "plots/compare_abm_ide_lct_ode/exponential/")

    else:
        # Path where simulation results are stored.
        result_dir = os.path.join(os.path.dirname(
            __file__), "../../..", "simulation_results/compare_abm_ide_lct_ode/different_dists/")
        # Path where plots will be stored.
        plot_dir = os.path.join(os.path.dirname(
            __file__), "../../..", "plots/compare_abm_ide_lct_ode/different_dists/")

    plot_model_comparison_all_compartments([os.path.join(result_dir, f"ode"),
                                            os.path.join(result_dir, f"lct"),
                                            os.path.join(result_dir, f"ide"),
                                            os.path.join(result_dir, f"ide")],  # last file needs to be changed to ABM later on
                                           save_dir=plot_dir)

    plot_model_comparison_one_compartment([os.path.join(result_dir, f"ode"),
                                           os.path.join(result_dir, f"lct"),
                                           os.path.join(result_dir, f"ide"),
                                           os.path.join(result_dir, f"ide")],  # last file needs to be changed to ABM later on
                                          5,
                                          exponential_scenario,
                                          save_dir=plot_dir)
