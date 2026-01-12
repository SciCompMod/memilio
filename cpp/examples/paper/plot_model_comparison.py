import h5py
import os
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from memilio.epidata import getDataIntoPandasDataFrame as gd

from plotting_settings import *

model_colors = [colors['Green'], colors['Rose'],
                colors['Blue'], colors['Teal']]

# Define compartments
secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
              5: 'ICU', 6: 'Recovered', 7: 'Dead'}


labels = ['ODE', 'LCT', 'IDE', 'ABM']
linestyles = ['-', '--', ':', '-']


def plot_one_compartment_impl(result_dir, ax, compartment_index,  percentiles_init=["05", "50", "95"],  plot_init=True, linewidth=3, num_age_groups=6):

    percentiles = copy.deepcopy(percentiles_init)
    start_time = 0
    if (plot_init):
        df = pd.read_csv(result_dir + f"comps.csv")
        values = df.iloc[:, 1 + compartment_index]
        for age in range(1, num_age_groups):
            values += df.iloc[:, 1 + compartment_index + age * len(secir_dict)]
        ax.plot(df["Time"], values, label=f"Initial ABM run",
                color="grey", linestyle=linestyles[3], linewidth=linewidth)
    while len(percentiles) > 0:
        if len(percentiles) == 1:
            df = pd.read_csv(result_dir + f"ABM_p{percentiles[0]}.csv")
            values = df.iloc[:, 1 + compartment_index]
            for age in range(1, num_age_groups):
                values += df.iloc[:, 1 +
                                  compartment_index + age * len(secir_dict)]
            ax.plot(df["Time"] + start_time, values, label=f"p{percentiles[0]} ABM",
                    color=model_colors[3], linestyle=linestyles[3], linewidth=linewidth)
            del percentiles[0]
        else:
            df_low = pd.read_csv(result_dir + f"ABM_p{percentiles[0]}.csv")
            df_high = pd.read_csv(result_dir + f"ABM_p{percentiles[-1]}.csv")
            values_low = df_low.iloc[:, 1 + compartment_index]
            values_high = df_high.iloc[:, 1 + compartment_index]
            for age in range(1, num_age_groups):
                values_low += df_low.iloc[:, 1 +
                                          compartment_index + age * len(secir_dict)]
                values_high += df_high.iloc[:, 1 +
                                            compartment_index + age * len(secir_dict)]
            ax.plot(df_low["Time"] + start_time, values_low,
                    color=colors['Teal'], alpha=0.3)
            ax.plot(df_high["Time"] + start_time, values_high,
                    color=colors['Teal'], alpha=0.3)
            ax.fill_between(df_low["Time"] + start_time, values_low, values_high, alpha=0.3,
                            color=model_colors[3], label=f"p{int(percentiles[0])} - p{int(percentiles[-1])} ABM")
            del percentiles[0]
            del percentiles[-1]

    # Add results to plot for IDE, LCT and ODE.
    filenames = [os.path.join("ode"),
                 os.path.join("lct"),
                 os.path.join("ide_no_init")]
    dates = []
    for file in range(len(filenames)):
        # Load data.
        h5file = h5py.File(result_dir + str(filenames[file]) + '.h5', 'r')

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

            ax.plot(dates,
                    total[:, compartment_index_ode], label=labels[file],  color=model_colors[file], linestyle=linestyles[file], linewidth=linewidth)

        # LCT and IDE
        elif file == 1 or file == 2:
            ax.plot(dates,
                    total[:, compartment_index], label=labels[file],  color=model_colors[file], linestyle=linestyles[file], linewidth=linewidth)

            h5file.close()

    return dates


def plot_model_comparison_one_compartment(result_dir, compartment_index, save_dir="", percentiles=["05", "50", "95"], plot_init=True, num_age_groups=6):
    """
    Plots the result of the simulation.

    @param[in] files Expects list of two files with ODE and IDE simulation results for compartments, respectively, in
        this order.
    @param[in] fileending Determines file ending of saved plot. Default is an empty string leading to no further
        specification.
    @param[in] save_dir Directory where plot will be stored. Default is an empty string leading to the plot not being
        saved.
    """

    # Define plot.
    fig, axs = plt.subplots(1, 1, sharex='all', figsize=(4, 3))
    linewidth = 3

    plot_one_compartment_impl(result_dir, axs, compartment_index,
                              percentiles, plot_init, linewidth, num_age_groups)

    axs.ticklabel_format(axis='y',
                         style='sci', scilimits=(0, 0))

    axs.set_xlabel('Simulation time [days]', labelpad=10)
    axs.set_ylabel(f'{secir_dict[compartment_index]} [#]')

    plt.tight_layout()

    # Save legend in seperate figure
    handles, labels = axs.get_legend_handles_labels()
    fig_leg = plt.figure(figsize=(5, 4))
    fig_leg.legend(
        handles,
        labels,
        loc="center",
        ncol=4,
        frameon=False
    )
    fig_leg.patch.set_visible(False)

    set_fontsize()

    # Save result.
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        fig.savefig(save_dir + f"model_comparison_{secir_dict[compartment_index]}.png",
                    bbox_inches='tight', dpi=dpi)
        fig_leg.savefig(save_dir + f"model_comparison_legend.png",
                        bbox_inches='tight', dpi=dpi)

    plt.close()


def plot_model_comparison_all_compartments(result_dir, percentiles, num_age_groups=6, plot_init=False, save_dir=""):
    """
    Plots the result of the simulation.

    @param[in] files Expects list of two files with ODE and IDE simulation results for compartments, respectively, in
        this order.
    @param[in] fileending Determines file ending of saved plot. Default is an empty string leading to no further
        specification.
    @param[in] save_dir Directory where plot will be stored. Default is an empty string leading to the plot not being
        saved.
    """

    # Define plot.
    fig, axs = plt.subplots(4, 2, sharex='all')
    num_plots = 8
    linewidth = 1

    for i in range(num_plots):

        plot_one_compartment_impl(result_dir, axs[int(i/2), i % 2], i,
                                  percentiles, plot_init, linewidth, num_age_groups)

    # Define some characteristics of the plot
    for i in range(num_plots):
        axs[int(i/2), i % 2].set_title(secir_dict[i], fontsize=8)
        # axs[int(i/2), i % 2].set_xlim(left=0, right=dates[-1])
        axs[int(i/2), i % 2].grid(True, linestyle='--', alpha=0.3)
        axs[int(i/2), i % 2].ticklabel_format(axis='y',
                                              style='sci', scilimits=(0, 0))
        axs[int(i/2), i % 2].tick_params(labelsize=8)

    # Create legend with all unique labels from subplots
    # Use dummy handle so that legend is sorted in a nicer way.
    if plot_init:
        all_handles = []
        all_labels = []
    else:
        dummy_handle = plt.Line2D([], [], color='none')
        all_handles = [dummy_handle]
        all_labels = [""]
    for ax in axs.flat:
        handles, labels_temp = ax.get_legend_handles_labels()
        for handle, label in zip(handles, labels_temp):
            if label not in all_labels:
                all_handles.append(handle)
                all_labels.append(label)

    fig.legend(all_handles, all_labels, loc="lower right",
               fancybox=False, shadow=False, ncol=2, fontsize=10)

    fig.supxlabel('Simulation time [days]', fontsize=12)
    fig.supylabel('Number of individuals', fontsize=12)
    plt.subplots_adjust(left=None, bottom=0.2, right=None,
                        top=None, wspace=None, hspace=0.6)

    set_fontsize()

    # Save result.
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f"all_compartments" + ".png",
                    bbox_inches='tight', dpi=500)

    plt.close()


def plot_model_comparison_all_scenarios_Exposed_ICU_Deaths(result_dirs_per_scenario, percentiles=["05", "50", "95"], num_age_groups=6, plot_init=True, save_dir=""):
    """
    Plots model comparison results across multiple scenarios for Exposed, ICU, and Dead compartments.

    @param[in] scenario_dirs List of directories containing simulation results for each scenario
    @param[in] scenario_names List of names for each scenario (for column labels)
    @param[in] percentiles List of percentiles to plot for ABM results
    @param[in] num_age_groups Number of age groups in the simulation
    @param[in] plot_init Whether to plot initial/deterministic ABM results
    @param[in] save_dir Directory where plot will be stored. Default is an empty string leading to the plot not being saved.
    """

    # Compartments to plot: Exposed (1), ICU (5), Dead (7)
    compartments_to_plot = [1, 5, 7]
    compartment_names = ['Exposed', 'ICU', 'Dead']

    # Define plot with 3 rows (compartments) and 4 columns (scenarios)
    fig, axs = plt.subplots(3, 4,
                            figsize=(16, 10))
    linewidth = 3

    # Iterate over scenarios (columns)
    for scenario_idx in range(4):
        # Iterate over compartments (rows)
        for comp_row_idx, compartment_index in enumerate(compartments_to_plot):
            ax = axs[comp_row_idx, scenario_idx]

            dates = plot_one_compartment_impl(result_dirs_per_scenario[scenario_idx], ax, compartment_index,
                                              percentiles, plot_init, linewidth, num_age_groups)

            # Set labels and formatting
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            ticks = np.arange(dates[0], dates[-1] + 1, 20)
            if 0 not in ticks:
                ticks = np.sort(np.unique(np.concatenate((ticks, [0]))))
            ax.set_xticks(ticks)
            # ax.grid(True, linestyle='--', alpha=0.3)

            # Row labels (compartment names on the left)
            if scenario_idx == 0:
                ax.set_ylabel(
                    f'{compartment_names[comp_row_idx]} [#]')

            # # X-axis labels at the bottom
            # if comp_row_idx == len(compartments_to_plot) - 1:
            #     ax.set_xlabel('Simulation time [days]')

    # X-axis label for all plots at the bottom
    fig.supxlabel('Simulation time [days]')

    # align y-labels of the left column so they share the same horizontal position
    fig.align_ylabels(axs[:, 0])

    # ensure enough left margin so labels are not clipped
    plt.subplots_adjust(left=0.12)

    # Create legend with all unique labels from subplots
    all_handles = []
    all_labels = []
    for ax_row in axs:
        for ax in ax_row if len(axs.shape) > 1 else [ax_row]:
            handles, labels_temp = ax.get_legend_handles_labels()
            for handle, label in zip(handles, labels_temp):
                if label not in all_labels:
                    all_handles.append(handle)
                    all_labels.append(label)

    fig.legend(all_handles, all_labels, loc="lower center", bbox_to_anchor=(0.5, -0.05),
               fancybox=False, shadow=False, ncol=len(all_labels))

    set_fontsize()
    plt.tight_layout()  # rect=[0, 0.08, 1, 1]

    # Save result
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        fig.savefig(save_dir + f"model_comparison_all_scenarios_E_U_D.png",
                    bbox_inches='tight', dpi=dpi)

    plt.close()


def plot_ABM_results_one_compartments(file, compartment_index, percentiles, num_age_groups, num_comps, save_dir="", plot_init=True):
    figsize = (5, 3.5)
    panel = [0.2, 0.18, 0.78, 0.8]
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(panel)
    secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                  5: 'ICU', 6: 'Recovered', 7: 'Dead'}
    start_time = 0
    if (plot_init):
        df = pd.read_csv(file + f"comps.csv")
        values = df.iloc[:, 1 + compartment_index]
        for age in range(1, num_age_groups):
            values += df.iloc[:, 1 + compartment_index + age * num_comps]
        ax.plot(df["Time"], values, color="grey")
    while len(percentiles) > 0:
        if len(percentiles) == 1:
            df = pd.read_csv(file + f"ABM_p{percentiles[0]}.csv")
            values = df.iloc[:, 1 + compartment_index]
            for age in range(1, num_age_groups):
                values += df.iloc[:, 1 + compartment_index + age * num_comps]
            ax.plot(df["Time"] + start_time, values, color=colors['Teal'])
            del percentiles[0]
        else:
            df_low = pd.read_csv(file + f"ABM_p{percentiles[0]}.csv")
            df_high = pd.read_csv(file + f"ABM_p{percentiles[-1]}.csv")
            values_low = df_low.iloc[:, 1 + compartment_index]
            values_high = df_high.iloc[:, 1 + compartment_index]
            for age in range(1, num_age_groups):
                values_low += df_low.iloc[:, 1 +
                                          compartment_index + age * num_comps]
                values_high += df_high.iloc[:, 1 +
                                            compartment_index + age * num_comps]
            # ax.plot(df_low["Time"] + start_time, values_low, color=colors['Teal'], alpha=0.5)
            # ax.plot(df_high["Time"]+ start_time, values_high, color=colors['Teal'], alpha=0.5)
            ax.fill_between(df_low["Time"] + start_time, values_low,
                            values_high, alpha=0.3, color=colors['Teal'])
            del percentiles[0]
            del percentiles[-1]
    ax.set_ylabel(f"{secir_dict[compartment_index]} [#]")
    ax.set_xlabel("Time [days]")
    fig.savefig(save_dir + f"ABM_results_{compartment_index}.png", dpi=dpi)


def plot_single_ABM_run(file, seeds, num_age_groups, num_comps, save_dir=""):
    figsize = (5, 3.5)
    panel = [0.18, 0.18, 0.8, 0.8]
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(panel)
    secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                  5: 'ICU', 6: 'Recovered', 7: 'Dead'}
    df = pd.read_csv(file + f"comps_{seeds}.csv")
    for comp in range(1, num_comps):
        values = df.iloc[:, 1 + comp]
        for age in range(1, num_age_groups):
            values += df.iloc[:, 1 + comp + age * num_comps]
        ax.plot(df["Time"], values, label=secir_dict[comp])
    ax.set_ylabel(f"Individuals [#]")
    ax.set_xlabel("Time [days]")
    ax.legend()
    fig.savefig(save_dir + f"ABM_results_{seeds}.png", dpi=dpi)


if __name__ == '__main__':

    root_dir = os.path.join(os.path.dirname(
        __file__),  f"simulation_results/compare_abm_ide_lct_ode")

    seed = "Seed1"
    num_runs = 100

    exponential_scenario = False
    one_location = False

    if exponential_scenario:
        distribution = "exponential"
    else:
        distribution = "different_dists"

    if one_location:
        location = "one_location"
    else:
        location = "multiple_locations"

    # Path where simulation results are stored.
    result_dir = f"{root_dir}/{num_runs}_runs/{distribution}/{location}/{seed}/"
    # Path where plots will be stored.
    plot_dir = f"{root_dir}/{num_runs}_runs/{distribution}/{location}/{seed}/"
    os.makedirs(plot_dir, exist_ok=True)

    num_age_groups = 6
    plot_init = True
    percentiles = ["05", "50", "95"]

    # for i in range(8):
    #     plot_model_comparison_one_compartment(result_dir,
    #                                           i,
    #                                           plot_dir, percentiles, plot_init, num_age_groups)

    # plot_model_comparison_all_compartments(
    #     result_dir, percentiles, num_age_groups, plot_init, plot_dir)

    result_dirs_per_scenario = [f"{root_dir}/{num_runs}_runs/exponential/one_location/{seed}/",
                                f"{root_dir}/{num_runs}_runs/exponential/multiple_locations/{seed}/",
                                f"{root_dir}/{num_runs}_runs/different_dists/one_location/{seed}/",
                                f"{root_dir}/{num_runs}_runs/different_dists/multiple_locations/{seed}/"]

    plot_model_comparison_all_scenarios_Exposed_ICU_Deaths(
        result_dirs_per_scenario, percentiles, num_age_groups, plot_init, plot_dir)

    # plot_ABM_results_one_compartments(
    #     result_dir, 0, ["05", "50", "95"], 6, 8, plot_dir)
    # plot_ABM_results_one_compartments(
    #     result_dir, 1, ["05", "50", "95"], 6, 8, plot_dir)
    # plot_ABM_results_one_compartments(
    #     result_dir, 3, ["05", "50", "95"], 6, 8, plot_dir)
    # plot_ABM_results_one_compartments(
    #     result_dir, 5, ["05", "50", "95"], 6, 8, plot_dir)
    # plot_ABM_results_one_compartments(
    #     result_dir, 6, ["05", "50", "95"], 6, 8, plot_dir)
    # plot_ABM_results_one_compartments(
    #     result_dir, 7, ["05", "50", "95"], 6, 8, plot_dir)

    # plot_single_ABM_run(
    #     result_dir, "518254265_179139074_1937166324_3882038653_1776323086_1261445406", 6, 8)
