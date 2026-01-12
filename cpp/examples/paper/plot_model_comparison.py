import h5py
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from memilio.epidata import getDataIntoPandasDataFrame as gd

from plotting_settings import *

model_colors = [colors['Green'], colors['Rose'],
                colors['Blue'], colors['Teal']]


def plot_model_comparison_all_compartments(result_dir,  percentiles, num_age_groups=6, file_suffix="", plot_init=False,  save_dir=""):
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

    # ABM
    start_time = 0
    if (plot_init):
        df = pd.read_csv(result_dir + f"comps.csv")
        for i in range(num_plots):
            values = df.iloc[:, 1 + i]
            # Start with index 1 as first age group is already contained in values
            for age in range(1, num_age_groups):
                values += df.iloc[:, 1 + i + age * len(secir_dict)]
            axs[int(i/2), i % 2].plot(df["Time"], values,
                                      color=model_colors[3], linestyle=linestyles[3], linewidth=linewidth)
    while len(percentiles) > 0:
        if len(percentiles) == 1:
            df = pd.read_csv(result_dir + f"ABM_p{percentiles[0]}.csv")
            for i in range(num_plots):
                values = df.iloc[:, 1 + i]
                # Start with index 1 as first age group is already contained in values
                for age in range(1, num_age_groups):
                    values += df.iloc[:, 1 +
                                      i + age * len(secir_dict)]
                axs[int(i/2), i % 2].plot(df["Time"] +
                                          start_time, values, label=f"p{percentiles[0]} ABM", color=model_colors[3], linestyle=linestyles[3], linewidth=linewidth)
            del percentiles[0]
        else:
            df_low = pd.read_csv(result_dir + f"ABM_p{percentiles[0]}.csv")
            df_high = pd.read_csv(result_dir + f"ABM_p{percentiles[-1]}.csv")
            for i in range(num_plots):
                values_low = df_low.iloc[:, 1 + i]
                values_high = df_high.iloc[:, 1 + i]
                # Start with index 1 as first age group is already contained in values_low and values_high
                for age in range(1, num_age_groups):
                    values_low += df_low.iloc[:, 1 +
                                              i + age * len(secir_dict)]
                    values_high += df_high.iloc[:, 1 +
                                                i + age * len(secir_dict)]
                axs[int(i/2), i % 2].plot(df_low["Time"] + start_time, values_low,
                                          color=colors['Teal'], alpha=0.3)
                axs[int(i/2), i % 2].plot(df_high["Time"] + start_time,
                                          values_high, color=colors['Teal'], alpha=0.3)
                axs[int(i/2), i % 2].fill_between(df_low["Time"] + start_time, values_low,
                                                  values_high, alpha=0.3, color=model_colors[3], label=f"{int(percentiles[-1])-int(percentiles[0])}% quantile ABM")
            del percentiles[0]
            del percentiles[-1]

    # Add results to plot for IDE, LCT and ODE.
    files = ["ode", "lct", "ide_no_init"]
    files = [file + file_suffix for file in files]
    dates = []
    for file in range(len(files)):

        # Load data.
        h5file = h5py.File(result_dir + str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")

        data = h5file[list(h5file.keys())[0]]
        # print(data)

        total = data['Total'][:, :]
        # print(total)

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

        h5file.close()

    # Define some characteristics of the plot
    for i in range(num_plots):
        axs[int(i/2), i % 2].set_title(secir_dict[i], fontsize=8)
        axs[int(i/2), i % 2].set_xlim(left=0, right=dates[-1])
        # axs[int(i/2), i % 2].set_ylim(bottom=np.min(total[:, i]*0.7),
        #                               top=np.max(total[:, i])*1.5)
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
    # fig.legend(labels, loc="lower right",
    #            fancybox=False, shadow=False, ncol=2, fontsize=10)  # bbox_to_anchor=(0.1, -0.73, 0.8, 0.8),

    fig.supxlabel('Simulation time [days]', fontsize=12)
    fig.supylabel('Number of individuals', fontsize=12)
    plt.subplots_adjust(left=None, bottom=0.2, right=None,
                        top=None, wspace=None, hspace=0.6)

    # plt.tight_layout()

    # Save result.
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f"all_compartments" + file_suffix + ".png",
                    bbox_inches='tight', dpi=500)

    plt.clf()


def plot_model_comparison_one_compartment(files, compartment_index, exponential_scenario=False, save_dir="", percentiles=["05", "50", "95"], plot_init=True):
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
    fig, axs = plt.subplots(1, 1, sharex='all', figsize=(4, 3))
    linewidth = 3
    labels = ['ODE', 'LCT', 'IDE', 'ABM']  # , 'IDE'
    linestyles = ['-', '--', ':', '-']

    start_time = 0
    if (plot_init):
        df = pd.read_csv(result_dir + f"comps.csv")
        values = df.iloc[:, 1 + compartment_index]
        for age in range(1, num_age_groups):
            values += df.iloc[:, 1 + compartment_index + age * len(secir_dict)]
        axs.plot(df["Time"], values,
                 color="grey", linestyle=linestyles[3], linewidth=linewidth)
    while len(percentiles) > 0:
        if len(percentiles) == 1:
            df = pd.read_csv(result_dir + f"ABM_p{percentiles[0]}.csv")
            values = df.iloc[:, 1 + compartment_index]
            for age in range(1, num_age_groups):
                values += df.iloc[:, 1 +
                                  compartment_index + age * len(secir_dict)]
            axs.plot(df["Time"] + start_time, values, label=f"p{percentiles[0]} ABM",
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
            axs.plot(df_low["Time"] + start_time, values_low,
                     color=colors['Teal'], alpha=0.3)
            axs.plot(df_high["Time"] + start_time, values_high,
                     color=colors['Teal'], alpha=0.3)
            axs.fill_between(df_low["Time"] + start_time, values_low, values_high, alpha=0.3,
                             color=model_colors[3], label=f"{int(percentiles[-1])-int(percentiles[0])}% quantile ABM")
            del percentiles[0]
            del percentiles[-1]

    # colors = ["C0", "limegreen"]
    # Add results to plot for IDE, LCT and ODE.
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
            if (file == 2):
                print("Compartment ", compartment_index)
                print(total[-50:, compartment_index])
            axs.plot(dates,
                     total[:, compartment_index], label=labels[file],  color=model_colors[file], linestyle=linestyles[file], linewidth=linewidth)

            h5file.close()

    # Define some characteristics of the plot
    # if exponential_scenario:
    #     axs.set_title("S1, one location")
    # else:
    #     axs.set_title("S2")

    # axs.set_xlim(left=dates[0], right=dates[-1])
    # axs.grid(True, linestyle='--', alpha=0.5)
    axs.ticklabel_format(axis='y',
                         style='sci', scilimits=(0, 0))

    axs.set_xlabel('Simulation time [days]', labelpad=10)
    axs.set_ylabel(f'{secir_dict[compartment_index]} [#]')
    # if compartment_index == 3:
    #     axs.set_ylabel('Mildly symptomatic individuals', labelpad=10)
    # if compartment_index == 5:
    #     axs.set_ylabel('Individuals on ICU', labelpad=10)

    # fig.legend(labels, bbox_to_anchor=(0.55, 0.85),
    #            fancybox=False, shadow=False, ncol=2)

    # axs.legend()

    # plt.subplots_adjust(left=0.5, bottom=None, right=0.7,
    #                     top=None, wspace=None, hspace=None)

    set_fontsize()

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

    # Save result.
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        fig.savefig(save_dir + f"model_comparison_{secir_dict[compartment_index]}.png",
                    bbox_inches='tight', dpi=dpi)
        fig_leg.savefig(save_dir + f"model_comparison_legend.png",
                        bbox_inches='tight', dpi=dpi)

    plt.clf()


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

    exponential_scenario = False
    one_location = False
    set_fontsize()

    if exponential_scenario:
        distribution = "exponential"
    else:
        distribution = "different_dists"

    if one_location:
        location = "one_location"
    else:
        location = "multiple_locations"

    folder = "Seed1"
    num_runs = 100

    # Path where simulation results are stored.
    result_dir = os.path.join(os.path.dirname(
        __file__),  f"simulation_results/compare_abm_ide_lct_ode/{num_runs}_runs/{distribution}/{location}/{folder}/")
    # Path where plots will be stored.
    plot_dir = os.path.join(os.path.dirname(
        __file__),  f"plots/compare_abm_ide_lct_ode/{num_runs}_runs/{distribution}/{location}/{folder}/")
    os.makedirs(plot_dir, exist_ok=True)

    num_age_groups = 6
    plot_init = True

    percentiles = ["05", "50", "95"]

    file_suffix = ""
    plot_model_comparison_all_compartments(
        result_dir, percentiles, num_age_groups, file_suffix, plot_init, save_dir=plot_dir)

    for i in range(8):
        plot_model_comparison_one_compartment([os.path.join(result_dir, f"ode"),
                                               os.path.join(
                                                   result_dir, f"lct"),
                                               os.path.join(result_dir, f"ide_no_init")],
                                              i,
                                              exponential_scenario,
                                              save_dir=plot_dir, percentiles=["05", "50", "95"])
        plt.close()

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
