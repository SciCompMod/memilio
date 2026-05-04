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


def smoothercos(time, damping, damping_time, smootherwindow, cont_freq):

    xleft = damping_time-smootherwindow
    xright = damping_time

    yleft = cont_freq
    yright = (1.-damping)*cont_freq

    if time <= xleft:
        return yleft

    if time >= xright:
        return yright

    else:
        result = 0.5*(yleft-yright) * np.cos(np.pi/(xright-xleft)
                                             * (time-xleft)) + 0.5 * (yleft + yright)
        return result


def smoothercos_deriv(time, damping, damping_time, smootherwindow, cont_freq):

    xleft = damping_time-smootherwindow
    xright = damping_time

    yleft = cont_freq
    yright = (1.-damping)*cont_freq

    if time <= xleft or time >= xright:
        return 0

    else:
        return -0.5 * (yleft-yright)/(xright-xleft) * np.pi * np.sin(np.pi / (xright-xleft)*(time-xleft))


def smoothstep(time, damping, damping_time, smootherwindow, cont_freq):

    xleft = damping_time-smootherwindow
    xright = damping_time

    yleft = cont_freq
    yright = (1.-damping)*cont_freq

    if time <= xleft:
        return yleft

    if time >= xright:
        return yright

    else:

        normalized_time = (time - xleft) / (xright - xleft)

        smoothed_value = yleft + (yright - yleft) * \
            (3. * normalized_time**2 - 2 * normalized_time**3)

        return smoothed_value


def smoothstep_deriv(time, damping, damping_time, smootherwindow, cont_freq):

    xleft = damping_time-smootherwindow
    xright = damping_time

    yleft = cont_freq
    yright = (1.-damping)*cont_freq

    if time <= xleft or time >= xright:
        return 0

    else:
        normalized_time = (time - xleft) / (xright - xleft)
        normalized_time_deriv = 1. / (xright - xleft)
        deriv = (yright - yleft) * (6. * normalized_time * normalized_time_deriv -
                                    6 * normalized_time ** 2 * normalized_time_deriv)
        return deriv


def plot_susceptibles(files, fileending, save_dir, damping, damping_time, smootherwindow, cont_freq, smoothercos_func=True):
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
    secir_dict = {0: 'Susceptible',  1: "S' analytical", 2: "S' numerical"}

    # Define plot.
    num_plots = 3
    fig, axs = plt.subplots(1, num_plots, sharex='all', num='Compare files')

    colors = ["C0", "limegreen"]
    linestyles = ['-', '--',]
    linewidth = 1
    labels = ["Groundtruth", "IDE result"]

    dates = []

    # Add results to plot.

    # Load data.
    h5file = h5py.File(str(files[0]) + '.h5', 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]

    # As there should be only one Group, total is the simulation result.
    total = data['Total'][:, :]

    dates = data['Time'][:]

    dates_groundtruth = np.linspace(
        data['Time'][:][0], data['Time'][:][-1], 1000)

    if smoothercos_func:
        axs[0].plot(dates_groundtruth,
                    [
                        smoothercos(timepoint, damping, damping_time, smootherwindow, cont_freq) for timepoint in dates_groundtruth], label=labels[0],  linestyle=linestyles[0], color=colors[0], linewidth=linewidth)
        axs[1].plot(dates_groundtruth,
                    [
                        smoothercos_deriv(timepoint, damping, damping_time, smootherwindow, cont_freq) for timepoint in dates_groundtruth], label=labels[0],  linestyle=linestyles[0], color=colors[0], linewidth=linewidth)
        axs[2].plot(dates_groundtruth,
                    [
                        smoothercos_deriv(timepoint, damping, damping_time, smootherwindow, cont_freq) for timepoint in dates_groundtruth], label=labels[0],  linestyle=linestyles[0], color=colors[0], linewidth=linewidth)

    else:
        axs[0].plot(dates_groundtruth,
                    [
                        smoothstep(timepoint, damping, damping_time, smootherwindow, cont_freq) for timepoint in dates_groundtruth], label=labels[0],  linestyle=linestyles[0], color=colors[0], linewidth=linewidth)
        axs[1].plot(dates_groundtruth,
                    [
                        smoothstep_deriv(timepoint, damping, damping_time, smootherwindow, cont_freq) for timepoint in dates_groundtruth], label=labels[0],  linestyle=linestyles[0], color=colors[0], linewidth=linewidth)
        axs[2].plot(dates_groundtruth,
                    [
                        smoothstep_deriv(timepoint, damping, damping_time, smootherwindow, cont_freq) for timepoint in dates_groundtruth], label=labels[0],  linestyle=linestyles[0], color=colors[0], linewidth=linewidth)

    # Plot data.
    for i in range(num_plots):
        axs[i].plot(dates,
                    total[:, i], label=labels[1],  linestyle=linestyles[1], color=colors[1], linewidth=linewidth)

    h5file.close()

    # Define some characteristics of the plot
    for i in range(num_plots):
        axs[i].set_title(secir_dict[i], fontsize=8)
        axs[i].set_xlim(left=dates[0], right=dates[-1])
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


def subfolders_scandir(path):
    # path = os.path.dirname(path)
    print(path)
    with os.scandir(path) as it:
        return [entry.name for entry in it if entry.is_dir()]


def get_damping_from_dir_name(dir_name):
    damping_string = [x for x in dir_name.split("_") if "damping=" in x]
    damping = float(damping_string[0].split("=")[-1])

    return damping


def get_dampingtime_from_dir_name(dir_name):
    string = [x for x in dir_name.split("_") if "dampingtime=" in x]
    value = int(string[0].split("=")[-1])

    return value


def get_contfreq_from_dir_name(dir_name):
    string = [x for x in dir_name.split("_") if "contfreq" in x]
    value = int(string[0].split("=")[-1])

    return value


if __name__ == '__main__':

    fdorder = 4
    smootherwindow = 10

    smoothercos_func = False

    # dir_name = "detailed_init_exponential_t0ide=50_tmax=51_finite_diff=1_tolexp=8"
    root_dir = os.path.join(os.path.dirname(
        __file__), "../simulation_results")
    main_dir = f"2026-04-30/smoothstep_fdordercontacts={fdorder}_smootherwindow={smootherwindow}"
    relevant_dir = os.path.join(root_dir, main_dir)

    sub_dirs = subfolders_scandir(relevant_dir)

    for dir_name in sub_dirs:
        print(main_dir + "/" + dir_name)

        # Path where simulation results are stored.
        result_dir = os.path.join(os.path.dirname(
            __file__),  f"../simulation_results/{main_dir}/{dir_name}/")
        # Path where plots will be stored.
        plot_dir = os.path.join(os.path.dirname(
            __file__),  f"../plots/{main_dir}/{dir_name}/")

        gregory_orders = [3]

        damping = get_damping_from_dir_name(dir_name)
        damping_time = get_dampingtime_from_dir_name(dir_name)
        cont_freq = get_contfreq_from_dir_name(dir_name)

        for ide_exponent in [0, 1, 2, 3, 4]:
            plot_susceptibles([
                os.path.join(result_dir, f"result_ide_dt=1e-{ide_exponent}")],
                f"dt=1e-{ide_exponent}", plot_dir, damping, damping_time, smootherwindow, cont_freq, smoothercos_func)
