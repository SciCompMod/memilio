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
import numpy as np
import matplotlib.pyplot as plt

import memilio.epidata.getDataIntoPandasDataFrame as gd

from matplotlib.markers import MarkerStyle
from matplotlib.transforms import Affine2D


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
        return 0.5*(yleft-yright) * np.cos(np.pi/(xright-xleft)*(time-xleft)) + 0.5 * (yleft + yright)


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


def define_groundtruth(timesteps, t0, tmax, damping, damping_time, smootherwindow, cont_freq, smoothercos_func=True):

    model = 'groundtruth'
    results = {model: []}

    for timestep in timesteps:
        groundtruth_at_timepoints = []
        for compartment in range(3):
            timepoints = np.linspace(t0, tmax, int((tmax-t0)/timestep)+1)

            if smoothercos_func:
                if compartment == 0:
                    groundtruth_at_timepoints.append([smoothercos(timepoint, damping, damping_time, smootherwindow, cont_freq)
                                                      for timepoint in timepoints])
                if compartment == 1:
                    groundtruth_at_timepoints.append([
                        smoothercos_deriv(timepoint, damping, damping_time, smootherwindow, cont_freq) for timepoint in timepoints])
                if compartment == 2:
                    groundtruth_at_timepoints.append([
                        smoothercos_deriv(timepoint, damping, damping_time, smootherwindow, cont_freq) for timepoint in timepoints])
            else:
                if compartment == 0:
                    groundtruth_at_timepoints.append([smoothstep(timepoint, damping, damping_time, smootherwindow, cont_freq)
                                                      for timepoint in timepoints])
                if compartment == 1:
                    groundtruth_at_timepoints.append([
                        smoothstep_deriv(timepoint, damping, damping_time, smootherwindow, cont_freq) for timepoint in timepoints])
                if compartment == 2:
                    groundtruth_at_timepoints.append([
                        smoothstep_deriv(timepoint, damping, damping_time, smootherwindow, cont_freq) for timepoint in timepoints])

        results[model].append(groundtruth_at_timepoints)

    return results


def read_groundtruth(data_dir, ide_exponents):
    """ Read data into a dict, where the keys correspond to the respective model.
    At the moment we are only storing results of the IDE model here. There, we have an array that contains all results
    obtained with the IDE model for all time points for each time step size that is investigated. The results can
    either be compartments or flows as indicated by the flag 'flows'.
    """
    model = 'groundtruth'
    results = {model: []}

    for exponent in ide_exponents:

        h5file = h5py.File(os.path.join(
            data_dir, f'groundtruth_dt=1e-{exponent:.0f}.h5'), 'r')

        data = h5file[list(h5file.keys())[0]]

        if len(data['Total'][0]) == 3:
            # As there should be only one Group, total is the simulation result.
            results[model].append(data['Total'][:, :])
        else:
            raise gd.DataError(
                "Expected a different size of vector in time series.")

        h5file.close()

    return results


def read_data(data_dir, ide_exponents):
    """ Read data into a dict, where the keys correspond to the respective model.
    At the moment we are only storing results of the IDE model here. There, we have an array that contains all results
    obtained with the IDE model for all time points for each time step size that is investigated. The results can
    either be compartments or flows as indicated by the flag 'flows'.
    """
    model = 'ide'
    results = {model: []}

    for exponent in ide_exponents:

        h5file = h5py.File(os.path.join(
            data_dir, f'result_{model}_dt=1e-{exponent:.0f}.h5'), 'r')

        data = h5file[list(h5file.keys())[0]]

        if len(data['Total'][0]) == 3:
            # As there should be only one Group, total is the simulation result.
            results[model].append(data['Total'][:, :])
        else:
            raise gd.DataError(
                "Expected a different size of vector in time series.")

        h5file.close()

    return results


def compute_l2_norm(timeseries, timestep):
    """ Computes L2 norm of a time series.

    @param[in] timeseries Considered timeseries.
    @param[in] timestep Time step size.
    @returns Norm.
    """
    norm = np.sqrt(timestep * np.sum(timeseries**2))
    return norm


def compute_errors_l2(groundtruth, results, timesteps_ide, t0_ide, cut_off=0):
    """ Computes L2 norm of the difference between time series from ODE and time series
    from IDE for all compartments/flows.

    @param[in] groundtruth Result obtained with ODE model.
    @param[in] results Results obtained with IDE model for different time step sizes.
    @param[in] save_exponent The results of the ODE model were saved using the step size 10^{-save_exponent}.
    @param[in] timesteps_ide List of time steps used in IDE simulations.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False.
    @param[in] Array that contains computed errors.
    """
    num_errors = 3

    errors = []

    # Compute error.
    for i in range(len(results['ide'])):
        errors.append([])
        for compartment in range(num_errors):
            timestep = timesteps_ide[i]

            result_ode = np.array(groundtruth['groundtruth'][i][:, compartment][int(
                t0_ide/timestep)+cut_off::])
            result_ide = np.array(results['ide'][i][int(
                t0_ide/timestep)+cut_off::][:, compartment])

            difference = result_ode - result_ide

            errors[i].append(compute_l2_norm(
                difference, timestep))

    return np.array(errors)


def compute_max_norm(timeseries):
    """ Computes maximum norm of a time series.

    @param[in] timeseries Considered timeseries.
    @returns Norm.
    """
    norm = np.max(np.abs(timeseries))
    return norm


def compute_errors_max(groundtruth, results,  timesteps_ide, t0_ide,  cut_off=0):
    """ Computes maximum norm of the difference between time series from ODE and time series
    from IDE for all compartments.
    """
    num_errors = 3

    errors = []

    # Compute error.
    for i in range(len(results['ide'])):
        errors.append([])
        for compartment in range(num_errors):
            timestep = timesteps_ide[i]

            difference = groundtruth['groundtruth'][i][:, compartment][int(
                t0_ide/timestep)+cut_off::]-results['ide'][i][int(t0_ide/timestep)+cut_off::][:, compartment]

            errors[i].append(compute_max_norm(
                difference))

    return np.array(errors)


def plot_convergence(errors_all_gregory_orders, timesteps_ide,
                     l2=True, maxnorm=False, save_dir=""):
    """ Plots errors against timesteps with a subplot for each compartment /flow.

    @param[in] errors Array that contains computed errors of IDE model compared to groundtruth.
    @param[in] timesteps_ide List of time steps used in IDE simulations.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False.
    @param[in] save_dir Directory where plot will be stored. Default is an empty string leading to the plot not being
        saved.
    """
    # Define subplots and labels.
    num_plots = 3
    num_plotted_results = 1

    fig, axs = plt.subplots(1, num_plots, sharex=True,
                            figsize=(10, 3))
    secir_dict = {0: 'Susceptible', 1:  'Infected', 2:  'Recovered'}
    # labels = [
    #     f"Gregory order {gregory_order}" for gregory_order in gregory_orders_simulation]
    # labels.insert(0, "")
    labels = ["Error"]
    labels.append(r"$\mathcal{O}(\Delta t)$")
    labels.append(r"$\mathcal{O}(\Delta t^2)$")
    labels.append(r"$\mathcal{O}(\Delta t^3)$")
    labels.append(r"$\mathcal{O}(\Delta t^4)$")

    handles = []

    # Define colors.
    colors_ = [plt.cm.viridis(x)
               for x in np.linspace(0, 1, 3)]
    colors = ["darkorange", colors_[1], "darkred"]

    # Plot comparison line for linear convergence as well as second, third and fourth order.
    for i in range(num_plots):
        plotted_timesteps = timesteps_ide[:6]
        comparison = [0.5*errors_all_gregory_orders[0]
                      [0, i]*dt for dt in timesteps_ide]
        first = axs[i].plot(plotted_timesteps, comparison,
                            '--', color='gray', linewidth=1.2, label=r"$\mathcal{O}(\Delta t)$")

        plotted_timesteps = timesteps_ide[:6]
        comparison = [0.5*errors_all_gregory_orders[0]
                      [0, i]*dt**2 for dt in plotted_timesteps]
        second = axs[i].plot(plotted_timesteps, comparison,
                             '--', color=colors[0], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^2)$")

        plotted_timesteps = timesteps_ide[:6]
        comparison = [0.5*errors_all_gregory_orders[0]
                      [0, i]*dt**3 for dt in plotted_timesteps]
        third = axs[i].plot(plotted_timesteps, comparison,
                            '--', color=colors[1], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^3)$")

        plotted_timesteps = timesteps_ide[:6]
        comparison = [0.5*errors_all_gregory_orders[0]
                      [0, i]*dt**4 for dt in plotted_timesteps]
        fourth = axs[i].plot(plotted_timesteps, comparison,
                             '--', color=colors[2], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^4)$")

        # Plot results.

        line = axs[i].plot(timesteps_ide,
                           errors_all_gregory_orders[0][:, i], '-o', color=colors[2], label=labels[0])
        if i == 0:
            handles.append(line[0])

        # Append lines to handles for legend.
        if i == 0:
            handles.append(first[0])
            handles.append(second[0])
            handles.append(third[0])
            handles.append(fourth[0])
            # Invert x axis only for one plot so that sharex=True and invert_xaxis work as intended.
            axs[i].invert_xaxis()
        titles = [r"S", r"S' analytical", r"S' numerical"]
        axs[i].set_title(titles[i], fontsize=10)

        # Adapt plots.
        axs[i].set_xscale("log", base=10)
        axs[i].set_yscale("log", base=10)

        axs[i].grid(True, linestyle='--', alpha=0.6)

    fig.supxlabel(r'Time step $\Delta t$', fontsize=12)

    ylabel = fig.supylabel(
        r"$err$", fontsize=12)

    # print(handles)

    legend = fig.legend(handles=handles, labels=labels, ncol=2,  loc='lower right',
                        fontsize=8, bbox_transform=fig.transFigure, bbox_to_anchor=(1., -0.1))  # bbox_to_anchor=(0.5, -0.06),
    fig.tight_layout()

    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)

        filename = ""

        if l2:
            filename = f'{save_dir}/convergence_all_compartments_l2'
        elif maxnorm:
            filename = f'{save_dir}/convergence_all_compartments_max'

        filename = filename + "_abs"

        plt.savefig(filename + ".png", format='png', bbox_extra_artists=(legend, ylabel), bbox_inches='tight',
                    dpi=500)

    plt.close()


def plot_difference_per_timestep(groundtruth, results, timesteps_ide, t0_ide, cut_off=0, save_dir="", damping_time=-1, smootherwindow=1):
    num_errors = 3

    errors = []

    compartments = [r"S", r"S' analytical", r"S' numerical"]

    # Compute error.
    difference = []

    num_plots = 3
    figsize_x = 12

    for i, timestep in enumerate(timesteps_ide):
        fig, axs = plt.subplots(1, num_plots, sharex=True,
                                figsize=(figsize_x, 3))
        if damping_time != -1:
            fig_zoom, axs_zoom = plt.subplots(1, num_plots, sharex=True,
                                              figsize=(figsize_x, 3))

        errors.append([])
        for compartment in range(num_errors):

            difference = groundtruth['groundtruth'][i][:, compartment][int(
                t0_ide/timestep)+cut_off::]-results['ide'][i][int(t0_ide/timestep)+cut_off::][:, compartment]

            indices = np.linspace(
                t0_ide, t0_ide+(len(difference)-1)*timestep, len(difference))

            axs[compartment].scatter(indices[0:], difference[0:], s=1)
            axs[compartment].set_title(f"{compartments[compartment]}")

            if damping_time != -1:

                difference_zoom = difference[int((damping_time-t0_ide-smootherwindow)/timestep):int(
                    (damping_time-t0_ide)/timestep)+1]
                indices_zoom = indices[int((damping_time-t0_ide-smootherwindow)/timestep):int(
                    (damping_time-t0_ide)/timestep)+1]

                axs_zoom[compartment].scatter(
                    indices_zoom, difference_zoom, s=1)
                axs_zoom[compartment].set_title(f"{compartments[compartment]}")

        fig.supxlabel("Time")
        fig_zoom.supxlabel("Time")

        # plt.show()

        if save_dir != "":
            if not os.path.isdir(f"{save_dir}/differences"):
                os.makedirs(f"{save_dir}/differences")

            filename = f"{save_dir}/differences/timestep={timestep}"

            fig.savefig(filename + ".png", format='png',
                        dpi=500)

            if damping_time != -1:
                if not os.path.isdir(f"{save_dir}/differences_zoom"):
                    os.makedirs(f"{save_dir}/differences_zoom")

                filename = f"{save_dir}/differences_zoom/timestep={timestep}"

                fig_zoom.savefig(filename + ".png", format='png',
                                 dpi=500)

        plt.close()
        plt.close()


def compute_order_of_convergence(errors, timesteps_ide):
    """ Compute order of convergence between two consecutive time step sizes.

    @param[in] errors Array that contains computed errors of IDE model compared to groundtruth.
    @param[in] timesteps_ide List of time steps used in IDE simulations.
    """
    num_orders = 3

    order = []
    for compartment in range(num_orders):
        order.append([])
        for i in range(len(errors)-1):
            order[compartment].append(np.log(errors[i+1][compartment]/errors[i][compartment]) /
                                      np.log(timesteps_ide[i+1]/timesteps_ide[i]))
    return np.array(order)


def subfolders_scandir(path):
    # path = os.path.dirname(path)
    print(path)
    with os.scandir(path) as it:
        return [entry.name for entry in it if entry.is_dir()]


def get_ide_exponents(data_dir):
    files = os.listdir(data_dir)

    max_exponent = 0
    for possible_exponent in range(5):
        if f'result_{"ide"}_dt=1e-{possible_exponent}.h5' in files:
            max_exponent = possible_exponent

    return range(max_exponent+1)


def get_t0_from_dir_name(dir_name):
    t0_string = [x for x in dir_name.split("_") if "t0" in x]
    t0 = int(t0_string[0].split("=")[-1])

    return t0


def get_tinit_from_dir_name(dir_name):
    tinit_string = [x for x in dir_name.split("_") if "tinit" in x]
    tinit = int(tinit_string[0].split("=")[-1])

    return tinit


def get_tmax_from_dir_name(dir_name):
    tmax_string = [x for x in dir_name.split("_") if "tmax" in x]
    tmax = int(tmax_string[0].split("=")[-1])

    return tmax


def get_damping_from_dir_name(dir_name):
    damping_string = [x for x in dir_name.split("_") if "damping=" in x]
    damping = float(damping_string[0].split("=")[-1])

    return damping


def get_dampingtime_from_dir_name(dir_name):
    string = [x for x in dir_name.split("_") if "dampingtime" in x]
    value = int(string[0].split("=")[-1])

    return value


def get_contfreq_from_dir_name(dir_name):
    string = [x for x in dir_name.split("_") if "contfreq" in x]
    value = int(string[0].split("=")[-1])

    return value


def main():

    fdorder = 4
    smootherwindow = 2

    # smoother_func_str = "smoothercos"

    root_dir = os.path.join(os.path.dirname(
        __file__), "../simulation_results")
    main_dir = f"2026-05-05/sigmoid_fdordercontacts=4_smootherwindow=5_sigmoidparam=5"
    relevant_dir = os.path.join(root_dir, main_dir)
    # print(relevant_dir)
    sub_dirs = subfolders_scandir(relevant_dir)
    # sub_dirs = [sub_dirs[-1]]

    cut_off = 0

    for dir_index, dir_name in enumerate(sub_dirs):
        print(dir_name)

        t0 = get_t0_from_dir_name(dir_name)
        tinit = get_tinit_from_dir_name(dir_name)
        tmax = get_tmax_from_dir_name(dir_name)

        damping = get_damping_from_dir_name(dir_name)
        damping_time = get_dampingtime_from_dir_name(dir_name)
        cont_freq = get_contfreq_from_dir_name(dir_name)

        # Path where simulation results are stored.
        result_dir = os.path.join(os.path.dirname(
            __file__),  f"{relevant_dir}/{dir_name}/")

        # Path where plots will be stored.
        plot_dir = os.path.join(os.path.dirname(
            __file__),  f"../plots/{main_dir}/{dir_name}/")

        # errors_all_gregory_orders_l2_rel = []
        errors_all_gregory_orders_l2_abs = []

        # errors_all_gregory_orders_max_rel = []
        errors_all_gregory_orders_max_abs = []

        # Get exponents for which IDE simulations have been computed for considered subdirectory.
        ide_exponents = get_ide_exponents(result_dir)

        # Calculate time steps resulting from ide_exponents.
        timesteps_ide = []
        for exp in ide_exponents:
            timesteps_ide.append(pow(10, -exp))

        # groundtruth = define_groundtruth(
        #     timesteps_ide, t0, tmax, damping, damping_time, smootherwindow, cont_freq, smoothercos_func)
        groundtruth = read_groundtruth(result_dir, ide_exponents)
        # Read results from IDE simulations.
        results = read_data(result_dir, ide_exponents)

        # Compute errors of IDE results compared to groundtruth.
        errors_l2_abs = compute_errors_l2(
            groundtruth, results, timesteps_ide, t0,  cut_off)
        errors_all_gregory_orders_l2_abs.append(errors_l2_abs)

        errors_max_abs = compute_errors_max(
            groundtruth, results, timesteps_ide, t0,  cut_off)
        errors_all_gregory_orders_max_abs.append(errors_max_abs)

        plot_difference_per_timestep(
            groundtruth, results, timesteps_ide, t0, cut_off, plot_dir, damping_time, smootherwindow)

        # L2 norm
        l2 = True
        maxnorm = False
        plot_convergence(errors_all_gregory_orders_l2_abs, timesteps_ide,
                         l2, maxnorm,  plot_dir)

        # max norm
        l2 = False
        maxnorm = True
        plot_convergence(errors_all_gregory_orders_max_abs, timesteps_ide,
                         l2, maxnorm,  plot_dir)


if __name__ == '__main__':
    main()
