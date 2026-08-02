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

from matplotlib.markers import MarkerStyle
from matplotlib.transforms import Affine2D


def read_groundtruth_ide(data_dir, groundtruth_exponent, gregory_order):
    """ Read groundtruth from data. We define the groundtruth as the results obtained by the ODE model with timestep dt=1e-6.

    @param[in] data_dir Directory where h5 files are stored.
    @param[in] ode_exponent Exponent that determines time step size via dt =10^{-ode_exponent}.
    @param[in] save_exponent The results of the ODE model were saved using the step size 10^{-save_exponent}.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False.
    @returns Dict with results of ODE model.
    """
    model = 'ide'
    results = []

    h5file = h5py.File(os.path.join(
        data_dir, f'result_{model}_dt=1e-{groundtruth_exponent:.0f}_gregoryorder={gregory_order}.h5'), 'r')

    if (len(list(h5file.keys())) > 1):
        print("File should contain one dataset.")
        return
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        print("Expected only one group.")
        return

    data = h5file[list(h5file.keys())[0]]

    if len(data['Total'][0]) == 3:
        # As there should be only one Group, total is the simulation result
        results.append(data['Total'][:, :])
    else:
        raise gd.DataError(
            'Expected a different size of vector in time series.')

    dates = data['Time'][:]

    h5file.close()

    return results


def read_groundtruth_ode(data_dir, groundtruth_exponent, groundtruth_save_exponent):
    """ Read groundtruth from data. We define the groundtruth as the results obtained by the ODE model with timestep dt=1e-6.

    @param[in] data_dir Directory where h5 files are stored.
    @param[in] ode_exponent Exponent that determines time step size via dt =10^{-ode_exponent}.
    @param[in] save_exponent The results of the ODE model were saved using the step size 10^{-save_exponent}.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False.
    @returns Dict with results of ODE model.
    """
    model = 'ode'
    results = []

    h5file = h5py.File(os.path.join(
        data_dir, f'result_{model}_dt=1e-{groundtruth_exponent:.0f}_savedt=1e-{groundtruth_save_exponent:.0f}.h5'), 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]

    if len(data['Total'][0]) == 3:
        # As there should be only one Group, total is the simulation result
        results.append(data['Total'][:, :])
    else:
        raise gd.DataError(
            'Expected a different size of vector in time series.')

    dates = data['Time'][:]

    h5file.close()

    return results


def read_groundtruth(result_dir, groundtruth_exponent, groundtruth_save_exponent,  gregory_order_groundtruth=3, groundtruth_ode=True):

    if not groundtruth_ode:
        return read_groundtruth_ide(
            result_dir, groundtruth_exponent, gregory_order_groundtruth)

    return read_groundtruth_ode(
        result_dir, groundtruth_exponent, groundtruth_save_exponent)


def read_data(data_dir, ide_exponents, gregory_order):
    """ Read data into a dict, where the keys correspond to the respective model.
    At the moment we are only storing results of the IDE model here. There, we have an array that contains all results
    obtained with the IDE model for all time points for each time step size that is investigated. The results can
    either be compartments or flows as indicated by the flag 'flows'.
    """
    model = 'ide'
    results = []

    for exponent in ide_exponents:

        h5file = h5py.File(os.path.join(
            data_dir, f'result_{model}_dt=1e-{exponent:.0f}_gregoryorder={gregory_order}.h5'), 'r')

        data = h5file[list(h5file.keys())[0]]

        if len(data['Total'][0]) == 3:
            # As there should be only one Group, total is the simulation result.
            results.append(data['Total'][:, :])
            # print("ide exponent: ", exponent)
            # print("value at t_init:", results[0][0])
        else:
            raise gd.DataError(
                "Expected a different size of vector in time series.")

        h5file.close()

    return results


def compute_errors_endpoint(groundtruth, results, relative_error=True):
    """ Computes relative L2 norm of the difference between time series from ODE and time series
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

    # Compute error. Here, we define the error by the absolute value of the difference at the last time point between
    # groundtruth and simulation results.
    for i in range(len(results)):
        errors.append([])

        for compartment in range(num_errors):

            # model = list(groundtruth.keys())[0]

            if relative_error:

                error = np.abs(
                    (groundtruth[0][-1, compartment]-results[i][-1, compartment])/np.abs(groundtruth[0][-1, compartment]))

            else:
                error = np.abs(groundtruth[0][-1, compartment] -
                               results[i][-1, compartment])

            errors[i].append(error)

    return np.array(errors)


def compute_l2_norm(timeseries, timestep):
    """ Computes L2 norm of a time series.

    @param[in] timeseries Considered timeseries.
    @param[in] timestep Time step size.
    @returns Norm.
    """
    norm = np.sqrt(timestep * np.sum(timeseries**2))
    return norm


def compute_errors_l2(groundtruth, results, groundtruth_save_exponent, timesteps_ide, t0_ide, t_init, error_type):
    """ Computes relative L2 norm of the difference between time series from ODE and time series
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
    for i in range(len(results)):
        errors.append([])
        for compartment in range(num_errors):
            timestep = timesteps_ide[i]
            scale_timesteps = timestep/pow(10, -groundtruth_save_exponent)
            num_timepoints = len(results[i])

            groundtruth_sim_interval = groundtruth[0][int(
                pow(10, groundtruth_save_exponent)*(t0_ide))::int(scale_timesteps)][:, compartment]

            results_sim_interval = results[i][int(
                (t0_ide-t_init)/timestep)::][:, compartment]

            difference = groundtruth_sim_interval-results_sim_interval

            if error_type == "abs":
                errors[i].append(compute_l2_norm(
                    difference, timestep))

            if error_type == "rel":
                relative_difference = difference/groundtruth_sim_interval
                errors[i].append(compute_l2_norm(
                    relative_difference, timestep))

            if error_type == "weighted":  # weighted
                norm_groundtruth = compute_l2_norm(
                    groundtruth_sim_interval, timestep)
                errors[i].append(compute_l2_norm(
                    difference, timestep)/norm_groundtruth)

    return np.array(errors)


def compute_max_norm(timeseries):
    """ Computes maximum norm of a time series.

    @param[in] timeseries Considered timeseries.
    @returns Norm.
    """
    # print(timeseries)
    norm = np.max(np.abs(timeseries))
    return norm


def compute_errors_max(groundtruth, results, groundtruth_save_exponent, timesteps_ide, t0_ide, t_init,  error_type):
    """ Computes relative maximum norm of the difference between time series from ODE and time series
    from IDE for all compartments.
    """
    # num_errors = 3

    errors = []

    compartments = ["S", "I", "R"]

    # Compute error.

    for i in range(len(results)):
        # if not relative_error:
        #     print(f"Timestep {timesteps_ide[i]}")
        errors.append([])
        for compartment in range(len(compartments)):
            timestep = timesteps_ide[i]
            scale_timesteps = timestep/pow(10, -groundtruth_save_exponent)
            # print("scale timesteps: ", scale_timesteps)

            groundtruth_sim_interval = groundtruth[0][int(
                pow(10, groundtruth_save_exponent)*(t0_ide))::int(scale_timesteps)][:, compartment]

            results_sim_interval = results[i][int(
                (t0_ide-t_init)/timestep)::][:, compartment]

            difference = groundtruth_sim_interval-results_sim_interval

            # print("timestep: ", timestep)
            # print("difference at t_init: ", difference[0])
            # print()

            # to debug:
            # groundtruth[0][int(pow(10, groundtruth_save_exponent)*(t_init))::int(scale_timesteps)][:, compartment] - results[i][0::][:, compartment]

            if error_type == "abs":
                errors[i].append(compute_max_norm(
                    difference))

            if error_type == "rel":
                relative_difference = difference/groundtruth_sim_interval
                errors[i].append(compute_max_norm(
                    relative_difference))

            if error_type == "weighted":  # weighted
                norm_groundtruth = compute_max_norm(groundtruth_sim_interval)
                errors[i].append(compute_max_norm(
                    difference)/norm_groundtruth)

    return np.array(errors)


def plot_difference_per_timestep(groundtruth, results, groundtruth_save_exponent, timesteps_ide, t0_ide, t_init,  gregory_order,  save_dir="", damping_time=float(-1)):
    num_errors = 3

    # errors = []
    # errors_zoom = []

    compartments = ["S", "I", "R"]

    # Compute error.
    difference = []

    num_plots = 3
    figsize_x = 12

    for i, timestep in enumerate(timesteps_ide):
        fig, axs = plt.subplots(1, num_plots, sharex=True,
                                figsize=(figsize_x, 3))

        # if damping_time != -1:
        fig_zoom, axs_zoom = plt.subplots(1, num_plots, sharex=True,
                                          figsize=(figsize_x, 3))

        for compartment in range(num_errors):
            scale_timesteps = timestep/pow(10, -groundtruth_save_exponent)

            difference = groundtruth[0][int(
                pow(10, groundtruth_save_exponent)*(t0_ide))::int(scale_timesteps)][:, compartment]-results[i][int((t0_ide-t_init)/timestep)::][:, compartment]

            indices = np.linspace(
                t0_ide, t0_ide+(len(difference)-1)*timestep, len(difference))

            axs[compartment].scatter(indices, difference, s=1)
            axs[compartment].set_title(f"{compartments[compartment]}")

            if damping_time != -1:

                padding_index = 2
                difference_zoom = difference[int((damping_time-t0_ide-1-padding_index)/timestep):int(
                    (damping_time-t0_ide+padding_index)/timestep)+1]
                indices_zoom = indices[int((damping_time-t0_ide-1-padding_index)/timestep):int(
                    (damping_time-t0_ide+padding_index)/timestep)+1]

                axs_zoom[compartment].scatter(
                    indices_zoom, difference_zoom, s=1)
                axs_zoom[compartment].set_title(f"{compartments[compartment]}")

        fig.supxlabel("Time")
        fig_zoom.supxlabel("Time")

        if save_dir != "":
            if not os.path.isdir(f"{save_dir}/differences"):
                os.makedirs(f"{save_dir}/differences")

            filename = f"{save_dir}/differences/gregoryorder={gregory_order}_timestep={timestep}"

            fig.savefig(filename + ".png", format='png',
                        dpi=500)

            if damping_time != -1:
                if not os.path.isdir(f"{save_dir}/differences_zoom"):
                    os.makedirs(f"{save_dir}/differences_zoom")

                filename = f"{save_dir}/differences_zoom/gregoryorder={gregory_order}_timestep={timestep}"

                fig_zoom.savefig(filename + ".png", format='png',
                                 dpi=500)

        plt.close()


def plot_convergence(errors_all_gregory_orders, timesteps_ide,
                     gregory_orders_simulation, norm, error_type, save_dir="", only_S=False):
    """ Plots errors against timesteps with a subplot for each compartment /flow.

    @param[in] errors Array that contains computed errors of IDE model compared to groundtruth.
    @param[in] timesteps_ide List of time steps used in IDE simulations.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False.
    @param[in] save_dir Directory where plot will be stored. Default is an empty string leading to the plot not being
        saved.
    """
    # Define subplots and labels.
    if only_S:
        num_plots = 1
        figsize_x = 5
    else:
        num_plots = 3
        figsize_x = 10

    num_plotted_results = len(gregory_orders_simulation)

    fig, axs = plt.subplots(1, num_plots, sharex=True,
                            figsize=(figsize_x, 3.5))
    secir_dict = {0: 'Susceptible', 1:  'Infected', 2:  'Recovered'}
    labels = [
        f"Gregory order {gregory_order}" for gregory_order in gregory_orders_simulation]
    # labels = ["Trapez. rule", "I Gregory rule", "II Gregory rule"]
    labels.insert(0, "")
    labels.append(r"$\mathcal{O}(\Delta t)$")
    labels.append(r"$\mathcal{O}(\Delta t^2)$")
    labels.append(r"$\mathcal{O}(\Delta t^3)$")
    labels.append(r"$\mathcal{O}(\Delta t^4)$")

    handles = [plt.Line2D([], [], color='none')]

    # # Define colors.
    # colors_ = [plt.cm.viridis(x)
    #            for x in np.linspace(0, 1, num_plotted_results)]
    # if len(colors_) > 2:
    #     colors = ["darkorange", colors_[1], "darkred"]
    # else:
    #     colors = ["darkred"]

    colors = ["#332288", "#44AA99", "#882255"]  # sand, teal, wine
    gray = "#888888"

    for i in range(num_plots):

        if only_S:
            ax_obj = axs
        else:
            ax_obj = axs[i]

        for j in range(len(gregory_orders_simulation)):
            # Plot results.
            if i == 0:
                line = ax_obj.plot(timesteps_ide,
                                   errors_all_gregory_orders[j][:, i], '-o', color=colors[j], label=labels[j])
                handles.append(line[0])
            else:
                line = ax_obj.plot(timesteps_ide,
                                   errors_all_gregory_orders[j][:, i], '-o', color=colors[j])
                # handles.append(line)

        # Plot comparison line for linear convergence as well as second, third and fourth order.
        plotted_timesteps = timesteps_ide[:6]
        first_timestep = timesteps_ide[0]
        comparison = [0.5*1/first_timestep*errors_all_gregory_orders[0]
                      [0, i]*dt for dt in timesteps_ide]
        first = ax_obj.plot(plotted_timesteps, comparison,
                            '--', color=gray, linewidth=1.2, label=r"$\mathcal{O}(\Delta t)$")

        plotted_timesteps = timesteps_ide[:6]
        comparison = [0.5*(1/first_timestep)**2*errors_all_gregory_orders[0]
                      [0, i]*dt**2 for dt in plotted_timesteps]
        second = ax_obj.plot(plotted_timesteps, comparison,
                             '--', color=colors[0], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^2)$")

        plotted_timesteps = timesteps_ide[:6]
        comparison = [0.5*(1/first_timestep)**3*errors_all_gregory_orders[1]
                      [0, i]*dt**3 for dt in plotted_timesteps]
        third = ax_obj.plot(plotted_timesteps, comparison,
                            '--', color=colors[1], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^3)$")

        plotted_timesteps = timesteps_ide[:6]
        comparison = [0.5*(1/first_timestep)**4*errors_all_gregory_orders[2]
                      [0, i]*dt**4 for dt in plotted_timesteps]
        fourth = ax_obj.plot(plotted_timesteps, comparison,
                             '--', color=colors[2], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^4)$")

        # Append lines to handles for legend.
        if i == 0:
            handles.append(first[0])
            handles.append(second[0])
            handles.append(third[0])
            handles.append(fourth[0])
            # Invert x axis only for one plot so that sharex=True and invert_xaxis work as intended.
            ax_obj.invert_xaxis()

        # Adapt plots.
        ax_obj.set_xscale("log", base=10)
        ax_obj.set_yscale("log", base=10)

        ax_obj.set_title(secir_dict[i], fontsize=10)
        ax_obj.grid(True, linestyle='--', alpha=0.6)

    ax_label = axs if only_S else axs[1]
    ax_label.set_xlabel(r'Time step $\Delta t$', fontsize=12, labelpad=15)

    if error_type == "abs":
        ylabel = fig.supylabel(
            r"$err_{abs}$", fontsize=12)

    if error_type == "rel":
        ylabel = fig.supylabel(
            r"$err_{rel}$", fontsize=12)

    if error_type == "weighted":
        ylabel = fig.supylabel(
            r"$err_{weighted}$", fontsize=12)

    # print(handles)

    labels_reordered = [labels[0], labels[4], labels[1],
                        labels[5], labels[2], labels[6], labels[3], labels[7]]

    handles_reordered = [handles[0], handles[4], handles[1],
                         handles[5], handles[2], handles[6], handles[3], handles[7]]

    legend = fig.legend(handles=handles_reordered, labels=labels_reordered, ncol=4,  loc='lower center',
                        fontsize=8, bbox_transform=fig.transFigure, bbox_to_anchor=(0.53, -0.1))  # bbox_to_anchor=(0.5, -0.2), # bbox_to_anchor=(1., -0.1)
    fig.tight_layout()

    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)

        filename = f'{save_dir}/convergence_all_compartments'
        if norm == "max":
            filename += "_max"
        elif norm == "l2":
            filename += "_l2"

        if error_type == "abs":
            filename += "_abs"
        elif error_type == "rel":
            filename += "_rel"
        elif error_type == "weighted":
            filename += "_weighted"

        plt.savefig(filename + ".png", format='png', bbox_extra_artists=(legend, ylabel), bbox_inches='tight',
                    dpi=500)

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


def get_total_pop_end(results):

    total_pop_end = []
    for timestep in range(len(results)):
        total_pop_end.append(results[timestep][-1].sum())

    return total_pop_end


def plot_total_pop_diff(gregory_orders_simulation, fd_orders, timesteps_ide, total_pop_all_fd_orders, total_pop_reference, save_dir=""):

    # Plot relative difference per time step for FD order for a fixed Gregory order.
    fig, axs = plt.subplots(1, 1, sharex=True, figsize=(9, 3))
    labels = [
        f"FD order {fd_order}" for fd_order in fd_orders]

    # Define colors.
    colors = ["purple", "green", "mediumblue"]

    gregory_index = 2

    for fd_index in range(len(fd_orders)):
        # Plot results.
        print("FD index: ", fd_orders[fd_index])

        for timestep in range(len(timesteps_ide)):
            print(np.abs(1 - total_pop_all_fd_orders[fd_index]
                  [gregory_index][timestep]/total_pop_reference))
        line = axs.plot(timesteps_ide,
                        np.abs(1 - total_pop_all_fd_orders[fd_index][gregory_index][:]/total_pop_reference), '-o', color=colors[fd_index], label=labels[fd_index])

    axs.invert_xaxis()
    axs.set_xscale("log", base=10)
    axs.set_yscale("log", base=10)

    fig.legend(bbox_to_anchor=(0.97, 0.9))

    fig.supxlabel(r'Time step $\Delta t$', fontsize=12)
    fig.supylabel(
        r"Relative deviation at $t_\max$", fontsize=12)

    axs.set_title(
        f"Gregory order {gregory_orders_simulation[gregory_index]}")

    plt.tight_layout()

    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/mass_conservation_diff_gregory={gregory_orders_simulation[gregory_index]}.png', format='png',  bbox_inches='tight',
                    dpi=500)

    plt.close()


def subfolders_scandir(path, prefix=None):
    # path = os.path.dirname(path)
    print(path)
    with os.scandir(path) as it:
        return [entry.name for entry in it
                if entry.is_dir() and (prefix is None or entry.name.startswith(prefix))]


def get_ide_exponents(data_dir):
    files = os.listdir(data_dir)

    max_exponent = 0
    for possible_exponent in range(5):
        if f'result_{"ide"}_dt=1e-{possible_exponent}_gregoryorder=3.h5' in files:
            max_exponent = possible_exponent

    return range(max_exponent+1)


def get_t0_ide_from_dir_name(dir_name):
    t0_string = [x for x in dir_name.split(
        "_") if ("t0" in x)]
    t0 = float(t0_string[0].split("=")[-1])

    return t0


def get_tinit_from_dir_name(dir_name):
    tinit_string = [x for x in dir_name.split(
        "_") if ("tinit" in x)]
    t_init = float(tinit_string[0].split("=")[-1])

    return t_init


def get_dampingtime_ide_from_dir_name(dir_name):
    dampingtime_string = [x for x in dir_name.split(
        "_") if ("dampingtime" in x)]

    if dampingtime_string==[]:
        return float(-1)

    else:
        return float(dampingtime_string[0].split("=")[-1])

 


def get_tmax_ide_from_dir_name(dir_name):
    tmax_string = [x for x in dir_name.split(
        "_") if ("tmax" in x)]
    tmax = float(tmax_string[0].split("=")[-1])

    return tmax


# def get_savedt_from_dir_name(dir_name):
#     savedt_string = [x for x in dir_name.split(
#         "_") if ("savedt" in x)]
#     savedt = float(savedt_string[0].split("=")[-1])

#     return savedt


def main():

    groundtruth_exponent = 6
    groundtruth_save_exponent = 2
    only_S = False

    main_dir = f"2026-07-30/convergence_lct_dtode=1e-6_t0ode=0_timeinf=2_contfreq=0.73"

    ##############################################

    root_dir = os.path.join(os.path.dirname(
        __file__), "../simulation_results")
    relevant_dir = os.path.join(root_dir, main_dir)
    # print(relevant_dir)
    sub_dirs = subfolders_scandir(relevant_dir)
    # sub_dirs = [sub_dirs[-1]]

    total_pop_reference = 0
    total_pop_all_fd_orders = []

    gregory_orders_simulation = [1, 2, 3]

    for dir_index, dir_name in enumerate(sub_dirs):
        print(dir_name)

        parent_dir = os.path.join(relevant_dir, dir_name)

        t0_ide = get_t0_ide_from_dir_name(dir_name)
        # savedt = get_savedt_from_dir_name(dir_name)

        # Read groundtruth from the top-level directory
        groundtruth = read_groundtruth(
            parent_dir, groundtruth_exponent, groundtruth_save_exponent)

        ide_sub_dirs = subfolders_scandir(parent_dir, prefix="tinit=")
        if not ide_sub_dirs:
            ide_sub_dirs = [""]

        for ide_sub_dir in ide_sub_dirs:
            if ide_sub_dir != "":
                t_init = get_tinit_from_dir_name(ide_sub_dir)
                ide_result_dir = os.path.join(parent_dir, ide_sub_dir)
                plot_dir = os.path.join(os.path.dirname(
                    __file__),  f"../plots/{main_dir}/{dir_name}/{ide_sub_dir}")
            else:
                t_init = get_tinit_from_dir_name(dir_name)
                ide_result_dir = parent_dir
                plot_dir = os.path.join(os.path.dirname(
                    __file__),  f"../plots/{main_dir}/{dir_name}")

            # errors_all_gregory_orders_l2_rel = []
            errors_all_gregory_orders_l2_abs = []

            errors_all_gregory_orders_max_abs = []
            errors_all_gregory_orders_max_rel = []
            errors_all_gregory_orders_max_weighted = []

            # Get exponents for which IDE simulations have been computed for considered directory.
            ide_exponents = get_ide_exponents(ide_result_dir)

            # Calculate time steps resulting from ide_exponents.
            timesteps_ide = []
            for exp in ide_exponents:
                timesteps_ide.append(pow(10, -exp))

            # Compute errors and total population at end of simulation.
            for gregory_order_simulation in gregory_orders_simulation:
                print(f"Gregory order {gregory_order_simulation}")
                # Read results from IDE simulations.
                results = read_data(ide_result_dir, ide_exponents,
                                    gregory_order_simulation)

                # Compute errors of IDE results compared to groundtruth.
                errors_l2_abs = compute_errors_l2(
                    groundtruth, results, groundtruth_save_exponent, timesteps_ide, t0_ide, t_init, error_type="abs")
                errors_all_gregory_orders_l2_abs.append(errors_l2_abs)

                errors_max_abs = compute_errors_max(
                    groundtruth, results, groundtruth_save_exponent, timesteps_ide, t0_ide, t_init, error_type="abs")
                errors_all_gregory_orders_max_abs.append(errors_max_abs)

                errors_max_abs_rel = compute_errors_max(
                    groundtruth, results, groundtruth_save_exponent, timesteps_ide, t0_ide, t_init, error_type="rel")
                errors_all_gregory_orders_max_rel.append(errors_max_abs_rel)

                errors_max_abs_weighted = compute_errors_max(
                    groundtruth, results, groundtruth_save_exponent, timesteps_ide, t0_ide, t_init, error_type="weighted")
                errors_all_gregory_orders_max_weighted.append(
                    errors_max_abs_weighted)

                damping_time = get_dampingtime_ide_from_dir_name(dir_name)
                plot_difference_per_timestep(
                    groundtruth, results, groundtruth_save_exponent, timesteps_ide, t0_ide, t_init, gregory_order_simulation, plot_dir, damping_time=damping_time)

                print()
                # print(f"Gregory order {gregory_order_simulation}")
                # print("Errors: ")
                # print(errors[:, :])

                # # Determine order of convergence
                # order = compute_order_of_convergence(
                #     errors, timesteps_ide)

                # print(
                #     f"Orders of convergence: ")
                # print(order.T)

            # Plot convergence of all compartments separately.

            # L2 norm
            norm = "l2"
            # absolute error
            error_type = "abs"
            plot_convergence(errors_all_gregory_orders_l2_abs, timesteps_ide,
                             gregory_orders_simulation, norm, error_type, plot_dir, only_S)

            # max norm
            norm = "max"
            # absolute error
            error_type = "abs"
            plot_convergence(errors_all_gregory_orders_max_abs, timesteps_ide,
                             gregory_orders_simulation, norm, error_type, plot_dir, only_S)
            # relative error
            error_type = "rel"
            plot_convergence(errors_all_gregory_orders_max_rel, timesteps_ide,
                             gregory_orders_simulation, norm, error_type, plot_dir, only_S)
            # weighted error
            error_type = "weighted"
            plot_convergence(errors_all_gregory_orders_max_weighted, timesteps_ide,
                             gregory_orders_simulation, norm, error_type, plot_dir, only_S)


if __name__ == '__main__':
    main()
