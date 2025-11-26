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


def read_groundtruth_ode(data_dir, groundtruth_exponent):
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
        data_dir, f'result_{model}_dt=1e-{groundtruth_exponent:.0f}.h5'), 'r')

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


def read_groundtruth(result_dir, groundtruth_exponent, gregory_order_groundtruth=3, groundtruth_ode=True):

    if not groundtruth_ode:
        return read_groundtruth_ide(
            result_dir, groundtruth_exponent, gregory_order_groundtruth)

    return read_groundtruth_ode(
        result_dir, groundtruth_exponent)


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


def compute_errors_l2(groundtruth, results, groundtruth_exponent, timesteps_ide, t0_ide, relative_error=True):
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
            scale_timesteps = timestep/pow(10, -groundtruth_exponent)
            num_timepoints = len(results[i])

            difference = groundtruth[0][int(
                pow(10, groundtruth_exponent)*(t0_ide))::int(scale_timesteps)][:, compartment]-results[i][int(t0_ide/timestep)::][:, compartment]

            if relative_error:
                norm_groundtruth = compute_l2_norm(groundtruth[0][int(
                    pow(10, groundtruth_exponent)*(t0_ide))::int(scale_timesteps)][:, compartment], timestep)
                errors[i].append(compute_l2_norm(
                    difference, timestep)/norm_groundtruth)

            else:
                errors[i].append(compute_l2_norm(
                    difference, timestep))

    return np.array(errors)


def compute_errors_l2_sumIR(groundtruth, results, groundtruth_exponent, timesteps_ide, t0_ide, relative_error=True):
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
            scale_timesteps = timestep/pow(10, -groundtruth_exponent)

            if compartment == 1:

                sum_ide = np.abs(results[i][int(
                    t0_ide/timestep)::][:, 1] + results[i][int(t0_ide/timestep)::][:, 2])

                sum_ode = np.abs(groundtruth[0][int(
                    pow(10, groundtruth_exponent)*(t0_ide))::int(scale_timesteps)][:, 1]) + np.abs(groundtruth[0][int(
                        pow(10, groundtruth_exponent)*(t0_ide))::int(scale_timesteps)][:, 2])

                difference = sum_ode - sum_ide

                if relative_error:

                    errors[i].append(compute_l2_norm(
                        difference, timestep)/compute_l2_norm(sum_ode, timestep))

                else:
                    errors[i].append(compute_l2_norm(
                        difference, timestep))
            else:
                errors[i].append(0)

    return np.array(errors)


def compute_max_norm(timeseries):
    """ Computes maximum norm of a time series.

    @param[in] timeseries Considered timeseries.
    @returns Norm.
    """
    # print(timeseries)
    norm = np.max(np.abs(timeseries))
    return norm


def compute_errors_max(groundtruth, results, groundtruth_exponent, timesteps_ide, t0_ide, relative_error=True):
    """ Computes relative maximum norm of the difference between time series from ODE and time series
    from IDE for all compartments.
    """
    num_errors = 3

    errors = []

    # Compute error.

    for i in range(len(results)):
        if not relative_error:
            print(f"Timestep {timesteps_ide[i]}")
        errors.append([])
        for compartment in range(num_errors):
            timestep = timesteps_ide[i]
            scale_timesteps = timestep/pow(10, -groundtruth_exponent)

            difference = groundtruth[0][int(
                pow(10, groundtruth_exponent)*(t0_ide))::int(scale_timesteps)][:, compartment]-results[i][int(t0_ide/timestep)::][:, compartment]

            if relative_error:
                norm_groundtruth = compute_max_norm(groundtruth[0][int(
                    pow(10, groundtruth_exponent)*(t0_ide))::int(scale_timesteps)][:, compartment])
                errors[i].append(compute_max_norm(
                    difference)/norm_groundtruth)

            else:
                print(
                    f"Error is highest at index {np.argmax(difference)} out of {len(difference)-1}")
                errors[i].append(compute_max_norm(
                    difference))

    return np.array(errors)


def plot_convergence(errors_all_gregory_orders, timesteps_ide,
                     gregory_orders_simulation, fd_order=1, l2=True, maxnorm=False, norm_of_sum=False, relative_error=True, save_dir="", only_S=False):
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
        figsize_x = 12

    num_plotted_results = len(gregory_orders_simulation)

    fig, axs = plt.subplots(1, num_plots, sharex=True,
                            figsize=(figsize_x, 3))
    secir_dict = {0: 'Susceptible', 1:  'Infected', 2:  'Recovered'}
    labels = [
        f"Gregory order {gregory_order}" for gregory_order in gregory_orders_simulation]
    labels.insert(0, "")
    labels.append(r"$\mathcal{O}(\Delta t)$")
    labels.append(r"$\mathcal{O}(\Delta t^2)$")
    labels.append(r"$\mathcal{O}(\Delta t^3)$")
    labels.append(r"$\mathcal{O}(\Delta t^4)$")

    handles = [plt.Line2D([], [], color='none')]

    # Define colors.
    colors_ = [plt.cm.viridis(x)
               for x in np.linspace(0, 1, num_plotted_results)]
    colors = ["darkorange", colors_[1], "darkred"]

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
        comparison = [0.5*10*errors_all_gregory_orders[0]
                      [0, i]*dt for dt in timesteps_ide]
        first = ax_obj.plot(plotted_timesteps, comparison,
                            '--', color='gray', linewidth=1.2, label=r"$\mathcal{O}(\Delta t)$")

        plotted_timesteps = timesteps_ide[:6]
        comparison = [0.5*100*errors_all_gregory_orders[0]
                      [0, i]*dt**2 for dt in plotted_timesteps]
        second = ax_obj.plot(plotted_timesteps, comparison,
                             '--', color=colors[0], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^2)$")

        plotted_timesteps = timesteps_ide[:6]
        comparison = [0.5*1000*errors_all_gregory_orders[1]
                      [0, i]*dt**3 for dt in plotted_timesteps]
        third = ax_obj.plot(plotted_timesteps, comparison,
                            '--', color=colors[1], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^3)$")

        plotted_timesteps = timesteps_ide[:6]
        comparison = [0.5*10000*errors_all_gregory_orders[2]
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

        if l2:
            filename = f'{save_dir}/convergence_all_compartments_l2'
        elif maxnorm:
            filename = f'{save_dir}/convergence_all_compartments_max'
        elif norm_of_sum:
            filename = f'{save_dir}/convergence_all_compartments_sumIR'
        else:
            filename = f'{save_dir}/convergence_all_compartments_endpoint'

        if relative_error:
            filename = filename + "_rel"
        else:
            filename = filename + "_abs"

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


def subfolders_scandir(path):
    # path = os.path.dirname(path)
    print(path)
    with os.scandir(path) as it:
        return [entry.name for entry in it if entry.is_dir()]


def get_ide_exponents(data_dir):
    files = os.listdir(data_dir)

    max_exponent = 0
    for possible_exponent in range(5):
        if f'result_{"ide"}_dt=1e-{possible_exponent}_gregoryorder=3.h5' in files:
            max_exponent = possible_exponent

    return range(max_exponent+1)


def get_t0_ide_from_dir_name(dir_name):
    t0_string = [x for x in dir_name.split("_") if "t0ide" in x]
    t0 = int(t0_string[0].split("=")[-1])

    return t0


def main():

    groundtruth_exponent = 6
    groundtruth_ode = False
    only_S = True

    main_dir = "2025-11-17/analytical_example"

    ##############################################

    root_dir = os.path.join(os.path.dirname(
        __file__), "../simulation_results")
    relevant_dir = os.path.join(root_dir, main_dir)
    print(relevant_dir)
    sub_dirs = subfolders_scandir(relevant_dir)

    total_pop_reference = 0
    total_pop_all_fd_orders = []

    gregory_orders_simulation = [1, 2, 3]

    for dir_index, dir_name in enumerate(sub_dirs):
        print(dir_name)

        t0_ide = get_t0_ide_from_dir_name(dir_name)

        # Path where simulation results are stored.
        result_dir = os.path.join(os.path.dirname(
            __file__),  f"{relevant_dir}/{dir_name}/")

        # Path where plots will be stored.
        plot_dir = os.path.join(os.path.dirname(
            __file__),  f"../plots/{main_dir}/{dir_name}/")

        # # Read groundtruth.
        groundtruth = read_groundtruth(
            result_dir, groundtruth_exponent, groundtruth_ode=groundtruth_ode)

        errors_all_gregory_orders_endpoint_rel = []
        errors_all_gregory_orders_endpoint_abs = []

        errors_all_gregory_orders_l2_rel = []
        errors_all_gregory_orders_l2_abs = []

        errors_all_gregory_orders_max_rel = []
        errors_all_gregory_orders_max_abs = []

        total_pop_end_all_gregory_orders = []

        # Get exponents for which IDE simulations have been computed for considered subdirectory.
        # ide_exponents = get_ide_exponents(result_dir)
        ide_exponents = [1, 2, 3, 4]

        # Calculate time steps resulting from ide_exponents.
        timesteps_ide = []
        for exp in ide_exponents:
            timesteps_ide.append(pow(10, -exp))

        # Compute errors and total population at end of simulation.
        for gregory_order_simulation in gregory_orders_simulation:
            print(f"Gregory order {gregory_order_simulation}")
            # Read results from IDE simulations.
            results = read_data(result_dir, ide_exponents,
                                gregory_order_simulation)

            # Compute errors of IDE results compared to groundtruth.
            errors_tmax_rel = compute_errors_endpoint(
                groundtruth, results, True)
            errors_tmax_abs = compute_errors_endpoint(
                groundtruth, results, False)
            errors_all_gregory_orders_endpoint_rel.append(errors_tmax_rel)
            errors_all_gregory_orders_endpoint_abs.append(errors_tmax_abs)

            errors_l2_rel = compute_errors_l2(
                groundtruth, results, groundtruth_exponent, timesteps_ide, t0_ide, True)
            errors_l2_abs = compute_errors_l2(
                groundtruth, results, groundtruth_exponent, timesteps_ide, t0_ide, False)
            errors_all_gregory_orders_l2_rel.append(errors_l2_rel)
            errors_all_gregory_orders_l2_abs.append(errors_l2_abs)

            errors_max_rel = compute_errors_max(
                groundtruth, results, groundtruth_exponent, timesteps_ide, t0_ide, True)
            errors_max_abs = compute_errors_max(
                groundtruth, results, groundtruth_exponent, timesteps_ide, t0_ide, False)
            errors_all_gregory_orders_max_rel.append(errors_max_rel)
            errors_all_gregory_orders_max_abs.append(errors_max_abs)

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

            # print(
            #     f"Total population at end for time step {timesteps_ide[-1]}: {results[-1][-1].sum()}")

            total_pop_end = get_total_pop_end(results)
            total_pop_end_all_gregory_orders.append(total_pop_end)

            total_pop_reference = results[-1][0].sum()

        print("Max norm")
        for i in range(len(errors_all_gregory_orders_max_abs)):
            print(errors_all_gregory_orders_max_abs[i])
        total_pop_all_fd_orders.append(total_pop_end_all_gregory_orders)

        # Plot convergence of all compartments separately.
        fd_order = 1  # dummy right now

        computed_errors_rel = [errors_all_gregory_orders_endpoint_rel,
                               errors_all_gregory_orders_l2_rel, errors_all_gregory_orders_max_rel]

        computed_errors_abs = [errors_all_gregory_orders_endpoint_abs,
                               errors_all_gregory_orders_l2_abs, errors_all_gregory_orders_max_abs]

        computed_errors = [computed_errors_rel,
                           computed_errors_abs]  # computed_errors_rel,

        for i in range(len(computed_errors)):
            if i == 0:
                relative_error = True
            else:
                relative_error = False

            # endpoint norm
            l2 = False
            maxnorm = False
            norm_of_sum = False
            plot_convergence(computed_errors[i][0], timesteps_ide,
                             gregory_orders_simulation, fd_order, l2, maxnorm, norm_of_sum, relative_error, plot_dir, only_S)

            # L2 norm
            l2 = True
            maxnorm = False
            norm_of_sum = False
            plot_convergence(computed_errors[i][1], timesteps_ide,
                             gregory_orders_simulation, fd_order, l2, maxnorm, norm_of_sum,  relative_error, plot_dir, only_S)

            # max norm
            l2 = False
            maxnorm = True
            norm_of_sum = False
            plot_convergence(computed_errors[i][2], timesteps_ide,
                             gregory_orders_simulation, fd_order, l2, maxnorm, norm_of_sum, relative_error, plot_dir, only_S)

    # Path where plots will be stored.
    # plot_dir = os.path.join(os.path.dirname(
    #     __file__),  f"../plots/{main_dir}/")
    # plot_total_pop_diff(gregory_orders_simulation, fd_orders, timesteps_ide,
    #                     total_pop_all_fd_orders, total_pop_reference, plot_dir)


if __name__ == '__main__':
    main()
