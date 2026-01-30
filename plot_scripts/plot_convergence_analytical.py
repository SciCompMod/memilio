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


def define_groundtruth(timesteps, t0, tmax):

    model = 'ode'
    results = {model: []}

    for timestep in timesteps:
        groundtruth_at_timepoints = []
        for compartment in range(3):
            timepoints = np.linspace(t0, tmax, int((tmax-t0)/timestep)+1)
            if compartment == 0:
                groundtruth_at_timepoints.append([np.cosh(timepoint)
                                                  for timepoint in timepoints])  # np.cosh(timepoint)
            if compartment == 1:
                groundtruth_at_timepoints.append([
                    np.sinh(timepoint) for timepoint in timepoints])  # np.sinh(timepoint)
            if compartment == 2:
                groundtruth_at_timepoints.append([
                    np.cosh(timepoint) for timepoint in timepoints])

        results[model].append(groundtruth_at_timepoints)

    return results


def read_data(data_dir, ide_exponents, gregory_order):
    """ Read data into a dict, where the keys correspond to the respective model.
    At the moment we are only storing results of the IDE model here. There, we have an array that contains all results
    obtained with the IDE model for all time points for each time step size that is investigated. The results can
    either be compartments or flows as indicated by the flag 'flows'.
    """
    model = 'ide'
    results = {model: []}

    for exponent in ide_exponents:

        h5file = h5py.File(os.path.join(
            data_dir, f'result_{model}_dt=1e-{exponent:.0f}_gregoryorder={gregory_order}.h5'), 'r')

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


def compute_errors_l2(groundtruth, results, groundtruth_exponent, timesteps_ide, t0_ide, relative_error=True, cut_off=0):
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
    for i in range(len(results['ide'])):
        errors.append([])
        for compartment in range(num_errors):
            timestep = timesteps_ide[i]
            scale_timesteps = timestep/pow(10, -groundtruth_exponent)
            num_timepoints = len(results['ide'][i])

            result_ode = np.array(groundtruth['ode'][i][compartment][int(
                t0_ide/timestep)+cut_off::])
            result_ide = np.array(results['ide'][i][int(
                t0_ide/timestep)+cut_off::][:, compartment])

            difference = result_ode - result_ide

            if relative_error:
                norm_groundtruth = compute_l2_norm(
                    result_ode, timestep)
                errors[i].append(compute_l2_norm(
                    difference, timestep)/norm_groundtruth)

            else:
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


def compute_errors_max(groundtruth, results, groundtruth_exponent, timesteps_ide, t0_ide, relative_error=True, cut_off=0):
    """ Computes relative maximum norm of the difference between time series from ODE and time series
    from IDE for all compartments.
    """
    num_errors = 3

    errors = []

    # Compute error.
    for i in range(len(results['ide'])):
        errors.append([])
        for compartment in range(num_errors):
            timestep = timesteps_ide[i]
            scale_timesteps = timestep/pow(10, -groundtruth_exponent)

            difference = groundtruth['ode'][i][compartment][int(
                t0_ide/timestep)+cut_off::]-results['ide'][i][int(t0_ide/timestep)+cut_off::][:, compartment]

            if relative_error:
                norm_groundtruth = compute_max_norm(
                    groundtruth['ode'][i][int(t0_ide/timestep)::][compartment])
                errors[i].append(compute_max_norm(
                    difference)/norm_groundtruth)

            else:
                errors[i].append(compute_max_norm(
                    difference))

    return np.array(errors)


def plot_convergence(errors_all_gregory_orders, timesteps_ide,
                     gregory_orders_simulation, l2=True, maxnorm=False, relative_error=True, save_dir=""):
    """ Plots errors against timesteps with a subplot for each compartment /flow.

    @param[in] errors Array that contains computed errors of IDE model compared to groundtruth.
    @param[in] timesteps_ide List of time steps used in IDE simulations.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False.
    @param[in] save_dir Directory where plot will be stored. Default is an empty string leading to the plot not being
        saved.
    """
    # Define subplots and labels.
    num_plots = 3
    num_plotted_results = len(gregory_orders_simulation)

    fig, axs = plt.subplots(1, num_plots, sharex=True,
                            figsize=(10, 3))
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
        comparison = [0.5*errors_all_gregory_orders[1]
                      [0, i]*dt**3 for dt in plotted_timesteps]
        third = axs[i].plot(plotted_timesteps, comparison,
                            '--', color=colors[1], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^3)$")

        plotted_timesteps = timesteps_ide[:6]
        comparison = [0.5*errors_all_gregory_orders[2]
                      [0, i]*dt**4 for dt in plotted_timesteps]
        fourth = axs[i].plot(plotted_timesteps, comparison,
                             '--', color=colors[2], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^4)$")

        # Plot results.
        for j in range(len(gregory_orders_simulation)):
            line = axs[i].plot(timesteps_ide,
                               errors_all_gregory_orders[j][:, i], '-o', color=colors[j], label=labels[j])
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
        titles = [r"S(t)", r"S'(t)", r"S(0)+\int_0^t S'(x)\,dx"]
        axs[i].set_title(titles[i], fontsize=10)

        # Adapt plots.
        axs[i].set_xscale("log", base=10)
        axs[i].set_yscale("log", base=10)

        axs[i].grid(True, linestyle='--', alpha=0.6)

    fig.supxlabel(r'Time step $\Delta t$', fontsize=12)
    if relative_error:
        ylabel = fig.supylabel(
            r"$err_\text{rel}$", fontsize=12)
    else:
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

        if relative_error:
            filename = filename + "_rel"
        else:
            filename = filename + "_abs"

        plt.savefig(filename + ".png", format='png', bbox_extra_artists=(legend, ylabel), bbox_inches='tight',
                    dpi=500)

    plt.close()


def plot_difference_per_timestep(groundtruth, results, groundtruth_exponent, timesteps_ide, t0_ide, gregory_order, cut_off=0, save_dir=""):
    num_errors = 3

    errors = []

    compartments = ["S", "I", "R"]

    # Compute error.
    difference = []

    num_plots = 3
    figsize_x = 12

    for i, timestep in enumerate(timesteps_ide):
        fig, axs = plt.subplots(1, num_plots, sharex=True,
                                figsize=(figsize_x, 3))

        errors.append([])
        for compartment in range(num_errors):
            scale_timesteps = timestep/pow(10, -groundtruth_exponent)

            # difference = groundtruth[0][int(
            #     pow(10, groundtruth_exponent)*(t0_ide))::int(scale_timesteps)][:, compartment]-results[i][int(t0_ide/timestep)::][:, compartment]

            difference = groundtruth['ode'][i][compartment][int(
                t0_ide/timestep)+cut_off::]-results['ide'][i][int(t0_ide/timestep)+cut_off::][:, compartment]

            indices = range(len(difference))

            axs[compartment].scatter(indices[0:], difference[0:], s=1)
            axs[compartment].set_title(f"{compartments[compartment]}")

        plt.show()

        if save_dir != "":
            if not os.path.isdir(f"{save_dir}/differences"):
                os.makedirs(f"{save_dir}/differences")

            filename = f"{save_dir}/differences/gregoryorder={gregory_order}_timestep={timestep}"

            plt.savefig(filename + ".png", format='png',
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


def get_tmax_from_dir_name(dir_name):
    tmax_string = [x for x in dir_name.split("_") if "tmax" in x]
    tmax = int(tmax_string[0].split("=")[-1])

    return tmax


def main():

    groundtruth_exponent = 6

    root_dir = os.path.join(os.path.dirname(
        __file__), "../simulation_results")
    main_dir = "2026-01-27/groundtruth=cosh_S_computed_no_init/"
    relevant_dir = os.path.join(root_dir, main_dir)
    # print(relevant_dir)
    sub_dirs = subfolders_scandir(relevant_dir)
    # sub_dirs = ["t0ide=0_tmax=5"]

    gregory_orders_simulation = [1, 2, 3]

    cut_off = 10

    for dir_index, dir_name in enumerate(sub_dirs):
        print(dir_name)

        t0 = get_t0_ide_from_dir_name(dir_name)
        tmax = get_tmax_from_dir_name(dir_name)

        # Path where simulation results are stored.
        result_dir = os.path.join(os.path.dirname(
            __file__),  f"{relevant_dir}/{dir_name}/")

        # Path where plots will be stored.
        plot_dir = os.path.join(os.path.dirname(
            __file__),  f"../plots/{main_dir}/{dir_name}/cutoff={cut_off}/")

        # errors_all_gregory_orders_l2_rel = []
        errors_all_gregory_orders_l2_abs = []

        # errors_all_gregory_orders_max_rel = []
        errors_all_gregory_orders_max_abs = []

        # Get exponents for which IDE simulations have been computed for considered subdirectory.
        ide_exponents = get_ide_exponents(result_dir)
        # ide_exponents = [0, 1, 2, 3, 4]

        # Calculate time steps resulting from ide_exponents.
        timesteps_ide = []
        for exp in ide_exponents:
            timesteps_ide.append(pow(10, -exp))

        groundtruth = define_groundtruth(timesteps_ide, t0, tmax)

        # Compute errors and total population at end of simulation.
        for gregory_order_simulation in gregory_orders_simulation:
            # Read results from IDE simulations.
            results = read_data(result_dir, ide_exponents,
                                gregory_order_simulation)

            # Compute errors of IDE results compared to groundtruth.
            # errors_l2_rel = compute_errors_l2(
            #     groundtruth, results, groundtruth_exponent, timesteps_ide, t0, True)
            errors_l2_abs = compute_errors_l2(
                groundtruth, results, groundtruth_exponent, timesteps_ide, t0, False, cut_off)
            # errors_all_gregory_orders_l2_rel.append(errors_l2_rel)
            errors_all_gregory_orders_l2_abs.append(errors_l2_abs)

            # errors_max_rel = compute_errors_max(
            #     groundtruth, results, groundtruth_exponent, timesteps_ide, t0, True)
            errors_max_abs = compute_errors_max(
                groundtruth, results, groundtruth_exponent, timesteps_ide, t0, False, cut_off)
            # errors_all_gregory_orders_max_rel.append(errors_max_rel)
            errors_all_gregory_orders_max_abs.append(errors_max_abs)

            plot_difference_per_timestep(
                groundtruth, results, groundtruth_exponent, timesteps_ide, t0, gregory_order_simulation, cut_off, plot_dir)

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
        # computed_errors_rel = [
        #     errors_all_gregory_orders_l2_rel, errors_all_gregory_orders_max_rel]
        # computed_errors_abs = [
        #     errors_all_gregory_orders_l2_abs, errors_all_gregory_orders_max_abs]

        # computed_errors = [computed_errors_rel, computed_errors_abs]

        # for i in range(len(computed_errors)):
        #     if i == 0:
        #         relative_error = True
        #     else:
        relative_error = False

        # L2 norm
        l2 = True
        maxnorm = False
        plot_convergence(errors_all_gregory_orders_l2_abs, timesteps_ide,
                         gregory_orders_simulation, l2, maxnorm,  relative_error, plot_dir)

        # max norm
        l2 = False
        maxnorm = True
        plot_convergence(errors_all_gregory_orders_max_abs, timesteps_ide,
                         gregory_orders_simulation,  l2, maxnorm, relative_error, plot_dir)


if __name__ == '__main__':
    main()
