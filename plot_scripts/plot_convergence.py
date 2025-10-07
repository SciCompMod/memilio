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
    results = {model: []}

    h5file = h5py.File(os.path.join(
        data_dir, f'result_{model}_dt=1e-{groundtruth_exponent:.0f}_gregoryorder={gregory_order}.h5'), 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]

    if len(data['Total'][0]) == 3:
        # As there should be only one Group, total is the simulation result
        results[model].append(data['Total'][:, :])
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
    results = {model: []}

    h5file = h5py.File(os.path.join(
        data_dir, f'result_{model}_dt=1e-{groundtruth_exponent:.0f}.h5'), 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]

    if len(data['Total'][0]) == 3:
        # As there should be only one Group, total is the simulation result
        results[model].append(data['Total'][:, :])
    else:
        raise gd.DataError(
            'Expected a different size of vector in time series.')

    dates = data['Time'][:]

    h5file.close()

    return results


def read_groundtruth(result_dir, groundtruth_exponent, gregory_order_groundtruth=3, groundtruth_ode=True):

    # if "ode" in result_dir.split("/")[-2]:
    #     return read_groundtruth_ode(result_dir, groundtruth_exponent)

    # if "ide" in result_dir.split("/")[-2]:
    if not groundtruth_ode:
        return read_groundtruth_ide(
            result_dir, groundtruth_exponent, gregory_order_groundtruth)

    return read_groundtruth_ode(
        result_dir, groundtruth_exponent)


def read_data(data_dir, exponents_ide, gregory_order):
    """ Read data into a dict, where the keys correspond to the respective model.
    At the moment we are only storing results of the IDE model here. There, we have an array that contains all results
    obtained with the IDE model for all time points for each time step size that is investigated. The results can
    either be compartments or flows as indicated by the flag 'flows'.

    @param[in] data_dir Directory where h5 files are stored.
    @param[in] ode_exponent Exponent that determines time step size of ODE simulation via dt =10^{-ode_exponent}.
    @param[in] exponents_ide List of considered exponents that determine time step size of IDE simulation via
    dt =10^{-exponent_ide}.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False.
    @returns Dict with results of ODE model per time step size.
    """
    models = ['ide']
    results = {models[0]: []}
    for model in models:
        for exponent in exponents_ide:

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


def compute_errors(groundtruth, results, groundtruth_exponent, timesteps_ide, gregory_order):
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
    for i in range(len(results['ide'])):
        errors.append([])

        for compartment in range(num_errors):

            model = list(groundtruth.keys())[0]

            difference = np.abs(
                groundtruth[model][0][-1, compartment]-results['ide'][i][-1, compartment])

            errors[i].append(difference)

    return np.array(errors)


def plot_convergence(errors_all_gregory_orders, timesteps_ide,
                     gregory_orders_simulation, save_dir=""):
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

    fig, axs = plt.subplots(1, num_plots, sharex=True, figsize=(9, 3))
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
    colors = [plt.cm.viridis(x)
              for x in np.linspace(0, 1, num_plotted_results)]

    for i in range(num_plots):
        for j in range(len(gregory_orders_simulation)):
            # Plot results.
            if i == 0:
                line = axs[i].plot(timesteps_ide,
                                   errors_all_gregory_orders[j][:, i], '-o', color=colors[j], label=labels[j])
                handles.append(line[0])
            else:
                line = axs[i].plot(timesteps_ide,
                                   errors_all_gregory_orders[j][:, i], '-o', color=colors[j])
                # handles.append(line)

        # Plot comparison line for linear convergence as well as second, third and fourth order.
        # Susceptibles
        if i == 0:
            num_timesteps_ide = len(timesteps_ide)

            plotted_timesteps = timesteps_ide[:4]
            comparison = [0.1*dt**2 for dt in plotted_timesteps]
            second = axs[i].plot(plotted_timesteps, comparison,
                                 '--', color=colors[0], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^2)$")
            handles.append(second[0])

            plotted_timesteps = timesteps_ide[:3]
            comparison = [0.04*dt**3 for dt in plotted_timesteps]
            third = axs[i].plot(plotted_timesteps, comparison,
                                '--', color=colors[1], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^3)$")
            handles.append(third[0])

            plotted_timesteps = timesteps_ide[:2]
            comparison = [0.007*dt**4 for dt in plotted_timesteps]
            fourth = axs[i].plot(plotted_timesteps, comparison,
                                 '--', color=colors[2], linewidth=1.2, alpha=0.5, label=r"$\mathcal{O}(\Delta t^4)$")
            handles.append(fourth[0])

        # Infected
        if i == 1:
            num_timesteps_ide = len(timesteps_ide)

            plotted_timesteps = timesteps_ide[:4]
            comparison = [0.1*dt**2 for dt in plotted_timesteps]
            axs[i].plot(plotted_timesteps, comparison,
                        '--', color=colors[0], linewidth=1.2, alpha=0.5)

            plotted_timesteps = timesteps_ide[:3]
            comparison = [0.03*dt**3 for dt in plotted_timesteps]
            axs[i].plot(plotted_timesteps, comparison,
                        '--', color=colors[1], linewidth=1.2, alpha=0.5)

            plotted_timesteps = timesteps_ide[:2]
            comparison = [0.005*dt**4 for dt in plotted_timesteps]
            axs[i].plot(plotted_timesteps, comparison,
                        '--', color=colors[2], linewidth=1.2, alpha=0.5)

        # Recovered
        if i == 2:
            comparison = [2*dt for dt in timesteps_ide]
            linear = axs[i].plot(timesteps_ide, comparison,
                                 '--', color='gray', linewidth=1.2, label=r"$\mathcal{O}(\Delta t)$")
            handles.insert(4, linear[0])

        # Adapt plots.
        axs[i].set_xscale("log", base=10)
        axs[i].set_yscale("log", base=10)

        axs[i].set_title(secir_dict[i], fontsize=10)
        axs[i].grid(True, linestyle='--', alpha=0.6)

        # if i != 2:
        #     axs[i].set_ylim(1e-8, 5*1e-1)

    fig.supxlabel(r'Time step $\Delta t$', fontsize=12)
    fig.supylabel(
        r"$\Vert \widehat{Z}_{\text{IDE}} - \widehat{Z}_{\text{ODE}}\Vert_{2}$", fontsize=12)

    # Invert x axis only for one plot so that sharex=True and invert_xaxis work as intended.
    axs[0].invert_xaxis()

    # print(handles)

    legend = fig.legend(handles=handles, labels=labels, ncol=2,  loc='lower right',
                        fontsize=8, bbox_transform=fig.transFigure, bbox_to_anchor=(1., -0.1))  # bbox_to_anchor=(0.5, -0.06),
    plt.tight_layout()
    # plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.1)
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/convergence_all_compartments.png', format='png', bbox_extra_artists=(legend,), bbox_inches='tight',
                    dpi=500)


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
    for timestep in range(len(results['ide'])):
        total_pop_end.append(results['ide'][timestep][-1].sum())

    return total_pop_end


def plot_total_pop_diff(gregory_orders_simulation, timesteps_ide, total_pop_end_all_gregory_orders, total_pop_reference, save_dir=""):

    # Plot relative difference per time step for all Gregory orders for a fixed finite_difference_order.

    fig, axs = plt.subplots(1, 1, sharex=True, figsize=(9, 3))
    labels = [
        f"Gregory order {gregory_order}" for gregory_order in gregory_orders_simulation]

    # Define colors.
    colors = [plt.cm.viridis(x)
              for x in np.linspace(0, 1, len(gregory_orders_simulation))]

    for j in range(len(gregory_orders_simulation)):
        # Plot results.
        line = axs.plot(timesteps_ide,
                        1 - total_pop_end_all_gregory_orders[j][:]/total_pop_reference, '-o', color=colors[j], label=labels[j])
        # handles.append(line[0])

    axs.invert_xaxis()
    axs.set_xscale("log", base=10)
    axs.set_yscale("log", base=10)

    fig.legend(bbox_to_anchor=(0.97, 0.9))

    fig.supxlabel(r'Time step $\Delta t$', fontsize=12)
    fig.supylabel(
        r"Relative deviation", fontsize=12)

    plt.tight_layout()

    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/mass_conservation_diff.png', format='png',  bbox_inches='tight',
                    dpi=500)


def main():

    groundtruth_exponent = 6
    gregory_order_groundtruth = 3

    finite_difference_order = 4

    dir_name = "time_infected=2/detailed_init_exponential_t0ide=20_tmax=30_finite_diff=2"
    print(dir_name)

    groundtruth_ode = True

    gregory_orders_simulation = [1, 2, 3]
    ide_exponents = [0, 1, 2, 3, 4]

    ###################################################################################################################

    # Path where simulation results (generated with ide_convergence_rate.cpp) are stored.
    result_dir = os.path.join(os.path.dirname(
        __file__),  f"../simulation_results/{dir_name}/")

    # Path where plots will be stored.
    plot_dir = os.path.join(os.path.dirname(
        __file__),  f"../plots/{dir_name}/")

    # The IDE model was simulated using a fixed step size dt=10^{-ide_exponent} for ide_exponent in ide_exponents.

    # Calculate time steps resulting from exponents_ide.
    timesteps_ide = []
    for exp in ide_exponents:
        timesteps_ide.append(pow(10, -exp))

    # # Read groundtruth.
    groundtruth = read_groundtruth(
        result_dir, groundtruth_exponent, gregory_order_groundtruth, groundtruth_ode)

    errors_all_gregory_orders = []

    total_pop_reference = 0
    total_pop_end_all_gregory_orders = []

    for gregory_order_simulation in gregory_orders_simulation:
        # Read results from IDE simulations.
        results = read_data(result_dir, ide_exponents,
                            gregory_order_simulation)

        # Compute errors of IDE results compared to groundtruth.
        errors = compute_errors(
            groundtruth, results, groundtruth_exponent, timesteps_ide, gregory_order_simulation)

        errors_all_gregory_orders.append(errors)

        print()
        print(f"Gregory order {gregory_order_simulation}")
        # print("Errors: ")
        # print(errors[:, :])

        # Determine order of convergence
        order = compute_order_of_convergence(
            errors, timesteps_ide)

        print(
            f"Orders of convergence: ")
        print(order.T)

        # print(
        #     f"Total population at end for time step {timesteps_ide[-1]}: {results['ide'][-1][-1].sum()}")

        total_pop_end = get_total_pop_end(results)
        total_pop_end_all_gregory_orders.append(total_pop_end)

        total_pop_reference = results['ide'][-1][0].sum()

    # Plot convergence of all compartments separately.
    plot_convergence(errors_all_gregory_orders, timesteps_ide,
                     gregory_orders_simulation, plot_dir)

    plot_total_pop_diff(gregory_orders_simulation, timesteps_ide,
                        total_pop_end_all_gregory_orders, total_pop_reference, plot_dir)


if __name__ == '__main__':
    main()
