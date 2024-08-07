#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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


def read_groundtruth(data_dir, dt_ode, setting, flows=False):
    """ Read groundtruth from data- We define the groundtruth as the results obtained by the ODE model with timestep dt=1e-6."""

    model = 'ode'
    results = {model: []}

    if flows:
        h5file = h5py.File(os.path.join(
            data_dir, f'result_{model}_flows_dt={dt_ode:.4f}_setting{setting}.h5'), 'r')
    else:
        h5file = h5py.File(os.path.join(
            data_dir, f'result_{model}_dt={dt_ode:.4f}_setting{setting}.h5'), 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]

    if flows:
        # For flows we divide the values by dt_ode to go from \tilde{\sigma} to \hat{\sigma} and make results comparable
        results[model].append(data['Total'][:, :]/dt_ode)
    else:
        if len(data['Total'][0]) == 8:
            # As there should be only one Group, total is the simulation result
            results[model].append(data['Total'][:, :])
        elif len(data['Total'][0]) == 10:
            # in ODE there are two compartments we don't use, throw these out
            results[model].append(
                data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]])
        else:
            raise gd.DataError(
                'Expected a different size of vector in time series.')

    dates = data['Time'][:]

    h5file.close()

    return results


def read_data(data_dir, dt_ode, timesteps_ide, setting, flows=False):
    """ Read data into a dict, where the keys correspond to the respective model.
    At the moment we are only storing results of the IDE model here. There
    we have an array that contains all results for SECIHURD for all time points
    for each time step size that is investigated."""

    models = ['ide']
    results = {models[0]: []}
    for model in models:
        for timestep in timesteps_ide:
            if flows:
                h5file = h5py.File(os.path.join(
                    data_dir, f'result_{model}_flows_dt={timestep:.4f}_init_dt_ode={dt_ode:.4f}_setting{setting}.h5'), 'r')
            else:
                h5file = h5py.File(os.path.join(
                    data_dir, f'result_{model}_dt={timestep:.4f}_init_dt_ode={dt_ode:.4f}_setting{setting}.h5'), 'r')

            data = h5file[list(h5file.keys())[0]]

            if flows:
                # For flows we divide the values by dt_ode to go from \tilde{\sigma} to \hat{\sigma} and make results comparable
                results[model].append(data['Total'][:, :]/timestep)
            else:
                if len(data['Total'][0]) == 8:
                    # As there should be only one Group, total is the simulation result
                    results[model].append(data['Total'][:, :])
                elif len(data['Total'][0]) == 10:
                    # in ODE there are two compartments we don't use, throw these out
                    results[model].append(
                        data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]])
                else:
                    raise gd.DataError(
                        "Expected a different size of vector in time series.")

            h5file.close()

    return results


def compute_l2_norm(timeseries, timestep):
    """ Compute norm of a time series."""

    norm = np.sqrt(timestep * np.sum(timeseries**2))

    return norm


def compute_error_norm_l2(groundtruth, results, dt_ode, timesteps_ide, flows=False):
    """ Compute norm of the difference between time series from ODE and time series
    from IDE for all compartments/flows."""

    if flows:
        num_errors = 10
    else:
        num_errors = 8
    errors = []

    # Compute error.
    for i in range(len(results['ide'])):
        errors.append([])
        for compartment in range(num_errors):
            timestep = timesteps_ide[i]
            scale_timesteps = timestep/dt_ode
            num_timepoints = len(results['ide'][i])
            # Only consider compartments of ODE model for t>=t0_IDE. Use that 2*t0_IDE = t_max.
            if flows:
                difference = groundtruth['ode'][0][int(scale_timesteps*(num_timepoints/2-1))::int(
                    scale_timesteps)][:, compartment]-results['ide'][i][int(num_timepoints/2-1):][:, compartment]
            # Only consider compartments of ODE model for t>=t0_IDE.
            else:
                difference = groundtruth['ode'][0][int(scale_timesteps*(num_timepoints -
                                                                        1))::int(scale_timesteps)][:, compartment]-results['ide'][i][:, compartment]
            errors[i].append(compute_l2_norm(difference, timestep))

    return np.array(errors)


def compute_error_norm_tmax(groundtruth, results, timesteps_ide, flows=False):
    """ Compute norm of the difference between ODE and IDE at tmax for all compartments/flows."""

    if flows:
        num_errors = 10
    else:
        num_errors = 8
    errors = []

    # Compute error in IDE simulation at tmax.
    for i in range(len(results['ide'])):
        errors.append([])
        # Compute difference for all compartments.
        for compartment in range(num_errors):
            timestep = timesteps_ide[i]
            num_timepoints = len(results['ide'][i])
            difference = groundtruth['ode'][0][-1][compartment] - \
                results['ide'][i][-1][compartment]
            errors[i].append(np.sqrt(difference ** 2))

    return np.array(errors)


def plot_convergence(errors, timesteps_ide, setting, flows=False, compartment=None, save=False):
    """ Plot errors against timesteps."""

    if flows:
        secir_dict = {0: r"$\sigma_S^E$", 1: r"$\sigma_E^C$", 2: r"$\sigma_C^I$", 3: r"$\sigma_C^R$", 4: r"$\sigma_I^H$",
                      5: r"$\sigma_I^R$", 6: r"$\sigma_H^U$", 7: r"$\sigma_H^R$", 8: r"$\sigma_U^D$", 9: r"$\sigma_U^R$"}

        compartments = [r"$\sigma_S^E$", r"$\sigma_E^C$", r"$\sigma_C^I$", r"$\sigma_C^R$", r"$\sigma_I^H$",
                        r"$\sigma_I^R$", r"$\sigma_H^U$", r"$\sigma_H^R$", r"$\sigma_U^D$", r"$\sigma_U^R$"]
    else:
        secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                      5: 'ICU', 6: 'Recovered', 7: 'Dead'}
        compartments = ['S', 'E', 'C', 'I', 'H', 'U', 'R', 'D']

    # helmholtzdarkblue, helmholtzclaim
    colors = [(0, 40/255, 100/255), (20/255, 200/255, 255/255)]

    # Plot only one compartment.
    if compartment != None:

        print('Need to uncomment')

        # # TODO: include check if compartment is in 0,...,8

        # fig, ax = plt.subplots()

        # ax.plot(timesteps, errors[:, compartment],
        #         '-o', color=colors[0], label='Results')
        # comparison = [1800 * dt for dt in timesteps]
        # ax.plot(timesteps, comparison, color='lightgray',
        #         label=r"$\mathcal{O}(\Delta t)$")

        # ax.set_xscale("log", base=10)
        # ax.set_yscale("log", base=10)
        # ax.invert_xaxis()

        # fig.supxlabel('Time step')
        # fig.supylabel(r"$\Vert {compartment}_{IDE}(t_{max}) - {compartment}_{ODE}(t_{max})\Vert$".replace(
        #     'compartment', compartments[compartment]))

        # plt.legend()

        # if save:
        #     if not os.path.isdir('plots'):
        #         os.makedirs('plots')
        #     plt.savefig('plots/convergence_{compartments[compartment]}_setting{}.png'.format(setting),
        #                 bbox_inches='tight', dpi=500)

        # else:
        #     plt.show()

    # Plot all compartments.
    else:
        if flows:
            fig, axs = plt.subplots(5, 2, sharex=True, figsize=(10, 10))
            num_plots = 10
        else:
            fig, axs = plt.subplots(4, 2, sharex=True, figsize=(10, 8))
            num_plots = 8

        for i in range(num_plots):

            # Plot results.
            axs[int(i/2), i % 2].plot(timesteps_ide,
                                      errors[:, i], '-o', color=colors[1], label='Results')

            # Plot comparison line for linear convergence.
            factor = 50 * errors[0, i]
            comparison = [factor * dt for dt in timesteps_ide]
            axs[int(i/2), i % 2].plot(timesteps_ide, comparison, color='gray',
                                      label=r"$\mathcal{O}(\Delta t)$")

            # Adapt plots.
            axs[int(i/2), i % 2].set_xscale("log", base=10)
            axs[int(i/2), i % 2].set_yscale("log", base=10)

            axs[int(i/2), i % 2].set_title(secir_dict[i], fontsize=10)

        fig.supxlabel('Time step')
        fig.supylabel(
            r"$\Vert {Z}_{IDE}(t_{max}) - {Z}_{ODE}(t_{max})\Vert$")

        # Invert x axis only for one plot so that sharex=True and invert_xaxis work as intended.
        axs[0, 0].invert_xaxis()

        labels = ['Results', r"$\mathcal{O}(\Delta t)$", ]
        fig.legend(labels, bbox_to_anchor=(0.1, -0.73, 0.8, 0.8),
                   fancybox=False, shadow=False, ncol=1)

        # plt.tight_layout(rect=[0, 0, 0.83, 1])

        if save:
            if flows:
                if not os.path.isdir('plots/flows'):
                    os.makedirs('plots/flows')
                plt.savefig(f'plots/flows/convergence_all_flows_setting={setting}.png', format='png',
                            dpi=500)
            else:
                if not os.path.isdir('plots/compartments'):
                    os.makedirs('plots/compartments')
                plt.savefig(f'plots/compartments/convergence_all_compartments_setting={setting}.png', format='png',
                            dpi=500)


def plot_convergence_oneplot(errors, timesteps_ide, setting, flows=False, save=False):
    """ Plot errors against timesteps. This function creates one plot to depict all compartments/flows, respectively."""

    fig, ax = plt.subplots()

    if flows:
        secir_dict = {0: r"$\sigma_S^E$", 1: r"$\sigma_E^C$", 2: r"$\sigma_C^I$", 3: r"$\sigma_C^R$", 4: r"$\sigma_I^H$",
                      5: r"$\sigma_I^R$", 6: r"$\sigma_H^U$", 7: r"$\sigma_H^R$", 8: r"$\sigma_U^D$", 9: r"$\sigma_U^R$"}
        fig.supylabel(
            r"$\Vert {\widehat{\sigma}}_{IDE} - {\widehat{\sigma}}_{ODE}\Vert_2$")
    else:
        secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                      5: 'ICU', 6: 'Recovered', 7: 'Dead'}
        fig.supylabel(
            r"$\Vert \widehat{Z}_{IDE} - \widehat{Z}_{ODE}\Vert_2$")

    if flows:
        # fig, axs = plt.subplots(5, 2, sharex=True, figsize=(10, 10))
        num_plots = 10
    else:
        # fig, axs = plt.subplots(4, 2, sharex=True, figsize=(10, 8))
        num_plots = 8

    colors = [plt.cm.viridis(x) for x in np.linspace(0, 1, num_plots)]

    angles = [0, 45, 0, 45, 0, 45, 0, 45, 0, 45]
    if not flows:
        angles[6] = 45
        angles[7] = 0

    for i in range(num_plots):

        # Plot results.
        rotation = Affine2D().rotate_deg(angles[i])
        ax.plot(timesteps_ide,
                errors[:, i], '-', marker=MarkerStyle('x', 'full', rotation), markersize=5, color=colors[i], label=secir_dict[i])

    # Plot comparison line for linear convergence.
    factor = 50 * min(errors[0, :])
    comparison = [factor * dt for dt in timesteps_ide]
    ax.plot(timesteps_ide, comparison, '--', color='gray',
            label=r"$\mathcal{O}(\Delta t)$")

    # Adapt plots.
    ax.set_xscale("log", base=10)
    ax.set_yscale("log", base=10)
    ax.invert_xaxis()

    fig.supxlabel('Time step')

    fig.legend()

    # plt.tight_layout(rect=[0, 0, 0.83, 1])

    if save:
        if flows:
            if not os.path.isdir('plots/flows'):
                os.makedirs('plots/flows')
            plt.savefig(f'plots/flows/convergence_flows_setting={setting}.png', format='png',
                        dpi=500)
        else:
            if not os.path.isdir('plots/compartments'):
                os.makedirs('plots/compartments')
            plt.savefig(f'plots/compartments/convergence_compartments_setting={setting}.png', format='png',
                        dpi=500)


def compute_order_of_convergence(errors, timesteps_ide, flows=False):
    """ Compute order of convergence between two consecutive time step sizes."""
    if flows:
        num_orders = 10
    else:
        num_orders = 8

    timesteps_ide
    order = []
    for compartment in range(num_orders):
        order.append([])
        for i in range(len(errors)-1):
            order[compartment].append(np.log(errors[i+1][compartment]/errors[i][compartment]) /
                                      np.log(timesteps_ide[i+1]/timesteps_ide[i]))
    return np.array(order)


def print_results(groundtruth, results, timesteps_ide):
    """ Print results of ODE and IDE simulations at tmax."""

    compartments = ['S', 'E', 'C', 'I', 'H', 'U', 'R', 'D']

    for compartment in range(8):
        print('\n')
        print(f'{compartments[compartment]}: ')
        print('Groundtruth (using ODE): ',
              groundtruth['ode'][0][-1][compartment])
        print('\n')
        print('IDE:')
        for i in range(len(timesteps_ide)):
            print('Timestep ', timesteps_ide[i], ':',
                  results['ide'][i][-1][compartment])


def print_initial_values(groundtruth, results, dt_ode, timesteps_ide):
    """ Print results of ODE and IDE simulations at t0_ide."""

    compartments = ['S', 'E', 'C', 'I', 'H', 'U', 'R', 'D']

    for compartment in range(8):
        print('\n')
        print(f'{compartments[compartment]}: ')
        # get values of ODE model  at t0_ide = 35
        print('Groundtruth (using ODE): ',
              groundtruth['ode'][0][int(35/dt_ode)][compartment])
        print('\n')
        print('IDE:')
        for i in range(len(timesteps_ide)):
            print('Timestep ', timesteps_ide[i], ':',
                  results['ide'][i][0][compartment])

    print(0)


def print_total_population(results, timesteps_ide):

    for i in range(len(timesteps_ide)):
        print('Timestep ', timesteps_ide[i], ':')
        total_population = 0
        for compartment in range(8):
            total_population += results['ide'][i][-1][compartment]

        print('\n')
        print('Total population (IDE):', total_population)


def print_errors(errors, timesteps_ide):

    compartments = ['S', 'E', 'C', 'I', 'H', 'U', 'R', 'D']

    for compartment in range(8):
        print('\n')
        print(f'{compartments[compartment]}: ')
        print('\n')
        print('Errors of IDE (compared to ODE):')
        for i in range(len(timesteps_ide)):
            print('Timestep ', timesteps_ide[i], ':',
                  errors[i][compartment])


def main():
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "results")

    setting = 2

    print('Setting: ', setting)

    dt_ode = 1e-4

    timesteps_ide = [1e-2, 1e-3, 1e-4]

    flows = True

    groundtruth = read_groundtruth(data_dir, dt_ode, setting, flows=flows)

    results = read_data(data_dir, dt_ode, timesteps_ide, setting, flows=flows)

    errors_l2 = compute_error_norm_l2(
        groundtruth, results, dt_ode, timesteps_ide, flows=flows)

    errors = compute_error_norm_tmax(
        groundtruth, results, timesteps_ide, flows=flows)

    # plot_convergence(errors, timesteps_ide, setting, flows=flows, save=True)

    plot_convergence_oneplot(
        errors_l2, timesteps_ide, setting, flows=flows, save=True)

    order = compute_order_of_convergence(errors_l2, timesteps_ide, flows=flows)

    print('Orders of convergence: ', order)

    # print_initial_values(groundtruth, results, dt_ode, timesteps_ide)

    # print_total_population(results, timesteps_ide)

    # errors = compute_error_norm_tmax(
    #     groundtruth, results, timesteps_ide, flows=flows)

    # print_results(groundtruth, results, timesteps_ide)

    # print_errors(errors, timesteps_ide)

    return


if __name__ == '__main__':
    main()
