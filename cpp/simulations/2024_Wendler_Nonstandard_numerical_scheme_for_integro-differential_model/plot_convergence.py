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


def read_groundtruth(data_dir, ode_exponent, save_exponent, flows=False):
    """ Read groundtruth from data. We define the groundtruth as the results obtained by the ODE model with timestep dt=1e-6.

    @param[in] data_dir Directory where h5 files are stored. 
    @param[in] ode_exponent Exponent that determines time step size via dt =10^{-ode_exponent}.
    @param[in] save_exponent The results of the ODE model were saved using the step size 10^{-save_exponent}.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False. 
    @returns Dict with results of ODE model.    
    """
    model = 'ode'
    results = {model: []}
    if flows:
        h5file = h5py.File(os.path.join(
            data_dir, f'result_{model}_flows_dt=1e-{ode_exponent:.0f}_savefrequency{save_exponent:.0f}.h5'), 'r')
    else:
        h5file = h5py.File(os.path.join(
            data_dir, f'result_{model}_dt=1e-{ode_exponent:.0f}_savefrequency{save_exponent:.0f}.h5'), 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]

    if flows:
        # Flows are already scaled to one day.
        results[model].append(data['Total'][:, :])
    else:
        if len(data['Total'][0]) == 8:
            # As there should be only one Group, total is the simulation result
            results[model].append(data['Total'][:, :])
        elif len(data['Total'][0]) == 10:
            # In ODE there are two compartments we don't use, throw these out
            results[model].append(
                data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]])
        else:
            raise gd.DataError(
                'Expected a different size of vector in time series.')

    dates = data['Time'][:]

    h5file.close()

    return results


def read_data(data_dir, ode_exponent, exponents_ide, flows=False):
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
            if flows:
                h5file = h5py.File(os.path.join(
                    data_dir, f'result_{model}_flows_dt=1e-{exponent:.0f}_init_dt_ode=1e-{ode_exponent:.0f}.h5'), 'r')
            else:
                h5file = h5py.File(os.path.join(
                    data_dir, f'result_{model}_dt=1e-{exponent:.0f}_init_dt_ode=1e-{ode_exponent:.0f}.h5'), 'r')

            data = h5file[list(h5file.keys())[0]]

            if flows:
                # Flows are already scaled to one day.
                results[model].append(data['Total'][:, :])
            else:
                if len(data['Total'][0]) == 8:
                    # As there should be only one Group, total is the simulation result.
                    results[model].append(data['Total'][:, :])
                elif len(data['Total'][0]) == 10:
                    # In ODE there are two compartments we don't use, throw these out.
                    results[model].append(
                        data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]])
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


def compute_relerror_norm_l2(groundtruth, results, save_exponent, timesteps_ide, flows=False):
    """ Computes relative L2 norm of the difference between time series from ODE and time series
    from IDE for all compartments/flows.

    @param[in] groundtruth Result obtained with ODE model.
    @param[in] results Results obtained with IDE model for different time step sizes. 
    @param[in] save_exponent The results of the ODE model were saved using the step size 10^{-save_exponent}.
    @param[in] timesteps_ide List of time steps used in IDE simulations.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False. 
    @param[in] Array that contains computed errors.
    """
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
            scale_timesteps = timestep/pow(10, -save_exponent)
            num_timepoints = len(results['ide'][i])
            # Only consider compartments of ODE model for t>=t0_IDE. Use that 2*t0_IDE = t_max.
            if flows:
                difference = groundtruth['ode'][0][int(scale_timesteps*(num_timepoints/2-1))::int(
                    scale_timesteps)][:, compartment]-results['ide'][i][int(num_timepoints/2-1):][:, compartment]
                norm_groundtruth = compute_l2_norm(groundtruth['ode'][0][int(scale_timesteps*(num_timepoints/2-1))::int(
                    scale_timesteps)][:, compartment], timestep)
                errors[i].append(compute_l2_norm(
                    difference, timestep)/norm_groundtruth)
            else:
                difference = groundtruth['ode'][0][int(scale_timesteps*(num_timepoints -
                                                                        1))::int(scale_timesteps)][:, compartment]-results['ide'][i][:, compartment]
                norm_groundtruth = compute_l2_norm(groundtruth['ode'][0][int(scale_timesteps*(num_timepoints -
                                                                                              1))::int(scale_timesteps)][:, compartment], timestep)
                errors[i].append(compute_l2_norm(
                    difference, timestep)/norm_groundtruth)

    return np.array(errors)


def plot_convergence(errors, timesteps_ide, flows=False, save_dir=""):
    """ Plots errors against timesteps with a subplot for each compartment /flow.

    @param[in] errors Array that contains computed errors of IDE model compared to groundtruth.
    @param[in] timesteps_ide List of time steps used in IDE simulations.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False. 
    @param[in] save_dir Directory where plot will be stored. Default is an empty string leading to the plot not being 
        saved.
    """
    # Define subplots and labels.
    if flows:
        fig, axs = plt.subplots(
            5, 2, sharex=True,  sharey=True, figsize=(10, 10))
        num_plots = 10
        secir_dict = {0: r"$\sigma_S^E$", 1: r"$\sigma_E^C$", 2: r"$\sigma_C^I$", 3: r"$\sigma_C^R$", 4: r"$\sigma_I^H$",
                      5: r"$\sigma_I^R$", 6: r"$\sigma_H^U$", 7: r"$\sigma_H^R$", 8: r"$\sigma_U^D$", 9: r"$\sigma_U^R$"}
        labels = [
            r"$\Vert {\widehat{\sigma}}_{\text{IDE}} - {\widehat{\sigma}}_{\text{ODE}}\Vert_{2,\text{rel}}$", r"$\mathcal{O}(\Delta t)$"]
    else:
        fig, axs = plt.subplots(4, 2, sharex=True, figsize=(10, 8))
        num_plots = 8
        secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                      5: 'ICU', 6: 'Recovered', 7: 'Dead'}
        labels = [
            r"$\Vert \widehat{Z}_{\text{IDE}} - \widehat{Z}_{\text{ODE}}\Vert_{2,\text{rel}}$", r"$\mathcal{O}(\Delta t)$"]

    # Define colors, we use helmholtzdarkblue and helmholtzclaim.
    colors = [(0, 40/255, 100/255), (20/255, 200/255, 255/255)]

    for i in range(num_plots):
        # Plot results.
        axs[int(i/2), i % 2].plot(timesteps_ide,
                                  errors[:, i], '-o', color=colors[1])

        # Plot comparison line for linear convergence.
        comparison = [dt for dt in timesteps_ide]
        axs[int(i/2), i % 2].plot(timesteps_ide, comparison,
                                  '--', color='gray', linewidth=1.2)

        # Adapt plots.
        axs[int(i/2), i % 2].set_xscale("log", base=10)
        axs[int(i/2), i % 2].set_yscale("log", base=10)

        axs[int(i/2), i % 2].set_title(secir_dict[i], fontsize=10)
        axs[int(i/2), i % 2].grid(True, linestyle='--', alpha=0.6)

    fig.supxlabel(r'Time step $\Delta t$', fontsize=12)

    # Invert x axis only for one plot so that sharex=True and invert_xaxis work as intended.
    axs[0, 0].invert_xaxis()

    lgd = fig.legend(labels, ncol=2,  loc='outside lower center',
                     fontsize=14, bbox_to_anchor=(0.5, -0.06), bbox_transform=fig.transFigure)
    plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.1)
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        if flows:
            plt.savefig(f'{save_dir}/convergence_all_flows.png', format='png', bbox_extra_artists=(lgd,), bbox_inches='tight',
                        dpi=500)
        else:
            plt.savefig(f'{save_dir}/convergence_all_compartments.png', format='png', bbox_extra_artists=(lgd,), bbox_inches='tight',
                        dpi=500)


def plot_convergence_oneplot(errors, timesteps_ide, flows=False, save_dir=""):
    """ Plot errors against timesteps. This function creates one plot to depict all compartments/flows, respectively.

    @param[in] errors Array that contains computed errors of IDE model compared to groundtruth.
    @param[in] timesteps_ide List of time steps used in IDE simulations.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False. 
    @param[in] save_dir Directory where plot will be stored. Default is an empty string leading to the plot not being 
        saved.
    """
    plt.figure()

    if flows:
        secir_dict = {0: r"$\sigma_S^E$", 1: r"$\sigma_E^C$", 2: r"$\sigma_C^I$", 3: r"$\sigma_C^R$", 4: r"$\sigma_I^H$",
                      5: r"$\sigma_I^R$", 6: r"$\sigma_H^U$", 7: r"$\sigma_H^R$", 8: r"$\sigma_U^D$", 9: r"$\sigma_U^R$"}
        plt.ylabel(
            r"$err_{\text{rel}}$", fontsize=10)
    else:
        secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                      5: 'ICU', 6: 'Recovered', 7: 'Dead'}
        plt.ylabel(
            r"$err_{\text{rel}}$", fontsize=10)

    if flows:
        num_lines = 10
    else:
        num_lines = 8

    colors = [plt.cm.viridis(x) for x in np.linspace(0, 1, num_lines)]

    # Angles to rotate markers of the plot.
    angles = [0, 45, 0, 45, 0, 45, 0, 45, 0, 45]

    for i in range(num_lines):
        # Plot results.
        rotation = Affine2D().rotate_deg(angles[i])
        plt.plot(timesteps_ide,
                 errors[:, i], '-', marker=MarkerStyle('x', 'full', rotation), markersize=5, color=colors[i], label=secir_dict[i])

    # Plot comparison line for linear convergence.
    comparison = [dt for dt in timesteps_ide]
    plt.plot(timesteps_ide, comparison, '--', color='gray',
             label=r"$\mathcal{O}(\Delta t)$")

    # Adapt plots.
    plt.xscale("log", base=10)
    plt.yscale("log", base=10)
    plt.gca().invert_xaxis()

    plt.xlabel(r'Time step $\Delta t$', fontsize=10)

    plt.legend(fontsize=10, framealpha=0.5, ncol=2)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        if flows:
            plt.savefig(f'{save_dir}/convergence_flows.png', format='png',
                        dpi=500)
        else:
            plt.savefig(f'{save_dir}/convergence_compartments.png', format='png',
                        dpi=500)


def compute_order_of_convergence(errors, timesteps_ide, flows=False):
    """ Compute order of convergence between two consecutive time step sizes.

    @param[in] errors Array that contains computed errors of IDE model compared to groundtruth.
    @param[in] timesteps_ide List of time steps used in IDE simulations.
    @param[in] flows Bool that determines whether we consider flows or compartments. Default is False. 
    """
    if flows:
        num_orders = 10
    else:
        num_orders = 8

    order = []
    for compartment in range(num_orders):
        order.append([])
        for i in range(len(errors)-1):
            order[compartment].append(np.log(errors[i+1][compartment]/errors[i][compartment]) /
                                      np.log(timesteps_ide[i+1]/timesteps_ide[i]))
    return np.array(order)


def main():
    # Paths are valid if script is executed e.g. in
    # memilio/cpp/simulations/2024_Wendler_Nonstandard_numerical_scheme_for_integro-differential_model.

    # Path where simulation results (generated with ide_convergence_rate.cpp) are stored.
    result_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/simulation_results/convergence/")

    # Path where plots will be stored.
    plot_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/plots/convergence/")

    # The ODE model was simulated using a fixed step size dt=10^{-ode_exponent}.
    ode_exponent = 6
    # The results of the ODE model were saved using the step size 10^{-save_exponent}
    # as for very small step sizes used for the simulation, the number of time points stored gets very big.
    save_exponent = 4
    # The IDE model was simulated using a fixed step size dt=10^{-ide_exponent} for ide_exponent in ide_exponents.
    ide_exponents = [1, 2, 3, 4]
    # Calculate time steps resulting from exponents_ide.
    timesteps_ide = []
    for exp in ide_exponents:
        timesteps_ide.append(pow(10, -exp))

    # Plot compartments and flows.
    flow_bools = [False, True]

    for flow_bool in flow_bools:
        # Read groundtruth (from ODE model).
        groundtruth = read_groundtruth(
            result_dir, ode_exponent, save_exponent, flow_bool)

        # Read results from IDE simulations.
        results = read_data(result_dir, ode_exponent, ide_exponents, flow_bool)

        # Compute relative L2 error norm of IDE results compared to groundtruth.
        relerrors_l2 = compute_relerror_norm_l2(
            groundtruth, results, save_exponent, timesteps_ide, flow_bool)

        # Plot convergence of all compartments/flows in one plot, respectively.
        plot_convergence_oneplot(
            relerrors_l2, timesteps_ide, flow_bool, plot_dir)
        # Plot convergence of all compartments/flows separately.
        plot_convergence(relerrors_l2, timesteps_ide,  flow_bool, plot_dir)

        # Determine order of convergence
        order = compute_order_of_convergence(
            relerrors_l2, timesteps_ide, flow_bool)

        print('Orders of convergence: ', order)


if __name__ == '__main__':
    main()
