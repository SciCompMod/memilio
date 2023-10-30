import h5py
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import memilio.epidata.getDataIntoPandasDataFrame as gd

# Read data into a dict, where the keys correspond to ODE and IDE models. There
# we have an array that contains all results for SECIHURD for all time points
# for each time step size that is investigated.


def read_data(data_dir, timesteps):

    models = ['ode', 'ide']
    results = {models[0]: [], models[1]: []}
    for model in models:
        for timestep in timesteps:
            h5file = h5py.File(os.path.join(data_dir, 'result_{}_dt={}_setting2'.format(
                model, timestep)) + '.h5', 'r')

            if (len(list(h5file.keys())) > 1):
                raise gd.DataError("File should contain one dataset.")
            if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
                raise gd.DataError("Expected only one group.")

            data = h5file[list(h5file.keys())[0]]

            if len(data['Total'][0]) == 8:
                # As there should be only one Group, total is the simulation result
                results[model].append(data['Total'][:, :])
            elif len(data['Total'][0]) == 10:
                # in ODE there are two compartments we don't use, throw these out
                results[model].append(
                    data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]])

            dates = data['Time'][:]

            if (results[model][-1].shape[1] != 8):
                raise gd.DataError(
                    "Expected a different number of compartments.")

            h5file.close()

    return results

# Compute norm of a one-dimensional timeseries


def compute_norm(timeseries, timestep):

    norm = np.sqrt(np.sum(timeseries**2 * timestep))

    return norm

# Compute norm of the difference between time series for S from ODE and time series for S from IDE
# TODO: Überlegen, ob wir hier den ersten Zeitschritt rausnehmen, damit
# besser vergleichbar zwischen verschiedenen Zeitschrittgrößen


def compute_error_norm_S_timeseries(results, timesteps):
    errors = []

    # Compute error for S for every time step
    for i in range(len(results['ode'])):
        timestep = timesteps[i]
        num_timepoints = len(results['ide'][i])
        # for now, compute only difference for S
        difference = results['ode'][i][num_timepoints -
                                       1:][:, 0]-results['ide'][i][:, 0]
        errors.append(compute_norm(difference, timestep))

    return errors

# Compute norm of the difference between S from ODE and IDE at tmax, respectively.


def compute_error_norm_S_endpoint(results, timesteps):
    errors = []

    # Compute error for S for every time step
    for i in range(len(results['ode'])):
        timestep = timesteps[i]
        num_timepoints = len(results['ide'][i])
        # for now, compute only difference for S
        difference = results['ode'][i][-1][0]-results['ide'][i][-1][0]
        errors.append(np.sqrt(difference ** 2))

    return errors

# Plot errors against timesteps.


def plot_convergence(errors, timesteps, save=False):
    fig, ax = plt.subplots()

    ax.plot(timesteps, errors, '-o', label='Results')
    comparison = [2400 * dt for dt in timesteps]
    ax.plot(timesteps, comparison, color='lightgray',
            label=r"$\mathcal{O}(\Delta t)$")
    # ax.scatter(timesteps, errors)
    # ax.plot(np.linspace(timesteps[0], timesteps[-1],1000), )
    ax.set_xscale("log", base=10)
    ax.set_yscale("log", base=10)
    ax.invert_xaxis()

    fig.supxlabel('Time step')
    fig.supylabel(r"$\Vert S_{IDE}(t_{max}) - S_{ODE}(t_{max})\Vert$")

    plt.legend()

    if save:
        if not os.path.isdir('Plots'):
            os.makedirs('Plots')
        plt.savefig('Plots/convergence.png',
                    bbox_inches='tight', dpi=500)

    plt.show()

# Compute order of convergence between two consecutive time step sizes


def compute_order_of_convergence(errors, timesteps):
    order = []
    for i in range(len(errors)-1):
        order.append(np.log(errors[i+1]/errors[i]) /
                     np.log(timesteps[i+1]/timesteps[i]))

    return order

# Print relatie and absolute results of ODE and IDE simulations


def print_results(results, timesteps):

    for i in range(len(timesteps)-1):
        print('ODE: ', results['ode'][i][-1][0]/results['ode'][i+1][-1][0])

    for i in range(len(timesteps)-1):
        print('IDE: ', results['ide'][i][-1][0]/results['ide'][i+1][-1][0])

    for i in range(len(timesteps)):
        print('ODE: ', results['ode'][i][-1][0])
        print('IDE: ', results['ide'][i][-1][0])


if __name__ == '__main__':
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "results")

    timesteps = ['1e-2', '1e-3', '1e-4']

    results = read_data(data_dir, timesteps)

    timesteps = [1e-2, 1e-3, 1e-4]
    errors = compute_error_norm_S_endpoint(results, timesteps)

    plot_convergence(errors, timesteps, save=True)

    order = compute_order_of_convergence(errors, timesteps)

    print('Orders of convergence: ', order)

    print_results(results, timesteps)
