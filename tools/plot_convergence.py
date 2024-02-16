import h5py
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

# matplotlib.use('tkagg')

# import memilio.epidata.getDataIntoPandasDataFrame as gd


# We define the groundtruth as the results obtained by the ODE model with timestep dt=1e-6.
def read_groundtruth(data_dir, setting):

    model = 'ode'
    results = {model: []}

    h5file = h5py.File(os.path.join(data_dir, 'result_{}_dt=1e-3_setting{}'.format(
        model, setting)) + '.h5', 'r')

    # if (len(list(h5file.keys())) > 1):
    #     raise gd.DataError("File should contain one dataset.")
    # if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
    #     raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]

    if len(data['Total'][0]) == 8:
        # As there should be only one Group, total is the simulation result
        results[model].append(data['Total'][:, :])
    elif len(data['Total'][0]) == 10:
        # in ODE there are two compartments we don't use, throw these out
        results[model].append(
            data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]])

    dates = data['Time'][:]

    # if (results[model][-1].shape[1] != 8):
    #     raise gd.DataError(
    #         "Expected a different number of compartments.")

    h5file.close()

    return results


# Read data into a dict, where the keys correspond to ODE and IDE models. There
# we have an array that contains all results for SECIHURD for all time points
# for each time step size that is investigated.
def read_data(data_dir, timesteps, setting):

    models = ['ode', 'ide']
    results = {models[0]: [], models[1]: []}
    for model in models:
        for timestep in timesteps:
            h5file = h5py.File(os.path.join(data_dir, 'result_{}_dt={}_setting{}'.format(
                model, timestep, setting)) + '.h5', 'r')

            # if (len(list(h5file.keys())) > 1):
            #     raise gd.DataError("File should contain one dataset.")
            # if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            #     raise gd.DataError("Expected only one group.")

            data = h5file[list(h5file.keys())[0]]

            if len(data['Total'][0]) == 8:
                # As there should be only one Group, total is the simulation result
                results[model].append(data['Total'][:, :])
            elif len(data['Total'][0]) == 10:
                # in ODE there are two compartments we don't use, throw these out
                results[model].append(
                    data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]])

            dates = data['Time'][:]

            # if (results[model][-1].shape[1] != 8):
            #     raise gd.DataError(
            #         "Expected a different number of compartments.")

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

# Compute norm of the difference between ODE and IDE at tmax for all compartments


def compute_error_norm_tmax(groundtruth, results, timesteps):
    errors = []

    # Compute error for S for every used time step in IDE simulation at tmax
    for i in range(len(results['ide'])):
        errors.append([])
        # compute difference for all compartments
        for compartment in range(8):
            timestep = timesteps[i]
            num_timepoints = len(results['ide'][i])
            # for now, compute only difference for S
            difference = groundtruth['ode'][0][-1][compartment] - \
                results['ide'][i][-1][compartment]
            errors[i].append(np.sqrt(difference ** 2))

    return np.array(errors)

# Plot errors against timesteps.


def plot_convergence(errors, timesteps, setting, compartment=None, save=False):

    secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                  5: 'ICU', 6: 'Recovered', 7: 'Dead'}
    compartments = ['S', 'E', 'C', 'I', 'H', 'U', 'R', 'D']

    # helmholtzdarkblue, helmholtzclaim
    colors = [(0, 40/255, 100/255), (20/255, 200/255, 255/255)]
    if compartment != None:

        # TODO: include check if compartment is in 0,...,8

        fig, ax = plt.subplots()

        ax.plot(timesteps, errors[:, compartment],
                '-o', color=colors[0], label='Results')
        comparison = [1800 * dt for dt in timesteps]
        ax.plot(timesteps, comparison, color='lightgray',
                label=r"$\mathcal{O}(\Delta t)$")

        ax.set_xscale("log", base=10)
        ax.set_yscale("log", base=10)
        ax.invert_xaxis()

        fig.supxlabel('Time step')
        fig.supylabel(r"$\Vert {compartment}_{IDE}(t_{max}) - {compartment}_{ODE}(t_{max})\Vert$".replace(
            'compartment', compartments[compartment]))

        plt.legend()

        if save:
            if not os.path.isdir('plots'):
                os.makedirs('plots')
            plt.savefig('plots/convergence_{compartments[compartment]}_setting{}.png'.format(setting),
                        bbox_inches='tight', dpi=500)

        else:
            plt.show()

    else:
        fig, ax = plt.subplots(4, 2, sharex=True)

        for i in range(8):

            # plot comparison line for linear convergence
            comparison = [500 * dt for dt in timesteps]
            ax[int(i/2), i % 2].plot(timesteps, comparison, color='lightgray',
                                     label=r"$\mathcal{O}(\Delta t)$")

            # plot results
            ax[int(i/2), i % 2].plot(timesteps,
                                     errors[:, i], '-o', color=colors[0], label='Results')

            # adapt plots
            ax[int(i/2), i % 2].set_xscale("log", base=10)
            ax[int(i/2), i % 2].set_yscale("log", base=10)

            ax[int(i/2), i % 2].set_title(secir_dict[i], fontsize=8)

            fig.supxlabel('Time step')
            fig.supylabel(
                r"$\Vert {K}_{IDE}(t_{max}) - {K}_{ODE}(t_{max})\Vert$")

        # invert x axis only for one plot so that sharex=True and invert_xaxis work as intended
        ax[0, 0].invert_xaxis()

        labels = [r"$\mathcal{O}(\Delta t)$", 'Results']
        fig.legend(labels, loc='center right',
                   fancybox=False, shadow=False, ncol=1)

        plt.tight_layout(rect=[0, 0, 0.83, 1])

        if save:
            if not os.path.isdir('plots'):
                os.makedirs('plots')
            plt.savefig(f'plots/convergence_all_setting={setting}.png', format='png',
                        dpi=500)  # bbox_inches='tight',

        # plt.show()


# Compute order of convergence between two consecutive time step sizes


def compute_order_of_convergence(errors, timesteps):
    order = []
    for compartment in range(8):
        order.append([])
        for i in range(len(errors)-1):
            order[compartment].append(np.log(errors[i+1][compartment]/errors[i][compartment]) /
                                      np.log(timesteps[i+1]/timesteps[i]))
    return np.array(order)

# Print relative and absolute results of ODE and IDE simulations


def print_results(groundtruth, results, timesteps):

    compartments = ['S', 'E', 'C', 'I', 'H', 'U', 'R', 'D']

    for compartment in range(8):
        print('\n')
        print(f'{compartments[compartment]}: ')
        print('Groundtruth (using ODE): ',
              groundtruth['ode'][0][-1][compartment])
        print('\n')
        print('IDE:')
        for i in range(len(timesteps)):
            print('Timestep ', timesteps[i], ':',
                  results['ide'][i][-1][compartment])


def print_errors(errors, timesteps):

    compartments = ['S', 'E', 'C', 'I', 'H', 'U', 'R', 'D']

    for compartment in range(8):
        print('\n')
        print(f'{compartments[compartment]}: ')
        print('\n')
        print('Errors of IDE (compared to ODE):')
        for i in range(len(timesteps)):
            print('Timestep ', timesteps[i], ':',
                  errors[i][compartment])


def main():
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "results")

    setting = 2

    timesteps = ['1e-2', '1e-3', '1e-4']  # , '1e-4'

    groundtruth = read_groundtruth(data_dir, setting)

    results = read_data(data_dir, timesteps, setting)

    timesteps = [1e-2, 1e-3, 1e-4]  # , 1e-4
    errors = compute_error_norm_tmax(groundtruth, results, timesteps)

    plot_convergence(errors, timesteps, setting, save=True)

    # order = compute_order_of_convergence(errors, timesteps)

    # print('Orders of convergence: ', order)

    # print_results(groundtruth, results, timesteps)

    # print_errors(errors, timesteps)

    return


if __name__ == '__main__':
    main()
