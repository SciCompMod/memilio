import h5py
import os
import numpy as np
import matplotlib.pyplot as plt

# matplotlib.use('tkagg')

# import memilio.epidata.getDataIntoPandasDataFrame as gd


# We define the groundtruth as the results obtained by the ODE model with timestep dt=1e-6.
def read_groundtruth(data_dir, dt_ode, setting, flows=False):

    model = 'ode'
    results = {model: []}

    if flows:
        h5file = h5py.File(os.path.join(data_dir, 'result_{}_flows_dt={}_setting{}'.format(
            model, dt_ode, setting)) + '.h5', 'r')
    else:
        h5file = h5py.File(os.path.join(data_dir, 'result_{}_dt={}_setting{}'.format(
            model, dt_ode, setting)) + '.h5', 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]

    if flows:
        # For flows we divide the values by dt_ode to go from \tilde{\sigma} to \hat{\sigma} and make results comparable
        results[model].append(data['Total'][:, :]/float(dt_ode))
    else:
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


# Read data into a dict, where the keys correspond to the respective model.
# At the moment we are only storing results of the IDE model here. There
# we have an array that contains all results for SECIHURD for all time points
# for each time step size that is investigated.
def read_data(data_dir, dt_ode, timesteps_ide, setting, flows=False):

    models = ['ide']
    results = {models[0]: []}
    for model in models:
        for timestep in timesteps_ide:
            if flows:
                h5file = h5py.File(os.path.join(data_dir, 'result_{}_flows_dt={}_init_dt_ode={}_setting{}'.format(
                    model, timestep, dt_ode, setting)) + '.h5', 'r')
            else:
                h5file = h5py.File(os.path.join(data_dir, 'result_{}_dt={}_init_dt_ode={}_setting{}'.format(
                    model, timestep, dt_ode, setting)) + '.h5', 'r')

            # if (len(list(h5file.keys())) > 1):
            #     raise gd.DataError("File should contain one dataset.")
            # if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            #     raise gd.DataError("Expected only one group.")

            data = h5file[list(h5file.keys())[0]]

            if flows:
                # For flows we divide the values by dt_ode to go from \tilde{\sigma} to \hat{\sigma} and make results comparable
                results[model].append(data['Total'][:, :]/float(timestep))
            else:
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


def compute_error_norm_S_timeseries(groundtruth, results, dt_ode, timesteps_ide):
    errors = []

    # Compute error for S for every time step
    for i in range(len(results['ide'])):
        timestep = timesteps_ide[i]
        scale_timesteps = timestep/float(dt_ode)
        num_timepoints = len(results['ide'][i])
        # for now, compute only difference for S
        difference = groundtruth['ode'][0][int(scale_timesteps*(num_timepoints -
                                                                1))::int(scale_timesteps)][:, 0]-results['ide'][i][:, 0]
        errors.append(compute_norm(difference, timestep))

    return errors

# Compute norm of the difference between ODE and IDE at tmax for all compartments


def compute_error_norm_tmax(groundtruth, results, timesteps_ide, flows=False):

    if flows:
        num_errors = 10
    else:
        num_errors = 8
    errors = []

    # Compute error of compartments in IDE simulation at tmax
    for i in range(len(results['ide'])):
        errors.append([])
        # compute difference for all compartments
        for compartment in range(num_errors):
            timestep = timesteps_ide[i]
            num_timepoints = len(results['ide'][i])
            difference = groundtruth['ode'][0][-1][compartment] - \
                results['ide'][i][-1][compartment]
            errors[i].append(np.sqrt(difference ** 2))

    return np.array(errors)

# Plot errors against timesteps.


def plot_convergence(errors, timesteps_ide, setting, flows=False, compartment=None, save=False):

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
            fig, ax = plt.subplots(5, 2, sharex=True)
            num_plots = 10
        else:
            fig, ax = plt.subplots(4, 2, sharex=True)
            num_plots = 8

        for i in range(num_plots):

            # plot comparison line for linear convergence
            comparison = [500 * dt for dt in timesteps_ide]
            ax[int(i/2), i % 2].plot(timesteps_ide, comparison, color='lightgray',
                                     label=r"$\mathcal{O}(\Delta t)$")

            # plot results
            ax[int(i/2), i % 2].plot(timesteps_ide,
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

        # plt.show()


# Compute order of convergence between two consecutive time step sizes


def compute_order_of_convergence(errors, timesteps_ide, flows=False):
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

# Print results of ODE and IDE simulations at tmax


def print_results(groundtruth, results, timesteps_ide):

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

# Print results of ODE and IDE simulations at t0_ide


def print_initial_values(groundtruth, results, dt_ode, timesteps_ide):

    compartments = ['S', 'E', 'C', 'I', 'H', 'U', 'R', 'D']

    for compartment in range(8):
        print('\n')
        print(f'{compartments[compartment]}: ')
        # get values of ODE model  at t0_ide = 35
        print('Groundtruth (using ODE): ',
              groundtruth['ode'][0][int(35/float(dt_ode))][compartment])
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

    setting = '2'

    print('Setting: ', setting)

    dt_ode = '1e-4'

    timesteps_ide = ['1e-2', '1e-3']  # , '1e-4'

    flows = True

    groundtruth = read_groundtruth(data_dir, dt_ode, setting, flows=flows)

    results = read_data(data_dir, dt_ode, timesteps_ide, setting, flows=flows)

    timesteps_ide = [1e-2, 1e-3]  # , 1e-4

    errors = compute_error_norm_tmax(
        groundtruth, results, timesteps_ide, flows=flows)

    plot_convergence(errors, timesteps_ide, setting, flows=False, save=True)

    # print_initial_values(groundtruth, results, dt_ode, timesteps_ide)

    # print_total_population(results, timesteps_ide)

    order = compute_order_of_convergence(errors, timesteps_ide, flows=flows)

    print('Orders of convergence: ', order)

    # print_results(groundtruth, results, timesteps_ide)

    # print_errors(errors, timesteps_ide)

    return


if __name__ == '__main__':
    main()
