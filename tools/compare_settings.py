import h5py
import os
import numpy as np
import matplotlib.pyplot as plt

# compare the results of two IDE simulations in different settings


def read_data_for_comparison(data_dir, timesteps, settings, ide=True):

    results = {settings[0]: [], settings[1]: []}
    for setting in settings:
        for timestep in timesteps:
            if ide:
                h5file = h5py.File(os.path.join(data_dir, 'result_ide_dt={}_init_dt_ode=1e-4_setting{}'.format(
                    timestep, setting)) + '.h5', 'r')
            else:
                h5file = h5py.File(os.path.join(data_dir, 'result_ide_dt={}_init_dt_ode=1e-4_setting{}_new'.format(
                    timestep, setting)) + '.h5', 'r')

            data = h5file[list(h5file.keys())[0]]

            if len(data['Total'][0]) == 8:
                # As there should be only one Group, total is the simulation result
                results[setting].append(data['Total'][:, :])
            elif len(data['Total'][0]) == 10:
                # in ODE there are two compartments we don't use, throw these out
                results[setting].append(
                    data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]])

            dates = data['Time'][:]

            # if (results[model][-1].shape[1] != 8):
            #     raise gd.DataError(
            #         "Expected a different number of compartments.")

            h5file.close()

    return results

# assume that the first setting is the ussual one and the second is R and D switched


def compare_data(results, settings, switch=False):
    compartments = ['S', 'E', 'C', 'I', 'H', 'U', 'R', 'D']

    if switch:
        # compare SECIHU
        for i in range(5):
            # compare initial values
            if (results[settings[0]][0][0][i] == results[settings[1]][0][0][i]):
                print('Results for ', compartments[i], 'are equal at t0.')
            else:
                print('Results for ', compartments[i], 'are not equal at t0.')
            # compare values at tmax
            if (results[settings[0]][0][-1][i] == results[settings[1]][0][-1][i]):
                print('Results for ', compartments[i], 'are equal at tmax.')
            else:
                print('Results for ',
                      compartments[i], 'are not equal at tmax.')

        # compare R and D
        # compare initial values
        if (results[settings[0]][0][0][6] == results[settings[1]][0][0][7]) and (results[settings[0]][0][0][7] == results[settings[1]][0][0][6]):
            print('Results for R and D are equal at t0.')
        else:
            print('Results for R and D are not equal at t0.')
        # compare values at tmax
        if (results[settings[0]][0][-1][6] == results[settings[1]][0][-1][7]) and (results[settings[0]][0][-1][7] == results[settings[1]][0][-1][6]):
            print('Results for R and D are equal at tmax.')
        else:
            print('Results for R and D are not equal at tmax.')

    else:
        # compare SECIHURD
        for i in range(len(compartments)):
            # compare initial values
            if (results[settings[0]][0][0][i] == results[settings[1]][0][0][i]):
                print('Results for ', compartments[i], 'are equal at t0.')
            else:
                print('Results for ', compartments[i], 'are not equal at t0.')
            # compare values at tmax
            if (results[settings[0]][0][-1][i] == results[settings[1]][0][-1][i]):
                print('Results for ', compartments[i], 'are equal at tmax.')
            else:
                print('Results for ',
                      compartments[i], 'are not equal at tmax.')


def main():
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "results")

    timesteps = ['1e-2']

    settings = ['2', '2']

    results = read_data_for_comparison(
        data_dir, timesteps, settings, ide=False)

    compare_data(results, settings)


if __name__ == '__main__':
    main()
