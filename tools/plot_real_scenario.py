import h5py
import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from memilio.epidata import getDataIntoPandasDataFrame as gd

# Define parameters used for simulation, used for plotting real data.
parameters = {
    'TimeExposed':  3.335,
    'TimeInfectedNoSymptoms':  3.31331,
    'TimeInfectedSymptoms': 6.94547,
    'TimeInfectedSevere': 11.634346,
    'TimeInfectedCritical': 17.476959,
    'RecoveredPerInfectedNoSymptoms':  0.206901,
    'start_date': pd.Timestamp('2020.10.01') - pd.DateOffset(days=20),
    'end_date': pd.Timestamp('2020.10.01') + pd.DateOffset(days=30),
    'scaleConfirmed': 2.
}


def load_data(file, start_date):
    """ Loads RKI data and computes 'InfectedSymptoms', 'Deaths' and 'NewInfectionsDay' using scales, dates etc from the dictionary parameters.
    Method matches the method for computing initial values for the LCT model. See also cpp/models/lct_secir/parameters_io.h.
    @param[in] file Path to the RKI data file for whole Germany. Can be downloaded eg via pycode/memilio-epidata/memilio/epidata/getCaseData.py.
    """

    parameters['start_date'] = start_date - pd.DateOffset(days=20)
    parameters['end_date'] = start_date + pd.DateOffset(days=30)
    # Read data.
    df = pd.read_json(file)
    df = df.drop(columns=['Recovered'])

    # Remove unnecessary dates.
    df = df[(df['Date'] >= parameters['start_date']+pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'])))
            & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.ceil(parameters['TimeExposed']+parameters['TimeInfectedNoSymptoms'])))]
    # Scale confirmed cases because of undetected infections.
    df['Confirmed'] = parameters['scaleConfirmed']*df['Confirmed']
    # df2 stores the result of the computation.
    df2 = df.copy()
    df2 = df2[(df['Date'] >= parameters['start_date'])
              & (df['Date'] <= parameters['end_date'])]
    df2 = df2.reset_index()
    df2 = df2.drop(columns=['index', 'Confirmed', 'Deaths'])
    # Calculate individuals in compartment InfectedSymptoms.
    help_I = df['Confirmed'][(df['Date'] >= parameters['start_date'])
                             & (df['Date'] <= parameters['end_date'])].to_numpy()
    help_I = help_I - (1-math.fmod(parameters['TimeInfectedSymptoms'], 1))*df['Confirmed'][(df['Date'] >= parameters['start_date']+pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptoms'])))
                                                                                           & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptoms'])))].to_numpy()
    help_I = help_I - math.fmod(parameters['TimeInfectedSymptoms'], 1) * df['Confirmed'][(df['Date'] >= parameters['start_date']+pd.DateOffset(days=-math.ceil(
        parameters['TimeInfectedSymptoms']))) & (df['Date'] <= parameters['end_date']+pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptoms'])))].to_numpy()
    df2['InfectedSymptoms'] = help_I
    # Calculate number of dead individuals.
    help_D = (1-(1-math.fmod(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'], 1)))*df['Deaths'][(df['Date'] >= parameters['start_date']+pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'])))
                                                                                                                                                       & (df['Date'] <= parameters['end_date']+pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'])))].to_numpy()
    help_D = help_D + (1-math.fmod(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'], 1))*df['Deaths'][(df['Date'] >= parameters['start_date']+pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'])))
                                                                                                                                                            & (df['Date'] <= parameters['end_date']+pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'])))].to_numpy()
    df2['Deaths'] = help_D
    # Calculate new infections per day.
    fmod = math.fmod(
        parameters['TimeInfectedNoSymptoms']+parameters['TimeExposed'], 1)
    help_newE = fmod*df['Confirmed'][(df['Date'] >= parameters['start_date']+pd.DateOffset(days=math.ceil(parameters['TimeInfectedNoSymptoms']+parameters['TimeExposed'])))
                                     & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.ceil(parameters['TimeInfectedNoSymptoms']+parameters['TimeExposed'])))].to_numpy()
    help_newE = help_newE+(1-2*fmod)*df['Confirmed'][(df['Date'] >= parameters['start_date']+pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptoms']+parameters['TimeExposed'])))
                                                     & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptoms']+parameters['TimeExposed'])))].to_numpy()
    help_newE = help_newE-(1-fmod)*df['Confirmed'][(df['Date'] >= parameters['start_date']+pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptoms']+parameters['TimeExposed']-1)))
                                                   & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptoms']+parameters['TimeExposed']-1)))].to_numpy()
    df2['NewInfectionsDay'] = help_newE / \
        (1-parameters['RecoveredPerInfectedNoSymptoms'])
    return df2


def plot_new_infections(files, start_date, legendplot, flows=True, fileending="", save=True, save_dir='plots/'):

    datafile = os.path.join(os.path.dirname(
        __file__), "..", "data", "pydata", "Germany", "cases_all_germany.json")
    data_rki = load_data(datafile, start_date)

    fig, ax = plt.subplots()

    # helmholtzdarkblue, helmholtzclaim
    colors = [(0, 40/255, 100/255), (20/255, 200/255, 255/255)]
    linestyles = ['-', '--']
    # add results to plot
    for file in range(len(files)):
        # load data
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]

        if flows:
            # As there should be only one Group, total is the simulation result
            total = data['Total'][:, :]
        else:
            if len(data['Total'][0]) == 8:
                # As there should be only one Group, total is the simulation result
                total = data['Total'][:, :]
            elif len(data['Total'][0]) == 10:
                # in ODE there are two compartments we don't use, throw these out
                total = data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]]

        dates = data['Time'][:]

        # get indices where dates are >=0
        # indices = np.where(dates >= 0)
        # plot data
        if flows:
            # ODE
            if file == 0:
                # transform cumulative flows to flows absolute flows
                # then transform from flows over time interval to flows at time points
                ax.plot(dates[1:], np.diff(total[:, 0])/np.diff(dates), label=legendplot[file],
                        color=colors[file], linestyle=linestyles[file])
            # IDE
            elif file == 1:
                # transform from flows over time interval to flows at time points
                ax.plot(dates[1:], total[1:, 0]/np.diff(dates), label=legendplot[file],
                        color=colors[file], linestyle=linestyles[file])

            ax.scatter(np.linspace(-20, 30, 51),
                       data_rki["NewInfectionsDay"], s=10)
        else:
            incidence = (total[:-1, 0]-total[1:, 0])/(dates[1:]-dates[:-1])
            ax.plot(dates, incidence, label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])

        h5file.close()

        # ax.set_title(secir_dict[i], fontsize=8)
        # axs[int(i/2), i % 2].set_ylim(bottom=0)
        ax.set_xlim(left=-20)
        ax.grid(True, linestyle='--')
        ax.legend(fontsize=8)

    fig.supxlabel(' Time')
    fig.supylabel('Number of new infections')
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=None, hspace=0.6)

    # save result
    if save:
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f"real_{fileending}.png",
                    bbox_inches='tight', dpi=500)


def plot_infectedsymptoms_deaths(files, start_date, legendplot, fileending="", save=True, save_dir='plots/'):

    datafile = os.path.join(os.path.dirname(
        __file__), "..", "data", "pydata", "Germany", "cases_all_germany.json")
    data_rki_ma = load_data(datafile, start_date)

    # helmholtzdarkblue, helmholtzclaim
    colors = [(0, 40/255, 100/255), (20/255, 200/255, 255/255)]
    linestyles = ['-', '--']
    # add results to plot

    compartments = [["InfectedSymptoms", 3], ["Deaths", 7]]

    for compartment in range(len(compartments)):

        fig, ax = plt.subplots()

        ax.scatter(np.linspace(-20, 30, 51),
                   data_rki_ma[compartments[compartment][0]], s=10)

        for file in range(len(files)):
            # load data
            h5file = h5py.File(str(files[file]) + '.h5', 'r')

            if (len(list(h5file.keys())) > 1):
                raise gd.DataError("File should contain one dataset.")
            if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
                raise gd.DataError("Expected only one group.")

            data = h5file[list(h5file.keys())[0]]

            if len(data['Total'][0]) == 8:
                # As there should be only one Group, total is the simulation result
                total = data['Total'][:, :]
            elif len(data['Total'][0]) == 10:
                # in ODE there are two compartments we don't use, throw these out
                total = data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]]

            dates = data['Time'][:]

            ax.plot(dates, total[:, compartments[compartment][1]], label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])

            h5file.close()

        # ax.set_title(secir_dict[i], fontsize=8)
        # axs[int(i/2), i % 2].set_ylim(bottom=0)
        ax.set_xlim(left=0)
        ax.grid(True, linestyle='--')
        ax.legend(fontsize=8)

        fig.supxlabel(' Time')
        fig.supylabel(
            f'Number of individuals in {compartments[compartment][0]}')
        plt.subplots_adjust(left=None, bottom=None, right=None,
                            top=None, wspace=None, hspace=0.6)

        # save result
        if save:
            if not os.path.isdir(save_dir):
                os.makedirs(save_dir)
            plt.savefig(save_dir + f"real_{compartments[compartment][0]}_{fileending}.png",
                        bbox_inches='tight', dpi=500)


if __name__ == '__main__':
    # Path to simulation results
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "results/real/")

    legendplot = list(["ODE", "IDE"])

    plot_new_infections([os.path.join(data_dir, f"ode_2020-10-01_30_0.1000_flows"),
                        os.path.join(data_dir, f"ide_2020-10-01_30_0.1000_flows")],
                        pd.Timestamp('2020.10.01'),
                        legendplot, flows=True, fileending="2020-10-01_30_0.1000_flows", save=True, save_dir='plots/real/newinfections/')

    plot_infectedsymptoms_deaths([os.path.join(data_dir, f"ode_2020-10-01_30_0.1000_compartments"),
                                  os.path.join(data_dir, f"ide_2020-10-01_30_0.1000_compartments")],
                                 pd.Timestamp('2020.10.01'),
                                 legendplot, fileending="2020-10-01_30_0.1000", save=True, save_dir='plots/real/newinfections/')

    # plot_real_scenario([os.path.join(data_dir, f"ode_2020-06-01_30_0.1000_flows"),
    #                     os.path.join(data_dir, f"ide_2020-06-01_30_0.1000_flows")],
    #                    pd.Timestamp('2020.06.01'),
    #                    legendplot, flows=True, fileending="2020-06-01_30_0.1000_flows", save=True, save_dir='plots/real/newinfections/')
