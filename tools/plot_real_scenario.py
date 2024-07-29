import h5py
import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from memilio.epidata import getDataIntoPandasDataFrame as gd

# Define parameters used for simulation, used for plotting real data.
parameters = {
    'TimeExposed': 4.5,
    'TimeInfectedNoSymptoms': 3.18163,
    'TimeInfectedSymptoms': 7.85313,
    'TimeInfectedSevere': 11.9713,
    'TimeInfectedCritical': 15.2303,
    'TimeInfectedNoSymptomsToInfectedSymptoms': 1.1,
    'TimeInfectedSymptomsToInfectedSevere': 6.60011,
    'TimeInfectedSevereToInfectedCritical': 1.52924,  # 1.52924
    'TimeInfectedCriticalToDead': 10.7,
    'InfectedSymptomsPerInfectedNoSymptoms': 0.698315,
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
    df = df[(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'])))
            & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.ceil(parameters['TimeExposed'] + parameters['TimeInfectedNoSymptomsToInfectedSymptoms'])))]
    # Scale confirmed cases because of undetected infections.
    df['Confirmed'] = parameters['scaleConfirmed'] * df['Confirmed']
    # df2 stores the result of the computation.
    df2 = df.copy()
    df2 = df2[(df['Date'] >= parameters['start_date'])
              & (df['Date'] <= parameters['end_date'])]
    df2 = df2.reset_index()
    df2 = df2.drop(columns=['index', 'Confirmed', 'Deaths'])
    # Calculate individuals in compartment InfectedSymptoms.
    help_I = df['Confirmed'][(df['Date'] >= parameters['start_date'])
                             & (df['Date'] <= parameters['end_date'])].to_numpy()
    help_I = help_I - (1 - math.fmod(parameters['TimeInfectedSymptoms'], 1)) * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptoms'])))
                                                                                               & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptoms'])))].to_numpy()
    help_I = help_I - math.fmod(parameters['TimeInfectedSymptoms'], 1) * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.ceil(
        parameters['TimeInfectedSymptoms']))) & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptoms'])))].to_numpy()
    df2['InfectedSymptoms'] = help_I
    # Calculate number of dead individuals.
    help_D = (1 - (1 - math.fmod(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'], 1))) * df['Deaths'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'])))
                                                                                                                                                                                                         & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'])))].to_numpy()
    help_D = help_D + (1 - math.fmod(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'], 1)) * df['Deaths'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'])))
                                                                                                                                                                                                            & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'])))].to_numpy()
    df2['Deaths'] = help_D
    # df2['Deaths'] = df['Deaths'][(df['Date'] >= parameters['start_date'])
    #                          & (df['Date'] <= parameters['end_date'])].to_numpy()
    # Calculate new infections per day.
    fmod = math.fmod(
        parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'], 1)
    help_newE = fmod * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=math.ceil(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'])))
                                       & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.ceil(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'])))].to_numpy()
    help_newE = help_newE + (1 - 2 * fmod) * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'])))
                                                             & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'])))].to_numpy()
    help_newE = help_newE - (1 - fmod) * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'] - 1)))
                                                         & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'] - 1)))].to_numpy()
    df2['NewInfectionsDay'] = help_newE / \
        parameters['InfectedSymptomsPerInfectedNoSymptoms']
    return df2


def plot_new_infections(files, start_date, simulation_time,
                        legendplot, fileending="", save=True, save_dir='plots/'):

    datafile = os.path.join(os.path.dirname(
        __file__), "..", "data", "pydata", "Germany", "cases_all_germany.json")
    data_rki = load_data(datafile, start_date)

    datafile_ma7 = os.path.join(os.path.dirname(
        __file__), "..", "data", "pydata", "Germany", "cases_all_germany_ma7.json")
    data_rki_ma7 = load_data(datafile, start_date)

    print("New infections from RKI (ma7)  on 1.10.2020: ",
          data_rki_ma7[data_rki_ma7["Date"] == "2020-10-01"]["NewInfectionsDay"].values[0])
    print("New infections from RKI (ma7)  on 2.10.2020: ",
          data_rki_ma7[data_rki_ma7["Date"] == "2020-10-02"]["NewInfectionsDay"].values[0])

    fig, ax = plt.subplots()

    # helmholtzdarkblue, helmholtzclaim
    colors = [(0, 40 / 255, 100 / 255), (20 / 255, 200 / 255, 255 / 255)]
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

        # As there should be only one Group, total is the simulation result
        total = data['Total'][:, :]

        dates = data['Time'][:]

        # get indices where dates are >=0
        # indices = np.where(dates >= 0)
        # plot data
        # ODE
        if file == 0:
            # transform cumulative flows to flows absolute flows
            # then transform from flows over time interval to flows at time
            # points
            ax.plot(dates[1:], np.diff(total[:, 0]) / np.diff(dates), label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])

            date_idx = 1
            print(f"ODE new infections at dates {dates[date_idx]}: ", (np.diff(
                total[:, 0]) / np.diff(dates))[date_idx - 1])
            date_idx = int(1 / np.diff(dates)[0])
            print(f"ODE new infections on {dates[date_idx]}: ", (np.diff(
                total[:, 0]) / np.diff(dates))[date_idx - 1])

        # IDE
        elif file == 1:
            # transform from flows over time interval to flows at time points
            ax.plot(dates[1:], total[1:, 0] / np.diff(dates), label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])

            timestep = np.diff(dates)[0]
            date_idx = int(-30 / timestep - 1)
            print(
                f"IDE new infections on {dates[date_idx]}: ", total[date_idx, 0] / timestep)
            date_idx = int(-30 / timestep)
            print(
                f"IDE new infections on {dates[date_idx]}: ", total[date_idx, 0] / timestep)
            date_idx = int(-29 / timestep - 1)
            print(
                f"IDE new infections at {dates[date_idx]}: ", total[date_idx, 0] / timestep)

        ax.scatter(np.linspace(-20, simulation_time, 51),
                   data_rki["NewInfectionsDay"], s=10)

        h5file.close()

        # ax.set_title(secir_dict[i], fontsize=8)
        # axs[int(i/2), i % 2].set_ylim(bottom=0)
        ax.set_xlim(left=-20, right=30)
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
        plt.savefig(save_dir + f"NewInfections_{fileending}.png",
                    bbox_inches='tight', dpi=500)


def plot_infectedsymptoms_deaths(
        files, start_date, simulation_time, legendplot, fileending="", save=True, save_dir='plots/'):

    datafile = os.path.join(os.path.dirname(
        __file__), "..", "data", "pydata", "Germany", "cases_all_germany.json")
    data_rki = load_data(datafile, start_date)

    datafile_ma7 = os.path.join(os.path.dirname(
        __file__), "..", "data", "pydata", "Germany", "cases_all_germany_ma7.json")
    data_rki_ma7 = load_data(datafile, start_date)

    # print("Infectedsymptoms from RKI (ma7)  on 1.10.2020: ", data_rki_ma7[data_rki_ma7["Date"]=="2020-10-01"]["InfectedSymptoms"].values[0])

    # helmholtzdarkblue, helmholtzclaim
    colors = [(0, 40 / 255, 100 / 255), (20 / 255, 200 / 255, 255 / 255)]
    linestyles = ['-', '--']
    # add results to plot

    compartments = [["InfectedSymptoms", 3], ["Deaths", 7]]

    for compartment in range(len(compartments)):

        print(f"{compartments[compartment][0]} from RKI (ma7)  on 1.10.2020: ",
              data_rki_ma7[data_rki_ma7["Date"] == "2020-10-01"][compartments[compartment][0]].values[0])

        fig, ax = plt.subplots()

        ax.scatter(np.linspace(-20, simulation_time, 51),
                   data_rki[compartments[compartment][0]], s=10)

        for file in range(len(files)):
            # load data
            h5file = h5py.File(str(files[file]) + '.h5', 'r')

            if (len(list(h5file.keys())) > 1):
                raise gd.DataError("File should contain one dataset.")
            if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
                raise gd.DataError("Expected only one group.")

            data = h5file[list(h5file.keys())[0]]

            if len(data['Total'][0]) == 8:
                # As there should be only one Group, total is the simulation
                # result
                total = data['Total'][:, :]
            elif len(data['Total'][0]) == 10:
                # in ODE there are two compartments we don't use, throw these
                # out
                total = data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]]

            dates = data['Time'][:]

            ax.plot(dates, total[:, compartments[compartment][1]], label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])

            if file == 0:
                print(f"{compartments[compartment][0]} in ODE on 01.10.2020: ",
                      total[:, compartments[compartment][1]][0])

            if file == 1:
                print(f"{compartments[compartment][0]} in IDE on 01.10.2020: ",
                      total[:, compartments[compartment][1]][0])

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
            plt.savefig(save_dir + f"{compartments[compartment][0]}_{fileending}.png",
                        bbox_inches='tight', dpi=500)


if __name__ == '__main__':
    # Path to simulation results
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "results/real/")

    legendplot = list(["ODE", "IDE"])

    start_date = '2020-10-01'
    simulation_time = 30
    timestep = "0.1000"

    plot_new_infections([os.path.join(data_dir, f"ode_{start_date}_{simulation_time}_{timestep}_flows"),
                        os.path.join(data_dir, f"ide_{start_date}_{simulation_time}_{timestep}_flows")],
                        pd.Timestamp(start_date), simulation_time,
                        legendplot, fileending=f"{start_date}_{simulation_time}_{timestep}_flows", save=True, save_dir=f"plots/real/{start_date}/")

    plot_infectedsymptoms_deaths([os.path.join(data_dir, f"ode_{start_date}_{simulation_time}_{timestep}_compartments"),
                                  os.path.join(data_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")],
                                 pd.Timestamp(start_date), simulation_time,
                                 legendplot, fileending=f"{start_date}_{simulation_time}_{timestep}", save=True, save_dir=f"plots/real/{start_date}/")

    # plot_real_scenario([os.path.join(data_dir, f"ode_2020-06-01_30_0.1000_flows"),
    #                     os.path.join(data_dir, f"ide_2020-06-01_30_0.1000_flows")],
    #                    pd.Timestamp('2020.06.01'),
    # legendplot, flows=True, fileending="2020-06-01_30_0.1000_flows",
    # save=True, save_dir='plots/real/newinfections/')
