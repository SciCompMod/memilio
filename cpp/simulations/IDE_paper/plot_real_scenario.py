#############################################################################
# Copyright (C) 2020-2024 German Aerospace Center (DLR-SC)
#
# Authors: Anna Wendler, Lena Ploetzke
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
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from memilio.epidata import getDataIntoPandasDataFrame as gd

# Define parameters used for simulation, also used for plotting reported data.
parameters = {
    'TimeExposed': 4.5,
    'TimeInfectedNoSymptoms':  2.527617,
    'TimeInfectedSymptoms': 7.889900,
    'TimeInfectedSevere': 15.225278,
    'TimeInfectedNoSymptomsToInfectedSymptoms': 1.1,
    'TimeInfectedSymptomsToInfectedSevere': 6.6,
    'TimeInfectedSymptomsToRecovered': 8.0,
    'TimeInfectedSevereToInfectedCritical': 1.5,
    'TimeInfectedCriticalToDead': 10.7,
    'InfectedSymptomsPerInfectedNoSymptoms': 0.793099,
    'SeverePerInfectedSymptoms': 0.078643,
    'start_date': pd.Timestamp('2020.10.01') ,
    'end_date': pd.Timestamp('2020.10.01') + pd.DateOffset(days=45),
    'scaleConfirmed': 1.
}

def load_data(file, start_date, simulation_time):
    """ Loads RKI data and computes the number of mildly symptomatic individuals (stored in column 'InfectedSymptoms'), 
    the number of dead individuals (stored in column 'Deaths') and the number of daily new transmissions (stored in 
    column 'NewInfectionsDay') using scales, dates etc. from the dictionary parameters.
    Method matches the method for computing initial values for the IDE model. 
    See also cpp/models/ide_secir/parameters_io.h.

    @param[in] file Path to the RKI data file for whole Germany. Can be downloaded eg via pycode/memilio-epidata/memilio/epidata/getCaseData.py.
    @param[in] start_date Start date of interest. 
    @param[in] simulation_time Number of days to be simulated.
    """
    # Set start date and end date according to input arguments.
    parameters['start_date'] = start_date
    parameters['end_date'] = start_date + pd.DateOffset(days=simulation_time)

    # Read data into df that will be used for computations below.
    df = pd.read_json(file)
    df = df.drop(columns=['Recovered'])

    # Define df_result where we will store the results of the computations. 
    df_result = df.copy()
    df_result = df_result[(df['Date'] >= parameters['start_date'])
              & (df['Date'] <= parameters['end_date'])]
    df_result = df_result.reset_index()
    df_result = df_result.drop(columns=['index', 'Confirmed', 'Deaths'])

    # Remove dates from df that are not necessary for any computations.
    df = df[(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'])))
            & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.ceil(parameters['TimeExposed'] + parameters['TimeInfectedNoSymptomsToInfectedSymptoms'])))]
    # Scale confirmed cases because of undetected infections.
    df['Confirmed'] = parameters['scaleConfirmed'] * df['Confirmed']

    # Calculate mildly symptomatic individuals, i.e. individuals in compartment InfectedSymptoms.
    help_I = df['Confirmed'][(df['Date'] >= parameters['start_date'])
                             & (df['Date'] <= parameters['end_date'])].to_numpy()
    # Shift according to TimeInfectedSymptomsToInfectedSevere and TimeInfectedSymptomsToRecovered.
    help_I = help_I - parameters["SeverePerInfectedSymptoms"]*((1 - math.fmod(parameters['TimeInfectedSymptomsToInfectedSevere'], 1)) * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptomsToInfectedSevere'])))
                                                                                                                                                        & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptomsToInfectedSevere'])))].to_numpy()) \
                    - (1-parameters["SeverePerInfectedSymptoms"])*((1 - math.fmod(parameters['TimeInfectedSymptomsToRecovered'], 1)) * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptomsToRecovered'])))
                                                                                                                                                       & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptomsToRecovered'])))].to_numpy())
    help_I = help_I - parameters["SeverePerInfectedSymptoms"]*(math.fmod(parameters['TimeInfectedSymptomsToInfectedSevere'], 1) * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.ceil(
        parameters['TimeInfectedSymptomsToInfectedSevere']))) & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptomsToInfectedSevere'])))].to_numpy()) \
        - (1-parameters["SeverePerInfectedSymptoms"])*(math.fmod(parameters['TimeInfectedSymptomsToRecovered'], 1) * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.ceil(
            parameters['TimeInfectedSymptomsToRecovered']))) & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptomsToRecovered'])))].to_numpy())
    df_result['InfectedSymptoms'] = help_I

    # Calculate number of dead individuals.
    help_D = (1 - (1 - math.fmod(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'], 1))) * df['Deaths'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'])))
                                                                                                                                                                                                         & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'])))].to_numpy()
    help_D = help_D + (1 - math.fmod(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'], 1)) * df['Deaths'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'])))
                                                                                                                                                                                                            & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptomsToInfectedSevere'] + parameters['TimeInfectedSevereToInfectedCritical'] + parameters['TimeInfectedCriticalToDead'])))].to_numpy()
    df_result['Deaths'] = help_D

    # Calculate daily new transmissions.
    fmod = math.fmod(
        parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'], 1)
    help_newE = fmod * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=math.ceil(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'])))
                                       & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.ceil(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'])))].to_numpy()
    help_newE = help_newE + (1 - 2 * fmod) * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'])))
                                                             & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'])))].to_numpy()
    help_newE = help_newE - (1 - fmod) * df['Confirmed'][(df['Date'] >= parameters['start_date'] + pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'] - 1)))
                                                         & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.floor(parameters['TimeInfectedNoSymptomsToInfectedSymptoms'] + parameters['TimeExposed'] - 1)))].to_numpy()
    df_result['NewInfectionsDay'] = help_newE / \
        parameters['InfectedSymptomsPerInfectedNoSymptoms']
    
    return df_result


def get_scale_contacts(files, reported_data_dir, start_date, simulation_time):
    """ Gets the scaling factor for the contacts so that the simulation results obtained with the IDE model match the
    reported number of daily new transmissions (as calculated in load_data()).

    @param[in] file Expects file with IDE simulation results for flows. 
    @param[in] data_dir Directory where RKI data is stored. 
    @param[in] start_date Start date of interest. 
    @param[in] simulation_time Number of days to be simulated.
    """

    datafile = os.path.join(reported_data_dir, "cases_all_germany.json")
    data_rki = load_data(datafile, start_date, simulation_time)

    # Load IDE data.
    for file in range(len(files)):
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]

        # As there should be only one Group, total is the simulation result
        total = data['Total'][:, :]

        dates = data['Time'][:]
        timestep = np.diff(dates)[0]

    # Date index referring to first time step of IDE simulation.
    date_idx = int(-simulation_time / timestep)

    # Get daily new transmissions from IDE simulation results. 
    new_transmissions_ide = total[date_idx, 0] / timestep
    print(f"IDE new infections on {dates[date_idx]}: ", new_transmissions_ide)

    # Get daily new transmissions from RKI data at first timestep.
    new_transmissions_rki = data_rki[data_rki["Date"] == start_date]["NewInfectionsDay"].values[0] + timestep * (data_rki[data_rki["Date"] == start_date + pd.DateOffset(
        days=1)]["NewInfectionsDay"].values[0] - data_rki[data_rki["Date"] == start_date]["NewInfectionsDay"].values[0])
    print(f"Expected new infections at {timestep}: {new_transmissions_rki}")

    # Compute scaling factor for contacts.
    scale_contacts = new_transmissions_rki/new_transmissions_ide

    return scale_contacts


def get_scale_confirmed_cases(files, reported_data_dir, start_date):
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

    icu_ide = total[:, 5][0]
    print(f"ICU in IDE on {start_date}: ",
          total[:, 5][0])

    datafile_icu = os.path.join(reported_data_dir, "germany_divi.json")

    df = pd.read_json(datafile_icu)

    icu_divi = df[df["Date"] == start_date]["ICU"].values[0]
    print(f"ICU from DIVI (ma7)  on {start_date}: ",
          df[df["Date"] == start_date]["ICU"].values[0])

    scale_confirmed_cases = icu_divi/icu_ide

    return scale_confirmed_cases


def plot_daily_new_transmissions(files, reported_data_dir, start_date, simulation_time, fileending="", save_dir=""):

    datafile = os.path.join(reported_data_dir, "cases_all_germany.json")
    data_rki = load_data(datafile, start_date, simulation_time)

    fig, ax = plt.subplots()

    ax.scatter(np.linspace(0, simulation_time, simulation_time + 1),
               data_rki["NewInfectionsDay"], marker="x",  s=20, color='gray', label="Extrapolated RKI data")

    legendplot = list(["ODE", "IDE"])
    # helmholtzdarkblue, helmholtzclaim
    colors = [(0, 40 / 255, 100 / 255), (20 / 255, 200 / 255, 255 / 255)]
    linestyles = ['-', '-']
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
        timestep = np.diff(dates)[0]

        # get indices where dates are >=0
        # indices = np.where(dates >= 0)
        # plot data
        # ODE
        if file == 0:
            print(f"New infections from RKI  on {start_date}: ",
                  data_rki[data_rki["Date"] == start_date]["NewInfectionsDay"].values[0])
            print(f"Expected new infections at {timestep}: ",
                  data_rki[data_rki["Date"] == start_date]["NewInfectionsDay"].values[0] + timestep * (data_rki[data_rki["Date"] == start_date + pd.DateOffset(days=1)]["NewInfectionsDay"].values[0] - data_rki[data_rki["Date"] == start_date]["NewInfectionsDay"].values[0]))

            print(f"New infections from RKI  on {start_date + pd.DateOffset(days=1)}: ",
                  data_rki[data_rki["Date"] == start_date + pd.DateOffset(days=1)]["NewInfectionsDay"].values[0])
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

            date_idx = int(-simulation_time / timestep - 1)
            print(
                f"IDE new infections on {dates[date_idx]}: ", total[date_idx, 0] / timestep)
            date_idx = int(-simulation_time / timestep)
            print(
                f"IDE new infections on {dates[date_idx]}: ", total[date_idx, 0] / timestep)
            date_idx = int(-(simulation_time - 1) / timestep - 1)
            print(
                f"IDE new infections at {dates[date_idx]}: ", total[date_idx, 0] / timestep)

        h5file.close()

    ax.set_xlim(left=0, right=simulation_time)
    ax.set_ylim(bottom=0)
    if start_date == pd.Timestamp("2020-6-1"):
        ax.set_ylim(bottom=0, top=1000)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(fontsize=12)

    # Define x-ticks.
    datelist = np.array(pd.date_range(parameters["start_date"].date(),
                                      periods=simulation_time+1, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int((simulation_time) / 5) + 1) * 5)
    plt.xticks(tick_range, datelist[tick_range],
               rotation=45, fontsize=12)
    plt.xticks(np.arange(simulation_time), minor=True)

    fig.supxlabel('Date')
    fig.supylabel(r'Daily new transmissions')
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=None, hspace=0.6)

    plt.tight_layout()

    # save result
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f"NewInfections_{fileending}.png",
                    bbox_inches='tight', dpi=500)


def plot_infectedsymptoms_deaths(
        files, reported_data_dir, start_date, simulation_time, fileending="", save_dir=""):

    datafile = os.path.join(reported_data_dir, "cases_all_germany.json")
    data_rki = load_data(datafile, start_date, simulation_time)

    datafile_ma7 = os.path.join(reported_data_dir,  "cases_all_germany_ma7.json")
    data_rki_ma7 = load_data(datafile_ma7, start_date,
                             simulation_time)

    legendplot = list(["ODE", "IDE"])
    # helmholtzdarkblue, helmholtzclaim
    colors = [(0, 40 / 255, 100 / 255), (20 / 255, 200 / 255, 255 / 255)]
    linestyles = ['-', '-']
    # add results to plot

    compartments = [["InfectedSymptoms", 3], ["Deaths", 7]]

    for compartment in range(len(compartments)):

        print(f"{compartments[compartment][0]} from RKI on {start_date}: ",
              data_rki_ma7[data_rki["Date"] == start_date][compartments[compartment][0]].values[0])

        fig, ax = plt.subplots()

        ax.scatter(np.linspace(0, simulation_time, simulation_time + 1),
                   data_rki[compartments[compartment][0]], marker="x",  s=20, color='gray', label="Extrapolated RKI data")

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
                print(f"{compartments[compartment][0]} in ODE on {start_date}: ",
                      total[:, compartments[compartment][1]][0])

            if file == 1:
                print(f"{compartments[compartment][0]} in IDE on {start_date}: ",
                      total[:, compartments[compartment][1]][0])

            h5file.close()

        ax.set_xlim(left=0, right=simulation_time)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(fontsize=12)

        # Define x-ticks.
        datelist = np.array(pd.date_range(parameters["start_date"].date(),
                                          periods=simulation_time+1, freq='D').strftime('%m-%d').tolist())
        tick_range = (np.arange(int((simulation_time) / 5) + 1) * 5)
        plt.xticks(tick_range, datelist[tick_range],
                   rotation=45, fontsize=12)
        plt.xticks(np.arange(simulation_time), minor=True)

        fig.supxlabel('Date')
        if compartment == 0:
            fig.supylabel(
                f'Mildly symptomatic individuals')
        if compartment == 1:
            fig.supylabel(
                f'Deaths')
        plt.subplots_adjust(left=None, bottom=None, right=None,
                            top=None, wspace=None, hspace=0.6)
        plt.tight_layout()

        # save result
        if save_dir != "":
            if not os.path.isdir(save_dir):
                os.makedirs(save_dir)
            plt.savefig(save_dir + f"{compartments[compartment][0]}_{fileending}.png",
                        bbox_inches='tight', dpi=500)


def plot_icu(files, reported_data_dir, start_date, simulation_time, fileending="", save_dir=""):

    # datafile_icu_ma7 = os.path.join(os.path.dirname(
    #     __file__), "..", "data", "pydata", "Germany", "germany_divi_ma7.json")

    datafile_icu = os.path.join(reported_data_dir,  "germany_divi.json")

    df = pd.read_json(datafile_icu)
    # data_icu_ma7 = load_data(datafile_icu_ma7, start_date, simulation_time)

    # print("Infectedsymptoms from RKI (ma7)  on 1.10.2020: ", data_rki_ma7[data_rki_ma7["Date"]=="2020-10-01"]["InfectedSymptoms"].values[0])

    # helmholtzdarkblue, helmholtzclaim
    legendplot = list(["ODE", "IDE"])
    colors = [(0, 40 / 255, 100 / 255), (20 / 255, 200 / 255, 255 / 255)]
    linestyles = ['-', '-']
    # add results to plot

    compartments = [["InfectedCritical", 5]]

    for compartment in range(len(compartments)):

        print(f"{compartments[compartment][0]} from DIVI (ma7)  on {start_date}: ",
              df[df["Date"] == start_date]["ICU"].values[0])

        fig, ax = plt.subplots()

        ax.scatter(np.linspace(0, simulation_time, simulation_time + 1),
                   df[(df["Date"] >= start_date) & (df["Date"] <= start_date + pd.DateOffset(days=simulation_time))]["ICU"], marker="x",  s=20, color='gray', label="DIVI data")

        for file in range(len(files)):
            # load data
            h5file = h5py.File(str(files[file]) + '.h5', 'r')

            if (len(list(h5file.keys())) > 1):
                raise gd.DataError("File should contain one dataset.")
            if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
                raise gd.DataError("Expected only one group.")

            data = h5file[list(h5file.keys())[0]]

            if len(data['Total'][0]) == 8:
                # As there should be only one Group, total is the simulation result.
                total = data['Total'][:, :]
            elif len(data['Total'][0]) == 10:
                # In ODE model there are two compartments we do not use, throw these out.
                total = data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]]

            dates = data['Time'][:]

            ax.plot(dates, total[:, compartments[compartment][1]], label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])

            if file == 0:
                print(f"{compartments[compartment][0]} in ODE on {start_date}: ",
                      total[:, compartments[compartment][1]][0])

            if file == 1:
                print(f"{compartments[compartment][0]} in IDE on {start_date}: ",
                      total[:, compartments[compartment][1]][0])

            h5file.close()

        ax.set_xlim(left=0, right=simulation_time)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(fontsize=12)

        # Define x-ticks.
        datelist = np.array(pd.date_range(parameters["start_date"].date(),
                                          periods=simulation_time+1, freq='D').strftime('%m-%d').tolist())
        tick_range = (np.arange(int((simulation_time) / 5) + 1) * 5)
        plt.xticks(tick_range, datelist[tick_range],
                   rotation=45, fontsize=12)
        plt.xticks(np.arange(simulation_time), minor=True)

        fig.supxlabel('Date')
        fig.supylabel(
            f'ICU patients')
        plt.subplots_adjust(left=None, bottom=None, right=None,
                            top=None, wspace=None, hspace=0.6)

        plt.tight_layout()

        # Save result.
        if save_dir != "":
            if not os.path.isdir(save_dir):
                os.makedirs(save_dir)
            plt.savefig(save_dir + f"{compartments[compartment][0]}_{fileending}.png",
                        bbox_inches='tight', dpi=500)


def main():
    # Paths are valid if file is executed e.g. in memilio/cpp/simulations/IDE_paper.
    # Path where simulation results (generated with ide_real_scenario.cpp) are stored. 
    result_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/simulation_results/covid_inspired_scenario/")
    
    # Path where RKI and DIVI data are stored.
    reported_data_dir =  os.path.join(os.path.dirname(
        __file__), "../../..", "data/pydata/Germany/")
    
    # Path where plots will be stored. 
    plot_dir =  os.path.join(os.path.dirname(
        __file__), "../../..", "data/plots/covid_inspired_scenario/")

    # Define start_date, simulation_time and timestep.
    start_date = '2020-10-1'
    simulation_time = 45
    timestep = "0.0100"

    # Plot daily new transmissions.
    plot_daily_new_transmissions([os.path.join(result_dir, f"ode_{start_date}_{simulation_time}_{timestep}_flows"),
                        os.path.join(result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_flows")], reported_data_dir,
                        pd.Timestamp(start_date), simulation_time,
                        fileending=f"{start_date}_{simulation_time}_{timestep}",save_dir=plot_dir)

    # Plot mildly symptomatic and dead individuals. 
    plot_infectedsymptoms_deaths([os.path.join(result_dir, f"ode_{start_date}_{simulation_time}_{timestep}_compartments"),
                                  os.path.join(result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")], reported_data_dir,
                                 pd.Timestamp(
                                     start_date), simulation_time,
                                 fileending=f"{start_date}_{simulation_time}_{timestep}", save_dir=plot_dir)

    # Plot ICU patients. 
    plot_icu([os.path.join(result_dir, f"ode_{start_date}_{simulation_time}_{timestep}_compartments"),
              os.path.join(result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")],
              reported_data_dir,
             pd.Timestamp(start_date), simulation_time, fileending=f"{start_date}_{simulation_time}_{timestep}", save_dir=plot_dir)

if __name__ == '__main__':
    main()
