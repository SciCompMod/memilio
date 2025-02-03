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
import subprocess
import os
import math
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
    'start_date': pd.Timestamp('2020.10.01'),
    'end_date': pd.Timestamp('2020.10.01') + pd.DateOffset(days=45),
    'scaleConfirmed': 1.
}


def run_covid_inspired_scenario(result_dir, data_dir, start_date, simulation_time, timestep, scale_contacts):
    """ Run the Covid inspired scenarias defined in ide_covid_inspired_scenario.cpp with the here defined paths and
    parameters.

    @param[in] result_dir Directory where simiulation results are stored.
    @param[in] data_dir Directory where folders with data on contacts and reported data are located. 
    @param[in] start_date Start date of the simulations.
    @param[in] simulation_time Duration of the simulation.
    @param[in] timestep Time step used for the simulations. 
    @param[in] scale_contacts Scaling factor that is applied to the used contact matrix.
    """
    year = start_date.split("-")[0]
    month = start_date.split("-")[1]
    day = start_date.split("-")[2]

    subprocess.call([f"./../../../build/bin/ide_covid_inspired_scenario", data_dir, result_dir,
                     f"{year}", f"{month}", f"{day}", f"{simulation_time}", f"{timestep}", f"{scale_contacts}"])


def get_scale_contacts(files, reported_data_dir, start_date, simulation_time):
    """ Gets the scaling factor for the contacts so that after the first time step of the simulation, the simulation 
    results obtained with the IDE model match the (interpolated) reported number of daily new transmissions (as calculated in 
    load_data() below).

    @param[in] file Expects file with IDE simulation results for flows. 
    @param[in] reported_data_dir Directory where RKI data is stored. 
    @param[in] start_date Start date of interest. 
    @param[in] simulation_time Number of days to be simulated.
    @returns scale_contacts Scaling factor for contacts to scale IDE results to RKI data.
    """

    datafile = os.path.join(
        reported_data_dir, "cases_all_germany_all_dates.json")
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
    date_idx = math.floor(-simulation_time / timestep)

    # Get daily new transmissions from IDE simulation results.
    new_transmissions_ide = total[date_idx, 0] / timestep
    print("")
    print("Comparing daily new transmissions after first time step:")
    print(
        f"IDE new infections at time {dates[date_idx]}: ", new_transmissions_ide)

    # Get daily new transmissions from RKI data at first timestep.
    new_transmissions_rki = data_rki[data_rki["Date"] == start_date]["NewInfectionsDay"].values[0] + timestep * (data_rki[data_rki["Date"] == start_date + pd.DateOffset(
        days=1)]["NewInfectionsDay"].values[0] - data_rki[data_rki["Date"] == start_date]["NewInfectionsDay"].values[0])
    print(
        f"Expected new infections (based on RKI reports) at time {timestep}: {new_transmissions_rki}")

    # Compute scaling factor for contacts.
    scale_contacts = new_transmissions_rki/new_transmissions_ide

    print(f"New scaling factor for contacts is {scale_contacts}.")
    print("")

    return scale_contacts


def load_data(file, start_date, simulation_time):
    """ Loads RKI data and computes the number of mildly symptomatic individuals (stored in column 'InfectedSymptoms'), 
    the number of dead individuals (stored in column 'Deaths') and the number of daily new transmissions (stored in 
    column 'NewInfectionsDay') using scales, dates etc. from the dictionary parameters.
    Method matches the method for computing initial values for the IDE model. 
    See also cpp/models/ide_secir/parameters_io.h.

    @param[in] file Path to the RKI data file for whole Germany. Can be downloaded eg via pycode/memilio-epidata/memilio/epidata/getCaseData.py.
    @param[in] start_date Start date of interest. 
    @param[in] simulation_time Number of days to be simulated.
    @returns df_result Dataframe containing processed RKI data.
    """
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


def plot_daily_new_transmissions(files, reported_data_dir, start_date, simulation_time, fileending="", save_dir=""):
    """ Plots the number of daily new transmissions obtained by simulations with an ODE and an IDE model as well as 
    the number of daily new transmissions extracted from RKI data.

    @param[in] files Expects list of two files with ODE and IDE simulation results for flows, respectively, in this order.
    @param[in] reported_data_dir Directory that contains files with RKI data.
    @param[in] start_date Start date of the simulations.
    @param[in] simulation_time Duration of the simulation.
    @param[in] fileending String that further specifies filename of plot. 
    @param[in] save_dir Directory where plot will be stored. If this is not set, the plot will not be stored. 
    """
    # Read RKI data.
    datafile = os.path.join(
        reported_data_dir, "cases_all_germany_all_dates.json")
    data_rki = load_data(datafile, start_date, simulation_time)

    # Add results to plot.
    fig, ax = plt.subplots()

    legendplot = list(["ODE", "IDE"])

    # Define colors, we use helmholtzdarkblue and helmholtzclaim.
    colors = [(0, 40 / 255, 100 / 255), (20 / 255, 200 / 255, 255 / 255)]
    linestyles = ['-', '-']

    # Plot RKI data.
    ax.scatter(np.linspace(0, simulation_time, simulation_time + 1),
               data_rki["NewInfectionsDay"], marker="x",  s=20, color='gray', label="Extrapolated RKI data")

    # Read simulation results.
    for file in range(len(files)):
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]

        # As there is only one Group, total is the simulation result.
        total = data['Total'][:, :]

        dates = data['Time'][:]
        timestep = np.diff(dates)[0]

        # Plot simulation results.
        # ODE
        if file == 0:
            # Print expected results after first timestep according to RKI (computed by linear extrapolation).
            print("")
            print(f"Expected new transmissions (based on RKI reports) at time {timestep}: ",
                  data_rki[data_rki["Date"] == start_date]["NewInfectionsDay"].values[0] + timestep * (data_rki[data_rki["Date"] == start_date + pd.DateOffset(days=1)]["NewInfectionsDay"].values[0] - data_rki[data_rki["Date"] == start_date]["NewInfectionsDay"].values[0]))

            # The ODE simulation yields cumulative flows. For comparison with IDE model, we first transform cumulative
            # flows to flows per time interval. Then we transform from flows over time interval to flows at time
            # points and add these to plot.
            ax.plot(dates[1:], np.diff(total[:, 0]) / np.diff(dates), label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])

            # Print simulation result after first timestep.
            date_idx = 1
            print(f"ODE new transmissions at time {dates[date_idx]}: ", (np.diff(
                total[:, 0]) / np.diff(dates))[date_idx - 1])

        # IDE
        elif file == 1:
            # The IDE simulation yields flows per time interval. We transform these flows per time interval to flows at
            # time points and add these to plot.
            ax.plot(dates[1:], total[1:, 0] / np.diff(dates), label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])

            # Print simulation result after first timestep.
            date_idx = math.floor(-simulation_time / timestep)
            print(
                f"IDE new transmissions at time {dates[date_idx]}: ", total[date_idx, 0] / timestep)
            print("")

        h5file.close()

    # Adjust plot.
    ax.set_xlim(left=0, right=simulation_time)
    ax.set_ylim(bottom=0)
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

    # Save result.
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f"NewInfections_{fileending}.png",
                    bbox_inches='tight', dpi=500)


def plot_infectedsymptoms_deaths(
        files, reported_data_dir, start_date, simulation_time, fileending="", save_dir=""):
    """ Plots the number of mildly symptomatic individuals and deaths obtained by simulations with an ODE and an IDE model as well as 
    the respective number extracted from RKI data.

    @param[in] files Expects list of two files with ODE and IDE simulation results for compartments, respectively, in 
        this order.
    @param[in] reported_data_dir Directory that contains files with RKI data.
    @param[in] start_date Start date of the simulations.
    @param[in] simulation_time Duration of the simulation.
    @param[in] fileending String that further specifies filename of plot. 
    @param[in] save_dir Directory where plot will be stored. If this is not set, the plot will not be stored. 
    """
    # Define compartment_name and compartment_indices of InfectedSymptoms and Deaths in simulation results.
    compartment_names = ["InfectedSymptoms", "Deaths"]
    compartment_indices = [3, 7]

    # Read RKI data.
    datafile = os.path.join(
        reported_data_dir, "cases_all_germany_all_dates.json")
    data_rki = load_data(datafile, start_date, simulation_time)

    # Add results to plot.

    legendplot = list(["ODE", "IDE"])
    # Define colors, we use helmholtzdarkblue and helmholtzclaim.
    colors = [(0, 40 / 255, 100 / 255), (20 / 255, 200 / 255, 255 / 255)]
    linestyles = ['-', '-']

    for i, compartment in enumerate(compartment_names):

        # Print number of individuals in respective compartment at start_date according to RKI data.
        print(f"{compartment} from RKI on {start_date.strftime('%Y-%m-%d')}: ",
              data_rki[data_rki["Date"] == start_date][compartment].values[0])

        fig, ax = plt.subplots()

        # Plot RKI data.
        ax.scatter(np.linspace(0, simulation_time, simulation_time + 1),
                   data_rki[compartment], marker="x",  s=20, color='gray', label="Extrapolated RKI data")

        # Read simulation results.
        for file in range(len(files)):
            h5file = h5py.File(str(files[file]) + '.h5', 'r')

            if (len(list(h5file.keys())) > 1):
                raise gd.DataError("File should contain one dataset.")
            if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
                raise gd.DataError("Expected only one group.")

            data = h5file[list(h5file.keys())[0]]

            if len(data['Total'][0]) == 8:
                # As there is one Group, total is the simulation
                # result
                total = data['Total'][:, :]
            elif len(data['Total'][0]) == 10:
                # In ODE there are two compartments we do not use, throw these
                # out.
                total = data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]]

            dates = data['Time'][:]

            # Plot simulation results.
            ax.plot(dates, total[:, compartment_indices[i]], label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])

            if file == 0:
                print(f"{compartment} in ODE on {start_date.strftime('%Y-%m-%d')}: ",
                      total[:, compartment_indices[i]][0])

            if file == 1:
                print(f"{compartment} in IDE on {start_date.strftime('%Y-%m-%d')}: ",
                      total[:, compartment_indices[i]][0])

            h5file.close()

        print("")

        # Adjust plot.
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

        # Save result.
        if save_dir != "":
            if not os.path.isdir(save_dir):
                os.makedirs(save_dir)
            plt.savefig(save_dir + f"{compartment}_{fileending}.png",
                        bbox_inches='tight', dpi=500)


def plot_icu(files, reported_data_dir, start_date, simulation_time, fileending="", save_dir=""):
    """ Plots the number of ICU patients obtained by simulations with an ODE and an IDE model as well as 
    the number of ICU patients as reported by DIVI.

    @param[in] files Expects list of two files with ODE and IDE simulation results for compartments, respectively, in 
    this order.
    @param[in] reported_data_dir Directory that contains files with DIVI data.
    @param[in] start_date Start date of the simulations.
    @param[in] simulation_time Duration of the simulation.
    @param[in] fileending String that further specifies filename of plot. 
    @param[in] save_dir Directory where plot will be stored. If this is not set, the plot will not be stored. 
    """

    # Define compartment_name and compartment_index of InfectedCritical in simulation results.
    compartment_name = "InfectedCritical"
    compartment_idx = 5

    # Read DIVI data.
    datafile_divi = os.path.join(
        reported_data_dir,  "germany_divi_all_dates.json")
    df = pd.read_json(datafile_divi)

    # Print number of ICU patients at start_date according to DIVI data.
    print(f"{compartment_name} from DIVI on {start_date.strftime('%Y-%m-%d')}: ",
          df[df["Date"] == start_date]["ICU"].values[0])

    # Add results to plot.
    fig, ax = plt.subplots()

    legendplot = list(["ODE", "IDE"])
    # Define colors, we use helmholtzdarkblue and helmholtzclaim.
    colors = [(0, 40 / 255, 100 / 255), (20 / 255, 200 / 255, 255 / 255)]
    linestyles = ['-', '-']

    # Plot DIVI data.
    ax.scatter(np.linspace(0, simulation_time, simulation_time + 1),
               df[(df["Date"] >= start_date) & (df["Date"] <= start_date + pd.DateOffset(days=simulation_time))]["ICU"], marker="x",  s=20, color='gray', label="DIVI data")

    # Read simulation results.
    for file in range(len(files)):
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

        # Plot simulation results.
        ax.plot(dates, total[:, compartment_idx], label=legendplot[file],
                color=colors[file], linestyle=linestyles[file])

        if file == 0:
            # Print number of ICU patients at start_date according to ODE simulation.
            print(f"{compartment_name} in ODE on {start_date.strftime('%Y-%m-%d')}: ",
                  total[:, compartment_idx][0])

        if file == 1:
            # Print number of ICU patients at start_date according to IDE simulation.
            print(f"{compartment_name} in IDE on {start_date.strftime('%Y-%m-%d')}: ",
                  total[:, compartment_idx][0])

        h5file.close()

    print("")

    # Adjust plot.
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
        plt.savefig(save_dir + f"{compartment_name}_{fileending}.png",
                    bbox_inches='tight', dpi=500)


def plot_covid_inspired_scenario(result_dir, data_dir, plot_dir, start_date, simulation_time, timestep):
    """ Plots daily new transmissions, the number of mildly symptomatic indiciduals and deaths as well as the number
    of ICU patients. See the respective functions for further details. 

    @param[in] result_dir Directory where simiulation results are stored.
    @param[in] data_dir Directory where folders with data on contacts and reported data are located. 
    @param[in] plot_dir Directory where plots will be stored.
    @param[in] start_date Start date of the simulations.
    @param[in] simulation_time Duration of the simulation.
    @param[in] timestep Time step used for the simulations. 
    """
    plot_daily_new_transmissions([os.path.join(result_dir, f"ode_{start_date}_{simulation_time}_{timestep}_flows"),
                                  os.path.join(result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_flows")], data_dir,
                                 pd.Timestamp(start_date), simulation_time,
                                 fileending=f"{start_date}_{simulation_time}_{timestep}", save_dir=plot_dir)

    plot_infectedsymptoms_deaths([os.path.join(result_dir, f"ode_{start_date}_{simulation_time}_{timestep}_compartments"),
                                  os.path.join(result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")],
                                 data_dir,
                                 pd.Timestamp(
                                     start_date), simulation_time,
                                 fileending=f"{start_date}_{simulation_time}_{timestep}",  save_dir=plot_dir)

    plot_icu([os.path.join(result_dir, f"ode_{start_date}_{simulation_time}_{timestep}_compartments"),
              os.path.join(result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")],
             data_dir,
             pd.Timestamp(start_date), simulation_time,  fileending=f"{start_date}_{simulation_time}_{timestep}", save_dir=plot_dir)


def run_and_plot_scenario_with_adjusted_contact_scaling(result_dir, data_dir, plot_dir, start_date, simulation_time, timestep):
    """ Here, we run the scenario with a contact scaling that is adjusted to the data as reported by RKI.
    For this, we first run the scenario with a contact sclaing of 1. We then compare the IDE simulation results with
    the reported data at the first time step. Based on this, we compute the contact scaling such that the IDE 
    simulation results match the reported data. Finally, we run the simulation again with the adjusted contact scaling 
    and plot the results. 

    @param[in] result_dir Directory where simiulation results are stored.
    @param[in] data_dir Directory where folders with data on contacts and reported data are located. 
    @param[in] plot_dir Directory where plots will be stored.
    @param[in] start_date Start date of the simulations.
    @param[in] simulation_time Duration of the simulation.
    @param[in] timestep Time step used for the simulations. 
    """
    # First run the simulation with a contact scaling of 1.
    scale_contacts = 1.
    run_covid_inspired_scenario(result_dir, data_dir, start_date, simulation_time,
                                timestep, scale_contacts)

    # Then determine contact scaling such that IDE results and RKI new transmissions match at first timestep.
    reported_data_dir = data_dir + "/pydata/Germany/"

    scale_contacts = get_scale_contacts([os.path.join(
        result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_flows")], reported_data_dir, pd.Timestamp(start_date), simulation_time)

    # Run simulation with new contact scaling.
    run_covid_inspired_scenario(result_dir, data_dir, start_date, simulation_time,
                                timestep, scale_contacts)

    plot_covid_inspired_scenario(result_dir, reported_data_dir, plot_dir, start_date,
                                 simulation_time, timestep)


def main():
    # Scenario starting on October 1st, 2020.
    start_date = '2020-10-1'
    simulation_time = 45
    timestep = "0.0100"

    # Paths are valid if script is executed e.g. in
    # memilio/cpp/simulations/2024_Wendler_Nonstandard_numerical_scheme_for_integro-differential_model.

    # Path where simulation results (generated with ide_covid_inspired_scenario.cpp) are stored.
    result_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/simulation_results/covid_inspired_scenario/")

    # Path where data for contacts and reported data is stored.
    data_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/")

    # Path where plots will be stored.
    plot_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/plots/covid_inspired_scenario/")

    run_and_plot_scenario_with_adjusted_contact_scaling(result_dir, data_dir, plot_dir, start_date, simulation_time,
                                                        timestep)


if __name__ == "__main__":
    main()
