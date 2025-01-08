#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Lena Ploetzke
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
"""@plot_covid19_inspired.py
Functions to create the plots for section 4.5 of the paper.

The simulation results to be plotted should be stored in a '../../data/simulation_lct_covid19' folder 
as .h5 files. Have a look at the README for an explanation of how to create the simulation result data.
Additionally, reported RKI case data and DIVI data have to be downloaded beforehand.
See the README for further instructions on how to download the datasets.

Please note that the memilio.epidata package needs to be installed beforehand. 
Have a look at the pycode and the pycode/memilio-epidata READMEs.
"""

import h5py
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

import memilio.epidata.getDataIntoPandasDataFrame as gd

# Define compartments and RKI age groups.
secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
              5: 'ICU', 6: 'Recovered', 7: 'Dead'}
Age_RKI_names = {0: 'A00-A04', 1: 'A05-A14',
                 2: 'A15-A34',  3: 'A35-A59', 4: 'A60-A79', 5: 'A80+'}

# Define parameters used for plotting real data.
TimeExposed = [3.335, 3.335, 3.335, 3.335, 3.335, 3.335]
TimeInfectedNoSymptoms = [2.74, 2.74, 2.565, 2.565, 2.565, 2.565]
TimeInfectedSymptoms = [7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775]
TimeInfectedSevere = [5., 5., 5.925, 7.55, 8.5, 11.]
TimeInfectedCritical = [6.95, 6.95, 6.86, 17.36, 17.1, 11.6]
RecoveredPerInfectedNoSymptoms = [
    1 - 0.75, 1 - 0.75, 1 - 0.8, 1 - 0.8, 1 - 0.8, 1 - 0.8]
age_group_sizes = [3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434]

# Define color to be used while plotting for different models to make plots consistent.
color_dict = {'ODE': '#1f77b4',
              'LCT3': '#2ca02c',
              'LCT10': '#ff7f0e',
              'LCT20': '#9467bd',
              'LCT50':  '#e377c2',
              'LCTvar':  '#bbf90f',
              }
fontsize_labels = 14
fontsize_legends = 11
plotfolder = 'Plots/Plots_covid19'


def load_rki_data(file, start_date, tmax, scale_confirmed_cases=1.):
    """ Loads RKI data and computes 'InfectedSymptoms', 'Deaths' and 'DailyNewTransmissions'. 
    The result is an age resolved data frame with the data from start_date to start_date+tmax. 
    Method matches the method for computing age resolved initial values for the LCT model. 
    See also cpp/models/lct_secir/parameters_io.h.

    @param[in] file: Path to the RKI data file.
    @param[in] start_date: The data will be extracted from start_date on.
    @param[in] tmax: The time frame the data should be extracted for in days. 
    @param[in] scale_confirmed_cases: Scaling factor for the number of confirmed cases 
        to incorporate a detection ratio.
    """
    # Read data.
    df = pd.read_json(file)
    df = df.drop(columns=['Recovered'])

    end_date = start_date+pd.DateOffset(days=tmax)

    # Remove unnecessary dates (some days before start_date and some days after end_date are required for the
    # computation of the disease state values).
    df = df[(df['Date'] >= start_date +
             pd.DateOffset(days=-math.ceil(max(TimeInfectedSymptoms) + max(TimeInfectedSevere) + max(TimeInfectedCritical))))
            & (df['Date'] <= end_date + pd.DateOffset(days=math.ceil(max(TimeExposed) + max(TimeInfectedNoSymptoms))))]
    # Scale confirmed cases to incorporate a detection ratio.
    df['Confirmed'] = scale_confirmed_cases * df['Confirmed']
    # As unknown age groups are omitted in the cpp epi_data used for initialization, we omit it here too.
    df = df[(df['Age_RKI'] != 'unknown')]

    # df2 stores the result of the computation. Define suitable columns and rows first.
    df2 = df.copy()
    df2 = df2[(df['Date'] >= start_date) & (df['Date'] <= end_date)]
    df2 = df2.reset_index()
    df2 = df2.drop(columns=['index', 'Confirmed', 'Deaths'])
    df2['InfectedSymptoms'] = pd.Series(dtype='double')
    df2['Deaths'] = pd.Series(dtype='double')
    df2['DailyNewTransmissions'] = pd.Series(dtype='double')

    for age_group in range(len(Age_RKI_names)):
        # Calculate individuals in compartment InfectedSymptoms using TimeInfectedSymptoms.
        df_age = df[(df['Age_RKI'] == Age_RKI_names[age_group])]
        help_I_age = df_age['Confirmed'][(df_age['Date'] >= start_date) & (
            df_age['Date'] <= end_date)].to_numpy()

        help_I_age = help_I_age - (1 - math.fmod(TimeInfectedSymptoms[age_group], 1)) * df_age['Confirmed'][(
            df_age['Date'] >= start_date + pd.DateOffset(days=-math.floor(TimeInfectedSymptoms[age_group]))) & (
            df_age['Date'] <= end_date + pd.DateOffset(days=-math.floor(TimeInfectedSymptoms[age_group])))].to_numpy()

        help_I_age = help_I_age - math.fmod(TimeInfectedSymptoms[age_group], 1) * df_age['Confirmed'][(
            df_age['Date'] >= start_date + pd.DateOffset(days=-math.ceil(TimeInfectedSymptoms[age_group]))) & (
            df_age['Date'] <= end_date + pd.DateOffset(days=-math.ceil(TimeInfectedSymptoms[age_group])))].to_numpy()

        df2.loc[(df2['Age_RKI'] == Age_RKI_names[age_group]),
                "InfectedSymptoms"] = help_I_age

        # Calculate number of dead individuals.
        # Compute Dead by shifting RKI data according to mean stay times.
        # This is done because the RKI reports death with the date of positive test instead of the date of deaths.
        timeInfSym_Sev_Crit_age = TimeInfectedSymptoms[age_group] + \
            TimeInfectedSevere[age_group] + TimeInfectedCritical[age_group]
        help_D_age = (1 - (1 - math.fmod(timeInfSym_Sev_Crit_age, 1))) * df_age['Deaths'][(
            df_age['Date'] >= start_date + pd.DateOffset(days=-math.ceil(timeInfSym_Sev_Crit_age)))
            & (df_age['Date'] <= end_date + pd.DateOffset(
                days=-math.ceil(timeInfSym_Sev_Crit_age)))].to_numpy()
        help_D_age = help_D_age + (1 - math.fmod(timeInfSym_Sev_Crit_age, 1)) * df_age['Deaths'][
            (df_age['Date'] >= start_date +
             pd.DateOffset(days=-math.floor(timeInfSym_Sev_Crit_age)))
            & (df_age['Date'] <= end_date + pd.DateOffset(days=-math.floor(timeInfSym_Sev_Crit_age)))].to_numpy()

        df2.loc[(df2['Age_RKI'] == Age_RKI_names[age_group]),
                'Deaths'] = help_D_age

        # Calculate new transmissions per day.
        timeE_INS_age = TimeInfectedNoSymptoms[age_group] + \
            TimeExposed[age_group]
        fmod = math.fmod(timeE_INS_age, 1)
        help_newTrans_age = fmod * df_age['Confirmed'][(df_age['Date'] >= start_date + pd.DateOffset(days=math.ceil(
            timeE_INS_age)))
            & (df_age['Date'] <= end_date + pd.DateOffset(days=math.ceil(
               timeE_INS_age)))].to_numpy()
        help_newTrans_age = help_newTrans_age + (1 - 2 * fmod) * df_age['Confirmed'][(
            df_age['Date'] >= start_date + pd.DateOffset(days=math.floor(timeE_INS_age)))
            & (df_age['Date'] <= end_date + pd.DateOffset(days=math.floor(timeE_INS_age)))].to_numpy()
        help_newTrans_age = help_newTrans_age-(1-fmod)*df_age['Confirmed'][(df_age['Date'] >= start_date + pd.DateOffset(days=math.floor(timeE_INS_age - 1)))
                                                                           & (df_age['Date'] <= end_date + pd.DateOffset(days=math.floor(timeE_INS_age - 1)))].to_numpy()

        df2.loc[(df2['Age_RKI'] == Age_RKI_names[age_group]),
                'DailyNewTransmissions'] = help_newTrans_age / (1-RecoveredPerInfectedNoSymptoms[age_group])
    return df2


def load_divi_data(file, start_date, tmax):
    """ Loads the DIVI data. This data set is not resolved by age.
    The result is a data frame with the DIVI data from start_date to start_date+tmax. 

    @param[in] file: Path to the DIVI data file.
    @param[in] start_date: The data will be extracted from start_date on.
    @param[in] tmax: The time frame the data should be extracted for in days. 
    """
    # Load data.
    df = pd.read_json(file)

    # Remove unnecessary dates.
    end_date = start_date+pd.DateOffset(days=tmax)
    df = df[(df['Date'] >= start_date) & (df['Date'] <= end_date)]
    return df


def plot_InfectedSymptoms_or_Deaths(files, datafile, start_date, tmax, scale_confirmed_cases, legend_labels, plot_deaths=False, file_name="", age_group=-1):
    """ Plots simulation results compared with real data for the compartments Deaths or InfectedSymptoms.
    The simulation results should consist of accumulated numbers for subcompartments in case of an LCT model.
    The function can plot results for one age group or results that are accumulated by age. 
    The simulation results should be suitable.

    @param[in] files: List of paths to files (without .h5 extension) containing simulation results to be plotted.
         Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). 
    @param[in] datafile: Path to the RKI data file.
    @param[in] start_date: Start date of the simulation.
    @param[in] tmax: The time frame to be plotted is [start_date, start_date+tmax]. 
    @param[in] scale_confirmed_cases: Scaling factor for the number of confirmed cases 
        to incorporate a detection ratio.
    @param[in] legend_labels: List of names for the results to be used for the plot legend.
    @param[in] plot_deaths: If False, InfectedSymptoms compartment is plotted, if True, deaths are plotted.
    @param[in] file_name: The name of the file where the plot will be saved.
            If an empty string is provided, the plot will not be saved.
    @param[in] age_group: Index of considered age group. 
            If the default value -1 is provided, non-age-resolved results are plotted.
    """
    plt.figure(file_name)

    data_rki = load_rki_data(datafile, start_date, tmax, scale_confirmed_cases)
    num_days = tmax + 1

    if (age_group == -1):
        # Non-age-resolved.
        data_rki = data_rki.drop(columns=['Age_RKI'])
        data_rki = data_rki.groupby(['Date']).sum()
        if (plot_deaths):
            plt.plot(range(num_days), data_rki['Deaths'],
                     linestyle='None', color='grey', marker='x', markersize=5)
            compartment_idx = 7
            plt.ylabel("Deaths", fontsize=fontsize_labels)
        else:
            plt.plot(range(num_days), data_rki['InfectedSymptoms'],
                     linestyle='None', color='grey', marker='x', markersize=5)
            compartment_idx = 3
            plt.ylabel("Mildly symptomatic individuals",
                       fontsize=fontsize_labels)
    else:
        # Age-resolved.
        plt.title(Age_RKI_names[age_group])
        if (plot_deaths):
            plt.plot(range(num_days), data_rki['Deaths'][(data_rki['Age_RKI'] == Age_RKI_names[age_group])],
                     linestyle='None', color='grey', marker='x', markersize=5)
            compartment_idx = 7
            plt.ylabel("Deaths", fontsize=fontsize_labels)
        else:
            plt.plot(range(num_days), data_rki['InfectedSymptoms'][(data_rki['Age_RKI'] == Age_RKI_names[age_group])],
                     linestyle='None', color='grey', marker='x', markersize=5)
            compartment_idx = 3
            plt.ylabel("Mildly symptomatic individuals",
                       fontsize=fontsize_labels)

    # Add simulation results to plot.
    for file in range(len(files)):
        # Load data.
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result.
        total = data['Total'][:, :]

        if age_group == -1:
            if (total.shape[1] != len(secir_dict)):
                raise gd.DataError(
                    "Expected a different number of compartments.")
            # Plot result.
            plt.plot(dates, total[:, compartment_idx],
                     linewidth=1.2, linestyle="solid", color=color_dict[legend_labels[1+file]])
        else:
            if (total.shape[1] != len(secir_dict)*len(Age_RKI_names)):
                raise gd.DataError(
                    "Expected a different number of compartments.")
            # Plot result.
            plt.plot(dates, total[:, len(secir_dict) * age_group + compartment_idx],
                     linewidth=1.2, linestyle="solid", color=color_dict[legend_labels[1+file]])

        h5file.close()

    plt.xlabel('Date', fontsize=fontsize_labels)
    plt.yticks(fontsize=9)
    plt.xlim(left=0, right=num_days-1)
    # Define x-ticks as dates.
    datelist = np.array(pd.date_range(start_date.date(),
                                      periods=num_days, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int((num_days-1) / 5) + 1) * 5)
    plt.xticks(tick_range, datelist[tick_range],
               rotation=45, fontsize=12)
    plt.xticks(np.arange(num_days), minor=True)

    plt.legend(legend_labels, fontsize=fontsize_legends, framealpha=0.5)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # Save result.
    if file_name:
        plt.savefig(plotfolder+'/'+file_name +
                    '.png', bbox_inches='tight', dpi=500)
    plt.close()


def plot_InfectedCritical(files, datafile, start_date, tmax, legend_labels, file_name=""):
    """ Plots simulation results compared with real DIVI data for number of patients in intensive care units.
    The simulation results should consist of accumulated numbers for subcompartments in case of an LCT model
    and accumulated numbers for age groups.

   @param[in] files: List of paths to files (without .h5 extension) containing (non-age-resolved) 
        simulation results to be plotted. Results should contain exactly 8 compartments 
        (so use accumulated numbers for the subcompartments of LCT models). 
    @param[in] datafile: Path to the DIVI data file.
    @param[in] start_date: Start date of the simulation.
    @param[in] tmax: The time frame to be plotted is [start_date, start_date+tmax]. 
    @param[in] scale_confirmed_cases: Scaling factor for the number of confirmed cases 
        to incorporate a detection ratio.
    @param[in] legend_labels: List of names for the results to be used for the plot legend.
    @param[in] file_name: The name of the file where the plot will be saved.
            If an empty string is provided, the plot will not be saved.
    """
    plt.figure(file_name)

    data_icu = load_divi_data(datafile, start_date, tmax)
    num_days = tmax + 1

    # Plot ICU data.
    plt.plot(range(num_days), data_icu['ICU'],
             linestyle='None', color='grey', marker='x', markersize=5)

    # Index of ICU compartment in simulation results.
    compartment_idx = 5

    # Add simulation results to plot.
    for file in range(len(files)):
        # Load data.
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result.
        total = data['Total'][:, :]

        if (total.shape[1] != len(secir_dict)):
            raise gd.DataError(
                "Expected a different number of compartments.")
        # Plot result.
        plt.plot(dates, total[:, compartment_idx],
                 linewidth=1.2, linestyle="solid", color=color_dict[legend_labels[1+file]])
        h5file.close()

    plt.xlabel("Date", fontsize=fontsize_labels)
    plt.ylabel("ICU patients", fontsize=fontsize_labels)
    plt.yticks(fontsize=9)
    plt.xlim(left=0, right=num_days-1)
    # Define x-ticks as dates.
    datelist = np.array(pd.date_range(start_date.date(),
                                      periods=num_days, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int((num_days-1) / 5) + 1) * 5)
    plt.xticks(tick_range, datelist[tick_range],
               rotation=45, fontsize=12)
    plt.xticks(np.arange(num_days), minor=True)

    plt.legend(legend_labels, fontsize=12, framealpha=0.5)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # Save result.
    if file_name:
        plt.savefig(plotfolder+'/'+file_name +
                    '.png', bbox_inches='tight', dpi=500)
    plt.close()


def plot_daily_new_transmissions(files, datafile, start_date, tmax, scale_confirmed_cases, legend_labels, file_name="", age_group=-1):
    """ Plots simulation results compared with real data for the daily new transmissions
    (equals the number of people transiting from Susceptible to Exposed within one day).
    The simulation results should consist of accumulated numbers for subcompartments in case of an LCT model.
    The function can plot results for one age group or results that are accumulated by age. 
    The simulation results should be suitable.

    @param[in] files: List of paths to files (without .h5 extension) containing simulation results to be plotted.
         Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). 
    @param[in] datafile: Path to the RKI data file.
    @param[in] start_date: Start date of the simulation.
    @param[in] tmax: The time frame to be plotted is [start_date, start_date+tmax]. 
    @param[in] scale_confirmed_cases: Scaling factor for the number of confirmed cases 
        to incorporate a detection ratio.
    @param[in] legend_labels: List of names for the results to be used for the plot legend.
    @param[in] file_name: The name of the file where the plot will be saved.
            If an empty string is provided, the plot will not be saved.
    @param[in] age_group: Index of considered age group. 
            If the default value -1 is provided, non-age-resolved results are plotted.
    """
    plt.figure(file_name)

    data_rki = load_rki_data(datafile, start_date, tmax, scale_confirmed_cases)
    num_days = tmax+1

    if (age_group == -1):
        # Non-age-resolved.
        data_rki = data_rki.drop(columns=['Age_RKI'])
        data_rki = data_rki.groupby(['Date']).sum()
        plt.plot(range(num_days), data_rki['DailyNewTransmissions'],
                 linestyle='None', color='grey', marker='x', markersize=5)
        print("Daily new transmissions at the first day of RKI data is: " +
              f"{data_rki.loc[start_date]['DailyNewTransmissions']}")
    else:
        # Age-resolved.
        plt.plot(range(num_days), data_rki['DailyNewTransmissions'][(data_rki['Age_RKI'] == Age_RKI_names[age_group])],
                 linestyle='None', color='grey', marker='x', markersize=5)
        plt.title(Age_RKI_names[age_group])
        print("Daily new transmissions at the first day of " +
              Age_RKI_names[age_group]+" of the RKI data is: "+f"{data_rki['DailyNewTransmissions'][(data_rki['Age_RKI'] == Age_RKI_names[age_group]) & ( data_rki['Date'] == start_date )].iloc[0]}")

    # Add simulation results to plot.
    for file in range(len(files)):
        # Load data.
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result.
        total = data['Total'][:, :]
        if age_group == -1:
            incidence = (total[:-1, 0] -
                         total[1:, 0])/(dates[1:]-dates[:-1])
        else:
            incidence = (total[:-1, len(secir_dict) * age_group] -
                         total[1:, len(secir_dict) * age_group])/(dates[1:]-dates[:-1])
        print("Daily new transmissions at the first day of the simulation " + legend_labels[file+1]+" is: "
              f"{incidence[0]}")

        # Plot result.
        if legend_labels[file] in color_dict:
            plt.plot(dates[1:], incidence, linewidth=1.2,
                     linestyle="solid", color=color_dict[legend_labels[1+file]])
        else:
            plt.plot(dates[1:], incidence, linewidth=1.2)
        h5file.close()

    plt.ylabel('Daily new transmissions', fontsize=fontsize_labels)
    plt.ylim(bottom=0)
    plt.xlabel('Date', fontsize=fontsize_labels)
    plt.xlim(left=0, right=num_days-1)
    # Define x-ticks as dates.
    datelist = np.array(pd.date_range(start_date.date(),
                                      periods=num_days, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int((num_days - 1) / 5) + 1) * 5)
    plt.xticks(tick_range, datelist[tick_range],
               rotation=45, fontsize=12)
    plt.xticks(np.arange(num_days), minor=True)

    plt.legend(legend_labels, fontsize=fontsize_legends, framealpha=0.5)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # Save result.
    if file_name:
        plt.savefig(plotfolder+'/'+file_name +
                    '.png', bbox_inches='tight', dpi=500)
    plt.close()
    print(" ")


def get_file_name(data_dir, start_date, num_subcompartments, boolageresolved=False):
    """ Gives a paths to a file with the simulation results for an LCT model with num_subcompartments subcompartments, 
    for an covid19 inspired scenario starting at start_date.
    This function uses standard defined naming convention of the lct covid19 inspired simulations.

    @param[in] data_dir: Data directory pointing to the folder where the simulation result file lies in. 
    @param[in] start_date: Start date of the simulation result.
    @param[in] num_subcompartments: Number of subcompartments for the LCT model used to obtain the simulation results.
    @param[in] boolageresolved: Specifies whether the result should be resolved by age (or accumulated results). 
    """
    filename = "lct_" + start_date+"_subcomp" + f"{ num_subcompartments}"+"_"
    if boolageresolved:
        filename += "ageresolved"
    else:
        filename += "accumulated"
    return os.path.join(data_dir, filename)


def main():
    if not os.path.isdir(plotfolder):
        os.makedirs(plotfolder)

    simulation_data_dir = os.path.join(os.path.dirname(
        __file__), "..", "..", "..", "data", "simulation_lct_covid19")
    path_to_rki_data = os.path.join(os.path.dirname(
        __file__), "..", "..", "..", "data", "pydata", "Germany", "cases_all_age_all_dates.json")
    path_to_divi_data = os.path.join(os.path.dirname(
        __file__), "..", "..", "..", "data", "pydata", "Germany", "germany_divi_all_dates.json")

    start_date = '2020-10-1'
    start_date_timestamp = pd.Timestamp(start_date)
    scale_confirmed_cases = 1.2
    tmax = 45

    # Define which subfigures of Figure 16 of the paper should be created.
    figures = ['topleft', 'topright', 'bottomleft', 'bottomright']
    if 'topleft' in figures:
        plot_daily_new_transmissions([get_file_name(simulation_data_dir, start_date, 1),
                                      get_file_name(
                                      simulation_data_dir, start_date, 3),
                                      get_file_name(
                                      simulation_data_dir, start_date, 10),
                                      get_file_name(
                                      simulation_data_dir, start_date, 50),
                                      get_file_name(simulation_data_dir, start_date, 0)],
                                     path_to_rki_data, start_date_timestamp, tmax, scale_confirmed_cases,
                                     legend_labels=list(
            ["Extrapolated RKI data", "ODE", "LCT3", "LCT10", "LCT50", "LCTvar"]),
            file_name="real_new_infections_"+start_date+"_allage")
    if 'topright' in figures:
        plot_InfectedSymptoms_or_Deaths([get_file_name(simulation_data_dir, start_date, 1),
                                         get_file_name(
            simulation_data_dir, start_date, 3),
            get_file_name(
            simulation_data_dir, start_date, 10),
            get_file_name(
            simulation_data_dir, start_date, 50),
            get_file_name(simulation_data_dir, start_date, 0)],
            path_to_rki_data, start_date_timestamp, tmax, scale_confirmed_cases, list(
            ["Extrapolated RKI data", "ODE", "LCT3", "LCT10", "LCT50", "LCTvar"]), plot_deaths=False, file_name="real_infected_"+start_date+"_allage")
    if 'bottomleft' in figures:
        plot_InfectedCritical([get_file_name(simulation_data_dir, start_date, 1),
                               get_file_name(simulation_data_dir,
                                             start_date, 3),
                               get_file_name(simulation_data_dir,
                                             start_date, 10),
                               get_file_name(simulation_data_dir,
                                             start_date, 50),
                               get_file_name(simulation_data_dir, start_date, 0)],
                              path_to_divi_data, start_date_timestamp, tmax,  list(
            ["Extrapolated RKI data", "ODE", "LCT3", "LCT10", "LCT50", "LCTvar"]),
            file_name="real_icu_"+start_date+"_allage")
    if 'bottomright' in figures:
        plot_InfectedSymptoms_or_Deaths([get_file_name(simulation_data_dir, start_date, 1),
                                         get_file_name(
                                             simulation_data_dir, start_date, 3),
                                         get_file_name(simulation_data_dir, start_date, 10), get_file_name(
                                             simulation_data_dir, start_date, 50),
                                         get_file_name(simulation_data_dir, start_date, 0)],
                                        path_to_rki_data, start_date_timestamp, tmax, scale_confirmed_cases,
                                        list(["Extrapolated RKI data", "ODE",
                                             "LCT3", "LCT10", "LCT50", "LCTvar"]),
                                        plot_deaths=True, file_name="real_deaths_"+start_date+"_allage")


if __name__ == "__main__":
    main()
