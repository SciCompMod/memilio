#############################################################################
# Copyright (C) 2020-2024 MEmilio
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
"""@plot_results_lct_secir_real.py
Functions to plot and compare results of simulations with different kind of models,
eg LCT, IDE or ODE SECIR models without division in agegroups with real data.

The data to be plotted should be stored in a '../data/simulation_lct/real' folder as .h5 files.
Data could be generated eg by executing the file ./cpp/examples/lct_secir_real_senarios.cpp.
"""

import h5py
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

import memilio.epidata.getDataIntoPandasDataFrame as gd

# Define compartments.
secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
              5: 'ICU', 6: 'Recovered', 7: 'Dead'}
Age_RKI_names = {0: 'A00-A04', 1: 'A05-A14',
                 2: 'A15-A34',  3: 'A35-A59', 4: 'A60-A79', 5: 'A80+'}

# Define parameters used for simulation, used for plotting real data.
TimeExposed = [3.335, 3.335, 3.335, 3.335, 3.335, 3.335]
TimeInfectedNoSymptoms = [2.74, 2.74, 2.565, 2.565, 2.565, 2.565]
TimeInfectedSymptoms = [7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775]
TimeInfectedSevere = [5., 5., 5.925, 7.55, 8.5, 11.]
TimeInfectedCritical = [6.95, 6.95, 6.86, 17.36, 17.1, 11.6]
RecoveredPerInfectedNoSymptoms = [
    1 - 0.75, 1 - 0.75, 1 - 0.8, 1 - 0.8, 1 - 0.8, 1 - 0.8]
age_group_sizes = [3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434]

# Define color and style to be used while plotting for different models to make plots consistent.
color_dict = {'ODE': '#1f77b4',
              'LCT3': '#2ca02c',
              'LCT10': '#ff7f0e',
              'LCT20': '#9467bd',
              'LCT50':  '#e377c2',
              }
linestyle_dict = {'ODE': 'solid',
                  'LCT3': 'solid',
                  'LCT10': 'solid',
                  'LCT20': 'solid',
                  'LCT50': 'solid'
                  }


def load_data(file, start_date, tmax, scaleConfirmed):
    """ Loads RKI data and computes 'InfectedSymptoms', 'Deaths' and 'DailyNewTransmissions' using scales, dates etc 
    from the dictionary parameters. Method matches the method for computing initial values for the LCT model. 
    See also cpp/models/lct_secir/parameters_io.h.
    @param[in] file Path to the RKI data file for whole Germany. 
        Can be downloaded eg via pycode/memilio-epidata/memilio/epidata/getCaseData.py.
    """

    # Read data.
    df = pd.read_json(file)
    df = df.drop(columns=['Recovered'])

    end_date = start_date+pd.DateOffset(days=tmax)

    # Remove unnecessary dates.
    df = df[(df['Date'] >= start_date +
             pd.DateOffset(days=-math.ceil(max(TimeInfectedSymptoms) + max(TimeInfectedSevere) + max(TimeInfectedCritical))))
            & (df['Date'] <= end_date + pd.DateOffset(days=math.ceil(max(TimeExposed) + max(TimeInfectedNoSymptoms))))]
    # Scale confirmed cases because of undetected infections.
    df['Confirmed'] = scaleConfirmed * df['Confirmed']
    # As unknown age groups are omitted in the cpp epi_data used for initialization, we omit it here too.
    df = df[(df['Age_RKI'] != 'unknown')]

    # df2 stores the result of the computation.
    df2 = df.copy()
    df2 = df2[(df['Date'] >= start_date) & (df['Date'] <= end_date)]
    df2 = df2.reset_index()
    df2 = df2.drop(columns=['index', 'Confirmed', 'Deaths'])
    df2['InfectedSymptoms'] = pd.Series(dtype='double')
    df2['Deaths'] = pd.Series(dtype='double')
    df2['DailyNewTransmissions'] = pd.Series(dtype='double')

    for age_group in range(len(Age_RKI_names)):
        # Calculate individuals in compartment InfectedSymptoms using TimeInfectedSymptoms..
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

        # Calculate new infections per day.
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


def compare_compartments_real(files, datafile, legendplot, deaths=False, filename_plot="compare_real"):
    """ Plots simulation results compared with real data for the compartments Deaths and InfectedSymptoms.
        The simulation results should consist of accumulated numbers for subcompartments in case of an LCT model.

    @param[in] files: Paths of the hdf5-files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] datafile: Path to the RKI data file for whole Germany. Can be downloaded eg via pycode/memilio-epidata/memilio/epidata/getCaseData.py.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] deaths: If False, InfectedSymptoms compartment is plotted, if True, deaths are plotted.
    @param[in] filename_plot: Name to use as the file name for the saved plot.
    """
    # Define plot.
    plt.figure(filename_plot)

    data_rki_ma = load_data(datafile)
    num_days = data_rki_ma.shape[0]

    if (deaths):
        plt.plot(range(num_days), data_rki_ma['Deaths'],
                 linestyle='dashed',  color='grey', linewidth=1.2)
        compartment_idx = 7
        labely = "Deaths"
    else:
        plt.plot(range(num_days), data_rki_ma['InfectedSymptoms'],
                 linestyle='dashed', color='grey', linewidth=1.2)
        compartment_idx = 3
        labely = "Number of people in I"

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
        if (total.shape[1] != 8):
            raise gd.DataError(
                "Expected a different number of compartments.")

        # Plot result.
        plt.plot(dates, total[:, compartment_idx],
                 linewidth=1.2, linestyle=linestyle_dict[legendplot[1+file]], color=color_dict[legendplot[1+file]])
        h5file.close()

    plt.xlabel('Date', fontsize=16)
    plt.ylabel(labely, fontsize=16)
    plt.yticks(fontsize=9)
    plt.xlim(left=0, right=num_days-1)
    # Define x-ticks.
    datelist = np.array(pd.date_range(start_date.date(),
                                      periods=num_days, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int((num_days-1) / 5) + 1) * 5)
    plt.xticks(tick_range, datelist[tick_range],
               rotation=45, fontsize=12)
    plt.xticks(np.arange(num_days), minor=True)

    plt.legend(legendplot, fontsize=12, framealpha=0.5)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # Save result.
    if not os.path.isdir('Plots_real'):
        os.makedirs('Plots_real')
    plt.savefig('Plots_real/'+filename_plot+'.png',
                bbox_inches='tight', dpi=500)


def plot_new_infections_real(files, age_group, datafile, start_date, tmax, scaleConfirmed, legendplot, filename_plot="compare_new_infections_real"):
    """ Plots simulation results compared with real data for new infections within one day.
        The simulation results should consist of accumulated numbers for subcompartments in case of an LCT model.

    @param[in] files: paths of the hdf5-files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] datafile: Path to the RKI data file for whole Germany. Can be downloaded eg via pycode/memilio-epidata/memilio/epidata/getCaseData.py.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """
    plt.figure(filename_plot)

    data_rki = load_data(datafile, start_date, tmax, scaleConfirmed)
    num_days = tmax+1

    if (age_group == -1):
        data_rki = data_rki.drop(columns=['Age_RKI'])
        data_rki = data_rki.groupby(['Date']).sum()
        plt.plot(range(num_days), data_rki['DailyNewTransmissions'],
                 linestyle='None', color='grey', marker='x', markersize=5)
        plt.title('All Age Groups')
        print("Daily new transmissions at the first day of RKI data is: " +
              f"{data_rki.loc[start_date]['DailyNewTransmissions']}")
    else:
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
        print("Daily new transmissions at the first day of the simulation " + legendplot[file+1]+" is: "
              f"{incidence[0]}")

        # Plot result.
        if legendplot[file] in linestyle_dict:
            plt.plot(dates[1:], incidence, linewidth=1.2,
                     linestyle=linestyle_dict[legendplot[1+file]], color=color_dict[legendplot[1+file]])
        else:
            plt.plot(dates[1:], incidence, linewidth=1.2)
        h5file.close()

    plt.ylabel('Daily new Transmissions', fontsize=14)
    plt.ylim(bottom=0)
    plt.xlabel('Date', fontsize=14)
    plt.xlim(left=0, right=num_days-1)
    # Define x-ticks.
    datelist = np.array(pd.date_range(start_date.date(),
                                      periods=num_days, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int((num_days - 1) / 5) + 1) * 5)
    plt.xticks(tick_range, datelist[tick_range],
               rotation=45, fontsize=12)
    plt.xticks(np.arange(num_days), minor=True)

    plt.legend(legendplot, fontsize=12, framealpha=0.5)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # Save result.
    if not os.path.isdir('Plots_real'):
        os.makedirs('Plots_real')
    plt.savefig('Plots_real/'+filename_plot+'.png',
                bbox_inches='tight', dpi=500)
    plt.close()
    print(" ")


def main():

    folder = "../data/simulation_lct_real/"
    start_date = pd.Timestamp('2020.10.01')

    load_data("../data/pydata/Germany/cases_all_age.json", start_date, 15, 1.)


if __name__ == "__main__":
    main()
