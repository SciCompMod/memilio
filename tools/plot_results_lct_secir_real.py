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
"""@plot_results_lct_secir.py
Functions to plot and compare results of simulations with different kind of models,
eg LCT, IDE or ODE SECIR models without division in agegroups.

The data to be plotted should be stored in a '../data/simulation_lct' folder as .h5 files.
Data could be generated eg by executing the file ./cpp/examples/lct_secir_initializations.cpp.
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

# Define parameters used for simulation, used for plotting real data.
parameters = {
    'TimeExposed':  3.335,
    'TimeInfectedNoSymptoms':  3.31331,
    'TimeInfectedSymptoms': 6.94547,
    'TimeInfectedSevere': 11.634346,
    'TimeInfectedCritical': 17.476959,
    'RecoveredPerInfectedNoSymptoms':  0.206901,
    'start_date': pd.Timestamp('2020.06.01'),
    'end_date': pd.Timestamp('2020.06.01')+pd.DateOffset(days=45),
    'scaleConfirmed': 1.
}

# Define color and style to be used while plotting for different models to make plots consistent.
color_dict = {'ODE': '#1f77b4',
              'LCT3': '#2ca02c',
              'LCT10': '#ff7f0e',
              'LCT20': '#9467bd',
              'LCTvar':   '#d62728',
              'IDE3': '#e377c2',
              'IDE10': '#17becf',
              }
linestyle_dict = {'ODE': 'solid',
                  'LCT3': 'solid',
                  'LCT10': 'solid',
                  'LCT20': 'solid',
                  'IDE3': 'dashed',
                  'IDE10': 'dashed',
                  }

def load_data(file):
    """ Loads RKI data and computes 'InfectedSymptoms', 'Deaths' and 'NewInfectionsDay' using scales, dates etc from the dictionary parameters.
    Method matches the method for computing initial values for the LCT model. See also cpp/models/lct_secir/parameters_io.h.
    @param[in] file Path to the RKI data file for whole Germany. Can be downloaded eg via pycode/memilio-epidata/memilio/epidata/getCaseData.py.
    """
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
        labely = "deaths"
    else:
        plt.plot(range(num_days), data_rki_ma['InfectedSymptoms'],
                 linestyle='dashed', color='grey', linewidth=1.2)
        compartment_idx = 3
        labely = "number of people in I"

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

    plt.xlabel('date', fontsize=16)
    plt.ylabel(labely, fontsize=16)
    plt.yticks(fontsize=9)
    plt.xlim(left=0, right=num_days-1)
    # Define x-ticks.
    datelist = np.array(pd.date_range(parameters["start_date"].date(),
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
    plt.savefig('Plots_real/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


def plot_new_infections_real(files, datafile, legendplot, filename_plot="compare_new_infections_real"):
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

    data_rki = load_data(datafile)
    num_days = data_rki.shape[0]

    plt.plot(range(num_days), data_rki['NewInfectionsDay'],
             linestyle='None', color='grey', marker='x', markersize=10)

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
        incidence = (total[:-1, 0]-total[1:, 0])/(dates[1:]-dates[:-1])

        # Plot result.
        plt.plot(dates[1:], incidence,
                 linewidth=1.2, linestyle=linestyle_dict[legendplot[1+file]], color=color_dict[legendplot[1+file]])
        h5file.close()

    plt.xlabel('Date', fontsize=16)
    plt.ylabel('New infections per day', fontsize=16)
    plt.xlim(left=0, right=num_days-1)
    # Define x-ticks.
    datelist = np.array(pd.date_range(parameters["start_date"].date(),
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
    plt.savefig('Plots_real/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


if __name__ == '__main__':
    # simulation results should be stored in folder "../data/simulation_lct".
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "data", "simulation_lct")

    # Defines which simulation case should be plotted.
    case = 1

    if case == 1:
        # Attention: also change values in the dictionary parameters!
        # Compare simulation results for a real situation starting on 01.06.2020.
        compare_compartments_real([os.path.join(data_dir, "real", "real_lct1_2020_6_1"), 
                                   os.path.join(data_dir, "real", "real_lct3_2020_6_1"), 
                                   os.path.join(data_dir, "real", "real_lct10_2020_6_1"),
                                   os.path.join(data_dir, "real", "real_lct20_2020_6_1")],
                                  os.path.join(os.path.dirname(
                                      __file__), "..", "data", "pydata", "Germany", "cases_all_germany_ma7.json"),
                                  legendplot=list(["real data", "ODE", "LCT3", "LCT10", "LCT20"]), 
                                  deaths=False, filename_plot="compare_infected_real_2020_6_1")
        
        plot_new_infections_real([os.path.join(data_dir, "real", "real_lct1_2020_6_1"), 
                                  os.path.join(data_dir, "real", "real_lct3_2020_6_1"), 
                                  os.path.join(data_dir, "real", "real_lct10_2020_6_1"),
                                  os.path.join(data_dir, "real", "real_lct20_2020_6_1")],
                                 os.path.join(os.path.dirname(
                                     __file__), "..", "data", "pydata", "Germany", "cases_all_germany_all_dates.json"),
                                 legendplot=list(["real data", "ODE", "LCT3", "LCT10", "LCT20"]), 
                                 filename_plot="compare_new_infections_real_2020_6_1")
        
    elif case == 2:
        # Attention: also change values in the dictionary parameters!
        # Compare simulation results for a real situation starting on 01.10.2020.
        compare_compartments_real([os.path.join(data_dir, "real", "real_lct1_2020_10_1"), 
                                   os.path.join(data_dir, "real", "real_lct3_2020_10_1"), 
                                   os.path.join(data_dir, "real", "real_lct10_2020_10_1")],
                                  os.path.join(os.path.dirname(
                                      __file__), "..", "data", "pydata", "Germany", "cases_all_germany_ma7.json"),
                                  legendplot=list(["real data", "ODE", "LCT3", "LCT10"]), 
                                  deaths=False, filename_plot="compare_infected_real_2020_10_1")
        
        plot_new_infections_real([os.path.join(data_dir, "real", "real_lct1_2020_10_1"),
                                  os.path.join(data_dir, "real", "real_lct3_2020_10_1"), 
                                  os.path.join(data_dir, "real", "real_lct10_2020_10_1")],
                                 os.path.join(os.path.dirname(
                                     __file__), "..", "data", "pydata", "Germany", "cases_all_germany_all_dates.json"),
                                 legendplot=list(["real data", "ODE", "LCT3", "LCT10"]), 
                                 filename_plot="compare_new_infections_real_2020_10_1")
