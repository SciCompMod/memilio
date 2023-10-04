#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
functions to plot results of a simulation with a LCT SECIR model with subcompartments.
There is also a method to compare different results of different models.

The data to be plotted should be stored in a '../data/simulation_lct' folder as .h5 files.
Data could be generated eg by executing the file ./cpp/examples/lct_secir_compare.cpp.
"""

import h5py
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

import memilio.epidata.getDataIntoPandasDataFrame as gd

# Define compartments
secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
              5: 'ICU', 6: 'Recovered', 7: 'Death'}

# Define parameters used for simulation, used for plotting real data.
parameters = {
    'TimeExposed':  3.335,
    'TimeInfectedNoSymptoms':  3.31331,
    'TimeInfectedSymptoms': 6.94547,
    'TimeInfectedSevere':  11.634346,
    'TimeInfectedCritical': 17.476959,
    'RecoveredPerInfectedNoSymptoms':  0.206901,
    'start_date': pd.Timestamp('2020.10.15'),
    'end_date': pd.Timestamp('2020.11.15'),
    'scaleConfirmed': 2.
}
color_dict = {'ODE': '#1f77b4',
              'RKI': 'grey',
              'LCT': '#ff7f0e',
              'LCT3': '#d62728',
              'LCT10': '#ff7f0e',
              'LCT20': '#9467bd',
              'LCTvar': '#2ca02c',
              'IDE10': '#17becf',
              'IDEvar': '#e377c2',
              }
linestyle_dict = {'ODE': 'solid',
                  'RKI': 'dashed',
                  'LCT': 'solid',
                  'LCT3': 'dashdot',
                  'LCT10': 'solid',
                  'LCT20': 'dashdot',
                  'LCTvar': 'solid',
                  'IDE10': 'dashed',
                  'IDEvar': 'dashed'}


def load_data(file):
    df = pd.read_json(file)
    df = df.drop(columns=['Recovered'])

    df = df[(df['Date'] >= parameters['start_date']+pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'])))
            & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=math.ceil(parameters['TimeExposed']+parameters['TimeInfectedNoSymptoms'])))]
    df['Confirmed'] = parameters['scaleConfirmed']*df['Confirmed']
    df2 = df.copy()
    df2 = df2[(df['Date'] >= parameters['start_date'])
              & (df['Date'] <= parameters['end_date'])]
    df2 = df2.reset_index()
    df2 = df2.drop(columns=['index', 'Confirmed', 'Deaths'])

    help_I = df['Confirmed'][(df['Date'] >= parameters['start_date'])
                             & (df['Date'] <= parameters['end_date'])].to_numpy()
    help_I = help_I - (1-math.fmod(parameters['TimeInfectedSymptoms'], 1))*df['Confirmed'][(df['Date'] >= parameters['start_date']+pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptoms'])))
                                                                                           & (df['Date'] <= parameters['end_date'] + pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptoms'])))].to_numpy()
    help_I = help_I - math.fmod(parameters['TimeInfectedSymptoms'], 1) * df['Confirmed'][(df['Date'] >= parameters['start_date']+pd.DateOffset(days=-math.ceil(
        parameters['TimeInfectedSymptoms']))) & (df['Date'] <= parameters['end_date']+pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptoms'])))].to_numpy()
    df2['InfectedSymptoms'] = help_I
    help_D = (1-(1-math.fmod(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'], 1)))*df['Deaths'][(df['Date'] >= parameters['start_date']+pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'])))
                                                                                                                                                       & (df['Date'] <= parameters['end_date']+pd.DateOffset(days=-math.ceil(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'])))].to_numpy()
    help_D = help_D + (1-math.fmod(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'], 1))*df['Deaths'][(df['Date'] >= parameters['start_date']+pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'])))
                                                                                                                                                            & (df['Date'] <= parameters['end_date']+pd.DateOffset(days=-math.floor(parameters['TimeInfectedSymptoms']+parameters['TimeInfectedSevere']+parameters['TimeInfectedCritical'])))].to_numpy()
    df2['Deaths'] = help_D
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
    """ Plots the result of a simulation with an LCT SECIR model compared with real data for the compartments deaths and infected.
        The simulation result should consist of accumulated numbers for subcompartments.

    @param[in] files: paths of the hdf5-files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] datafile: Path to the RKI data file for whole Germany. Can be downloaded eg via pycode/memilio-epidata/memilio/epidata/getCaseData.py.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """
    # define plot
    plt.figure(filename_plot)

    data_rki = load_data(datafile)
    num_days = data_rki.shape[0]

    if (deaths):
        plt.plot(range(num_days), data_rki['Deaths'],
                 linestyle=linestyle_dict['RKI'], color=color_dict['RKI'], linewidth=1.2)
        compartment_idx = 7
        labely = "Todesfälle"
    else:
        plt.plot(range(num_days), data_rki['InfectedSymptoms'],
                 linestyle=linestyle_dict['RKI'], color=color_dict['RKI'], linewidth=1.2)
        compartment_idx = 3
        labely = "Anzahl der Individuen in I"

    # add simulation results to plot
    for file in range(len(files)):
        # load data
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result
        total = data['Total'][:, :]
        if (total.shape[1] != 8):
            raise gd.DataError(
                "Expected a different number of compartments.")

        # plot result
        plt.plot(dates, total[:, compartment_idx],
                 linewidth=1.2, linestyle=linestyle_dict[legendplot[1+file]], color=color_dict[legendplot[1+file]])
        h5file.close()

    plt.xlabel('Datum', fontsize=14)
    plt.ylabel(labely, fontsize=14)
    plt.xlim(left=0, right=num_days-1)
    # define x-ticks
    datelist = np.array(pd.date_range(parameters["start_date"].date(
    ), periods=num_days, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int((num_days-1) / 3) + 1) * 3)
    plt.xticks(tick_range, datelist[tick_range],
               rotation=45, fontsize=14)
    plt.xticks(np.arange(num_days), minor=True)

    plt.legend(legendplot, fontsize=10)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # save result
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


def plot_new_infections_real(files, datafile, legendplot, filename_plot="compare_new_infections_real"):
    """ Plots the result of a simulation with an LCT SECIR model compared with real data for the compartments deaths and infected.
        The simulation result should consist of accumulated numbers for subcompartments.

    @param[in] files: paths of the hdf5-files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] datafile: Path to the RKI data file for whole Germany. Can be downloaded eg via pycode/memilio-epidata/memilio/epidata/getCaseData.py.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """
    # define plot
    plt.figure(filename_plot)

    data_rki = load_data(datafile)
    num_days = data_rki.shape[0]

    plt.plot(range(num_days), data_rki['NewInfectionsDay'],
             linestyle=linestyle_dict['RKI'], color=color_dict['RKI'], linewidth=1.2)

    # add simulation results to plot
    for file in range(len(files)):
        # load data
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result
        total = data['Total'][:, :]
        if (total.shape[1] != 8):
            raise gd.DataError(
                "Expected a different number of compartments.")
        incidence = (total[:-1, 0]-total[1:, 0])/(dates[1:]-dates[:-1])

        # plot result
        plt.plot(dates[1:], incidence,
                 linewidth=1.2, linestyle=linestyle_dict[legendplot[1+file]], color=color_dict[legendplot[1+file]])
        h5file.close()

    plt.xlabel('Datum', fontsize=14)
    plt.ylabel('Neuansteckungen pro Tag', fontsize=14)
    plt.xlim(left=0, right=num_days-1)
    # define x-ticks
    datelist = np.array(pd.date_range(parameters["start_date"].date(
    ), periods=num_days, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int((num_days - 1) / 3) + 1) * 3)
    plt.xticks(tick_range, datelist[tick_range],
               rotation=45, fontsize=14)
    plt.xticks(np.arange(num_days), minor=True)

    plt.legend(legendplot, fontsize=10)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # save result
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


def plot_lct_subcompartments(file, vec_subcompartments=[1, 20, 20, 20, 20, 20, 1, 1], filename_plot="LCT_subcompartments"):
    """ Plots the result of a simulation with an LCT SECIR model in an 4x2 Plot.
        In each subplot, one compartment with one line per subcompartment is plotted.

    @param[in] file: path of the file (without file extension .h5) with the simulation result of an lct model with subcompartments.
        The number of subcompartments should be fitting to the ones specified in get_subcompartments().
    @param[in] vec_subcompartments: vector with the number of subcompartments per compartment used for simulation.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """
    secir_dict_short = {0: 'S', 1: 'E', 2: 'C', 3: 'I', 4: 'H',
                        5: 'U', 6: 'R', 7: 'D'}
    # load data
    h5file = h5py.File(str(file) + '.h5', 'r')

    if (len(list(h5file.keys())) > 1):
        raise gd.DataError("File should contain one dataset.")
    if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
        raise gd.DataError("Expected only one group.")

    data = h5file[list(h5file.keys())[0]]
    dates = data['Time'][:]
    # As there should be only one Group, total is the simulation result
    total = data['Total'][:, :]
    if (total.shape[1] != sum(vec_subcompartments)):
        raise gd.DataError("Expected a different number of subcompartments.")
    h5file.close()

    # Plot result.
    plt.style.use('seaborn-v0_8-colorblind')
    fig, axs = plt.subplots(
        4, 2, sharex='all', num=filename_plot)
    for i in range(len(vec_subcompartments)):
        idx_start = sum(vec_subcompartments[0:i])
        axs[int(i/2), i % 2].plot(dates,
                                  total[:, idx_start:idx_start+vec_subcompartments[i]], linewidth=1.2)
        axs[int(i/2), i % 2].grid(True, linestyle='--')
        legendplot = []
        if (vec_subcompartments[i] == 1):
            legendplot.append(secir_dict_short[i])
        elif (vec_subcompartments[i] < 5):
            for j in range(vec_subcompartments[i]):
                legendplot.append(secir_dict_short[i]+str(j+1))

        axs[int(i/2), i % 2].legend(legendplot, fontsize=10)
        axs[int(i/2), i % 2].set_xlim(left=0, right=dates[-1])
    fig.supxlabel('Time')
    fig.supylabel('Number of persons')
    plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.1)
    plt.subplots_adjust(bottom=0.09)

    # save plot
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    fig.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


def compare_all_compartments(files, legendplot, filename_plot="compare_compartments"):
    """ Creates a 4x2 Plot with one subplot per compartment and one line per result one wants to compare.

    @param[in] files: paths of the files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """

    fig, axs = plt.subplots(
        4, 2, sharex='all', num=filename_plot, tight_layout=True)

    # add results to plot
    for file in range(len(files)):
        # load data
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result
        total = data['Total'][:, :]
        if (total.shape[1] != 8):
            raise gd.DataError(
                "Expected a different number of compartments.")
        # plot data
        if legendplot[file] in linestyle_dict:
            for i in range(8):
                axs[int(i/2), i % 2].plot(dates,
                                          total[:, i], label=legendplot[file], linewidth=1.2, linestyle=linestyle_dict[legendplot[file]], color=color_dict[legendplot[file]])
        else:
            for i in range(8):
                axs[int(i/2), i % 2].plot(dates,
                                          total[:, i], label=legendplot[file], linewidth=1.2)
        h5file.close()

    # define some characteristics of the plot
    for i in range(8):
        axs[int(i/2), i % 2].set_title(secir_dict[i], fontsize=10)
        axs[int(i/2), i % 2].set_xlim(left=0, right=dates[-1])
        axs[int(i/2), i % 2].grid(True, linestyle='--')
        axs[int(i/2), i % 2].tick_params(axis='y', labelsize=8)
        axs[int(i/2), i % 2].tick_params(axis='x', labelsize=8)

    fig.supxlabel('Zeit', fontsize=12)
    fig.supylabel('Anzahl an Personen', fontsize=12)

    lines, labels = axs[0, 0].get_legend_handles_labels()
    lgd = fig.legend(lines, labels, ncol=len(legendplot),  loc='outside lower center',
                     fontsize=10, bbox_to_anchor=(0.5, - 0.06), bbox_transform=fig.transFigure)

    plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.1)
    plt.subplots_adjust(bottom=0.09)

    # save result
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    fig.savefig('Plots/'+filename_plot+'.png',
                bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=500)


def plot_new_infections(files, legendplot, filename_plot="compare_new_infections"):
    """ Single plot to compare the incidence of different results. Incidence means the number of people leaving the susceptible class per day.

    @param[in] files: paths of the files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
     @param[in] filename_plot: name to use as the file name for the saved plot.
    """
    # define plot
    plt.figure(filename_plot)

    # add results to plot
    for file in range(len(files)):
        # load data
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result
        total = data['Total'][:, :]
        if (total.shape[1] != 8):
            raise gd.DataError(
                "Expected a different number of compartments.")
        incidence = (total[:-1, 0]-total[1:, 0])/(dates[1:]-dates[:-1])

        # plot result
        if legendplot[file] in linestyle_dict:
            plt.plot(dates[1:], incidence, linewidth=1.2,
                     linestyle=linestyle_dict[legendplot[file]], color=color_dict[legendplot[file]])
        else:
            plt.plot(dates[1:], incidence, linewidth=1.2)

        h5file.close()

    plt.xlabel('Zeit', fontsize=14)
    plt.ylabel('Neuansteckungen pro Tag', fontsize=14)
    plt.ylim(bottom=0)
    plt.xlim(left=0, right=dates[-1])
    plt.legend(legendplot, fontsize=10)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # save result
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


if __name__ == '__main__':
    # simulation results should be stored in folder "../data/simulation_lct"
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "data", "simulation_lct")

    # Defines which simulation case should be plotted
    case = 5

    if case == 1:
        plot_new_infections([os.path.join(data_dir, "init", "lct_init_transitions20"), os.path.join(data_dir, "init", "lct_init_mean20"), os.path.join(data_dir, "init", "lct_init_first20")],
                            legendplot=list(["Übergänge", "Mittelwerte", "Erstes Subkompartiment"]), filename_plot="compare_new_infections_initialization20")
        compare_all_compartments([os.path.join(data_dir, "init", "lct_init_transitions20"), os.path.join(data_dir, "init", "lct_init_mean20"), os.path.join(data_dir, "init", "lct_init_first20")],
                                 legendplot=list(["Übergänge", "Mittelwerte", "Erstes Subkompartiment"]), filename_plot="compare_compartments_initialization20")
        plot_new_infections([os.path.join(data_dir, "init", "lct_init_transitions_long20"), os.path.join(data_dir, "init", "lct_init_mean_long20"), os.path.join(data_dir, "init", "lct_init_first_long20")],
                            legendplot=list(["Übergänge", "Mittelwerte", "Erstes Subkompartiment"]), filename_plot="compare_new_infections_initialization_long20")
        compare_all_compartments([os.path.join(data_dir, "init", "lct_init_transitions_long20"), os.path.join(data_dir, "init", "lct_init_mean_long20"), os.path.join(data_dir, "init", "lct_init_first_long20")],
                                 legendplot=list(["Übergänge", "Mittelwerte", "Erstes Subkompartiment"]), filename_plot="compare_compartments_initialization_long20")
        plot_new_infections([os.path.join(data_dir, "init", "lct_init_transitions3"), os.path.join(data_dir, "init", "lct_init_transitions20"), os.path.join(data_dir, "init", "lct_init_first3"), os.path.join(data_dir, "init", "lct_init_first10"), os.path.join(data_dir, "init", "lct_init_first20")],
                            legendplot=list(["Übergänge3", "Übergänge20", "Erstes Subkompartiment3", "Erstes Subkompartiment10", "Erstes Subkompartiment20"]), filename_plot="compare_new_infections_initialization_diffnums")
        plot_lct_subcompartments(os.path.join(data_dir, "init",  "lct_init_transitions_subcompartments3"), vec_subcompartments=[
            1, 3, 3, 3, 3, 3, 1, 1], filename_plot="lct_init_transitions_subcompartments3")

    elif case == 2:
        plot_new_infections([os.path.join(data_dir, "dropR0", "fictional_drop_lct_1"), os.path.join(data_dir, "dropR0", "fictional_drop_lct_3"),  os.path.join(data_dir, "dropR0", "fictional_drop_lct_10"), os.path.join(data_dir, "dropR0", "fictional_drop_lct_20"),  os.path.join(data_dir, "dropR0", "fictional_drop_lct_0")],
                            legendplot=list(["ODE", "LCT3", "LCT10", "LCT20", "LCTvar"]),  filename_plot="compare_new_infections_lct_drop")
        compare_all_compartments([os.path.join(data_dir, "dropR0", "fictional_drop_lct_1"), os.path.join(data_dir, "dropR0", "fictional_drop_lct_3"),  os.path.join(data_dir, "dropR0", "fictional_drop_lct_10"), os.path.join(data_dir, "dropR0", "fictional_drop_lct_20"),  os.path.join(data_dir, "dropR0", "fictional_drop_lct_0")],
                                 legendplot=list(["ODE", "LCT3", "LCT10", "LCT20", "LCTvar"]),  filename_plot="compare_compartments_lct_drop")
        plot_new_infections([os.path.join(data_dir, "dropR0", "fictional_drop_lct_1"), os.path.join(data_dir, "dropR0", "fictional_drop_ide_10"),  os.path.join(data_dir, "dropR0", "fictional_drop_ide_0"), os.path.join(data_dir, "dropR0", "fictional_drop_lct_10"),  os.path.join(data_dir, "dropR0", "fictional_drop_lct_0")],
                            legendplot=list(["ODE", "IDE10", "IDEvar", "LCT10", "LCTvar"]),  filename_plot="compare_new_infections_ide_drop")
        compare_all_compartments([os.path.join(data_dir, "dropR0", "fictional_drop_lct_1"), os.path.join(data_dir, "dropR0", "fictional_drop_ide_10"),  os.path.join(data_dir, "dropR0", "fictional_drop_ide_0"), os.path.join(data_dir, "dropR0", "fictional_drop_lct_10"),  os.path.join(data_dir, "dropR0", "fictional_drop_lct_0")],
                                 legendplot=list(["ODE", "IDE10", "IDEvar", "LCT10", "LCTvar"]),  filename_plot="compare_compartments_ide_drop")

    if case == 3:
        plot_new_infections([os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_1"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_3"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_10"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_20"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_0")],
                            legendplot=list(["ODE", "LCT3", "LCT10", "LCT20", "LCTvar"]),  filename_plot="compare_new_infections_lct_rise4")
        compare_all_compartments([os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_1"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_3"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_10"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_20"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_0")],
                                 legendplot=list(["ODE", "LCT3", "LCT10", "LCT20", "LCTvar"]),  filename_plot="compare_compartments_lct_rise4")
        plot_new_infections([os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_1"), os.path.join(data_dir, "riseR0short", "fictional_rise_ide_4_10"),  os.path.join(data_dir, "riseR0short", "fictional_rise_ide_4_0"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_10"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_0")],
                            legendplot=list(["ODE", "IDE10", "IDEvar", "LCT10", "LCTvar"]),  filename_plot="compare_new_infections_ide_rise4")
        compare_all_compartments([os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_1"), os.path.join(data_dir, "riseR0short", "fictional_rise_ide_4_10"),  os.path.join(data_dir, "riseR0short", "fictional_rise_ide_4_0"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_10"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_4_0")],
                                 legendplot=list(["ODE", "IDE10", "IDEvar", "LCT10", "LCTvar"]),  filename_plot="compare_compartments_ide_rise4")
        plot_new_infections([os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_1"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_3"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_10"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_20"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_0")],
                            legendplot=list(["ODE", "LCT3", "LCT10", "LCT20", "LCTvar"]),  filename_plot="compare_new_infections_lct_rise2")
        compare_all_compartments([os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_1"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_3"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_10"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_20"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_0")],
                                 legendplot=list(["ODE", "LCT3", "LCT10", "LCT20", "LCTvar"]),  filename_plot="compare_compartments_lct_rise2")
        plot_new_infections([os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_1"), os.path.join(data_dir, "riseR0short", "fictional_rise_ide_2_10"),  os.path.join(data_dir, "riseR0short", "fictional_rise_ide_2_0"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_10"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_0")],
                            legendplot=list(["ODE", "IDE10", "IDEvar", "LCT10", "LCTvar"]),  filename_plot="compare_new_infections_ide_rise2")
        compare_all_compartments([os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_1"), os.path.join(data_dir, "riseR0short", "fictional_rise_ide_2_10"),  os.path.join(data_dir, "riseR0short", "fictional_rise_ide_2_0"), os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_10"),  os.path.join(data_dir, "riseR0short", "fictional_rise_lct_2_0")],
                                 legendplot=list(["ODE", "IDE10", "IDEvar", "LCT10", "LCTvar"]),  filename_plot="compare_compartments_ide_rise2")

    if case == 4:
        plot_new_infections([os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long1"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long3"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long10"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long20"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long0")],
                            legendplot=list(["ODE", "LCT3", "LCT10", "LCT20", "LCTvar"]),  filename_plot="compare_new_infections_lct_rise4long")
        compare_all_compartments([os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long1"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long3"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long10"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long20"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long0")],
                                 legendplot=list(["ODE", "LCT3", "LCT10", "LCT20", "LCTvar"]),  filename_plot="compare_compartments_lct_rise4long")
        plot_new_infections([os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long1"), os.path.join(data_dir, "riseR0long", "fictional_rise_ide_4long10"),  os.path.join(data_dir, "riseR0long", "fictional_rise_ide_4long0"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long10"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long0")],
                            legendplot=list(["ODE", "IDE10", "IDEvar", "LCT10", "LCTvar"]),  filename_plot="compare_new_infections_ide_rise4long")
        compare_all_compartments([os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long1"), os.path.join(data_dir, "riseR0long", "fictional_rise_ide_4long10"),  os.path.join(data_dir, "riseR0long", "fictional_rise_ide_4long0"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long10"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_4long0")],
                                 legendplot=list(["ODE", "IDE10", "IDEvar", "LCT10", "LCTvar"]),  filename_plot="compare_compartments_ide_rise4long")
        plot_new_infections([os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long1"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long3"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long10"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long20"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long0")],
                            legendplot=list(["ODE", "LCT3", "LCT10", "LCT20", "LCTvar"]),  filename_plot="compare_new_infections_lct_rise2long")
        compare_all_compartments([os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long1"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long3"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long10"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long20"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long0")],
                                 legendplot=list(["ODE", "LCT3", "LCT10", "LCT20", "LCTvar"]),  filename_plot="compare_compartments_lct_rise2long")
        plot_new_infections([os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long1"), os.path.join(data_dir, "riseR0long", "fictional_rise_ide_2long10"),  os.path.join(data_dir, "riseR0long", "fictional_rise_ide_2long0"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long10"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long0")],
                            legendplot=list(["ODE", "IDE10", "IDEvar", "LCT10", "LCTvar"]),  filename_plot="compare_new_infections_ide_rise2long")
        compare_all_compartments([os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long1"), os.path.join(data_dir, "riseR0long", "fictional_rise_ide_2long10"),  os.path.join(data_dir, "riseR0long", "fictional_rise_ide_2long0"), os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long10"),  os.path.join(data_dir, "riseR0long", "fictional_rise_lct_2long0")],
                                 legendplot=list(["ODE", "IDE10", "IDEvar", "LCT10", "LCTvar"]),  filename_plot="compare_compartments_ide_rise2long")

    if case == 5:
        compare_compartments_real([os.path.join(data_dir, "real", "real_ode_2020_10_15"), os.path.join(data_dir, "real", "real_lct_2020_10_15")],
                                  os.path.join(os.path.dirname(
                                      __file__), "..", "data", "pydata", "Germany", "cases_all_germany_ma7.json"),
                                  legendplot=list(["RKI", "ODE", "LCT"]), deaths=True, filename_plot="compare_deaths_real_2020_10_15")
        compare_compartments_real([os.path.join(data_dir, "real", "real_ode_2020_10_15"), os.path.join(data_dir, "real", "real_lct_2020_10_15")],
                                  os.path.join(os.path.dirname(
                                      __file__), "..", "data", "pydata", "Germany", "cases_all_germany_ma7.json"),
                                  legendplot=list(["RKI", "ODE", "LCT"]), deaths=False, filename_plot="compare_infected_real_2020_10_15")
        plot_new_infections_real([os.path.join(data_dir, "real", "real_ode_2020_10_15"), os.path.join(data_dir, "real", "real_lct_2020_10_15")],
                                 os.path.join(os.path.dirname(
                                     __file__), "..", "data", "pydata", "Germany", "cases_all_germany_ma7.json"),
                                 legendplot=list(["RKI", "ODE", "LCT"]), filename_plot="compare_new_infections_real_2020_10_15")
