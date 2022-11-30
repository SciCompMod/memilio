import pandas as pd
import numpy as np
import h5py
import os
import datetime
import copy
import matplotlib.pyplot as plt
from matplotlib import ticker
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['axes.facecolor'] = 'w'


# Define compartments
# Use this dict for new model outputs
secir_dict = {0: 'SusceptibleNaive', 1: 'SusceptiblePartialImmunity', 2: 'ExposedNaive', 3: 'ExposedPartialImmunity',
              4: 'ExposedImprovedImmunity', 5: 'InfectedNoSymptomsNaive', 6: 'InfectedNoSymptomsPartialImmunity',
              7: 'InfectedNoSymptomsImprovedImmunity ', 8: 'InfectedNoSymptomsNaiveConfirmed', 9: 'InfectedNoSymptomsPartialImmunityConfirmed',
              10: 'InfectedNoSymptomsImprovedImmunityConfirmed', 11: 'InfectedSymptomsNaive', 12: 'InfectedSymptomsPartialImmunity',
              13: 'InfectedSymptomsImprovedImmunity', 14: 'InfectedSymptomsNaiveConfirmed', 15: 'InfectedSymptomsPartialImmunityConfirmed',
              16: 'InfectedSymptomsImprovedImmunityConfirmed', 17: 'InfectedSevereNaive', 18: 'InfectedSeverePartialImmunity',
              19: 'InfectedSevereImprovedImmunity', 20: 'InfectedCriticalNaive', 21: 'InfectedCriticalPartialImmunity',
              22: 'InfectedCriticalImprovedImmunity', 23: 'SusceptibleImprovedImmunity', 24: 'DeadNaive', 25: 'DeadPartialImmunity',
              26: 'DeadImprovedImmunity', 27: 'TotalInfections', 28: "TemporaryImmunity1", 29: "TemporaryImmunity2"}

# Define aggregated infection states to plot and create 'concat_comps' dictionary for aggregation
secir_dict_aggregated = {
    0: 'Susceptible', 1: 'Partial Immunity', 2: 'Exposed', 3:
    'InfectedNoSymptoms', 4: 'InfectedSymptoms', 5: 'InfectedSevere', 6:
    'InfectedCritical', 7: 'Improved Immunity', 8: 'Dead', 9: 'Infected Total',
    10: 'TemporaryImmunity'}
concat_comps = {0: [0], 1: [1], 7: [23], 8: [24, 25, 26],
                9: [27], 10: [28, 29]}  # new, three dead compartments
for i in range(2, 7):
    concat_comps[i] = []
    for k, infstate in secir_dict.items():
        if secir_dict_aggregated[i] in infstate:
            if (infstate != 'TotalInfections') and (infstate != 'Infected Total'):
                concat_comps[i].append(k)


# Define age groups
age_groups = ['0-4 Years', '5-14 Years', '15-34 Years',
              '35-59 Years', '60-79 Years', '80+ Years']

# Define population data for incidence values and relative plots
base = 100000
age_group_sizes = np.array(
    [3961376, 7429883, 19117865, 28919134, 18057318, 5681135])

relative_dict = {}
for i in range(len(age_group_sizes)):
    relative_dict['Group' + str(i+1)] = age_group_sizes[i]/base

relative_dict['Total'] = np.sum(age_group_sizes)/base

# Define start day and simulation period
year, month, day = '2022', '6', '1'
start_date = pd.Timestamp(year + '.' + month.zfill(2) + '.' + day.zfill(2))
tmax = '10'
daysPlot = 45

# Define different scenario folders that will be read and plotted
date_str = '_' + str(year) + '_' + str(month) + '_' + str(day) + '_' + str(tmax)
path_sim = 'data/'
path_rki = 'data/'
scenario_list = [
    '']

# Provide a list of labels for corresponding plots


# Opens files from folder
# @param path_sim Path where simulation files have been written
# @param path_rki Path where extrapolated real data have been written
# @param spec_str_sim Specified string after results (e.g. date) that points to a specific set of scenario folders
# @param spec_str_rki1 Specified string in results folder (e.g. date) that points to a specific RKI data folder
# @param spec_str_rki2 Specified string in results file that points to a specific RKI data file
# @param scenario_list List of string indicators for scenarios to be plotted
# @param percentiles List of percentiles to be printed (sublist from ['p50','p25','p75','p05','p95'])
# @param read_reports_extrapolation Defines if extrapolated reporting data (from RKI) will be loaded
def open_files(
        path_sim=path_sim, spec_str_sim=date_str, path_rki=path_rki,
        spec_str_rki1=date_str, spec_str_rki2='', scenario_list=scenario_list,
        percentiles=['p50', 'p25', 'p75', 'p05', 'p95'],
        read_casereports_extrapolation=False):

    files = {}

    for scenario in scenario_list:
        files[scenario] = {}

        path = "results_test2020_0_szenario_1_vacc_camp_1_vacc_eff_1"

        for p in percentiles:
            files[scenario][p] = h5py.File(
                path + '/' + p + '/Results_sum.h5', 'r')

        if read_casereports_extrapolation:
            files[scenario]['RKI'] = h5py.File(
                path_rki + spec_str_rki1 + '/Results_rki_sum' + spec_str_rki2 +
                '.h5', 'r')

    return files

# Closes file handles in @files
# @param files File handles of open HDF5 files


def close_files(files):
    for group in files:
        for file in files[group]:
            files[group][file].close()


# define colors for age groups
def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


plt_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

colors = {}
colors['Total'] = plt_colors[0]
for i in range(len(age_groups)):
    colors['Group' + str(i+1)] = plt_colors[i+1]

plotRKI = True           # Plots RKI Data if true
plotRelative = False     # Plots incidence values if true
plotPercentiles = True  # Plots 25 and 75 percentiles if true
plotConfidence = False   # Plots 05 and 95 percentiles if true

savePlot = True          # saves plot file if true
if savePlot:
    try:
        os.mkdir('Plots')
    except:
        print('Directory "Plots" already exists')

opacity = 0.15
lineWidth = 3.5
fontsize = 28
figsize = (13, 10)

# define x-ticks for plots
datelist = np.array(pd.date_range(
    start_date.date(),
    periods=daysPlot, freq='D').strftime('%m-%d').tolist())
tick_range = (np.arange(int(daysPlot / 10) + 1) * 10)
tick_range[-1] -= 1

# Plots one compartment for an individual scenario


def plot_results(
        files, comp_idx, title, ylim=None, filename='', key='Total',
        plotLegend=True, addVal=0, regionid='0'):
    fig, ax = plt.subplots(figsize=figsize)

    if plotRelative:
        factor = relative_dict[key]
    else:
        factor = 1

    X = np.arange(daysPlot)

    ax.plot(
        X, (addVal + files['p50'][regionid][key][: daysPlot, comp_idx]) /
        factor, label='p50', color=colors[key],
        linewidth=lineWidth)
    if plotPercentiles:
        ax.plot(
            X, (addVal + files['p25'][regionid][key][: daysPlot, comp_idx]) /
            factor, '--', label='p25', color=colors[key],
            linewidth=lineWidth)
        ax.plot(
            X, (addVal + files['p75'][regionid][key][: daysPlot, comp_idx]) /
            factor, '--', label='p75', color=colors[key],
            linewidth=lineWidth)
        ax.fill_between(
            X, (addVal + files['p25'][regionid][key][: daysPlot, comp_idx]) /
            factor,
            (addVal + files['p75'][regionid][key][: daysPlot, comp_idx]) /
            factor, color=colors[key],
            alpha=opacity)
    if plotConfidence:
        ax.plot(
            X, (addVal + files['p05'][regionid][key][: daysPlot, comp_idx]) /
            factor, '--', label='p05', color=colors[key],
            linewidth=lineWidth)
        ax.plot(
            X, (addVal + files['p95'][regionid][key][: daysPlot, comp_idx]) /
            factor, '--', label='p95', color=colors[key],
            linewidth=lineWidth)
        ax.fill_between(
            X, (addVal + files['p05'][regionid][key][: daysPlot, comp_idx]) /
            factor,
            (addVal + files['p95'][regionid][key][: daysPlot, comp_idx]) /
            factor, color=colors[key],
            alpha=opacity)

    if plotRKI:
        if 'RKI' in files.keys():
            ax.plot(
                X, files['RKI'][regionid][key][: daysPlot, comp_idx] /
                factor, '--', label='extrapolated real data', color='gray',
                linewidth=lineWidth)
        else:
            print('Error: Plotting extrapolated real data demanded but not read in.')

    ax.set_title(title, fontsize=fontsize)
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45, fontsize=fontsize)
    if plotRelative:
        ax.set_ylabel('individuals relative per 100.000', fontsize=fontsize)
    else:
        ax.set_ylabel('number of individuals', fontsize=fontsize)
    if plotLegend:
        ax.legend(fontsize=fontsize, loc='upper left')
    plt.yticks(fontsize=fontsize)
    ax.grid(linestyle='dotted')

    # if str(ylim) != 'None':
    #     if '_high' in filename:
    #         ylim[filename.replace('_high', '')] = ax.get_ylim()[1]
    #     else:
    #         ax.set_ylim(top=ylim[filename])

    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    ax.yaxis.set_major_formatter(formatter)
    ax.yaxis.offsetText.set_fontsize(fontsize)

    fig.tight_layout()

    if savePlot:
        if plotRKI:
            fig.savefig(
                'Plots/RKI_' + title.replace(' ', '_') + filename + '.png')
        else:
            fig.savefig('Plots/' + title.replace(' ', '_') + filename + '.png')

    plt.clf()

    return ylim


# Plots one compartment for all scenarios


def plot_all_results(
        all_files, comp_idx, title, filename='', key='Total', show_perc=False,
        regionid='0'):
    fig, ax = plt.subplots(figsize=figsize)

    for scenario, color in zip(all_files, list(colors.values())[1:]):
        files = all_files[scenario]
        if plotRelative:
            factor = relative_dict[key]
        else:
            factor = 1

        X = np.arange(daysPlot)

        ax.plot(
            X, files['p50'][regionid][key][: daysPlot, comp_idx] / factor,
            label=scenario_label[scenario],
            color=color, linewidth=lineWidth)

        if show_perc:
            ax.plot(
                X, files['p25'][regionid][key][: daysPlot, comp_idx] /
                factor, '--', color=color, linewidth=lineWidth)
            ax.plot(
                X, files['p75'][regionid][key][: daysPlot, comp_idx] /
                factor, '--', color=color, linewidth=lineWidth)
            ax.fill_between(
                X, files['p25'][regionid][key][: daysPlot, comp_idx] /
                factor, files['p75'][regionid][key][: daysPlot, comp_idx] /
                factor, color=color, alpha=opacity)

    if plotRKI:
        if 'RKI' in files.keys():
            ax.plot(
                X, files['RKI'][regionid][key][: daysPlot, comp_idx] /
                factor, '--', label='extrapolated real data', color='gray',
                linewidth=lineWidth)
        else:
            print('Error: Plotting extrapolated real data demanded but not read in.')

    ax.set_title(title, fontsize=fontsize)
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45, fontsize=fontsize)
    if plotRelative:
        ax.set_ylabel('individuals relative per 100.000', fontsize=fontsize)
    else:
        ax.set_ylabel('number of individuals', fontsize=fontsize)
    ax.legend(fontsize=20)
    plt.yticks(fontsize=fontsize)
    ax.grid(linestyle='dotted')

    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    ax.yaxis.set_major_formatter(formatter)
    ax.yaxis.offsetText.set_fontsize(fontsize)

    fig.tight_layout()

    if savePlot:
        fig.savefig('Plots/' + title.replace(' ', '_') + filename + '.png')


fs = 20


def plot_bars(
        show_perc, name, files, columns, rows, lim=100, rel=False, compart=9,
        regionid='0'):
    color = plt.rcParams['axes.prop_cycle'].by_key()['color']

    n_rows = len(rows)
    num_groups = 6
    bar_width = 1/len(files)
    index = np.arange(n_rows)*bar_width
    scen_width = 8/len(files)

    keys = list(files.keys())

    cell_text = []
    fig, ax = plt.subplots(figsize=(20, 10))
    #ax = fig.add_axes([0,0,1,1])
    for i in range(n_rows):
        key = 'Group' + str(i+1)
        for j in range(len(files)):
            factor = relative_dict[key]
            ax.bar(
                index[i] + j * scen_width,
                (files[keys[j]]['p50'][regionid][key][-1, compart] -
                 files[keys[j]]['p50'][regionid][key][0, compart]) / factor,
                color=color[i],
                width=bar_width, edgecolor='black')
            if show_perc and not (j == 0 or j == 4):
                ax.bar(
                    index[i] + j * scen_width,
                    (files[keys[j]]['p75'][regionid][key][-1, compart] -
                     files[keys[j]]['p75'][regionid][key][0, compart]) /
                    factor, color=color[i],
                    width=bar_width, edgecolor='black', alpha=0.6)

        cell_text.append(
            ['%1.1f' %
             ((files[keys[x]]['p50'][regionid][key][-1, compart] -
               files[keys[x]]['p50'][regionid][key][0, compart]) / factor)
             for x in range(len(files))])

    if len(files) == 8:
        ax.set_xlim(-0.2, 7.8)
    else:
        ax.set_xlim(-0.5, 7.5)
    '''_, top = ax.get_ylim()
    if top > lim:
        ax.set_ylim(top=lim)'''

    # Add a table at the bottom of the axes
    the_table = plt.table(cellText=cell_text,
                          rowLabels=rows,
                          rowColours=color,
                          colLabels=columns,
                          fontsize=fs+10,
                          loc='bottom',
                          cellLoc='center')
    the_table.auto_set_font_size(False)
    the_table.scale(1, 2)
    the_table.set_fontsize(fs-2)
    # Adjust layout to make room for the table:
    #plt.subplots_adjust(left=0.2, bottom=0.8)

    if rel:
        ax.set_ylabel('age distributed infections [%]', fontsize=fs)
        plt.title(
            'Age distributed ratios of infected on September 3', fontsize=fs +
            10)
    else:
        ax.set_ylabel('individuals relative per 100.000', fontsize=fs)
        plt.title('Cumulative number of infections per 100.000', fontsize=fs+10)
    #plt.yticks(values * value_increment, ['%d' % val for val in values])
    plt.xticks([])
    plt.yticks(fontsize=fs)
    plt.tight_layout()
    plt.savefig(name + '.png')

# Combines all compartments of a type into one (e.g. H_s, H_pv, H_v  ->  H)


def concat_comparts(files, comparts, scenario_list, regionid='0'):
    new_files = {}
    for scenario in scenario_list:
        new_files[scenario] = {}
        percentile_list = ['p50', 'p25', 'p75', 'p05', 'p95']
        if plotRKI:
            if 'RKI' in files[scenario].keys():
                percentile_list += ['RKI']
            else:
                print(
                    'Error in concat_comparts(). Extrapolated real data demanded but not read in.')
        for p in percentile_list:
            new_files[scenario][p] = {'0': {}}
            for key in [
                    'Group' + str(group + 1) for group in range(6)] + ['Total']:
                new_files[scenario][p][regionid][key] = np.zeros(
                    (len(files[scenario][p][regionid]['Time']), len(comparts)))
                for new_comp in range(len(comparts)):
                    for old_comp in comparts[new_comp]:
                        new_files[scenario][p][regionid][key][:, new_comp] += \
                            files[scenario][p][regionid][key][:, old_comp]

    return new_files


# This Cell plots all Scenarios with combined compartments
year, month, day = '2022', '6', '1'
start_date = pd.Timestamp(year + '.' + month.zfill(2) + '.' + day.zfill(2))
tmax = '10'
daysPlot = 10

datelist = np.array(pd.date_range(
    start_date.date(),
    periods=daysPlot, freq='D').strftime('%m-%d').tolist())
tick_range = (np.arange(int(daysPlot / 10) + 1) * 10)
tick_range[-1] -= 1
plotRKI = True
plotLegend = True
ylim = {}

secir_dict_inv = dict([(val, key) for key, val in secir_dict.items()])
# List of integers corresponding to the compartments to plot
# e.g. comparts = [4, 6, 8] would only plot Infected, ICU and Dead
for compart, compart_label in secir_dict_aggregated.items():
    # if compart_label == 'Infected':
    for high in ['_high', '']:
        scenario_list = [""]

        files = open_files(
            spec_str_sim='_rev2', spec_str_rki1='',
            scenario_list=scenario_list, read_casereports_extrapolation=False)
        new_files = concat_comparts(files, concat_comps, scenario_list)
        #plot_all_results(new_files, compart, secir_dict_aggregated[compart], filename=high + '_all', key='Total')
        for scenario in scenario_list:
            if True:  # compart_label in ['Infected','ICU','Dead']:
                if compart_label == 'Infected':
                    plotLegend = True
                else:
                    plotLegend = False
                if compart_label == 'Dead' and plotRKI and 'RKI' in files.keys():
                    addBase = files[scenario]['RKI'][regionid]['Total'][
                        :, secir_dict_inv['Dead']][0]
                else:
                    addBase = 0
                ylim = plot_results(
                    new_files[scenario],
                    compart, secir_dict_aggregated[compart],
                    ylim, filename=scenario, key='Total',
                    plotLegend=plotLegend, addVal=addBase)

        close_files(files)


# # Plots for Future
# year, month, day = '2021', '10', '15'
# start_date = pd.Timestamp(year + '.' + month.zfill(2) + '.' + day.zfill(2))
# tmax = '90'
# daysPlot = 90

# datelist = np.array(pd.date_range(
#     start_date.date(),
#     periods=daysPlot, freq='D').strftime('%m-%d').tolist())
# tick_range = (np.arange(int(daysPlot / 10) + 1) * 10)
# tick_range[-1] -= 1
# plotRKI = False
# plotPercentiles = True
# ylim = {}
# comparts = [4]
# scenario_list = ['']
# for compart in comparts:

#     files = open_files(spec_str_sim='_rev2_', spec_str_rki1='',
#                        spec_str_rki2='_future', scenario_list=scenario_list,
#                        read_casereports_extrapolation=True)
#     new_files = concat_comparts(files, concat_comps, scenario_list)
#     plot_all_results(
#         new_files, compart, secir_dict_aggregated[compart],
#         filename=high + '_all', key='Total', show_perc=True)
#     for scenario in scenario_list:
#         ylim = plot_results(
#             new_files[scenario],
#             compart, secir_dict_aggregated[compart],
#             ylim, filename=scenario, key='Total')

#     close_files(files)


# # Barplot for all Scenarios except future
# scenario_list = ['',  '_late', '_mask_test', '_late_mask_test',
#                  '_high',  '_high_late', '_high_mask_test', '_high_late_mask_test']
# plotRKI = False
# files = open_files(spec_str_sim='_rev2', spec_str_rki1='',
#                    scenario_list=scenario_list,
#                    read_casereports_extrapolation=True)
# new_files = concat_comparts(files, concat_comps, scenario_list)


# columns = ('S1 - 40%', 'S2 - 40%', 'S3 - 40%', 'S4 - 40%',
#            'S1 - 60%', 'S2 - 60%', 'S3 - 60%', 'S4 - 60%')
# rows = ['0-4 Years', '5-14 Years', '15-34 Years',
#         '35-59 Years', '60-79 Years', '80+ Years']

# plot_bars(True, 'age_incidence', new_files, columns, rows, 200)
# close_files(files)

# # Barplot for Future Scenarios
# scenario_list = ['s1f_future_high_mask_test',
#                  's2f_future_high_mask_test', 's3f_future_high_mask_test']
# files = open_files(spec_str_sim='_rev2_', spec_str_rki1='',
#                    spec_str_rki2='_future', scenario_list=scenario_list,
#                    read_casereports_extrapolation=True)
# new_files = concat_comparts(files, concat_comps, scenario_list)

# fs = 24

# columns = ('S1F: contact reduc. 27 % [21-32 %]',
#            'S2F: contact reduc. 37% [32-42 %]',
#            'S3F: contact reduc. 42% [37-46]')
# rows = ['0-4 Years', '5-14 Years', '15-34 Years',
#         '35-59 Years', '60-79 Years', '80+ Years']

# plot_bars(False, 'age_incidence', new_files, columns, rows, 200)
# close_files(files)


# # Plot for Partially Vaccinated, Vaccinated and Cumulative number of Infections
# year, month, day = '2021', '6', '6'
# start_date = pd.Timestamp(year + '.' + month.zfill(2) + '.' + day.zfill(2))
# tmax = '90'
# daysPlot = 91
# regionid = '0'

# datelist = np.array(pd.date_range(
#     start_date.date(),
#     periods=daysPlot, freq='D').strftime('%m-%d').tolist())
# tick_range = (np.arange(int(daysPlot / 10) + 1) * 10)
# tick_range[-1] -= 1
# plotRKI = False
# plotPercentiles = False
# savePlot = True
# ylim = {}
# comparts = [1, 7, 9]
# scenario_list = ['_high_mask_test']

# files = files = open_files(
#     spec_str_sim='_rev2', spec_str_rki1='', scenario_list=scenario_list)
# new_files = concat_comparts(files, concat_comps, scenario_list)
# new_files[scenario_list[0]]['p50'][regionid]['Total'][:,
#                                                       7] -= new_files[scenario_list[0]]['p50'][regionid]['Total'][:, 9]
# temp_secir_dict = {
#     1: 'Partially Vaccinated',
#     7: 'Fully Vaccinated',
#     9: 'Cumulative number of Infections'
# }
# for compart in comparts:

#     for scenario in scenario_list:
#         plot_results(
#             new_files[scenario],
#             compart, temp_secir_dict[compart],
#             None, filename='vacc' + scenario, key='Total', plotLegend=False)

#     close_files(files)

# # This Cell Plots all Scenarios, where the Y limits for "high" and
# # not "high" Scenarios are equal (Without "concat_comparts")
# year, month, day = '2021', '6', '6'
# start_date = pd.Timestamp(year + '.' + month.zfill(2) + '.' + day.zfill(2))
# tmax = '90'
# daysPlot = 91

# datelist = np.array(pd.date_range(
#     start_date.date(),
#     periods=daysPlot, freq='D').strftime('%m-%d').tolist())
# tick_range = (np.arange(int(daysPlot / 10) + 1) * 10)
# tick_range[-1] -= 1
# plotRKI = False
# ylim = {}
# comparts = secir_dict
# for compart in comparts:
#     for high in ['_high', '']:
#         scenario_list = [high, high + '_late', high + '_mask_test',
#                          high + '_late_mask_test', '_long' + high + '_mask_test']

#         files = open_files('_init_test')
#         new_files = concat_comparts(files, concat_comps)
#         #plot_all_results(new_files, compart, new_secir_dict[compart], filename=high + '_all', key='Total')
#         for scenario in scenario_list:
#             ylim = plot_results(
#                 files[scenario],
#                 compart, secir_dict[compart],
#                 ylim, filename=scenario, key='Total')

#         close_files(files)


# # This Cell only plots one Scenario (Without "concat_comparts")
# year, month, day = '2021', '6', '6'
# start_date = pd.Timestamp(year + '.' + month.zfill(2) + '.' + day.zfill(2))
# tmax = '90'
# daysPlot = 91

# datelist = np.array(pd.date_range(
#     start_date.date(),
#     periods=daysPlot, freq='D').strftime('%m-%d').tolist())
# tick_range = (np.arange(int(daysPlot / 10) + 1) * 10)
# tick_range[-1] -= 1

# comparts = secir_dict
# for compart in comparts:
#     # List of strings corresponding to Scenarios
#     scenario_list = ['_high_mask_test']

#     files = open_files()
#     new_files = concat_comparts(files, concat_comps)
#     #plot_all_results(new_files, compart, new_secir_dict[compart], filename=high + '_all', key='Total')
#     ylim = None
#     for scenario in scenario_list:
#         ylim = plot_results(
#             files[scenario],
#             compart, secir_dict[compart],
#             ylim, filename=scenario, key='Total')

#     close_files(files)


# # This Cell plots all Scenarios with combined compartments

# year, month, day = '2021', '6', '6'
# start_date = pd.Timestamp(year + '.' + month.zfill(2) + '.' + day.zfill(2))
# tmax = '90'
# daysPlot = 91

# datelist = np.array(pd.date_range(
#     start_date.date(),
#     periods=daysPlot, freq='D').strftime('%m-%d').tolist())
# tick_range = (np.arange(int(daysPlot / 10) + 1) * 10)
# tick_range[-1] -= 1
# plotRKI = False
# ylim = {}
# # List of integers corresponding to the compartments to plot
# # e.g. comparts = [4, 6, 8] would only plot Infected, ICU and Dead
# comparts = new_secir_dict
# for compart in comparts:
#     # for high in ['_high', '']:
#     #scenario_list = [high, high + '_late', high + '_mask_test', high + '_late_mask_test','_long' +  high + '_mask_test']
#     scenario_list = ['_high_mask_test']
#     files = open_files()
#     new_files = concat_comparts(files, concat_comps)
#     #plot_all_results(new_files, compart, new_secir_dict[compart], filename=high + '_all', key='Total')
#     for scenario in scenario_list:
#         ylim = plot_results(
#             new_files[scenario],
#             compart, new_secir_dict[compart],
#             ylim, filename=scenario, key='Total')

#     close_files(files)
