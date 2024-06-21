import datetime as dt
import os.path
import imageio

import numpy as np
import pandas as pd

from memilio.epidata import modifyDataframeSeries as mdfs
import memilio.epidata.getPopulationData as gpd
import memilio.plot.plotMap as pm
from memilio.epidata import geoModificationGermany as geoger
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import SymLogNorm, LinearSegmentedColormap
from tqdm.auto import tqdm
from datetime import datetime, timedelta
import h5py

import seaborn as sns

sns.set_style("darkgrid")

start_date = "2020-12-01"
total_pop = 83278910.0
opacity = 0.15
lineWidth = 2
fontsize = 28
legendsize = 15
ticks = 15
# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors = sns.color_palette("tab20")


def read_total_results_h5(path, comp, group_key='Total'):
    f = h5py.File(path, 'r')
    group = f['0']
    total = group[group_key][()]
    comp_simulated = np.sum(total[:, comp], axis=1)
    f.close()

    return comp_simulated


def get_state_id(county_id):
    county_id_str = str(county_id)
    state_id = None
    if len(county_id_str) == 4:
        # four digit numbers: first digit is the State-ID
        state_id = int(county_id_str[0])
    elif len(county_id_str) == 5:
        # five digit numbers: first two digits are the State-ID
        state_id = int(county_id_str[:2])
    else:
        raise ValueError(
            "Invalid county_id length. Must be 4 or 5 digits long.")
    return state_id


def read_state_results_h5(path, comp, group_key='Total'):
    f = h5py.File(path, 'r')
    # create array with shape 16, county filled with zeros
    arr = np.zeros((16, f['1001']['Total'].shape[0]))
    for key in f.keys():
        group = f[key]
        total = group[group_key][()]
        comp_simulated = np.sum(total[:, comp], axis=1)
        state_id = get_state_id(key)
        # state_id - 1 because state_id starts at 1
        arr[state_id - 1] += comp_simulated
    f.close()
    return arr


def read_county_results_h5(path, comp, group_key='Total'):
    data = {}
    with h5py.File(path, 'r') as f:
        for key in f.keys():
            group = f[key]
            total = group[group_key][()]
            comp_simulated = np.sum(total[:, comp], axis=1)
            data[key] = comp_simulated
    df = pd.DataFrame(data)
    return df


def plot(ys, labels, path_plots, title="", log_scale=False, ylabel="Number Individuals", plot_percentiles=True):
    num_data = len(ys)

    # Set days for x-axis
    num_days = len(ys[0]) - 1
    if isinstance(ys[0], dict):
        num_days = len(ys[0]['p25']) - 1
    start_date_datetime = datetime.strptime(start_date, "%Y-%m-%d")
    end_date_datetime = start_date_datetime + timedelta(days=num_days)
    days = pd.date_range(start_date_datetime, end_date_datetime)
    months = pd.date_range(start=start_date_datetime,
                           end=end_date_datetime, freq='MS')

    # Creating subplots based on the number of data series
    if num_data > 8:
        fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    else:
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        axes = [ax]

    for ax in axes:
        # Set ticks and labels
        ax.set_xticks(months)
        ax.set_xticklabels(months.strftime('%B %Y'),
                           fontsize=ticks, rotation=45)
        ax.tick_params(axis='y', labelsize=ticks)  # Set y-axis tick size

        if title:
            ax.set_title(title, fontsize=fontsize)

        ax.set_xlabel("Time [days]", fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.grid(True)
        if log_scale:
            ax.set_yscale('log')

    # Plotting the data
    for indx_data in range(num_data):
        y = ys[indx_data]
        ax = axes[indx_data % len(axes)]
        if isinstance(y, dict):
            ax.plot(days, y["p50"], label=labels[indx_data], linewidth=lineWidth,
                    linestyle='-', color=colors[indx_data])
            if plot_percentiles:
                ax.plot(days, y["p25"], linewidth=lineWidth,
                        linestyle='--', color=colors[indx_data])
                ax.plot(days, y["p75"], linewidth=lineWidth,
                        linestyle='--', color=colors[indx_data])
                ax.fill_between(
                    days, y["p25"], y["p75"], color=colors[indx_data], alpha=opacity)
        else:
            ax.plot(days, y, label=labels[indx_data], linewidth=lineWidth,
                    linestyle='-', color=colors[indx_data])

    for ax in axes:
        ax.legend(fontsize=legendsize, loc='center left',
                  bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig(os.path.join(path_plots, f"plot_{title}.png"))
    plt.clf()
    return 0


def get_results(path_results, indx_comp, group_key='Total', results="total", percentiles=["p25", "p50", "p75"]):

    y = {}
    for px in percentiles:
        if results == "total":
            path = os.path.join(path_results, px, "Results_sum.h5")
            y[px] = read_total_results_h5(path, indx_comp, group_key)
        elif results == "state":
            path = os.path.join(path_results, px, "Results.h5")
            y[px] = read_state_results_h5(path, indx_comp, group_key)

    # change shape of data if results are for states
    if results == "state":
        y = [{px: y[px][state] for px in percentiles}
             for state in range(16)]

    return y


def plot_risk(path_results, path_plots, log_scale=False, plot_percentiles=False):
    create_folder_if_not_exists(path_plots)
    num_counties = 400
    path_risk_results = os.path.join(path_results, "FeedbackDamping", "risk")
    # total risk
    plot_data = []
    plot_data.append(get_results(path_risk_results, [
                     0], 'Group1', results="total"))
    # divide every value by num_counties
    for key in plot_data[0].keys():
        for indx in range(len(plot_data[0][key])):
            plot_data[0][key][indx] = plot_data[0][key][indx] / num_counties
    # get total and local population and scale the risk
    labels = ["Total Risk"]
    plot(plot_data, labels, path_plots, title="Total_Risk",
         log_scale=log_scale, ylabel="Perceived Risk")

    # Risk per State
    counties_per_state = geoger.get_stateid_to_countyids_map()
    plot_data = get_results(
        path_risk_results, [0], 'Group1', results="state")
    # divide every value by num_counties in the state
    for state in range(16):
        for key in plot_data[state].keys():
            for indx in range(len(plot_data[state][key])):
                plot_data[state][key][indx] = plot_data[state][key][indx] / \
                    len(counties_per_state[state+1])
    labels = geoger.get_state_names()
    plot(plot_data, labels, path_plots, title="State_Risk",
         log_scale=log_scale, ylabel="Perceived Risk", plot_percentiles=plot_percentiles)


def plot_compartments(path_results, path_plots, modes, compartments, labels, title, log_scale=False):
    create_folder_if_not_exists(path_plots)
    plot_data = []
    labels_modes = []
    for mode in modes:
        for index, compartment in enumerate(compartments):
            labels_modes.append(labels[index] + f" {mode}")
            path_results_mode = os.path.join(path_results, mode)
            plot_data.append(get_results(
                path_results_mode, compartment, results="total"))
    plot(plot_data, labels_modes, path_plots, title=title,
         log_scale=log_scale, ylabel="Number Individuals")
    return 0


def plot_flows(path_results, path_plots, modes, flow_indx, labels, title, log_scale=False):
    create_folder_if_not_exists(path_plots)
    plot_data = []
    labels_modes = []
    for mode in modes:
        for index, compartment in enumerate(flow_indx):
            labels_modes.append(labels[index] + f" {mode}")
            path_results_mode = os.path.join(path_results, mode, "flows")
            data = get_results(
                path_results_mode, compartment, results="total")
            data = {key: np.diff(data[key]) for key in data.keys()}
            plot_data.append(data)
    plot(plot_data, labels_modes, path_plots, title=title,
         log_scale=log_scale, ylabel="Number Individuals")
    return 0


def plot_r0(path_results, path_plots, modes, percentile="p50"):
    population = get_pop()
    total_pop = population['Population'].sum()

    r0_nums = []
    for mode in modes:
        path_results_mode = os.path.join(path_results, mode, "r0")
        path = os.path.join(path_results_mode, percentile,  "Results.h5")
        with h5py.File(path, 'r') as f:
            num_days = f['1001']['Group1'].shape[0]
            r0_mode = np.zeros(num_days)
            indx_dict = {
                key: population.index[population['ID_County'] == int(key)] for key in f.keys()}
            for key, group in f.items():
                indx_key = indx_dict[key]
                r0_mode += group['Group1'][()][:, 0] * \
                    population['Population'][indx_key].values[0] / total_pop
            r0_nums.append(r0_mode)

    plot(r0_nums, modes, path_plots, title="R0", ylabel="R0")


def plot_r0_county_level(path_results, path_plots, modes, percentile="p50"):
    df_modes = []
    for mode in modes:
        path_results_mode = os.path.join(path_results, mode, "r0")
        path = os.path.join(path_results_mode, percentile,  "Results.h5")
        data = {}
        with h5py.File(path, 'r') as f:
            for key in f.keys():
                group = f[key]
                total = group['Group1'][()]
                data[key] = np.sum(total[:], axis=1)
        df = pd.DataFrame(data)
        df_modes.append(df)

    # Box Plot
    plt.figure(figsize=(12 * len(modes), 6))
    y_min, y_max = float('inf'), float('-inf')

    for i, mode in enumerate(modes):
        df = df_modes[i]
        # every second day
        df = df.iloc[::10]
        df = df.T
        plt.subplot(1, len(modes), i + 1)
        sns.boxplot(data=df, orient="v")
        plt.xlabel('Time', fontsize=fontsize)
        if i == 0:  # Nur im ersten Plot das y-Label setzen
            plt.ylabel('Reproduction Number', fontsize=fontsize)
        tick_positions = np.arange(0, len(df.columns), 5)
        tick_labels = df.columns[tick_positions]
        plt.xticks(ticks=tick_positions, labels=tick_labels, fontsize=ticks)
        plt.title(f'Distribution of Reproduction Numbers Over Time {mode}')
        y_min = min(y_min, df.min().min())
        y_max = max(y_max, df.max().max())

    for i in range(len(modes)):  # Setzen der gleichen Skala für alle Plots
        plt.subplot(1, len(modes), i + 1)
        plt.ylim(0.4, 1.7)

    plt.savefig(os.path.join(path_plots, 'box_plots.png'))

    # Violin Plot
    plt.figure(figsize=(12 * len(modes), 6))
    y_min, y_max = float('inf'), float('-inf')

    for i, mode in enumerate(modes):
        df = df_modes[i]
        df = df.T
        plt.subplot(1, len(modes), i + 1)
        sns.violinplot(data=df, orient="v")
        plt.xlabel('Time', fontsize=fontsize)
        if i == 0:  # Nur im ersten Plot das y-Label setzen
            plt.ylabel('Reproduction Number', fontsize=fontsize)
        tick_positions = np.arange(0, len(df.columns), 5)
        tick_labels = df.columns[tick_positions]
        plt.xticks(ticks=tick_positions, labels=tick_labels, fontsize=ticks)
        plt.title(f'Distribution of Reproduction Numbers Over Time {mode}')
        y_min = min(y_min, df.min().min())
        y_max = max(y_max, df.max().max())

    for i in range(len(modes)):  # Setzen der gleichen Skala für alle Plots
        plt.subplot(1, len(modes), i + 1)
        plt.ylim(y_min, y_max)

    plt.savefig(os.path.join(path_plots, 'violin_plots.png'))


def plot_contacts(path_results, path_plots, modes, percentile="p50"):
    population = get_pop()
    total_pop = population['Population'].sum()
    # adjust age groups
    old_ages = [age.split()[0] for age in population.columns[2:]]
    new_ages = ["0-4", "5-14", "15-34", "35-59", "60-79", ">80"]
    population.rename(
        dict(zip(population.columns[2:].tolist(), old_ages)), axis=1, inplace=True)
    for county in population.ID_County.unique():
        population.loc[population.ID_County == county, new_ages] = mdfs.fit_age_group_intervals(
            population.loc[population.ID_County == county, old_ages], new_ages)
    population.drop(old_ages, axis=1, inplace=True)

    contacts_nums = []
    contacts_avg = []
    for mode in modes:
        path_results_mode = os.path.join(path_results, mode, "contacts")
        path = os.path.join(path_results_mode, percentile,  "Results.h5")
        with h5py.File(path, 'r') as f:
            num_days = f['1001']['Group1'].shape[0]
            contacts_mode = np.zeros(num_days)
            # Berechne die Indizes einmal und speichere sie in einem Wörterbuch
            indx_dict = {
                key: population.index[population['ID_County'] == int(key)] for key in f.keys()}
            for key, group in f.items():
                indx_key = indx_dict[key]
                avg_contacts = (group['Group1'][()][:, 0]*population[new_ages[0]][indx_key].values[0]/population['Population'][indx_key].values[0] +
                                group['Group2'][()][:, 0]*population[new_ages[1]][indx_key].values[0]/population['Population'][indx_key].values[0] +
                                group['Group3'][()][:, 0]*population[new_ages[2]][indx_key].values[0]/population['Population'][indx_key].values[0] +
                                group['Group4'][()][:, 0]*population[new_ages[3]][indx_key].values[0]/population['Population'][indx_key].values[0] +
                                group['Group5'][()][:, 0]*population[new_ages[4]][indx_key].values[0]/population['Population'][indx_key].values[0] +
                                group['Group6'][()][:, 0]*population[new_ages[5]][indx_key].values[0]/population['Population'][indx_key].values[0])
                contacts_mode += avg_contacts * \
                    population['Population'][indx_key].values[0] / total_pop

            # contacts_mode is accumulated. Get the diffs
            contacts_mode = np.diff(contacts_mode)
            contacts_nums.append(contacts_mode)

        if mode != "ClassicDamping":
            # get average of contacts_mode
            avg_contacts_mode = np.mean(contacts_mode)
            contacts_avg.append([avg_contacts_mode for _ in range(num_days-1)])

    labels = modes + ["Average FeedbackDamping"]
    contacts_nums += contacts_avg

    # calculate contacts needed to match average feedback
    total_group0 = population[new_ages[0]].sum()
    total_group1 = population[new_ages[1]].sum()
    total_group2 = population[new_ages[2]].sum()
    total_group3 = population[new_ages[3]].sum()
    total_group4 = population[new_ages[4]].sum()
    total_group5 = population[new_ages[5]].sum()

    contacts_group1 = 10.4699
    contacts_group2 = 7.3629
    contacts_group3 = 10.1193
    contacts_group4 = 9.3645
    contacts_group5 = 3.7438
    contacts_group6 = 2.5016

    contacts_avg = contacts_avg[0][0]

    x = contacts_avg / ((contacts_group1*total_group0/total_pop +
                         contacts_group2*total_group1/total_pop +
                         contacts_group3*total_group2/total_pop +
                         contacts_group4*total_group3/total_pop +
                         contacts_group5*total_group4/total_pop +
                         contacts_group6*total_group5/total_pop))

    print(f"Contact Reduction needed to match average feedback: {1- x:.10f}")

    plot(contacts_nums, labels, path_plots,
         title="Contacts", ylabel="Number of Contacts")


def plot_icu_comp(path_results, path_plots, modes, path_icu_data, log_scale=False):
    create_folder_if_not_exists(path_plots)
    icu_comp = [7]
    label = "ICU Occupancy"
    labels = []
    title = "ICU Occupancy per 100_000 inhabitants"

    plot_data = []
    for mode in modes:
        path_results_mode = os.path.join(path_results, mode)
        labels.append(label + f" {mode}")
        plot_data.append(get_results(
            path_results_mode, icu_comp, results="total"))
    # calculate ICU occupancy per 100_000 inhabitants
    for data in plot_data:
        for key in data.keys():
            for indx in range(len(data[key])):
                data[key][indx] = data[key][indx] / total_pop * 100_000

    # create dict with same shape and set constant value for ICU capacity
    icu_capacity_val = 12
    icu_capacity = {key: [icu_capacity_val for _ in range(
        len(plot_data[0][key]))] for key in plot_data[0].keys()}
    plot_data.append(icu_capacity)
    labels.append("ICU Capacity")

    num_days = len(plot_data[0]['p25']) - 1
    start_date_datetime = datetime.strptime(start_date, "%Y-%m-%d")
    end_date_datetime = start_date_datetime + timedelta(days=num_days)
    end_date = end_date_datetime.strftime("%Y-%m-%d")

    df = pd.read_json(path_icu_data)

    filtered_df = df.loc[(df['Date'] >= start_date) &
                         (df['Date'] <= end_date)]
    icu_divi = filtered_df['ICU'].to_numpy() / total_pop * 100_000
    plot_data.append(icu_divi)
    labels.append("ICU Divi Data")

    plot(plot_data, labels, path_plots, title=title,
         log_scale=log_scale, ylabel="ICU Occupancy per 100_000")
    return 0


def plot_peaks(path_results, path_plots, modes, target_indx, percentile="p50", flows=True, title="Peaks for each County"):
    create_folder_if_not_exists(path_plots)
    peaks_modes = []

    for mode in modes:
        for index, compartment in enumerate(target_indx):
            path_results_mode = os.path.join(path_results, mode)
            if flows:
                path_results_mode = os.path.join(path_results_mode, "flows")
            df_data = read_county_results_h5(
                os.path.join(path_results_mode, percentile, "Results.h5"), compartment, group_key='Total')
            if flows:
                # if flows, the data is accumulated.
                df_data = df_data.diff(axis=0)
                df_data = df_data.iloc[1:]
                df_data.reset_index(drop=True, inplace=True)
                # set all first values to 0
                # df_data.iloc[0] = 0
            # get the index of the max value in each column and append it to the list
            peaks = df_data.idxmax().values
            # count the number of peaks for each day
            peaks_modes.append(np.bincount(peaks, minlength=df_data.shape[0]))

    plot(peaks_modes, modes, path_plots,
         title=title, ylabel="Number of Peaks")


def plot_peaks_single(path_results, path_plots, target_indx, percentile="p50", flows=True):
    # list all dirs beginning with kmin
    path_plots_peaks = os.path.join(path_plots, "peaks")
    kmin_dirs = [d for d in os.listdir(path_results) if os.path.isdir(
        os.path.join(path_results, d)) and d.startswith("kmin")]
    for runs in kmin_dirs:
        path_results_kmin = os.path.join(path_results, runs)
        plot_peaks(path_results_kmin, path_plots_peaks, modes, target_indx,
                   percentile, flows, title=f"Peaks_{runs}")


def get_pop(file_format='h5'):
    try:
        population = pd.read_json(
            'data/pydata/Germany/county_current_population.json')
    # pandas>1.5 raise FileNotFoundError instead of ValueError
    except (ValueError, FileNotFoundError):
        print(
            "Population data was not found. Download it from the internet.")
        population = gpd.get_population_data(
            read_data=False, file_format=file_format,
            out_folder='data/pydata/Germany/', no_raw=True,
            split_gender=False, merge_eisenach=True)

    return population


def create_folder_if_not_exists(path):
    if not os.path.exists(path):
        os.makedirs(path)


def get_kmin_dirs(path_results):
    return [d for d in os.listdir(path_results) if os.path.isdir(os.path.join(path_results, d)) and d.startswith("kmin")]


def read_data(path, target_indx, percentile, flows):
    df = read_county_results_h5(os.path.join(
        path, percentile, "Results.h5"), target_indx[0], group_key='Total')
    if flows:
        df = df.diff().iloc[1:].reset_index(drop=True)
    return df


def plot_and_save(fig, path, filename):
    plt.tight_layout()
    fig.savefig(os.path.join(path, filename))
    plt.close(fig)


def plot_peak_values(path_results, path_plots, modes, target_indx, percentile="p50", flows=True, plot_type='kmin', title='Peak Value'):
    if len(modes) > 1:
        print("Only one mode is allowed for the grid peak plot.")
        return

    create_folder_if_not_exists(os.path.join(path_plots, "peaks"))

    kmin_dirs = get_kmin_dirs(path_results)
    all_kmins = sorted(set(float(d.split("_")[1]) for d in kmin_dirs))

    num_days = 0

    for kmin in all_kmins:
        fig, axes = plt.subplots(
            1, round(10 - kmin * 10), figsize=(15, 3), squeeze=False)
        all_peaks = []

        for run in kmin_dirs:
            if kmin != float(run.split("_")[1]):
                continue

            path = os.path.join(path_results, run,
                                modes[0], "flows" if flows else "")
            df = read_data(path, target_indx, percentile, flows)

            if num_days == 0:
                num_days = df.shape[0]

            if plot_type == 'kmin':
                peaks = np.bincount(df.idxmax().values, minlength=df.shape[0])
                all_peaks.append(peaks)
                ax = axes[0, len(all_peaks) - 1]
                ax.plot(peaks)
                if len(all_peaks) == 1:
                    ax.set_ylabel('Number of Peaks')
            else:
                peak_indices = df.idxmax().values
                peak_max = df.max().values / \
                    get_pop()['Population'].values * 100_000
                peak_values = {i: [] for i in range(df.shape[0])}
                for idx, val in zip(peak_indices, peak_max):
                    peak_values[idx].append(val)
                df_long = pd.DataFrame(
                    [{'Day': day, 'Peak Value': val} for day, vals in peak_values.items() for val in vals])
                all_peaks.extend(peak_values.values())
                # ensure that all days are in the dataframe
                all_days = pd.DataFrame({'Day': range(num_days)})
                df_long = pd.merge(all_days, df_long, on='Day', how='left')
                ax = axes[0, len(all_peaks) // len(df) - 1]
                sns.boxplot(x='Day', y='Peak Value', data=df_long, ax=ax)
                ax.set_ylabel(title + ' per 100,000k')

            kmax = float(run.split("_")[3])
            ax.set_title(f'kmin: {kmin}\nkmax: {kmax:.1f}', fontsize=12)

        global_min = min(min(peaks) for peaks in all_peaks if peaks)
        global_max = max(max(peaks) for peaks in all_peaks if peaks) + 10
        for ax in axes.flat:
            ax.set_ylim(global_min, global_max)
            ax.set_xticks(range(0, num_days, 100))
            ax.set_xlabel('')
            if ax != axes[0, 0]:
                ax.set_yticklabels([])
                ax.set_ylabel('')
            if ax == axes[0, len(axes[0]) // 2]:
                ax.set_xlabel('Days', fontsize=16)

        plot_and_save(fig, os.path.join(path_plots, "peaks"),
                      f'peaks_grid_{plot_type}_kmin_{kmin}.png')

        for ax in axes.flat:
            ax.set_yscale('symlog')
            ax.set_xticks(range(0, num_days, 100))
            ax.set_xlabel('')
            if ax != axes[0, 0]:
                ax.set_yticklabels([])
                ax.set_ylabel('')
            if ax == axes[0, len(axes[0]) // 2]:
                ax.set_xlabel('Days', fontsize=16)
        plot_and_save(fig, os.path.join(path_plots, "peaks"),
                      f'peaks_grid_{plot_type}_kmin_{kmin}_log.png')


if __name__ == '__main__':
    path_cwd = os.getcwd()
    path_results = os.path.join(path_cwd, "results")
    path_plots = os.path.join(path_cwd, "plots")
    path_icu_data = os.path.join(
        path_cwd, "data/pydata/Germany/germany_divi_ma7.json")

    modes = ["ClassicDamping", "FeedbackDamping"]

    icu_compartment = [[7]]
    infected_compartment = [[1, 2, 3, 4, 5, 6, 7]]
    dead_compartment = [[9]]
    flow_se = [[0]]

    # plot_r0_county_level(path_results, path_plots, modes)
    # plot_contacts(path_results, path_plots, modes)
    # plot_risk(path_results,
    #           path_plots, plot_percentiles=False)
    # plot_compartments(path_results, path_plots, modes,
    #                   icu_compartment, ["ICU Occupancy"], "ICU Occupancy")
    # plot_compartments(path_results, path_plots, modes,
    #                   infected_compartment, [""], "Total Infected")
    # plot_compartments(path_results, path_plots, modes,
    #                   dead_compartment, [""], "Total Deaths")
    # plot_flows(path_results, path_plots, modes,
    #            flow_se, [""], "Daily Infections")
    # plot_icu_comp(path_results, path_plots, modes, path_icu_data)
    # plot_r0(path_results, path_plots, modes)
    # plot_peaks(path_results, path_plots, modes, flow_se)

    # plot_peaks_single(path_results, path_plots, ["FeedbackDamping"], flow_se)
    # plot_peak_values(os.path.join(path_cwd, "results", "save"), path_plots, [
    #                  "FeedbackDamping"], flow_se, plot_type='kmin')
    plot_peak_values(os.path.join(path_cwd, "results", "save"), path_plots, [
                     "FeedbackDamping"], icu_compartment, plot_type='val', title='ICU Occupancy', flows=False)
