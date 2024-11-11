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
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks
import h5py

import seaborn as sns

# sns.set_style("darkgrid")

start_date = "2020-10-01"
total_pop = 83278910.0
opacity = 0.15
lineWidth = 3
fontsize = 25
legendsize = 15
ticks = 18
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors_state_risk = sns.color_palette("tab20")
color_map = 'viridis'


def read_total_results_h5(path, comp, group_key='Total'):
    f = h5py.File(path, 'r')
    group = f['0']
    total = group[group_key][()]
    comp_simulated = np.sum(total[:, comp], axis=1)
    f.close()

    return comp_simulated


def read_results_h5(path, comp, group_key='Total'):
    with h5py.File(path, 'r') as f:
        keys = list(f.keys())
        res = np.zeros((len(keys), f[keys[0]][group_key].shape[0]))
        for i, key in enumerate(keys):
            group = f[key]
            total = group[group_key][()]
            if total.shape[1] > 1:
                comp_simulated = np.sum(total[:, comp], axis=1)
            else:
                comp_simulated = total[:, comp]
            res[i] = comp_simulated

    return res


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
    create_folder_if_not_exists(path_plots)

    # Set days for x-axis
    num_days = len(ys[0]) - 1
    if isinstance(ys[0], dict):
        num_days = len(ys[0]['p50']) - 1
    start_date_datetime = datetime.strptime(start_date, "%Y-%m-%d")
    end_date_datetime = start_date_datetime + timedelta(days=num_days)
    days = pd.date_range(start_date_datetime, end_date_datetime)
    months = pd.date_range(start=start_date_datetime,
                           end=end_date_datetime, freq='MS')

    # Creating subplots based on the number of data series
    colors_plot = colors
    if num_data > 8:
        fig, axes = plt.subplots(1, 2, figsize=(20, 10))
        colors_plot = colors_state_risk
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

        cl = colors_plot[indx_data]

        if isinstance(y, dict):
            ax.plot(days, y["p50"], label=labels[indx_data], linewidth=lineWidth,
                    linestyle='-', color=cl)
            if plot_percentiles:
                ax.plot(days, y["p25"], linewidth=lineWidth,
                        linestyle='--', color=cl)
                ax.plot(days, y["p75"], linewidth=lineWidth,
                        linestyle='--', color=cl)
                ax.fill_between(
                    days, y["p25"], y["p75"], color=cl, alpha=opacity)
        else:
            ax.plot(days, y, label=labels[indx_data], linewidth=lineWidth,
                    linestyle='-', color=cl)

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


def plot_flows(path_results, path_plots, modes, flow_indx, labels, title, log_scale=False, plot_percentiles=True):

    case_data_path = "/localdata1/code_2024/memilio/data/pydata/Germany/cases_all_germany_ma7.json"
    # start_date is format "YYYY-MM-DD". how to get start_date - 1
    start_date_m1 = datetime.strptime(
        start_date, "%Y-%m-%d") - timedelta(days=1)
    df = pd.read_json(case_data_path)
    df = df.loc[(df['Date'] >= start_date_m1)].reset_index(drop=True)
    # drop all cols except 'Date' and 'Confirmed'
    df = df[['Date', 'Confirmed']]
    # data is aggregated. get the diff
    df['Confirmed'] = df['Confirmed'].diff()
    df = df.iloc[1:]

    create_folder_if_not_exists(path_plots)
    plot_data = []
    labels_modes = []
    for mode in modes:
        for index, compartment in enumerate(flow_indx):
            labels_modes.append(labels[index] + f" {mode}")
            path_results_mode = os.path.join(
                path_results, mode, "flows")
            data = get_results(
                path_results_mode, compartment, results="total")
            data = {key: np.diff(data[key]) for key in data.keys()}
            plot_data.append(data)

    num_days = len(plot_data[0]['p25'])
    df_cut = df.iloc[:num_days]
    plot_data.append(df_cut['Confirmed'].to_numpy())
    labels_modes.append("Reported Confirmed Cases")
    plot(plot_data, labels_modes, path_plots, title=title,
         log_scale=log_scale, ylabel="Number Individuals", plot_percentiles=plot_percentiles)
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


def plot_r0_all_scenarios(path_results, path_plots, modes, percentile="p50"):
    population = get_pop()
    total_pop = population['Population'].sum()

    # alle dirs in path_results
    dirs = [d for d in os.listdir(path_results) if os.path.isdir(
        os.path.join(path_results, d))]

    r0_nums = []
    for dir_name in dirs:
        mode = "FeedbackDamping"
        if dir_name.startswith("fixed"):
            mode = 'ClassicDamping'

        path_results_mode = os.path.join(
            path_results, dir_name,  mode, "r0")
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

    plot(r0_nums, dirs, path_plots, title="R0_all_scenarios", ylabel="R0")


def plot_icu_all_scenarios(path_results, path_plots, modes, percentile="p50", log_scale=False):

    icu_compartment = [[7]]

    # alle dirs in path_results
    dirs = [d for d in os.listdir(path_results) if os.path.isdir(
        os.path.join(path_results, d))]

    icu_nums = []
    for dir_name in dirs:
        mode = "FeedbackDamping"
        if dir_name.startswith("fixed"):
            mode = 'ClassicDamping'
        path_results_mode = os.path.join(
            path_results, dir_name,  mode)
        res = get_results(path_results_mode, icu_compartment,
                          results="total", percentiles=["p50"])
        # res only has one entry in a dict. transform to list
        icu_nums.append(res['p50'])

    plot(icu_nums, dirs, path_plots, title="ICU",
         ylabel="ICU Occupancy", log_scale=log_scale)


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
        plt.ylim(y_min - 0.5, y_max + 0.5)

    plt.savefig(os.path.join(path_plots, 'plot_r0_box_plots.png'))

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
        plt.ylim(y_min - 0.5, y_max + 0.5)

    plt.savefig(os.path.join(path_plots, 'violin_plots.png'))


def plot_risk_county_level(path_results, path_plots, modes, percentile="p50"):
    df_modes = []
    for mode in modes:
        path_results_mode = os.path.join(path_results, mode, "risk")
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
        if i == 0:
            plt.ylabel('Perceived Risk', fontsize=fontsize)
        tick_positions = np.arange(0, len(df.columns), 5)
        tick_labels = df.columns[tick_positions]
        plt.xticks(ticks=tick_positions, labels=tick_labels, fontsize=ticks)
        plt.title(f'Distribution of Perceived Risk Over Time {mode}')
        y_min = min(y_min, df.min().min())
        y_max = max(y_max, df.max().max())

    for i in range(len(modes)):
        plt.subplot(1, len(modes), i + 1)
        plt.ylim(y_min - 0.5, y_max + 0.5)

    plt.savefig(os.path.join(path_plots, 'plot_risk_box_plots.png'))

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
        plt.title(f'Distribution of Perceived Risk Over Time {mode}')
        y_min = min(y_min, df.min().min())
        y_max = max(y_max, df.max().max())

    for i in range(len(modes)):
        plt.subplot(1, len(modes), i + 1)
        plt.ylim(y_min - 0.5, y_max + 0.5)

    plt.savefig(os.path.join(path_plots, 'risk_violin_plots.png'))


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


def plot_icu_comp(path_results, path_plots, modes, path_icu_data, log_scale=False, icu_capacity_val=9, plot_percentiles=True, kmax=0.34):
    create_folder_if_not_exists(path_plots)
    icu_comp = [7]
    labels = ["bfact 0.0", "bfact 0.36", "bfact 0.7",
              "bfact 0.0 (Max)", "bfact 0.36 (Max)", "bfact 0.7 (Max)"]
    km_formatted = f"{kmax:.2f}"
    title = "ICU Occupancy per 100_000 (kmax = " + \
        km_formatted + ")"

    plot_data = []
    for path in path_results:
        for mode in modes:
            path_results_mode = os.path.join(path, mode)
            plot_data.append(get_results(
                path_results_mode, icu_comp, results="total", percentiles=["p50"]))
    # calculate ICU occupancy per 100_000 inhabitants
    for data in plot_data:
        for key in data.keys():
            for indx in range(len(data[key])):
                data[key][indx] = data[key][indx] / total_pop * 100_000

    # get curve of the max ICU occupancy in each day
    for path in path_results:
        data_path = os.path.join(path, "FeedbackDamping", "p50", "Results.h5")
        max_vals = np.zeros(len(data['p50']))
        with h5py.File(data_path, 'r') as f:
            keys = list(f.keys())
            for i, key in enumerate(keys):
                group = f[key]
                total = group['Total'][()][:, icu_comp]
                # scale per 100_000 inhabitants
                population = pd.read_json(
                    'data/pydata/Germany/county_current_population.json')
                pop_key = population.loc[population['ID_County'] == int(
                    key)]['Population'].values[0]
                total = total / pop_key * 100_000
                # get the max value in each day
                for indx in range(len(total)):
                    max_vals[indx] = max(max_vals[indx], total[indx])
        plot_data.append(max_vals)

    # plot
    num_data = len(plot_data)
    create_folder_if_not_exists(path_plots)

    # Set days for x-axis
    num_days = len(plot_data[-1]) - 1
    start_date_datetime = datetime.strptime(start_date, "%Y-%m-%d")
    end_date_datetime = start_date_datetime + timedelta(days=num_days)
    days = pd.date_range(start_date_datetime, end_date_datetime)
    months = pd.date_range(start=start_date_datetime,
                           end=end_date_datetime, freq='MS')

    # Creating subplots based on the number of data series
    colors_plot = colors
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
        ax.set_ylabel('ICU OCC per 100k', fontsize=fontsize)
        ax.grid(True)
        if log_scale:
            ax.set_yscale('log')

    # Plotting the data
    for indx_data in range(num_data):
        y = plot_data[indx_data]
        ax = axes[indx_data % len(axes)]

        cl = colors_plot[indx_data]
        linest = '-'
        if indx_data > 2:
            cl = colors_plot[indx_data - 3]
            linest = '--'
        if isinstance(y, dict):
            ax.plot(days, y["p50"], label=labels[indx_data], linewidth=lineWidth,
                    linestyle=linest, color=cl)
        else:
            ax.plot(days, y, label=labels[indx_data], linewidth=lineWidth,
                    linestyle=linest, color=cl)

    for ax in axes:
        ax.legend(fontsize=legendsize, loc='center left',
                  bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig(os.path.join(path_plots, f"plot_{title}.png"))
    plt.clf()
    return 0


def plot_icu_comp_all_dirs(path_results, path_plots, modes, path_icu_data, log_scale=False, icu_capacity_val=9, plot_percentiles=True):
    create_folder_if_not_exists(path_plots)

    # list all files in path_results and filter for kmin and fixed
    dirs_in_results = [entry for entry in os.listdir(
        path_results) if entry.startswith("kmin") or entry.startswith("fixed")]
    plots_dir = os.path.join(path_plots, "ICU_Plots")
    create_folder_if_not_exists(plots_dir)

    for dir_name in dirs_in_results:

        icu_comp = [7]
        label = "ICU Occupancy"
        labels = []
        title = "ICU Occupancy per 100_000 inhabitants"

        plot_data = []
        for mode in modes:
            path_results_mode = os.path.join(path_results, dir_name, mode)
            # check if dir exists, othewrise continue
            if not os.path.exists(path_results_mode):
                continue
            labels.append(label + f" {mode}")
            plot_data.append(get_results(
                path_results_mode, icu_comp, results="total"))

        # calculate ICU occupancy per 100_000 inhabitants
        for data in plot_data:
            for key in data.keys():
                for indx in range(len(data[key])):
                    data[key][indx] = data[key][indx] / total_pop * 100_000

        # create dict with same shape and set constant value for ICU capacity
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

        plot(plot_data, labels, plots_dir, title=dir_name,
             log_scale=log_scale, ylabel="ICU Occupancy per 100_000", plot_percentiles=plot_percentiles)
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


def get_dirs_starting_with_x(x, path_results):
    def get_kmax_value(dir_name):
        return float(dir_name.split('_')[3])
    dirs = [d for d in os.listdir(path_results) if os.path.isdir(
        os.path.join(path_results, d)) and d.startswith(x)]
    return sorted(dirs, key=get_kmax_value)


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


def get_risk_total(path):
    data_risk = read_results_h5(path, 0, 'Group1')
    pop = get_pop()
    total_pop = pop['Population'].sum()
    # iterate over all rows from data_risk and scale the risk with the population
    risk_total = np.zeros(data_risk.shape[1])
    for i in range(data_risk.shape[0]):
        risk_total += data_risk[i, :] * \
            pop['Population'].values[i] / total_pop
    return risk_total


def plot_peak_values_splitted(path_results, path_plots, modes, target_indx, percentile="p50", flows=True, plot_type='kmin', title='Peak Value', vertical=False, dir_value='kmin'):
    if len(modes) > 1:
        print("Only one mode is allowed for the grid peak plot.")
        return

    plot_dir = os.path.join(path_plots, "peaks", title)
    create_folder_if_not_exists(plot_dir)

    allowed_kmax = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    kmin_dirs = get_dirs_starting_with_x(dir_value, path_results)
    all_kmins = 0
    if dir_value == 'kmin':
        all_kmins = sorted(set(float(d.split("_")[1]) for d in kmin_dirs))
    else:
        all_kmins = sorted(set(float(d.split("_")[-1]) for d in kmin_dirs))

    num_days = 0

    for index, kmin in enumerate(all_kmins):
        fig_index = 0
        if dir_value == 'fixed':
            if kmin > all_kmins[0]:
                continue

        num_subplots = round(11 - index)
        if dir_value == 'fixed':
            num_subplots = 11

        # Adjusting the width to 88mm (3.46 inches)
        subplot_width = 3.46  # ~88mm
        subplot_height = 3  # Keeping height fixed

        num_rows = 5
        num_cols = 2

        fig_width = num_cols * subplot_width
        fig_height = num_rows * subplot_height
        fig = plt.figure(figsize=(fig_width, fig_height))
        gs = GridSpec(num_rows, num_cols, figure=fig)
        tiks_size = 8
        all_peaks = []
        all_risks = []

        for run in kmin_dirs:
            if dir_value == 'kmin':
                if kmin != float(run.split("_")[1]):
                    continue
            kmax = float(run.split("_")[3])
            if kmax not in allowed_kmax:
                continue

            path = os.path.join(path_results, run,
                                modes[0], "flows" if flows else "")
            df = read_data(path, target_indx, percentile, flows)

            if num_days == 0:
                num_days = df.shape[0]

            # smooth the data from the df and find the peaks
            def get_peak(df):
                smoothed = df.rolling(window=7, min_periods=1).mean()
                # iterate over each column and find the peaks
                num_counties_no_peak = 0
                peaks = []
                for col in smoothed.columns:
                    peak_indices = find_peaks(
                        smoothed[col].values, prominence=1, distance=21)[0]
                    if len(peak_indices) == 0:
                        peak_indices = find_peaks(
                            smoothed[col].values, prominence=0, distance=21)[0]
                    peak_values = smoothed[col].values[peak_indices]
                    if len(peak_indices) == 0:
                        num_counties_no_peak += 1
                        continue

                    peaks.append((peak_indices, peak_values))
                if num_counties_no_peak > 0:
                    print("Number of counties without peak: ",
                          num_counties_no_peak)
                return peaks

            peaks = get_peak(df)

            if plot_type == 'kmin':
                peaks_each_day = np.zeros(num_days)
                for peak_indices, _ in peaks:
                    for peak_index in peak_indices:
                        peaks_each_day[peak_index] += 1
                all_peaks.append(peaks_each_day)
                row = fig_index % num_rows
                col = fig_index // num_rows
                ax = fig.add_subplot(gs[row, col])
                ax.plot(peaks_each_day, color=colors[0], linewidth=lineWidth)
                if row == 2 and col == 0:
                    ax.set_ylabel('Number of Peaks', fontsize=fontsize)

            else:
                peak_indices = df.idxmax().values
                for _, p_values in peaks:
                    for i in range(len(p_values)):
                        p_values[i] = p_values[i] / \
                            get_pop()['Population'].values[i] * 100_000

                peak_values = {i: [] for i in range(df.shape[0])}
                for idx, val in peaks:
                    for i, v in zip(idx, val):
                        peak_values[i].append(v)

                df_long = pd.DataFrame(
                    [{'Day': day, 'Peak Value': val} for day, vals in peak_values.items() for val in vals])
                all_peaks.extend(peak_values.values())
                all_days = pd.DataFrame({'Day': range(num_days)})
                df_long = pd.merge(all_days, df_long, on='Day', how='left')
                row = fig_index % num_rows
                col = fig_index // num_rows
                ax = fig.add_subplot(gs[row, col])
                sns.boxplot(x='Day', y='Peak Value', data=df_long, ax=ax)
                ax.set_xlabel('')
                ax.set_ylabel('')
                if row == 2 and col == 0:
                    ax.set_ylabel(title + ' per 100,000k', fontsize=fontsize)

                if row == 4:
                    ax.set_xlabel('Days', fontsize=fontsize)

            fig_index += 1

            if modes[0] == 'FeedbackDamping':
                risk = get_risk_total(os.path.join(
                    path_results, run, "FeedbackDamping", "risk", percentile, "Results.h5"))
                all_risks.append(risk)

            # ax.set_title(f'kmin: {float(run.split("_")[-1])}', fontsize=8)
            kmax = float(run.split("_")[3])
            ax.set_title(
                rf'$\psi_{{max}}$: {kmax:.1f}', fontsize=fontsize)

        global_min = 0
        global_max = 0
        if plot_type == 'kmin':
            global_min = min(min(peaks) for peaks in all_peaks if peaks.size)
            global_max = max(max(peaks)
                             for peaks in all_peaks if peaks.size) + 10
        else:
            global_min = min(min(peaks) for peaks in all_peaks if peaks)
            global_max = max(max(peaks)
                             for peaks in all_peaks if peaks) + 10

        count_ax = 0
        for ax in fig.get_axes():
            if count_ax > num_subplots:
                break
            ax.grid()
            ax.set_yscale('symlog')
            ax.set_ylim(global_min, global_max)
            ax.set_xticks(range(0, num_days, 50))
            ax.tick_params(axis='x', labelsize=tiks_size)
            ax.tick_params(axis='y', labelsize=tiks_size)
            if vertical:
                if ax == fig.get_axes()[4] or ax == fig.get_axes()[9]:
                    ax.set_xlabel('Days', fontsize=tiks_size)

            count_ax += 1

        plot_and_save(
            fig, plot_dir, f'splitted_peaks_grid_{plot_type}_{dir_value}_{kmin}.pdf')

    # if len(modes) > 1:
    #     print("Only one mode is allowed for the grid peak plot.")
    #     return

    # plot_dir = os.path.join(path_plots, "peaks", title)
    # create_folder_if_not_exists(plot_dir)

    # allowed_kmax = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    # kmin_dirs = get_dirs_starting_with_x(dir_value, path_results)
    # all_kmins = 0
    # if dir_value == 'kmin':
    #     all_kmins = sorted(set(float(d.split("_")[1]) for d in kmin_dirs))
    # else:
    #     all_kmins = sorted(set(float(d.split("_")[-1]) for d in kmin_dirs))

    # num_days = 0

    # # Set the grid for 2 rows and 5 columns, handling the case where there are 10 plots.
    # num_rows = 5
    # num_cols = 2
    # subplot_width = 8
    # subplot_height = 3

    # fig_width = num_cols * subplot_width
    # fig_height = num_rows * subplot_height
    # fig = plt.figure(figsize=(fig_width, fig_height))
    # gs = GridSpec(num_rows, num_cols, figure=fig)

    # all_peaks = []
    # all_risks = []

    # fig_index = 0

    # for index, kmin in enumerate(all_kmins):
    #     if dir_value == 'fixed' and kmin > all_kmins[0]:
    #         continue

    #     for run in kmin_dirs:
    #         if dir_value == 'kmin' and kmin != float(run.split("_")[1]):
    #             continue
    #         kmax = float(run.split("_")[3])
    #         if kmax not in allowed_kmax:
    #             continue

    #         path = os.path.join(path_results, run,
    #                             modes[0], "flows" if flows else "")
    #         df = read_data(path, target_indx, percentile, flows)

    #         if num_days == 0:
    #             num_days = df.shape[0]

    #         # smooth the data from the df and find the peaks
    #         def get_peak(df):
    #             smoothed = df.rolling(window=7, min_periods=1).mean()
    #             # iterate over each column and find the peaks
    #             num_counties_no_peak = 0
    #             peaks = []
    #             for col in smoothed.columns:
    #                 peak_indices = find_peaks(
    #                     smoothed[col].values, prominence=1, distance=21)[0]
    #                 if len(peak_indices) == 0:
    #                     peak_indices = find_peaks(
    #                         smoothed[col].values, prominence=0, distance=21)[0]
    #                 peak_values = smoothed[col].values[peak_indices]
    #                 if len(peak_indices) == 0:
    #                     num_counties_no_peak += 1
    #                     continue
    #                 peaks.append((peak_indices, peak_values))
    #             if num_counties_no_peak > 0:
    #                 print("Number of counties without peak: ",
    #                       num_counties_no_peak)
    #             return peaks

    #         peaks = get_peak(df)

    #         if plot_type == 'kmin':
    #             # Calculate peaks per day
    #             peaks_each_day = np.zeros(num_days)
    #             for peak_indices, _ in peaks:
    #                 for peak_index in peak_indices:
    #                     peaks_each_day[peak_index] += 1
    #             all_peaks.append(peaks_each_day)

    #             # Plot in the grid
    #             row = fig_index % num_rows
    #             col = fig_index // num_rows
    #             ax = fig.add_subplot(gs[row, col])
    #             ax.plot(peaks_each_day, color=colors[0], linewidth=lineWidth)
    #             ax.set_title(
    #                 rf'$\psi_{{min}}$: {kmin}, $\psi_{{max}}$: {kmax:.1f}', fontsize=fontsize)

    #             if row == 2 and col == 0:
    #                 ax.set_ylabel('Number of Peaks', fontsize=fontsize)
    #         else:
    #             # Process peak values for each day
    #             for _, p_values in peaks:
    #                 for i in range(len(p_values)):
    #                     p_values[i] = p_values[i] / \
    #                         get_pop()['Population'].values[i] * 100_000

    #             # create an empty dict with the peak values for each day
    #             peak_values = {i: [] for i in range(df.shape[0])}
    #             for idx, val in peaks:
    #                 for i, v in zip(idx, val):
    #                     peak_values[i].append(v)

    #             # convert the dict to a format for sns.boxplot
    #             df_long = pd.DataFrame([{'Day': day, 'Peak Value': val}
    #                                    for day, vals in peak_values.items() for val in vals])
    #             all_days = pd.DataFrame({'Day': range(num_days)})
    #             df_long = pd.merge(all_days, df_long, on='Day', how='left')

    #             # Plot in the grid
    #             row = fig_index % num_rows
    #             col = fig_index // num_rows
    #             ax = fig.add_subplot(gs[row, col])
    #             sns.boxplot(x='Day', y='Peak Value', data=df_long, ax=ax)
    #             ax.set_ylabel(title + ' per 100,000k')

    #             if row == 2 and col == 0:
    #                 ax.set_ylabel('Peak Value per 100,000k', fontsize=fontsize)

    #         fig_index += 1

    # # Set the same axis limits for all plots
    # global_min = min(min(peaks) for peaks in all_peaks if len(peaks) > 0)
    # global_max = max(max(peaks) for peaks in all_peaks if len(peaks) > 0) + 10

    # for ax in fig.get_axes():
    #     ax.grid()
    #     ax.set_yscale('symlog')
    #     ax.set_ylim(global_min, global_max)
    #     ax.set_xticks(range(0, num_days, 50))
    #     ax.tick_params(axis='x', labelsize=ticks)
    #     ax.tick_params(axis='y', labelsize=ticks)

    # # Save the figure
    # plot_and_save(
    #     fig, plot_dir, f'splitted_peaks_grid_{plot_type}_{dir_value}_{kmin}.pdf')


def plot_peak_values(path_results, path_plots, modes, target_indx, percentile="p50", flows=True, plot_type='kmin', title='Peak Value', vertical=False, dir_value='kmin'):
    if len(modes) > 1:
        print("Only one mode is allowed for the grid peak plot.")
        return

    plot_dir = os.path.join(path_plots, "peaks", title)
    create_folder_if_not_exists(plot_dir)

    allowed_kmax = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    kmin_dirs = get_dirs_starting_with_x(dir_value, path_results)
    all_kmins = 0
    if dir_value == 'kmin':
        all_kmins = sorted(set(float(d.split("_")[1]) for d in kmin_dirs))
    else:
        all_kmins = sorted(set(float(d.split("_")[-1]) for d in kmin_dirs))

    num_days = 0

    for index, kmin in enumerate(all_kmins):
        if dir_value == 'fixed':
            if kmin > all_kmins[0]:
                continue

        num_subplots = round(11 - index)
        if dir_value == 'fixed':
            num_subplots = 11
        subplot_width = 8  # Breite eines einzelnen Subplots
        subplot_height = 3  # Höhe eines einzelnen Subplots

        if vertical:
            fig_height = num_subplots * subplot_height
            fig_width = subplot_width
            fig = plt.figure(figsize=(fig_width, fig_height))
            gs = GridSpec(num_subplots, 1, figure=fig)
        else:
            fig_width = num_subplots * subplot_width
            fig_height = subplot_height
            fig = plt.figure(figsize=(fig_width, fig_height))
            gs = GridSpec(1, num_subplots, figure=fig)

        all_peaks = []
        all_risks = []

        for run in kmin_dirs:
            if dir_value == 'kmin':
                if kmin != float(run.split("_")[1]):
                    continue
            kmax = float(run.split("_")[3])
            if kmax not in allowed_kmax:
                continue

            path = os.path.join(path_results, run,
                                modes[0], "flows" if flows else "")
            df = read_data(path, target_indx, percentile, flows)

            if num_days == 0:
                num_days = df.shape[0]

            # smooth the data from the df and find the peaks
            def get_peak(df):
                smoothed = df.rolling(window=7, min_periods=1).mean()
                # iterate over each column and find the peaks
                num_counties_no_peak = 0
                peaks = []
                for col in smoothed.columns:
                    peak_indices = find_peaks(
                        smoothed[col].values, prominence=1, distance=21)[0]
                    if len(peak_indices) == 0:
                        peak_indices = find_peaks(
                            smoothed[col].values, prominence=0, distance=21)[0]
                    peak_values = smoothed[col].values[peak_indices]
                    if len(peak_indices) == 0:
                        num_counties_no_peak += 1
                        continue

                    peaks.append((peak_indices, peak_values))
                if num_counties_no_peak > 0:
                    print("Number of counties without peak: ",
                          num_counties_no_peak)
                return peaks

            # peaks contains the peak indices and the peak values for each county
            peaks = get_peak(df)

            if plot_type == 'kmin':
                # peaks = np.bincount(df.idxmax().values, minlength=df.shape[0])
                # calculate the peaks per day out of peaks
                peaks_each_day = np.zeros(num_days)
                for peak_indices, _ in peaks:
                    for peak_index in peak_indices:
                        peaks_each_day[peak_index] += 1
                all_peaks.append(peaks_each_day)
                if vertical:
                    ax = fig.add_subplot(gs[len(all_peaks) - 1, 0])
                else:
                    ax = fig.add_subplot(gs[0, len(all_peaks) - 1])
                ax.plot(peaks_each_day, color=colors[0], linewidth=lineWidth)
                if len(all_peaks) == 1:
                    ax.set_ylabel('Number of Peaks')
            else:
                peak_indices = df.idxmax().values
                # iterate over peaks and scale the peak_value
                for _, p_values in peaks:
                    for i in range(len(p_values)):
                        p_values[i] = p_values[i] / \
                            get_pop()['Population'].values[i] * 100_000

                # create a  empty dict with the peak values for each day
                peak_values = {i: [] for i in range(df.shape[0])}
                for idx, val in peaks:
                    for i, v in zip(idx, val):
                        peak_values[i].append(v)

                # convert the dict to a dict to use sns.boxplot
                df_long = pd.DataFrame(
                    [{'Day': day, 'Peak Value': val} for day, vals in peak_values.items() for val in vals])
                all_peaks.extend(peak_values.values())
                all_days = pd.DataFrame({'Day': range(num_days)})
                df_long = pd.merge(all_days, df_long, on='Day', how='left')
                if vertical:
                    ax = fig.add_subplot(gs[len(all_peaks) // len(df) - 1, 0])
                else:
                    ax = fig.add_subplot(gs[0, len(all_peaks) // len(df) - 1])
                sns.boxplot(x='Day', y='Peak Value', data=df_long, ax=ax)
                ax.set_ylabel(title + ' per 100,000k')

            if modes[0] == 'FeedbackDamping':
                risk = get_risk_total(os.path.join(
                    path_results, run, "FeedbackDamping", "risk", percentile, "Results.h5"))
                all_risks.append(risk)

            kmax = float(run.split("_")[3])
            ax.set_title(
                rf'$\psi_{{min}}$: {kmin}, $\psi_{{max}}$: {kmax:.1f}', fontsize=fontsize)

            if dir_value == 'fixed':
                ax.set_title(
                    f'kmin: { float(run.split("_")[-1])}', fontsize=fontsize)
        global_min = 0
        global_max = 0
        if plot_type == 'kmin':
            global_min = min(min(peaks) for peaks in all_peaks if peaks.size)
            global_max = max(max(peaks)
                             for peaks in all_peaks if peaks.size) + 10
        else:
            global_min = min(min(peaks) for peaks in all_peaks if peaks)
            global_max = max(max(peaks)
                             for peaks in all_peaks if peaks) + 10

        count_ax = 0
        for ax in fig.get_axes():
            if count_ax > num_subplots:
                break
            # grid
            ax.grid()
            # log scale
            ax.set_yscale('symlog')
            ax.set_ylim(global_min, global_max)
            ax.set_xticks(range(0, num_days, 50))
            ax.tick_params(axis='x', labelsize=ticks)
            ax.tick_params(axis='y', labelsize=ticks)
            if not vertical:
                ax.set_xlabel('')
                if ax != fig.get_axes()[0]:
                    ax.set_yticklabels([])
                    ax.set_ylabel('')
                if ax == fig.get_axes()[len(fig.get_axes()) // 2]:
                    ax.set_xlabel('Days', fontsize=fontsize)
            else:
                ax.set_xlabel('')
                ax.set_ylabel('')
                if ax == fig.get_axes()[len(fig.get_axes()) // 2]:
                    ax.set_ylabel(title + ' per 100,000k', fontsize=fontsize)
                if ax == fig.get_axes()[len(fig.get_axes()) - 1]:
                    ax.set_xlabel('Days', fontsize=fontsize)

            count_ax += 1

        plot_and_save(
            fig, plot_dir, f'peaks_grid_{plot_type}_{dir_value}_{kmin}.pdf')


def get_peak_data(path_results, blending_fact, modes, target_indx, percentile="p50", flows=True, title='Peak Value', dir_type='peak', filter=None):
    if len(modes) > 1:
        print("Only one mode is allowed for the grid peak plot.")
        return

    population = get_pop()

    data = {'b_fact': [], 'kmax': [], 'day': [], 'peak_value': []}
    for b_fact in blending_fact:
        path_results_blended = os.path.join(
            path_results, f"BlendingFactorRegional_{b_fact:.2f}" + "0000")

        kmin_dirs = get_dirs_starting_with_x('kmin', path_results_blended)
        all_kmax = sorted(set(float(d.split("_")[-1]) for d in kmin_dirs))

        num_days = 0

        for index, kmax in enumerate(all_kmax):
            path = os.path.join(path_results_blended, kmin_dirs[index],
                                modes[0], "flows" if flows else "")
            df = read_data(path, target_indx, percentile, flows)

            # columns as integers
            df.columns = df.columns.astype(int)
            # sort the columns
            df = df.reindex(sorted(df.columns), axis=1)

            # apply filter if not None
            if filter is not None:
                # Flatten the filter list to match the column names
                filter_flat = [str(item[0]) for item in filter]
                df = df[filter_flat]

            if num_days == 0:
                num_days = df.shape[0]

            # smooth the data from the df and find the peaks
            def get_peak(df, population=population):
                smoothed = df.rolling(window=7, min_periods=1).mean()
                # iterate over each column and find the peaks
                num_counties_no_peak = 0
                peaks = []
                for col in smoothed.columns:
                    peak_indices = find_peaks(
                        smoothed[col].values, prominence=1, distance=21)[0]
                    if len(peak_indices) == 0:
                        peak_indices = find_peaks(
                            smoothed[col].values, prominence=0, distance=21)[0]
                    peak_values = smoothed[col].values[peak_indices]
                    if len(peak_indices) == 0:
                        num_counties_no_peak += 1
                        continue

                    # scale the peak values per 100_000 inhabitants
                    # search which index from population['ID_County] is equal to col
                    idx = population[population['ID_County'] == col].index[0]
                    pop = population['Population'].values[idx]
                    peak_values = peak_values / pop * 100_000

                    peaks.append((peak_indices, peak_values))
                if num_counties_no_peak > 0:
                    print("Number of counties without peak: ",
                          num_counties_no_peak)
                return peaks

            # peaks contains the peak indices and the peak values for each county
            if dir_type == 'peak':
                peaks = get_peak(df)

                peak_values = {i: [] for i in range(df.shape[0])}
                for idx, val in peaks:
                    for i, v in zip(idx, val):
                        peak_values[i].append(v)

                # delete all empty entries
                peak_values = {k: v for k,
                               v in peak_values.items() if len(v) > 0}

            elif dir_type == 'max_val':
                # get max value and index for each column
                max_vals_df = df.max()
                max_indices_df = df.idxmax()

                # scale max_vals_df per 100_000
                max_vals_df = max_vals_df / \
                    population['Population'].values * 100_000

                peak_values = {i: [] for i in range(df.shape[0])}
                for i in range(len(max_vals_df)):
                    peak_values[max_indices_df.values[i]].append(
                        max_vals_df.values[i])

                # delete all empty entries
                peak_values = {k: v for k,
                               v in peak_values.items() if len(v) > 0}

            for day, values in peak_values.items():
                data['b_fact'].append(b_fact)
                data['kmax'].append(kmax)
                data['day'].append(day)
                data['peak_value'].append(values)

    # Convert dictionary to DataFrame
    df_peaks = pd.DataFrame(data)

    return df_peaks


def plot_peak_time_kmax(path_results, path_plots, blending_fact, modes, target_indx, percentile="p50", flows=True, title='Peak Value', dir_type='kmin', log_scale=True):
    data = get_peak_data(path_results, blending_fact,
                         modes, target_indx, percentile, flows, title, dir_type)

    for b_fact in blending_fact:
        peak_data = data[data['b_fact'] == b_fact]
        # Delete all entries with peak_value == 0
        # peak_data = peak_data[peak_data['peak_value'] > 0].reset_index(drop=True)
        data_variance = {'b_fact': [], 'kmax': [], 'mean': [], 'variance': []}

        # Calculate the mean peak values for each (kmax, day) combination
        mean_peaks = peak_data.groupby(['kmax', 'day'])[
            'peak_value'].mean().reset_index()

        # delete all rows with peak_value == 0
        mean_peaks = mean_peaks[mean_peaks['peak_value'] > 0].reset_index(
            drop=True)
        all_kmax = mean_peaks['kmax'].unique()
        # create df with columns kmax
        for kmax in all_kmax:
            peak_data_kmax = mean_peaks[mean_peaks['kmax'] == kmax]
            peak_days = peak_data_kmax['day'].values
            # calculate the mean and the variance for peak_days
            mean_peak = peak_days.mean()
            var_peak = peak_days.var()
            data_variance['b_fact'].append(b_fact)
            data_variance['kmax'].append(kmax)
            data_variance['mean'].append(mean_peak)
            data_variance['variance'].append(var_peak)

    statistical_measures = ['mean', 'variance']
    for measure in statistical_measures:
        results_df = pd.DataFrame(data_variance)
        heatmap_data = results_df.pivot_table(
            values=measure, index='kmax', columns='b_fact', aggfunc=lambda x: x)

        # Create the heatmap
        plt.figure(figsize=(20, 10))
        if log_scale:
            sns.heatmap(heatmap_data, cmap="viridis",
                        norm=SymLogNorm(linthresh=1), fmt=".2f")
        else:
            sns.heatmap(heatmap_data, cmap="viridis", fmt=".2f")
        plt.xlabel('kmax', fontsize=fontsize)
        plt.ylabel('Day', fontsize=fontsize)
        plt.xticks(fontsize=ticks)
        tick_positions = np.arange(0,  len(heatmap_data.index), 20)
        tick_labels = [heatmap_data.index[i] for i in tick_positions]
        plt.yticks(tick_positions + 0.5, tick_labels, rotation=0,
                   fontsize=ticks)
        # plt.ylim(0, len(heatmap_data.index))
        fn = 'peak_time_heatmap_' + dir_type + '_b_fact_' + str(b_fact)
        if log_scale:
            fn += '_log'
        plt.savefig(os.path.join(path_plots, fn + '.png'))
        plt.clf()
        plt.close()


def plot_peak_integral(path_results, path_plots, blending_fact, modes, target_indx, rhos, percentile="p50", flows=True, title='Peak Value', log_scale=True):

    # filter_counties = pd.read_csv(
    #     '/localdata1/code_2024/memilio/tools/border_counties.txt')
    # filter_counties = filter_counties.values.tolist()
    filter_counties = None
    population = get_pop()

    for rho in rhos:
        path_rho = os.path.join(
            path_results, f"rho_{rho:.6f}")
        data = {'b_fact': [], 'kmax': [], 'day': [], 'mean': [],
                'std': [], 'mean_global': [], 'std_global': []}
        for b_fact in blending_fact:
            path_results_blended = os.path.join(
                path_rho, f"BlendingFactorRegional_{b_fact:.2f}" + "0000")

            kmin_dirs = get_dirs_starting_with_x('kmin', path_results_blended)
            all_kmax = sorted(set(float(d.split("_")[-1]) for d in kmin_dirs))

            num_days = 0

            for index, kmax in enumerate(all_kmax):
                path = os.path.join(path_results_blended, kmin_dirs[index],
                                    modes[0], "flows" if flows else "")
                df = read_data(path, target_indx, percentile, flows)

                # apply filter if not None
                if filter_counties is not None:
                    # Flatten the filter list to match the column names
                    filter_flat = [str(item[0]) for item in filter]
                    df = df[filter_flat]

                if num_days == 0:
                    num_days = df.shape[0]

                # get mean for each col
                mean_icu = df.mean(axis=0)

                # because each county has a different population, we need to normalize the peak values per 100,000 individuals
                mean_icu = mean_icu / population['Population'].values * 100_000
                # calculate the mean and std of mean_icu
                mean = mean_icu.mean()
                std = mean_icu.std()

                # global values
                sum_icu = df.sum(axis=1)

                data['b_fact'].append(b_fact)
                data['kmax'].append(kmax)
                data['day'].append(num_days)
                data['mean'].append(mean)
                data['std'].append(std)
                data['mean_global'].append(sum_icu.mean())
                data['std_global'].append(sum_icu.std())

        # Plot and save each heatmap separately
        size_ticks = 18

        xlabel_title = 'Impact Regional Proximity'  # 'Impact Federal State'  #
        statistical_measures = ['mean', 'std', 'mean_global', 'std_global']
        for measure in statistical_measures:
            # Plot for 'kmin'
            fig_kmin, ax_kmin = plt.subplots(
                figsize=(10, 8))  # Increase figure size
            results_df_icu = pd.DataFrame(data)
            heatmap_data_kmin = results_df_icu.pivot_table(
                values=measure, index='kmax', columns='b_fact', aggfunc=lambda x: x)
            # invert the y-axis
            heatmap_data_kmin = heatmap_data_kmin.iloc[::-1]
            sns.heatmap(heatmap_data_kmin, cmap='viridis', ax=ax_kmin, norm=SymLogNorm(
                linthresh=1) if log_scale else None, cbar=True)

            ax_kmin.set_xlabel(xlabel_title, fontsize=fontsize)
            ax_kmin.set_ylabel(rf'$\psi_{{max}}$', fontsize=fontsize)
            ax_kmin.tick_params(axis='x', labelsize=size_ticks)
            ax_kmin.tick_params(axis='y', labelsize=size_ticks)

            # increase ticks size colorbar
            cbar = ax_kmin.collections[0].colorbar
            cbar.ax.tick_params(labelsize=size_ticks)

            # To specify the number of ticks on both or any single axes
            allowed_labels = [0.0, 0.1, 0.2, 0.3,
                              0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            for n, label in enumerate(ax_kmin.xaxis.get_ticklabels()):
                if float(label.get_text()) not in allowed_labels:
                    label.set_visible(False)

            for n, label in enumerate(ax_kmin.yaxis.get_ticklabels()):
                if float(label.get_text()) not in allowed_labels:
                    label.set_visible(False)

            plt.yticks(rotation=0)  # Keep the labels horizontal
            plt.xticks(rotation=0)

            plt.savefig(os.path.join(
                path_plots, f'icu_{measure}_heatmap_kmin_{rho}.pdf'), dpi=600)
            plt.clf()
            plt.close(fig_kmin)

        return 0


def plot_peak_time_kmax_variance(path_results, path_plots, blending_fact, modes, target_indx, rhos, percentile="p50", flows=True, title='Peak Value', log_scale=True):

    # filter_counties = pd.read_csv(
    #     '/localdata1/code_2024/memilio/tools/border_counties.txt')
    # filter_counties = filter_counties.values.tolist()
    filter_counties = None

    for rho in rhos:
        path_rho = os.path.join(
            path_results, f"rho_{rho:.6f}")

        data = []
        if title == 'Peak Value':
            data = get_peak_data(path_rho, blending_fact,
                                 modes, target_indx, percentile, flows, title, 'peak', filter=filter_counties)
        elif title == 'Max Value':
            data = get_peak_data(path_rho, blending_fact,
                                 modes, target_indx, percentile, flows, title, 'max_val', filter=filter_counties)
        else:
            print('Title not supported')
            return

        # Process data for 'kmin'
        all_b_fact_kmin = data['b_fact'].unique()
        all_kmax_kmin = data['kmax'].unique()
        data_variance_kmin = {'b_fact': [], 'kmax': [],
                              'mean': [], 'std': [], 'variance': []}

        for b_fact in all_b_fact_kmin:
            for kmax in all_kmax_kmin:
                peak_data = data[(data['b_fact'] == b_fact) & (
                    data['kmax'] == kmax)]
                # i want a list of all peak_days. So  just the col 'day'. However, if there are n peak_values, the day should be n times in the list
                peak_day_list = []
                for i in range(len(peak_data)):
                    entry = peak_data.iloc[i]
                    for _ in range(len(entry['peak_value'])):
                        peak_day_list.append(entry['day'])
                # to np array
                peak_day_list = np.array(peak_day_list)
                mean_peak = peak_day_list.mean()
                std_peak = peak_day_list.std()
                var_peak = peak_day_list.var()

                data_variance_kmin['b_fact'].append(b_fact)
                data_variance_kmin['kmax'].append(kmax)
                data_variance_kmin['mean'].append(mean_peak)
                data_variance_kmin['std'].append(std_peak)
                data_variance_kmin['variance'].append(var_peak)

        # Process data for 'val'
        # TODO: Noch prüfen, obs jetzt korrekt ist!!
        data_variance_val = {'b_fact': [], 'kmax': [],
                             'mean': [], 'std': [], 'variance': []}

        for b_fact in all_b_fact_kmin:
            for kmax in all_kmax_kmin:
                peak_data = data[(data['b_fact'] == b_fact)
                                 & (data['kmax'] == kmax)]
                peak_data = peak_data[peak_data['peak_value'].apply(
                    lambda x: len(x) > 0)].reset_index(drop=True)
                peak_data = peak_data.explode('peak_value')
                peak_data['peak_value'] = peak_data['peak_value'].astype(float)

                mean_peak = peak_data['peak_value'].values.mean()
                std_peak = peak_data['peak_value'].values.std()
                var_peak = peak_data['peak_value'].values.var()

                data_variance_val['b_fact'].append(b_fact)
                data_variance_val['kmax'].append(kmax)
                data_variance_val['mean'].append(mean_peak)
                data_variance_val['std'].append(var_peak)
                data_variance_val['variance'].append(var_peak)

        # Plot and save each heatmap separately
        size_ticks = 18

        xlabel_title = 'Impact Federal State'  # 'Impact Regional Proximity'  #
        statistical_measures = ['mean', 'std', 'variance']
        for measure in statistical_measures:
            # Plot for 'kmin'
            fig_kmin, ax_kmin = plt.subplots(
                figsize=(10, 8))  # Increase figure size
            results_df_kmin = pd.DataFrame(data_variance_kmin)
            heatmap_data_kmin = results_df_kmin.pivot_table(
                values=measure, index='kmax', columns='b_fact', aggfunc=lambda x: x)
            # invert the y-axis
            heatmap_data_kmin = heatmap_data_kmin.iloc[::-1]
            sns.heatmap(heatmap_data_kmin, cmap='viridis', ax=ax_kmin, norm=SymLogNorm(
                linthresh=1) if log_scale else None, cbar=True)
            ax_kmin.set_xlabel(xlabel_title, fontsize=fontsize)
            ax_kmin.set_ylabel(rf'$\psi_{{max}}$', fontsize=fontsize)
            ax_kmin.tick_params(axis='x', labelsize=size_ticks)
            ax_kmin.tick_params(axis='y', labelsize=size_ticks)

            # increase ticks size colorbar
            cbar = ax_kmin.collections[0].colorbar
            cbar.ax.tick_params(labelsize=size_ticks)

            # To specify the number of ticks on both or any single axes
            allowed_labels = [0.0, 0.1, 0.2, 0.3,
                              0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            for n, label in enumerate(ax_kmin.xaxis.get_ticklabels()):
                if float(label.get_text()) not in allowed_labels:
                    label.set_visible(False)

            for n, label in enumerate(ax_kmin.yaxis.get_ticklabels()):
                if float(label.get_text()) not in allowed_labels:
                    label.set_visible(False)

            # Set specific y-ticks
            # ax_kmin.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])

            plt.yticks(rotation=0)  # Keep the labels horizontal
            plt.xticks(rotation=0)

            plt.savefig(os.path.join(
                path_plots, f'{title}_{measure}_heatmap_kmin_{rho}.pdf'), dpi=600)
            plt.clf()
            plt.close(fig_kmin)

            # Plot for 'val'
            fig_val, ax_val = plt.subplots(
                figsize=(10, 8))  # Increase figure size
            results_df_val = pd.DataFrame(data_variance_val)
            heatmap_data_val = results_df_val.pivot_table(
                values=measure, index='kmax', columns='b_fact', aggfunc=lambda x: x)
            # invert the y-axis
            heatmap_data_val = heatmap_data_val.iloc[::-1]
            sns.heatmap(heatmap_data_val, cmap='plasma', ax=ax_val, norm=SymLogNorm(
                linthresh=1) if log_scale else None, cbar=True)
            ax_val.set_xlabel(xlabel_title, fontsize=fontsize)
            ax_val.set_ylabel(rf'$\psi_{{max}}$', fontsize=fontsize)
            ax_val.tick_params(axis='x', labelsize=size_ticks)
            ax_val.tick_params(axis='y', labelsize=size_ticks)

            # increase ticks size colorbar
            cbar = ax_val.collections[0].colorbar
            cbar.ax.tick_params(labelsize=size_ticks)

            allowed_labels = [0.0, 0.1, 0.2, 0.3,
                              0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            for n, label in enumerate(ax_val.xaxis.get_ticklabels()):
                if float(label.get_text()) not in allowed_labels:
                    label.set_visible(False)

            for n, label in enumerate(ax_val.yaxis.get_ticklabels()):
                if float(label.get_text()) not in allowed_labels:
                    label.set_visible(False)

            plt.yticks(rotation=0)  # Keep the labels horizontal
            plt.xticks(rotation=0)

            plt.savefig(os.path.join(
                path_plots, f'{title}_{measure}_heatmap_val_{rho}.pdf'), dpi=600)
            plt.clf()
            plt.close(fig_val)

    return 0


def plot_peaks_all_blending_factors(path_results, path_plots, blending_fact):
    for b_fact in blending_fact:
        b_fact_add = str(b_fact) + "00000"
        if len(b_fact_add) == 9:
            b_fact_add = b_fact_add[:8]
        path_blended = os.path.join(
            path_results, "BlendingFactorRegional_" + b_fact_add)
        path_plot_b_fact = os.path.join(
            path_plots, "BlendingFactorRegional_" + b_fact_add)
        create_folder_if_not_exists(path_plot_b_fact)
        plot_peak_values_splitted(path_blended, path_plot_b_fact, [
            "FeedbackDamping"], icu_compartment, plot_type='kmin', title='ICU Occupancy', flows=False, vertical=True)
        plot_peak_values_splitted(path_blended, path_plot_b_fact, [
            "FeedbackDamping"], icu_compartment, plot_type='val', title='ICU Occupancy', flows=False, vertical=True)
        plot_peak_values_splitted(path_blended, path_plot_b_fact, [
            "FeedbackDamping"], flow_se, plot_type='kmin', title='Daily Infections', vertical=True)
        plot_peak_values_splitted(path_blended, path_plot_b_fact, [
            "FeedbackDamping"], flow_se, plot_type='val', title='Daily Infections', flows=True, vertical=True)


def get_all_blending_factors(path_results):
    blending_factors = []
    for entry in os.listdir(path_results):
        if entry.startswith("BlendingFactorRegional"):
            blending_factors.append(float(entry.split("_")[-1]))
    # return sorted list
    blending_factors.sort()
    return blending_factors


def get_all_kmaxs(path_results):
    kmaxs = []
    for entry in os.listdir(path_results):
        # entry is of type 'kmin_0.000000_kmax_0.040000'. how to get the last value?
        if entry.startswith("kmin"):
            kmaxs.append(float(entry.split("_")[-1]))
    # return sorted list
    kmaxs.sort()
    return kmaxs


def violin_plot_all_kmax(path_results, path_plots, modes, blending_fact, kmaxs, target_indx, percentile="p50", flows=False, title='Peak Value', dir_type='kmin', log_scale=True, rho=1.0):
    # Retrieve peak data using the get_peak_data function
    data = get_peak_data(path_results, blending_fact,
                         modes, target_indx, percentile, flows, title, dir_type)

    # Filter data for specific blending factors and kmax values
    data = data[data['b_fact'].isin(blending_fact)]
    data = data[data['kmax'].isin(kmaxs)]

    for b_fact in blending_fact:
        peak_data = data[data['b_fact'] == b_fact]

        # Reduce the dimension of peak_value
        peak_data = peak_data[peak_data['peak_value']
                              > 0.0].reset_index(drop=True)
        peak_data = peak_data.loc[peak_data.index.repeat(
            peak_data['peak_value'])].copy().reset_index(drop=True)
        peak_data.drop(columns=['peak_value'], inplace=True)

        # gehe durch alle Einträge und gucke ob ein eintrag mit day kleiner 10 existiert, derren peak_value > 0 ist. Setze diesen auf 0
        # for i in range(len(peak_data)):
        #     if peak_data.loc[i, 'day'] < 10:
        #         peak_data.loc[i, 'day'] = 0

        plt.figure(figsize=(10, 6))
        sns.violinplot(x='kmax', y='day', data=peak_data,
                       inner="box", linewidth=2, palette='Blues')

        # Increase tick size
        plt.xticks(fontsize=24)
        plt.yticks(fontsize=24)

        # Set labels and title
        plt.xlabel(r'$\psi_{\max}$', fontsize=24)  # Using raw string for latex
        plt.ylabel('Day', fontsize=24)
        # title = 'Impact Regional Proximity'  # 'Impact Federal State'  #
        # plt.title(title + f': {b_fact}', fontsize=18)

        # grid
        plt.grid(axis='y', linestyle='--', alpha=0.6)

        # Adjust layout to prevent cutting off
        plt.tight_layout()

        # set ylim between 50 und 190
        plt.ylim(0, 230)
        # y axis log scale
        # plt.yscale('log')

        # Save the plot as a PDF with tight bounding box
        pdf_filename = os.path.join(
            path_plots, f'violin_plot_{b_fact}_rho_{rho}.pdf')
        plt.savefig(pdf_filename, format='pdf', bbox_inches='tight')

        # Clear the plot after saving to avoid overlap when looping
        plt.close()


def plot_kmax_diff_rhos(path_results, path_plots, modes, blending_fact, kmaxs, target_indx, percentile="p50", flows=False, title='Peak Value', dir_type='kmin', log_scale=True):
    # list all dirs in path_results
    rho_dirs = [d for d in os.listdir(path_results) if os.path.isdir(
        os.path.join(path_results, d)) and d.startswith("rho")]

    # dirs sind in der Form "rho_0.60000"
    rhos = [float(r.split("_")[1]) for r in rho_dirs]

    kmax_peak = []
    indx_dir = 0
    for r in rhos:
        path_rho = os.path.join(path_results, rho_dirs[indx_dir])
        data = get_peak_data(path_rho, blending_fact,
                             modes, target_indx, percentile, flows, title, dir_type)
        # delete all rows where b_fact is not in blending_fact
        data = data[data['b_fact'].isin(blending_fact)]
        # same for kmax
        data = data[data['kmax'].isin(kmaxs)]

        for b_fact in blending_fact:
            peak_data = data[data['b_fact'] == b_fact]
            # Reduce dimension peak_value
            # i) delete all rows with peak_value = 0
            peak_data = peak_data[peak_data['peak_value']
                                  > 0.0].reset_index(drop=True)
            # get the mean peak value for each kmax
            all_kmax = peak_data['kmax'].unique()
            median_peaks = []
            for kmax in all_kmax:
                peak_data_kmax = peak_data[peak_data['kmax'] == kmax]
                expanded_days = np.repeat(
                    peak_data_kmax['day'], peak_data_kmax['peak_value']).values
                # Calculate the median of the expanded days
                median_peak_day = np.median(expanded_days)
                mean_peak_day = np.mean(expanded_days)
                median_peaks.append(mean_peak_day)
            indx_max = np.argmax(median_peaks)
            kmax_peak.append(all_kmax[indx_max])
        indx_dir += 1

    # Plot the kmax values
    plt.figure(figsize=(10, 6))
    plt.bar(rhos, kmax_peak, width=0.3)
    plt.xlabel('rho', fontsize=12)
    plt.ylabel('kmax', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.savefig(os.path.join(path_plots, 'kmax_diff_rho_bar.png'))
    plt.show()


if __name__ == '__main__':
    path_cwd = os.getcwd()
    icu_cap = [6, 9, 12, 15]
    cap_indx = 1
    path_results = os.path.join(
        path_cwd, "results/", "ICUCap_" + str(icu_cap[cap_indx]) + ".000000/")
    # path_results = os.path.join("/localdata2/feedback_paper/rho1_0/")

    #
    path_plots = os.path.join(
        path_cwd, "plots", "ICUCap_" + str(icu_cap[cap_indx]) + ".000000")
    path_icu_data = os.path.join(
        path_cwd, "data/pydata/Germany/germany_divi_ma7.json")

    #
    # blending_fact = get_all_blending_factors(path_results)
    # blending_fact = [0.3]
    blending_fact = [x / 100 for x in range(0, 71, 2)]
    # [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]  # [0.3]  #
    # blending_fact = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    # all_kmax = get_all_kmaxs(os.path.join(
    #     path_results, "BlendingFactorRegional_0.000000"))
    # blending_fact = [0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70]
    all_kmax = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    modes = ["FeedbackDamping"]  # "ClassicDamping",

    icu_compartment = [[7]]
    infected_compartment = [[1, 2, 3, 4, 5, 6, 7]]
    dead_compartment = [[9]]
    flow_se = [[0]]
    flow_ci = [[3]]

    # plot_r0_county_level(path_results, path_plots, modes)
    # plot_risk_county_level(path_results, path_plots, modes)
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
    #            flow_se, [""], "Daily Infections", plot_percentiles=True)
    # list with entry from 0.0 to 1 in 0.05 steps
    regional_proximity = False
    kmax = [x / 100 for x in range(0, 101, 2)]
    loc = 'ICU_Proximity' if regional_proximity else "ICU_Federal_State"
    path_plots = os.path.join(
        path_plots, loc)
    # for km in kmax:
    #     km_formatted = f"{km:.2f}"
    #     res_dir = "results"
    #     if regional_proximity:
    #         res_dir += "/new_regional_def"
    #     results_files = [
    #         os.path.join(
    #             path_cwd, res_dir, "ICUCap_9.000000/rho_1.000000/BlendingFactorRegional_0.000000", f'kmin_0.000000_kmax_{km_formatted}0000'),
    #         os.path.join(
    #             path_cwd, res_dir, "ICUCap_9.000000/rho_1.000000/BlendingFactorRegional_0.420000", f'kmin_0.000000_kmax_{km_formatted}0000'),
    #         os.path.join(
    #             path_cwd, res_dir, "ICUCap_9.000000/rho_1.000000/BlendingFactorRegional_0.700000", f'kmin_0.000000_kmax_{km_formatted}0000')
    #     ]
    #     plot_icu_comp(results_files, path_plots, modes,
    #                   path_icu_data, plot_percentiles=False, log_scale=False, kmax=km)
    # plot_r0(path_results, path_plots, modes)
    # plot_r0_all_scenarios(os.path.join(
    #     path_cwd, "results", "ICUCap_" + str(icu_cap[cap_indx]) + ".000000"), path_plots, modes)

    # rhos = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
    #         1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2]
    # rhos = [0.9,  0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
    #         1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1]
    # rhos = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9,  0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
    #         1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.2,
    #         1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2]
    rhos = [0.92, 1.0, 1.08]
    # for rho in rhos:
    #     # rho mit 6 Nachkommastellen
    #     path_rho = os.path.join(
    #         path_results, f"rho_{rho:.6f}")
    #     violin_plot_all_kmax(path_rho, path_plots, modes, blending_fact, all_kmax, icu_compartment,
    #                          flows=False, title='Peak Value', dir_type='kmin', log_scale=True, rho=rho)

    # plot_peak_time_kmax(path_results, path_plots, blending_fact, modes,
    #                     icu_compartment, flows=False, dir_type='kmin', log_scale=True)

    # plot_peaks_all_blending_factors(path_results, path_plots, blending_fact)

    # plot_peak_time_kmax_variance(path_results, path_plots, blending_fact, modes,
    #                              icu_compartment, rhos, flows=False, log_scale=False)

    plot_peak_time_kmax_variance(path_results, path_plots, blending_fact, modes,
                                 icu_compartment, rhos, flows=False, log_scale=False, title='Max Value')

    # plot_peak_integral(path_results, path_plots, blending_fact,
    #                    modes, icu_compartment, rhos, flows=False, log_scale=False)

    # plot_icu_all_scenarios(path_results, path_plots, modes, log_scale=True)
    # plot_peaks(path_results, path_plots, modes, flow_se)

    # plot_peaks_single(path_results, path_plots, ["FeedbackDamping"], flow_se)

    # plot_icu_comp_all_dirs(path_results, path_plots, modes,
    #                        path_icu_data, plot_percentiles=True)

    # plot_peak_values(path_results, path_plots, [
    #     "ClassicDamping"], icu_compartment, plot_type='val', title='ICU Occupancy', flows=False, vertical=True, dir_value='fixed')
    # plot_peak_values(path_results, path_plots, [
    #     "ClassicDamping"], flow_se, plot_type='val', title='Daily Infections', vertical=True, dir_value='fixed')

    # plot_kmax_diff_rhos(path_results, path_plots, modes, blending_fact, all_kmax, icu_compartment,
    #                     flows=False, title='Peak Value', dir_type='kmin', log_scale=True)
