import datetime as dt
import os.path
import imageio

import numpy as np
import pandas as pd

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

start_date = "2020-10-15"
total_pop = 83278910.0
opacity = 0.15
lineWidth = 2
fontsize = 28
legendsize = 15
ticks = 15
# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors = sns.color_palette("colorblind", 16)


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


def plot(ys, labels, path_plots, title="", log_scale=False, ylabel="Number Individuals"):
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
            ax.plot(days, y["p25"], linewidth=lineWidth,
                    linestyle='--', color=colors[indx_data % len(colors)])
            ax.plot(days, y["p75"], linewidth=lineWidth,
                    linestyle='--', color=colors[indx_data % len(colors)])
            ax.plot(days, y["p50"], label=labels[indx_data], linewidth=lineWidth,
                    linestyle='-', color=colors[indx_data % len(colors)])
            ax.fill_between(
                days, y["p25"], y["p75"], color=colors[indx_data % len(colors)], alpha=opacity)
        else:
            ax.plot(days, y, label=labels[indx_data], linewidth=lineWidth,
                    linestyle='-', color=colors[indx_data % len(colors)])

    for ax in axes:
        ax.legend(fontsize=legendsize, loc='center left',
                  bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig(os.path.join(path_plots, f"plot_{title}.png"))
    plt.clf()
    return 0


def get_results(path_results, indx_comp, group_key='Total', results="total"):
    percentiles = ["p25", "p50", "p75"]
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


def plot_risk(path_results, path_plots, log_scale=False):
    if not os.path.exists(path_plots):
        os.makedirs(path_plots)
    num_counties = 400
    path_risk_results = os.path.join(path_results, "FeedbackSim", "risk")
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
         log_scale=log_scale, ylabel="Perceived Risk")


def plot_compartments(path_results, path_plots, modes, compartments, labels, title, log_scale=False):
    if not os.path.exists(path_plots):
        os.makedirs(path_plots)
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
    if not os.path.exists(path_plots):
        os.makedirs(path_plots)
    plot_data = []
    labels_modes = []
    for mode in modes:
        for index, compartment in enumerate(flow_indx):
            labels_modes.append(labels[index] + f" {mode}")
            path_results_mode = os.path.join(path_results, mode, "flows")
            plot_data.append(get_results(
                path_results_mode, compartment, results="total"))
    plot(plot_data, labels_modes, path_plots, title=title,
         log_scale=log_scale, ylabel="Number Individuals")
    return 0


def plot_r0(path_results, path_plots, modes, percentile="p50"):
    file_format = 'h5'
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

    total_pop = population['Population'].sum()

    r0_nums = []

    for mode in modes:
        r0_mode = []
        path_results_mode = os.path.join(path_results, mode, "r0")
        path = os.path.join(path_results_mode, percentile,  "Results.h5")
        f = h5py.File(path, 'r')
        num_days = f['1001']['Group1'].shape[0]
        for day in range(num_days):
            r0_val = 0
            for key in f.keys():
                group = f[key]
                # search for key in population['ID_County]
                indx_key = population.index[population['ID_County'] == int(
                    key)]
                r0_val += group['Group1'][()][day][0] * \
                    population['Population'][indx_key].values[0] / total_pop
            r0_mode.append(r0_val)
        f.close()
        r0_nums.append(r0_mode)

    plot(r0_nums, modes, path_plots, title="R0", ylabel="R0")


def plot_icu_comp(path_results, path_plots, modes, path_icu_data, log_scale=False):
    if not os.path.exists(path_plots):
        os.makedirs(path_plots)
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

    # plot_risk(path_results, path_plots)
    # plot_compartments(path_results, path_plots, modes,
    #                   icu_compartment, ["ICU Occupancy"], "ICU Occupancy")
    # plot_compartments(path_results, path_plots, modes,
    #                   infected_compartment, [""], "Total Infected")
    # plot_compartments(path_results, path_plots, modes,
    #                   dead_compartment, [""], "Total Deaths")
    # plot_flows(path_results, path_plots, modes,
    #            flow_se, [""], "Daily Infections")
    # plot_icu_comp(path_results, path_plots, modes, path_icu_data)

    plot_r0(path_results, path_plots, modes)
