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
import h5py

import seaborn as sns

sns.set_style("darkgrid")

opacity = 0.15
lineWidth = 2
fontsize = 28
legendsize = 24
ticks = 20
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


def read_total_results_h5(path, comp):
    f = h5py.File(path, 'r')
    group = f['0']

    # This assumes group[some_key_inside_the_group] is a dataset,
    # and returns a np.array:
    # time = group['Time'][()]
    total = group['Total'][()]
    # group1 = group['Group1'][()]

    comp_simulated = np.sum(total[:, comp], axis=1)
    f.close()

    return comp_simulated


def get_num_more_than_x(path, comp, x):
    f = h5py.File(path, 'r')
    count = np.zeros(f['1001']['Total'].shape[0])
    # iteriere über alle keys in f
    for key in f.keys():
        group = f[key]
        total = group['Total'][()]
        comp_simulated = np.sum(total[:, comp], axis=1)
        # prüfe jeden Eintrag von comp_simulated, ob er größer als x ist.
        count += np.where(comp_simulated > x, 1, 0)
    f.close()

    return count


def moving_average_smoothing(data, window_size=3):
    padded_data = np.pad(
        data, ((0, 0), (window_size // 2, window_size // 2)), mode='edge')
    smoothed_data = np.zeros_like(data)

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            smoothed_data[i, j] = np.mean(padded_data[i, j:j + window_size])

    return smoothed_data


def plot_icu_occupancy_per_scenario(path_results, indx_comps, interventions, flows,  num_runs):
    for intervention in interventions:
        results = []
        for index, comp in enumerate(indx_comps):
            data_p25 = []
            data_p50 = []
            data_p75 = []
            for run in range(num_runs):
                flows_add = ""
                if flows[index]:
                    flows_add = "flows"
                path = os.path.join(
                    path_results, intervention, "run_" + str(run), flows_add)

                data_p25.append(read_total_results_h5(os.path.join(
                    path,  "p25", "Results_sum.h5"), comp))
                data_p50.append(read_total_results_h5(os.path.join(
                    path, "p50", "Results_sum.h5"), comp))
                data_p75.append(read_total_results_h5(os.path.join(
                    path, "p75", "Results_sum.h5"), comp))
            data_p25 = np.sort(data_p25, axis=0)
            data_p50 = np.sort(data_p50, axis=0)
            data_p75 = np.sort(data_p75, axis=0)

            if flows[index]:
                data_p25 = np.diff(data_p25, axis=1)
                data_p50 = np.diff(data_p50, axis=1)
                data_p75 = np.diff(data_p75, axis=1)

                # smooth the data
                # window = 3
                # data_p25 = moving_average_smoothing(data_p25, window)
                # data_p50 = moving_average_smoothing(data_p50, window)
                # data_p75 = moving_average_smoothing(data_p75, window)

            comp_label = "Total in Hosptialized"
            if comp == indx_comps[1]:
                comp_label = "Total in ICU"
            elif comp == indx_comps[2]:
                comp_label = "Daily Transmissions"

            # median of num_runs
            median_runs = 0  # int(np.median(range(num_runs)))
            entry = {
                "intervention": intervention,
                "compartment": comp_label,
                "data_p25": data_p25[median_runs],
                "data_p50": data_p50[median_runs],
                "data_p75": data_p75[median_runs]
            }
            results.append(entry)
        # Plotting
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))

        plt.xticks(fontsize=ticks)
        plt.yticks(fontsize=ticks)

        # ax.set_title(title)
        ax.set_xlabel("Time [days]", fontsize=fontsize)
        ax.set_ylabel("Number of Individuals", fontsize=fontsize)
        ax.grid(True)
        color_indx = 0
        for entry in results:
            linestyle = '--'
            linewidth = 3
            color = colors[color_indx]

            if entry["intervention"] == "No_intervention":
                label = "No intervention"
            elif entry["intervention"] == "20p_reduc":
                label = "20% reduction"
            elif entry["intervention"] == "40p_reduc":
                label = "40% reduction"
            elif entry["intervention"] == "60p_reduc":
                label = "60% reduction"
            else:
                print("Intervention not found")

            # cut x-cut first days for all plot data
            x_cut = 14

            ax.plot(
                entry["data_p25"][x_cut:],
                linewidth=linewidth, linestyle=linestyle, color=color)
            ax.plot(
                entry["data_p75"][x_cut:],
                linewidth=linewidth, linestyle='--', color=color)
            ax.plot(
                entry["data_p50"][x_cut:], label=entry["compartment"],
                linewidth=linewidth, linestyle='-', color=color)
            ax.fill_between(
                np.arange(0, len(entry["data_p25"][x_cut:])), entry["data_p25"][x_cut:], entry["data_p75"][x_cut:], color=color, alpha=opacity)
            color_indx += 1
        ax.legend(fontsize=legendsize)
        plt.tight_layout()
        plt.yscale('log')
        # plt.show()
        plt.savefig(
            os.path.join(plot_path, f"{intervention}.png"))
        plt.clf()


def get_num_more_than_x(path, comp, x):
    f = h5py.File(path, 'r')
    count = np.zeros(f['1001']['Total'].shape[0])
    # iteriere über alle keys in f
    for key in f.keys():
        group = f[key]
        total = group['Total'][()]
        comp_simulated = np.sum(total[:, comp], axis=1)
        # prüfe jeden Eintrag von comp_simulated, ob er größer als x ist.
        count += np.where(comp_simulated > x, 1, 0)
    f.close()

    return count[1:]


def plot_num_counties_more_than_x(x, path_results, indx_comp, interventions, regions,  tnt_factors, plot_flows, title, num_runs):
    for intervention in interventions:
        for region in regions:
            results = []
            for tnt_factor in tnt_factors:
                data_p50 = []
                for run in range(num_runs):
                    flows_add = ""
                    if plot_flows:
                        flows_add = "flows"
                    path = os.path.join(
                        path_results, intervention, region, tnt_factor + str(run), flows_add, "p50", "Results.h5")
                    data_p50.append(get_num_more_than_x(path, indx_comp, x))

                data_p50 = np.sort(data_p50, axis=0)

                if plot_flows:
                    data_p50 = np.diff(data_p50, axis=1)

                entry = {
                    "intervention": intervention,
                    "region": region,
                    "tnt_factor": tnt_factor,
                    "min_data": data_p50[0],
                    "max_data": data_p50[num_runs - 1]
                }
                results.append(entry)

            # Plotting
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))

            # ax.set_title(title)
            ax.set_xlabel("Time [days]", fontsize=fontsize)
            ax.set_ylabel("Number counties", fontsize=fontsize)

            plt.xticks(fontsize=ticks)
            plt.yticks(fontsize=ticks)

            # ax.set_ylim(0, 100)
            # ax.set_yscale('log')
            ax.grid(True)
            count_entry = 0
            for entry in results:
                linestyle = '--'
                linewidth = 3
                color = colors[count_entry]

                num_counties_inf = 40 if entry['region'] == regions[0] else 80
                # transform tnt factor to float
                tnt_fact_str = ''.join(
                    filter(lambda ch: ch.isdigit() or ch == '.', entry['tnt_factor']))
                tnt_fact = round(float(tnt_fact_str), 1)

                label = "No testing capacities"
                if tnt_fact > 0.0:
                    label = "Testing capacities"

                ax.plot(
                    entry["max_data"], label=label + " max",
                    linewidth=linewidth, linestyle=linestyle, color=color)
                ax.plot(
                    entry["min_data"], label=label + " min",
                    linewidth=linewidth, linestyle='--', color=color)
                ax.fill_between(
                    np.arange(0, len(entry["max_data"])), entry["min_data"], entry["max_data"], color=color, alpha=opacity)

                count_entry += 1

            ax.legend(fontsize=legendsize, loc='upper left')
            plt.tight_layout()
            plt.savefig(
                os.path.join(plot_path, f"week2_num_infected_{region}.png"))
            plt.clf()


if __name__ == '__main__':
    path_results = "/localdata1/test/memilio/test"

    plot_flows = False

    # num samples
    num_runs = 1

    # create dir plots in path_results
    plot_path = os.path.join(path_results, "plots", "2d")
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)

    interventions = ["No_intervention", "20p_reduc", "40p_reduc", "60p_reduc"]

    #  choose compartments, we want to plot
    infected_compartments = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                             14, 15, 16, 17, 18, 19, 20, 21, 22]

    susceptible = [0, 1, 23]
    hosp = [17, 18, 19]
    icu = [20, 21, 22]

    inf_naive = [5, 11]

    flows_se = [0, 17, 33]

    flows_hu = [11, 28, 44]

    all_comp = np.arange(0, 29, 1)

    indx_comps = [hosp, icu, flows_se]
    flows = [False, False, True]

    plot_icu_occupancy_per_scenario(path_results, indx_comps,
                                    interventions, flows, num_runs)
