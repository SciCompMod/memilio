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


def get_ids(path_results):
    path_metro = os.path.join(
        path_results, "Metropolis", num_infected[0])
    ids_metro = os.listdir(path_metro)
    ids_metro = [int(i) for i in ids_metro]

    path_rural = os.path.join(
        path_results, "Rural area", num_infected[0])
    ids_rural = os.listdir(path_rural)
    ids_rural = [int(i) for i in ids_rural]
    return ids_metro, ids_rural


def plot_total_transmissions(path_results, percentile, indx_comp, areas, num_infected, interventions,  plot_flows):

    id_metro, id_rural = get_ids(os.path.join(path_results, interventions[0]))
    # progress bar
    total_iterations = len(interventions) * len(areas) * \
        len(num_infected) * max(len(id_metro), len(id_rural))

    progress_bar = tqdm(total=total_iterations)

    for area in areas:
        ids = id_metro if area == "Metropolis" else id_rural
        for id in ids:
            results = []
            for infected in num_infected:
                for intervention in interventions:
                    path = os.path.join(
                        path_results, intervention, area, infected, str(id), percentile, "Results_sum.h5")
                    if plot_flows:
                        path = os.path.join(path_results, intervention, area, infected, str(
                            id), 'flows', percentile, "Results_sum.h5")

                    data = read_total_results_h5(path, indx_comp)
                    entry = {
                        "intervention": intervention,
                        "infected": infected,
                        "data": data
                    }
                    results.append(entry)
                    progress_bar.update(1)

            # Plotting
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            ids_all = geoger.get_county_names_and_ids(
                merge_berlin=True, merge_eisenach=True, zfill=True)
            county_name = ""
            county_name = next(
                (entry[0] for entry in ids_all if int(entry[1]) == id), None)
            ax.set_title(f"ID: {id}, {area}")
            ax.set_xlabel("Time [days]")
            ax.set_ylabel("Daily Transmissions")
            ax.set_ylim(0, 1200)
            # ax.set_yscale('log')
            ax.grid(True)
            ax.set_title(f"County: {county_name}, {area}")
            for entry in results:
                linestyle = '-'
                linewidth = 3
                color = ''
                if entry['intervention'] == 'No_intervention':
                    color = 'blue'
                elif entry['intervention'] == 'Fixed_Damping':
                    color = 'green'
                elif entry['intervention'] == 'Dynamic_NPI':
                    color = 'red'
                if entry['infected'] == '1_Infected':
                    linestyle = '--'
                    linewidth = 3
                ax.plot(
                    entry["data"], label=f"{entry['intervention']}_{entry['infected']}",
                    linewidth=linewidth, linestyle=linestyle, color=color)
            ax.legend()
            plt.savefig(
                os.path.join(plot_path, f"{area}_{id}.png"))
            plt.clf()


def plot_transmissions_per_scenario(path_results, percentile, indx_comp,
                                    areas, num_infected, interventions, plot_flows):
    id_metro, id_rural = get_ids(os.path.join(path_results, interventions[0]))
    # progress bar
    total_iterations = len(interventions) * len(areas) * \
        len(num_infected) * max(len(id_metro), len(id_rural))

    for infected in num_infected:
        results = []
        for area in areas:
            ids = id_metro if area == "Metropolis" else id_rural
            for intervention in interventions:
                for id in ids:
                    path = os.path.join(
                        path_results, intervention, area, infected, str(id), percentile, "Results_sum.h5")
                    if plot_flows:
                        path = os.path.join(path_results, intervention, area, infected, str(
                            id), 'flows', percentile, "Results_sum.h5")

                    data = read_total_results_h5(path, indx_comp)
                    entry = {
                        "area": area,
                        "intervention": intervention,
                        "id": id,
                        "data": data
                    }

                    results.append(entry)

            # erstelle eine neue results_minmax Liste, wo die min und max Werte für jedes Szenarion gespeichert werden. EIn szenario ist ein Tupel aus area und intervention
            results_minmax = []
            for intervention in interventions:
                datas = []
                for entry in results:
                    if entry['area'] == area and entry['intervention'] == intervention:
                        datas.append(entry['data'])

                min_data = np.min(datas, axis=0)
                max_data = np.max(datas, axis=0)
                results_minmax.append({
                    "area": area,
                    "intervention": intervention,
                    "min_data": min_data,
                    "max_data": max_data
                })

            # Plotting
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))

            ax.set_title(f"A: {area}, I: {infected}")
            ax.set_xlabel("Time [days]")
            ax.set_ylabel("Daily Transmissions")
            if infected == "1_Infected":
                ax.set_ylim(0, 150)
            else:
                ax.set_ylim(0, 1200)
            # ax.set_yscale('log')
            ax.grid(True)
            ax.set_title(f"{area}")
            for entry in results_minmax:
                linestyle = '--'
                linewidth = 3
                color = ''
                label = ""
                if entry['intervention'] == 'No_intervention':
                    color = colors[0]
                    label = "No Intervention"
                elif entry['intervention'] == 'Fixed_Damping':
                    color = colors[1]
                    label = "Fixed Damping"
                elif entry['intervention'] == 'Dynamic_NPI':
                    color = colors[2]
                    label = "Dynamic NPI"

                ax.plot(
                    entry["max_data"], label=f"{entry['area']} {label} max",
                    linewidth=linewidth, linestyle=linestyle, color=color)
                ax.plot(
                    entry["min_data"], label=f"{entry['area']} {label} min",
                    linewidth=linewidth, linestyle='--', color=color)
                ax.fill_between(
                    np.arange(0, len(entry["max_data"])), entry["min_data"], entry["max_data"], color=color, alpha=opacity)
            ax.legend()
            plt.tight_layout()
            plt.savefig(
                os.path.join(plot_path, f"{infected}_quartils_{area}.png"))
            plt.clf()


def plot_num_counties_more_than_x(x, path_results, percentile, indx_comp, areas, num_infected, interventions, plot_flows):
    id_metro, id_rural = get_ids(os.path.join(path_results, interventions[0]))
    # progress bar
    total_iterations = len(interventions) * len(areas) * \
        len(num_infected) * max(len(id_metro), len(id_rural))

    progress_bar = tqdm(total=total_iterations)

    for area in areas:
        ids = id_metro if area == "Metropolis" else id_rural
        for infected in num_infected:
            results = []
            for intervention in interventions:
                for id in ids:
                    path = os.path.join(
                        path_results, intervention, area, infected, str(id), percentile, "Results.h5")
                    if plot_flows:
                        path = os.path.join(path_results, intervention, area, infected, str(
                            id), 'flows', percentile, "Results.h5")

                    data = get_num_more_than_x(path, indx_comp, x)
                    entry = {
                        "area": area,
                        "intervention": intervention,
                        "id": id,
                        "data": data
                    }
                    results.append(entry)
                    progress_bar.update(1)

            # erstelle eine neue results_minmax Liste, wo die min und max Werte für jedes Szenarion gespeichert werden. EIn szenario ist ein Tupel aus area und intervention
            results_minmax = []
            for intervention in interventions:
                datas = []
                for entry in results:
                    if entry['area'] == area and entry['intervention'] == intervention:
                        datas.append(entry['data'])

                min_data = np.min(datas, axis=0)
                max_data = np.max(datas, axis=0)
                results_minmax.append({
                    "area": area,
                    "intervention": intervention,
                    "min_data": min_data,
                    "max_data": max_data
                })

            # Plotting
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))

            ax.set_title(f"A: {area}, I: {infected}")
            ax.set_xlabel("Time [days]")
            ax.set_ylabel("Counties with more than " +
                          str(int(x)) + " transmissions")
            ax.set_ylim(0, 100)
            # ax.set_yscale('log')
            ax.grid(True)
            ax.set_title(f"{area}")
            for entry in results_minmax:
                linestyle = '--'
                linewidth = 3
                color = ''
                label = ""
                if entry['intervention'] == 'No_intervention':
                    color = colors[0]
                    label = "No Intervention"
                elif entry['intervention'] == 'Fixed_Damping':
                    color = colors[1]
                    label = "Fixed Damping"
                elif entry['intervention'] == 'Dynamic_NPI':
                    color = colors[2]
                    label = "Dynamic NPI"

                ax.plot(
                    entry["max_data"], label=f"{entry['area']} {label} max",
                    linewidth=linewidth, linestyle=linestyle, color=color)
                ax.plot(
                    entry["min_data"], label=f"{entry['area']} {label} min",
                    linewidth=linewidth, linestyle='--', color=color)
                ax.fill_between(
                    np.arange(0, len(entry["max_data"])), entry["min_data"], entry["max_data"], color=color, alpha=opacity)
            ax.legend()
            plt.tight_layout()
            plt.savefig(
                os.path.join(plot_path, f"{infected}_quartils_{area}_count.png"))
            plt.clf()


if __name__ == '__main__':
    path_results = "/localdata1/test/memilio/test"

    plot_flows = True
    percentile = "p50"

    # create dir plots in path_results
    plot_path = os.path.join(path_results, "plots", "2d")
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)

    interventions = ["No_intervention", "Fixed_Damping", "Dynamic_NPI"]
    areas = ["Metropolis", "Rural area"]
    num_infected = ["1_Infected", "10_Infected"]

    #  choose compartments, we want to plot
    # TODO: Missing for flows
    timms = [27, 28]
    deads = [24, 25, 26]
    susceptible = [0, 1, 23]
    infected_compartments = [5, 6, 7, 8, 9, 10, 11, 12, 13,
                             14, 15, 16, 17, 18, 19, 20, 21, 22]
    infected_symp = [11, 12, 13, 14, 15, 16]
    symptomatic = [11, 12, 13, 14, 15, 16]  # 2, 3, 4,
    icu = [20, 21, 22]
    exposed = [2, 3, 4]

    flows_se = [0, 17, 33]  # , 53, 70, 86, 106, 123, 139,
    # 159, 176, 192, 212, 229, 245, 265, 282, 298]

    # plot_total_transmissions(path_results, percentile, flows_se,
    #                          areas, num_infected, interventions, plot_flows)

    plot_transmissions_per_scenario(path_results, percentile, flows_se,
                                    areas, num_infected, interventions, plot_flows)

    plot_num_counties_more_than_x(1., path_results, percentile, flows_se,
                                  areas, num_infected, interventions, plot_flows)
