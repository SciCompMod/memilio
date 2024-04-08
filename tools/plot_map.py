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


if __name__ == '__main__':
    path_results = "/localdata1/test/memilio/test"

    plot_flows = False
    percentile = "p50"

    # create dir plots in path_results
    plot_path = os.path.join(path_results, "plots")
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)

    interventions = ["No_intervention", "Fixed_Damping", "Dynamic_NPI"]
    areas = ["Metropolis", "Rural area"]
    num_infected = ["1_Infected", "10_Infected"]
    days = [0, 10, 20, 30]

    #  choose compartments, we want to plot
    # TODO: Missing for flows
    timms = [27, 28]
    deads = [24, 25, 26]
    susceptible = [0, 1, 23]
    infected = [5, 6, 7, 8, 9, 10, 11, 12, 13,
                14, 15, 16, 17, 18, 19, 20, 21, 22]
    infected_symp = [11, 12, 13, 14, 15, 16]
    symptomatic = [11, 12, 13, 14, 15, 16]  # 2, 3, 4,
    icu = [20, 21, 22]
    exposed = [2, 3, 4]

    id_metro, id_rural = get_ids(os.path.join(path_results, interventions[0]))

    get_max_val = 0

    # progress bar
    total_iterations = len(interventions) * len(areas) * \
        len(num_infected) * max(len(id_metro), len(id_rural)) * len(days)

    progress_bar = tqdm(total=total_iterations)

    for intervention in interventions:
        for area in areas:
            for infected in num_infected:
                ids = id_metro if area == "Metropolis" else id_rural
                for id in ids:
                    path = os.path.join(
                        path_results, intervention, area, infected, str(id), percentile, "Results")

                    files_input = {
                        'Data set 1': path}
                    file_format = 'h5'
                    age_groups = {0: '0-4', 1: '5-14', 2: '15-34',
                                  3: '35-59', 4: '60-79', 5: '80+'}
                    if len(age_groups) == 6:
                        filter_age = None
                    else:
                        if file_format == 'json':
                            filter_age = [val for val in age_groups.values()]
                        else:
                            filter_age = ['Group' + str(key)
                                          for key in age_groups.keys()]

                    relative = False
                    num_days = 31

                    min_val = 0
                    max_val = 100

                    norm = SymLogNorm(linthresh=1, linscale=0.7,
                                      vmin=min_val, vmax=max_val)

                    # Own colorbar
                    colors = ["green", "yellow", "red", "purple"]
                    cmap = LinearSegmentedColormap.from_list(
                        "my_colormap", colors)
                    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                    sm.set_array([])

                    # Create a new figure for the colorbar
                    cbar_fig, ax = plt.subplots(figsize=(8, 1))
                    cbar = plt.colorbar(sm, orientation='horizontal', cax=ax)
                    cbar.ax.tick_params(labelsize=10)
                    plt.tight_layout()
                    plt.savefig(os.path.join(
                        plot_path, 'colorbar.png'), dpi=300)
                    plt.clf()

                    for day in days:

                        i = 0
                        for file in files_input.values():
                            # MEmilio backend hdf5 example

                            df = pm.extract_data(
                                file, region_spec=None, column=None, date=day,
                                filters={'Group': filter_age,
                                         'InfectionState': infected_symp},
                                file_format=file_format)

                            # all entries to numeric
                            df = df.apply(pd.to_numeric, errors='coerce')
                            dfs_all = df

                            if relative:

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

                                # For fitting of different age groups we need format ">X".
                                age_group_values = list(age_groups.values())
                                age_group_values[-1] = age_group_values[-1].replace(
                                    '80+', '>79')
                                # scale data
                                df = pm.scale_dataframe_relative(
                                    df, age_group_values, population)

                        filename = area + "_" + infected + "_" + \
                            str(id) + "_" + intervention + "_day_" + str(day)

                        dfs_all = dfs_all.apply(pd.to_numeric, errors='coerce')

                        dfs_all_sorted = dfs_all.sort_values(by='Region')
                        dfs_all_sorted = dfs_all_sorted.reset_index(drop=True)

                        max_val = dfs_all_sorted['Count'].max()
                        if max_val > get_max_val:
                            get_max_val = max_val

                        pm.plot_map(norm,
                                    dfs_all_sorted, scale_colors=[
                                        min_val, max_val],
                                    legend=['', ''],
                                    title='Synthetic data (relative) day ' +
                                    str(day),
                                    plot_colorbar=False,
                                    output_path=plot_path,
                                    fig_name=filename, dpi=300,
                                    outercolor='white')

                        progress_bar.update()
    progress_bar.close()
    print("max value: ", get_max_val)

    # # create gif
    # frames = []
    # for day in range(0, num_days):
    #     if new_scheme == True:
    #         filename = 'newscheme_day_' + str(day)
    #     else:
    #         filename = 'oldscheme_day_' + str(day)
    #     filename = filename + '.png'
    #     output_path = "/localdata1/code/memilio/results_paper/plot_gif"
    #     image = imageio.v2.imread(
    #         os.path.join(output_path, filename))
    #     frames.append(image)

    # imageio.mimsave("/localdata1/code/memilio/results_paper/plot_gif/filename.gif",  # output gif
    #                 frames,          # array of input frames
    #                 duration=10)         # optional: frames per second
