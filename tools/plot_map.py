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

if __name__ == '__main__':
    new_scheme = True
    if new_scheme == True:
        path = "/localdata1/code/memilio/results_paper/mask_0_cologne/p75/Results"
    else:
        path = "/localdata1/code/memilio/results_paper/mask_0_cologne_old/p50/Results"

    files_input = {
        'Data set 1': path}
    file_format = 'h5'
    # Define age groups which will be considered through filtering
    # Keep keys and values as well as its assignment constant, remove entries
    # if only part of the population should be plotted or considered, e.g., by
    # setting:
    # age_groups = {1: '5-14', 2: '15-34'}
    age_groups = {0: '0-4', 1: '5-14', 2: '15-34',
                  3: '35-59', 4: '60-79', 5: '80+'}
    if len(age_groups) == 6:
        filter_age = None
    else:
        if file_format == 'json':
            filter_age = [val for val in age_groups.values()]
        else:
            filter_age = ['Group' + str(key) for key in age_groups.keys()]

    relative = False
    num_days = 91

    # if new_scheme:
    #     min_val = 0
    #     max_val = 90695.5836895524
    # else:
    #     min_val = 0
    #     max_val = 7783.476726023636

    # für t = 90
    min_val = 0
    max_val = 62000

    # für t = 20
    # min_val = 0
    # max_val = 5114.970085160735

    # für t = 30
    # min_val = 0
    # max_val = 49250

    # für t = 40
    # min_val = 0
    # max_val = 89824.62877071957

    # make plot for colorbar
    # norm = mcolors.LogNorm(
    #     vmin=min_val if min_val > 0 else 1e-13, vmax=max_val if max_val > 0 else 1e-13)

    norm = SymLogNorm(linthresh=1, linscale=0.7,
                      vmin=min_val, vmax=max_val)

    # Create a ScalarMappable with the same Norm and Colormap as your data
    colors = ["green", "yellow", "red", "purple"]
    cmap = LinearSegmentedColormap.from_list("my_colormap", colors)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Create a new figure for the colorbar
    cbar_fig, ax = plt.subplots(figsize=(8, 1))
    cbar = plt.colorbar(sm, orientation='horizontal', cax=ax)
    cbar.ax.tick_params(labelsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(
        "/localdata1/code/memilio/results_paper/plot_gif", 'colorbar.png'), dpi=300)
    plt.clf()

    for day in range(0, num_days):

        i = 0
        for file in files_input.values():
            # MEmilio backend hdf5 example

            timms = [27, 28]
            deads = [24, 25, 26]
            susceptible = [0, 1, 23]
            infected = [5, 6, 7, 8, 9, 10, 11,
                        12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
            symptomatic = [11, 12, 13, 14, 15, 16]  # 2, 3, 4,
            icu = [20, 21, 22]
            exposed = [2, 3, 4]

            df = pm.extract_data(
                file, region_spec=None, column=None, date=day,
                filters={'Group': filter_age, 'InfectionState': infected},
                file_format=file_format)

            # all entries to numeric
            df = df.apply(pd.to_numeric, errors='coerce')
            dfs_all = df

            # # MEmilio epidata json example
            # df = pm.extract_data(
            #     file, region_spec='ID_County',
            #     column='Synthetic data',
            #     date=dt.date(2021, 11, 18),
            #     filters={'ID_State': None,
            #              'Age_RKI': filter_age},
            #     file_format=file_format)

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

            # if i == 0:
            #     dfs_all = pd.DataFrame(df.iloc[:, 0])

            # dfs_all[df.columns[-1] + ' ' + str(i)] = df[df.columns[-1]]
            # i += 1

        if new_scheme == True:
            filename = 'newscheme_day_' + str(day)
        else:
            filename = 'oldscheme_day_' + str(day)

        dfs_all = dfs_all.apply(pd.to_numeric, errors='coerce')

        dfs_all_sorted = dfs_all.sort_values(by='Region')
        dfs_all_sorted = dfs_all_sorted.reset_index(drop=True)

        # set all values zero of count (rel) 0
        # dfs_all_sorted[dfs_all_sorted.columns[1]] = 0

        # delete cologne from dfs_all_sorted
        dfs_all_sorted = dfs_all_sorted[dfs_all_sorted.Region != 5315]

        # if (max_val < dfs_all_sorted['Count'].max()):
        #     max_val = dfs_all_sorted['Count'].max()
        # print("max_val: " + str(max_val))

        print("day: " + str(day))

        # # manipulate the data to plot the connected counties. Also important to comment the delete cologne block above
        # county_ids = geoger.get_county_ids()
        # # get index of entty 5315
        # id_cologne = county_ids.index(5315)
        # list_connected_counties = [337, 21, 324, 58, 108, 180, 30, 36, 53, 325, 293, 126, 28, 295, 48, 124, 145, 27, 329, 330, 139, 29, 296, 284, 26, 146, 333, 44, 55, 334, 335, 147, 323, 144, 123, 138, 140, 61, 62, 339, 148, 64, 65, 66, 322, 137, 67, 122, 149, 68, 300, 69, 70, 71, 72, 150, 73, 74, 136, 301, 75, 127, 76, 321, 151, 77, 78, 2, 16, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 102,
        #                            104, 105, 106, 107, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 125, 128, 129, 130, 131, 132, 133, 134, 143, 152, 153, 154, 155, 156, 157, 158, 160, 162, 163, 176, 177, 179, 181, 183, 186, 193, 194, 196, 197, 199, 204, 224, 233, 239, 259, 319, 357, 362, 304, 135, 234, 305, 221, 366, 220, 367, 215, 17, 214, 211, 207, 205, 379, 15, 101, 32, 386, 141, 184, 178, 11, 175, 33, 173, 172, 170, 9, 169, 168, 167, 164]

        # county_ids_connected = []
        # for i in list_connected_counties:
        #     county_ids_connected.append(county_ids[i])

        # commuter_path = "/localdata1/code/memilio/test/commuter_migration_with_locals.txt"
        # commuter_data = np.loadtxt(commuter_path, delimiter=' ')

        # commuting_to = commuter_data[:, id_cologne]
        # commuting_from = commuter_data[id_cologne, :]

        # dfs_all_sorted['Count'] = 0

        # # iteriere über alle Zeilen in dfs_all_sorted
        # for index, row in dfs_all_sorted.iterrows():
        #     region_row = row['Region']
        #     if region_row in county_ids_connected:
        #         dfs_all_sorted.at[index, 'Count'] = 1
        #     if region_row == 5315:
        #         dfs_all_sorted.at[index, 'Count'] = 2

        # Setze dann 'Count' auf 1 für alle Zeilen, die in list_connected_counties sind
        # dfs_all_sorted.loc[dfs_all_sorted.index.isin(
        #     list_connected_counties), 'Count'] = 50

        # # Setze 'Count' auf 2 für die Zeile mit dem Index 79
        # dfs_all_sorted.at[79, 'Count'] = 100

        pm.plot_map(norm,
                    dfs_all_sorted, scale_colors=[min_val, max_val],
                    legend=['', ''],
                    title='Synthetic data (relative) day ' + str(day),
                    plot_colorbar=False,
                    output_path="/localdata1/code/memilio/results_paper/plot_gif",
                    fig_name=filename, dpi=300,
                    outercolor='white')

    # create gif
    frames = []
    for day in range(0, num_days):
        if new_scheme == True:
            filename = 'newscheme_day_' + str(day)
        else:
            filename = 'oldscheme_day_' + str(day)
        filename = filename + '.png'
        output_path = "/localdata1/code/memilio/results_paper/plot_gif"
        image = imageio.v2.imread(os.path.join(output_path, filename))
        frames.append(image)

    imageio.mimsave("/localdata1/code/memilio/results_paper/plot_gif/filename.gif",  # output gif
                    frames,
                    duration=10)
