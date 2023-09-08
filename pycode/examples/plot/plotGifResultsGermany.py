import datetime as dt
import os.path
import imageio

import numpy as np
import pandas as pd

import memilio.epidata.getPopulationData as gpd
import memilio.plot.plotMap as pm
from memilio.epidata import geoModificationGermany as geoger


def create_plot_map(day, age_groups, relative, filename, file_format, files_input, output_path):

    i = 0
    for file in files_input.values():
        # MEmilio backend hdf5 example

        deads = [24, 25, 26]
        susceptible = [0, 1, 23]
        infected = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

        df = pm.extract_data(
            file, region_spec=None, column=None, date=day,
            filters={'Group': filter_age, 'InfectionState': infected},
            file_format=file_format)

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
                    out_folder='data/pydata/Germany/', no_raw=True, merge_eisenach=True)

            # For fitting of different age groups we need format ">X".
            age_group_values = list(age_groups.values())
            age_group_values[-1] = age_group_values[-1].replace(
                '80+', '>79')
            # scale data
            df = pm.scale_dataframe_relative(
                df, age_group_values, population)

        if i == 0:
            dfs_all = pd.DataFrame(df.iloc[:, 0])

        dfs_all[df.columns[-1] + ' ' + str(i)] = df[df.columns[-1]]
        i += 1

    dfs_all = dfs_all.apply(pd.to_numeric, errors='coerce')

    dfs_all_sorted = dfs_all.sort_values(by='Region')
    dfs_all_sorted = dfs_all_sorted.reset_index(drop=True)

    min_val = dfs_all_sorted[dfs_all_sorted.columns[1:]].min().min()
    max_val = dfs_all_sorted[dfs_all_sorted.columns[1:]].max().max()

    pm.plot_map(
        dfs_all_sorted, scale_colors=np.array([min_val, max_val]),
        legend=['', ''],
        title='Synthetic data (relative) day ' + str(day), plot_colorbar=True,
        output_path=output_path,
        fig_name=filename, dpi=300,
        outercolor=[205 / 255, 238 / 255, 251 / 255])


if __name__ == '__main__':

    files_input = {'Data set 1': 'p75/Results'}
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

    relative = True
    time_steps = pm.extract_time_steps(
        files_input[list(files_input.keys())[0]], file_format=file_format)
    filename = "file"

    # create gif
    frames = []
    output_path = "plot_gif"

    for day in range(0, time_steps):
        create_plot_map(day, age_groups, relative, filename,
                        file_format, files_input, output_path)
        image = imageio.v2.imread(os.path.join(output_path, filename + ".png"))
        frames.append(image)

    imageio.mimsave(os.path.join(output_path, 'sim.gif'),  # output gif
                    frames,          # array of input frames
                    duration=10,
                    loop=0)         # optional: frames per second
