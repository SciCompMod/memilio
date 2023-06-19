import datetime as dt
import os.path
import imageio

import numpy as np
import pandas as pd

import memilio.epidata.getPopulationData as gpd
import memilio.plot.plotMap as pm

if __name__ == '__main__':

    # files_input = {'Data set 1': '/localdata1/test/memilio/test2/p75/Results'}
    # file_format = 'h5'
    # # Define age groups which will be considered through filtering
    # # Keep keys and values as well as its assignment constant, remove entries
    # # if only part of the population should be plotted or considered, e.g., by
    # # setting:
    # # age_groups = {1: '5-14', 2: '15-34'}
    # age_groups = {0: '0-4', 1: '5-14', 2: '15-34',
    #               3: '35-59', 4: '60-79', 5: '80+'}
    # if len(age_groups) == 6:
    #     filter_age = None
    # else:
    #     if file_format == 'json':
    #         filter_age = [val for val in age_groups.values()]
    #     else:
    #         filter_age = ['Group' + str(key) for key in age_groups.keys()]

    # relative = True
    num_days = 15
    # for day in range(0, num_days):

    #     i = 0
    #     for file in files_input.values():
    #         # MEmilio backend hdf5 example
    #         df = pm.extract_data(
    #             file, region_spec=None, column=None, date=day,
    #             filters={'Group': filter_age, 'InfectionState': [3, 4]},
    #             file_format=file_format)
    #         # # MEmilio epidata json example
    #         # df = pm.extract_data(
    #         #     file, region_spec='ID_County',
    #         #     column='Synthetic data',
    #         #     date=dt.date(2021, 11, 18),
    #         #     filters={'ID_State': None,
    #         #              'Age_RKI': filter_age},
    #         #     file_format=file_format)

    #         if relative:

    #             try:
    #                 population = pd.read_json(
    #                     'data/pydata/Germany/county_current_population.json')
    #             # pandas>1.5 raise FileNotFoundError instead of ValueError
    #             except (ValueError, FileNotFoundError):
    #                 print(
    #                     "Population data was not found. Download it from the internet.")
    #                 population = gpd.get_population_data(
    #                     read_data=False, file_format=file_format,
    #                     out_folder='data/pydata/Germany/', no_raw=True,
    #                     split_gender=False, merge_eisenach=True)

    #             # For fitting of different age groups we need format ">X".
    #             age_group_values = list(age_groups.values())
    #             age_group_values[-1] = age_group_values[-1].replace(
    #                 '80+', '>79')
    #             # scale data
    #             df = pm.scale_dataframe_relative(
    #                 df, age_group_values, population)

    #         if i == 0:
    #             dfs_all = pd.DataFrame(df.iloc[:, 0])

    #         dfs_all[df.columns[-1] + ' ' + str(i)] = df[df.columns[-1]]
    #         i += 1

    #     min_val = dfs_all[dfs_all.columns[1:]].min().min()
    #     max_val = dfs_all[dfs_all.columns[1:]].max().max()

    #     filename = 'customPlot_day_' + str(day)

    #     pm.plot_map(
    #         dfs_all, scale_colors=np.array([min_val, max_val]),
    #         legend=['', ''],
    #         title='Synthetic data (relative)', plot_colorbar=True,
    #         output_path="/localdata1/test/memilio/pycode/examples/plot/plot_gif",
    #         fig_name=filename, dpi=300,
    #         outercolor=[205 / 255, 238 / 255, 251 / 255])

    # create gif
    frames = []
    for day in range(0, num_days):
        filename = 'customPlot_day_' + str(day) + '.png'
        output_path = "/localdata1/test/memilio/pycode/examples/plot/plot_gif"
        image = imageio.v2.imread(os.path.join(output_path, filename))
        frames.append(image)

    imageio.mimsave(os.path.join(output_path, 'sim.gif'),  # output gif
                    frames,          # array of input frames
                    duration=1000)         # optional: frames per second
