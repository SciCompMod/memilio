#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Henrik Zunker, Maximilian Betz
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################

import datetime as dt
import os.path
import imageio
import tempfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import memilio.epidata.getPopulationData as gpd
import memilio.plot.plotMap as pm
from memilio.epidata import geoModificationGermany as geoger
import memilio.epidata.progress_indicator as progind
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def create_plot_map(day, filename, files_input, output_path, compartments,  file_format='h5', relative=False,  age_groups={0: '0-4', 1: '5-14', 2: '15-34', 3: '35-59', 4: '60-79', 5: '80+'}):
    """ Plots region-specific information for a single day of the simulation.

    :param day: Day of the simulation.
    :param filename: Name of the file to be created.
    :param files_input: Dictionary of input files.
    :param output_path: Output path for the figure.
    :param compartments: List of compartments to be plotted.
    :param file_format: Format of the file to be created. Either 'h5' or 'json'. (Default value = 'h5')
    :param relative: Defines if data should be scaled relative to population. (Default value = False)
    :param age_groups: Dictionary of age groups to be considered. (Default value = {0: '0-4', 1: '5-14', 2: '15-34', 3: '35-59', 4: '60-79', 5: '80+'}) 

    """

    if len(age_groups) == 6:
        filter_age = None
    else:
        if file_format == 'json':
            filter_age = [val for val in age_groups.values()]
        else:
            filter_age = ['Group' + str(key) for key in age_groups.keys()]

    # In file_input there can be two different files. When we enter two files,
    # both files are plotted side by side in the same figure.
    file_index = 0
    for file in files_input.values():

        df = pm.extract_data(
            file, region_spec=None, column=None, date=day,
            filters={'Group': filter_age, 'InfectionState': compartments},
            file_format=file_format)

        if relative:

            try:
                population = pd.read_json(
                    'data/pydata/Germany/county_current_population.json')
            # pandas>1.5 raise FileNotFoundError instead of ValueError
            except (ValueError, FileNotFoundError):
                print(
                    "Population data was not found. Downloading it from the internet.")
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

        if file_index == 0:
            dfs_all = pd.DataFrame(df.iloc[:, 0])

        dfs_all[df.columns[-1] + ' ' + str(file_index)] = df[df.columns[-1]]
        file_index += 1

    dfs_all = dfs_all.apply(pd.to_numeric, errors='coerce')

    dfs_all_sorted = dfs_all.sort_values(by='Region')
    dfs_all_sorted = dfs_all_sorted.reset_index(drop=True)

    min_val = dfs_all_sorted[dfs_all_sorted.columns[1:]].min().min()
    max_val = dfs_all_sorted[dfs_all_sorted.columns[1:]].max().max()

    pm.plot_map(
        dfs_all_sorted, scale_colors=np.array([min_val, max_val]),
        legend=['', ''],
        title='Synthetic data (relative) day ' + f'{day:2d}', plot_colorbar=True,
        output_path=output_path,
        fig_name=filename, dpi=300,
        outercolor=[205 / 255, 238 / 255, 251 / 255])


def create_gif_map_plot(input_data, output_dir, compartments, filename="simulation", relative=True, age_groups={0: '0-4', 1: '5-14', 2: '15-34',
                                                                                                                3: '35-59', 4: '60-79', 5: '80+'}):
    """ Creates a gif of the simulation results by calling create_plot_map for each day of the simulation and then
    storing the single plots in a temporary directory. Currently only works for the results created by the parameter study.

    :param input_data: Path to the input data. The Path should contain a file called 'Results' which contains
        the simulation results. This is the default output folder of the parameter study.
    :param output_dir: Path where the gif should be stored.
    :param filename: Name of the temporary file. (Default value = "simulation")
    :param relative: Defines if data should be scaled relative to population. (Default value = True)
    :param age_groups: Dictionary of age groups to be considered. (Default value = {0: '0-4', 1: '5-14', 2: '15-34', 3: '35-59', 4: '60-79', 5: '80+'})

    """

    files_input = {'Data set':  input_data + '/Results'}
    file_format = 'h5'

    if len(age_groups) == 6:
        filter_age = None
    else:
        filter_age = ['Group' + str(key) for key in age_groups.keys()]

    num_days = pm.extract_time_steps(
        files_input[list(files_input.keys())[0]], file_format=file_format)

    # create gif
    frames = []
    with progind.Percentage() as indicator:
        with tempfile.TemporaryDirectory() as tmpdirname:
            for day in range(0, num_days):
                create_plot_map(day, filename, files_input, tmpdirname,
                                compartments, file_format, relative, age_groups)

                image = imageio.v2.imread(
                    os.path.join(tmpdirname, filename + ".png"))
                frames.append(image)

                # Close the current figure to free up memory
                plt.close('all')
                indicator.set_progress((day+1)/num_days)

    imageio.mimsave(os.path.join(output_dir, filename + '.gif'),
                    frames,         # array of input frames
                    duration=0.2,    # duration of each frame in seconds
                    loop=0)         # optional: frames per second
