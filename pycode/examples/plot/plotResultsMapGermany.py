#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Martin J. Kuehn
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

import numpy as np
import pandas as pd

import memilio.epidata.getPopulationData as gpd
import memilio.plot.plotMap as pm

import matplotlib

import os

matplotlib.use('Agg')

if __name__ == '__main__':

    # Get the current script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    files_input = {
        'Data set 1': os.path.join(script_dir, 'Synthetic_data_counties_sh'),
        'Data set 2 (plot set 1 again...)': os.path.join(script_dir, 'Synthetic_data_counties_sh')
    }

    file_format = 'json'
    age_groups = {0: '0-4', 1: '5-14', 2: '15-34', 3: '35-59', 4: '60-79', 5: '80+'}
    
    # Set filter_age depending on which age groups you're interested in
    if len(age_groups) == 6:
        filter_age = None
    else:
        if file_format == 'json':
            filter_age = [val for val in age_groups.values()]
        else:
            filter_age = ['Group' + str(key) for key in age_groups.keys()]

    relative = True

    i = 0
    for file in files_input.values():
        # Extract the data using the corrected file path
        df = pm.extract_data(
            file, region_spec='ID_County',
            column='Synthetic data',
            date=dt.date(2021, 11, 18),
            filters={'ID_State': None,
                     'Age_RKI': filter_age},
            file_format=file_format)

        if relative:
            try:
                population = pd.read_json(
                    os.path.join(script_dir, 'data/pydata/Germany/county_current_population.json'))
            except (ValueError, FileNotFoundError):
                print("Population data was not found. Download it from the internet.")
                population = gpd.get_population_data(
                    read_data=False, file_format=file_format,
                    out_folder=os.path.join(script_dir, 'data/pydata/Germany/'), no_raw=True,
                    merge_eisenach=True)

            age_group_values = list(age_groups.values())
            age_group_values[-1] = age_group_values[-1].replace('80+', '>79')
            # Scale data
            df = pm.scale_dataframe_relative(df, age_group_values, population)

        if i == 0:
            dfs_all = pd.DataFrame(df.iloc[:, 0])

        dfs_all[df.columns[-1] + ' ' + str(i)] = df[df.columns[-1]]
        i += 1

    min_val = dfs_all[dfs_all.columns[1:]].min().min()
    max_val = dfs_all[dfs_all.columns[1:]].max().max()

    # Plot the map using the corrected paths
    pm.plot_map(
        dfs_all, scale_colors=np.array([min_val, max_val]),
        legend=['', ''],
        title='Synthetic data (relative)', plot_colorbar=True,
        output_path=script_dir,
        fig_name='customPlot', dpi=300,
        outercolor=[205 / 255, 238 / 255, 251 / 255])
