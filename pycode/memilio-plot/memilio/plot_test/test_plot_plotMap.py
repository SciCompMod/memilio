######################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors:
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
######################################################################
import datetime as dt
import os
import os.path
import unittest
from unittest.mock import MagicMock, patch

import geopandas
import numpy as np
import pandas as pd
from pyfakefs import fake_filesystem_unittest

import memilio.epidata.getPopulationData as gpd
import memilio.plot.plotMap as pm


class TestPlotMap(fake_filesystem_unittest.TestCase):

    path = '/home/plot_data'

    # strings for read, download and update data
    # be careful: not completely realistic data
    # In the json file there is one data set for each county (ensuring that the
    # completeness check evaluates to true) that is combined with some addional
    # test data sets defined below

    # Get a file object with write permission.
    here = os.path.dirname(os.path.abspath(__file__))

    files_input = {
        'Data set 1': os.path.join(
            here, 'test_data', 'Synthetic_data_counties_sh'),
        'Data set 2 (plot set 1 again...)': os.path.join(
            here, 'test_data', 'Synthetic_data_counties_sh')}
    file_format = 'json'
    # Define age groups which will be considered through filtering
    # Keep keys and values as well as its assignment constant, remove entries
    # if only part of the population should be plotted or considered, e.g., by
    # setting:
    # age_groups = {1: '5-14', 2: '15-34'}
    age_groups = {0: '0-4', 1: '5-14', 2: '15-34',
                  3: '35-59', 4: '60-79', 5: '80+'}

    def test_extract_data(self):
        filter_age = None
        i = 0
        for file in self.files_input.values():
            df = pm.extract_data(
                file, region_spec='ID_County',
                column='Synthetic data',
                date=dt.date(2021, 11, 18),
                filters={'ID_State': None,
                         'Age_RKI': filter_age},
                file_format=self.file_format)

            if i == 0:
                dfs_all = pd.DataFrame(df.iloc[:, 0])

            dfs_all[df.columns[-1] + ' ' + str(i)] = df[df.columns[-1]]
            i += 1
        data_columns = dfs_all.columns[1:]
        assert len(data_columns) == 2

    @patch('memilio.plot.plotMap.geopandas.GeoDataFrame.plot',
           return_value=(MagicMock()))
    @patch('memilio.plot.plotMap.plt', return_value=(MagicMock()))
    @patch('memilio.plot.plotMap.geopandas.read_file',
           return_value=(geopandas.GeoDataFrame(
               columns=['ARS', 'GEN', 'NUTS', 'geometry'])))
    @patch('memilio.plot.plotMap.save_interactive', return_value=0)
    def test_plot_list(
            self, mock_save_interactive, mock_read_file, mock_plt,
            mock_geopandas_plot):

        dfs_all = pd.DataFrame(columns=['ID_County', 'Count 0', 'Count 1'], data=[
                               [1001, 2, 3], [1002, 3, 4]])
        data_columns = dfs_all.columns[1:]
        pm.plot_map(
            dfs_all, scale_colors=np.array([0, 1]),
            legend=['', ''],
            title='Synthetic data (relative)', plot_colorbar=True,
            output_path=self.path,
            fig_name='customPlot', dpi=300,
            outercolor=[205 / 255, 238 / 255, 251 / 255])
        data_columns = dfs_all.columns[1:]
        assert mock_save_interactive.call_count == len(data_columns)
        assert mock_geopandas_plot.call_count == len(data_columns)


if __name__ == '__main__':
    unittest.main()
