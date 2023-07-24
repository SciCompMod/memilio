#############################################################################
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
#############################################################################
import os
import unittest
import json
import pandas as pd

from unittest.mock import patch
from pyfakefs import fake_filesystem_unittest

from memilio.epidata import defaultDict as dd
from memilio.epidata import getPopulationData as gpd
from memilio.epidata import progress_indicator


class Test_getPopulationData(fake_filesystem_unittest.TestCase):

    path = '/home/Population_Data'

    here = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(
        here, 'test_data', 'TestSetPopulationExport.json')
    with open(filename) as file_object:
        df_pop = pd.DataFrame(json.load(file_object))

    filename = os.path.join(
        here, 'test_data', 'TestSetPopulationRaw.json')
    with open(filename) as file_object:
        df_pop_raw = pd.DataFrame(json.load(file_object))

    def setUp(self):
        self.setUpPyfakefs()
        progress_indicator.ProgressIndicator.disable_indicators(True)

    def test_export_population_data(self):

        result_df = gpd.export_population_dataframe(
            self.df_pop, self.path, 'json', True)
        # check if one file is written
        self.assertEqual(len(os.listdir(self.path)), 1)

        # test result
        self.assertEqual(result_df[result_df.ID_County == 16056].empty, True)
        self.assertEqual(
            result_df[result_df.ID_County == 16063]['Population'].values, 81222)
        self.assertEqual(result_df.shape, (400, 13))
        self.assertListEqual(list(result_df.columns), ['ID_County', 'Population', '<3 years', '3-5 years', '6-14 years', '15-17 years',
                                                       '18-24 years', '25-29 years', '30-39 years', '40-49 years',
                                                       '50-64 years', '65-74 years', '>74 years'])

    def test_read_population_data(self):

        directory = os.path.join(self.path, 'Germany/')

        # test file not found
        with self.assertRaises(FileNotFoundError) as error:
            df = gpd.read_population_data(
                username='', password='', read_data=True, directory=directory)

    @patch('memilio.epidata.getPopulationData.read_population_data',
           return_value=df_pop_raw)
    @patch('memilio.epidata.getPopulationData.assign_population_data', return_value=df_pop)
    @patch('memilio.epidata.getPopulationData.test_total_population')
    def test_get_population_data_full(self, mock_test, mock_export, mock_download):
        # should not raise any errors
        gpd.get_population_data(out_folder=self.path)


if __name__ == '__main__':
    unittest.main()
