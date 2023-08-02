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
import unittest
<<<<<<< HEAD
=======
import json
import pandas as pd

from unittest.mock import patch
>>>>>>> e13e4c705c931ad490aea943f701f918bdc8803b
from pyfakefs import fake_filesystem_unittest

import os
import pandas as pd
import numpy as np

from memilio.epidata import getPopulationData as gpd
<<<<<<< HEAD
from memilio.epidata import defaultDict as dd
from unittest.mock import patch
=======
from memilio.epidata import progress_indicator
>>>>>>> e13e4c705c931ad490aea943f701f918bdc8803b


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
<<<<<<< HEAD
    def test_get_new_counties(self):
        test = gpd.get_new_counties(self.test_old_counties)
        self.assertTrue(np.array_equal(test, self.test_new_counties))
    
    @patch('memilio.epidata.getPopulationData.load_population_data',
           return_value=(test_counties, test_zensus, test_reg_key))
    def test_get_population(self, mock_data):
=======
        progress_indicator.ProgressIndicator.disable_indicators(True)

    def test_export_population_data(self):

        result_df = gpd.export_population_dataframe(
            self.df_pop, self.path, 'json', True)
        # check if one file is written
        self.assertEqual(len(os.listdir(self.path)), 1)
>>>>>>> e13e4c705c931ad490aea943f701f918bdc8803b

        # test result
        self.assertEqual(result_df[result_df.ID_County == 16056].empty, True)
        self.assertEqual(
            result_df[result_df.ID_County == 16063]['Population'].values, 81222)
        self.assertEqual(result_df.shape, (400, 13))
        self.assertListEqual(list(result_df.columns), ['ID_County', 'Population', '<3 years', '3-5 years', '6-14 years', '15-17 years',
                                                       '18-24 years', '25-29 years', '30-39 years', '40-49 years',
                                                       '50-64 years', '65-74 years', '>74 years'])

<<<<<<< HEAD
        # add seven to the number of test_counties as this is the current workaround to add counties which are
        # not mentioned in old data sets
        test_df = pd.read_json(os.path.join(
            self.path, 'Germany/', 'county_current_population_dim' + str(len(self.test_counties)+7) +'.json'))
        test_df = test_df.drop(
            test_df[test_df[dd.EngEng['population']] == 0].index)
        pd.testing.assert_frame_equal(
            test_df.astype('int64'), self.test_current_population_result)

    @patch('memilio.epidata.getPopulationData.load_population_data',
           return_value=(test_counties, test_zensus, test_reg_key))
    def test_popul_split_gender(self, mock_data):

        test_df = gpd.get_population_data(
            read_data=False, file_format='json', out_folder=self.path,
            no_raw=False, split_gender=True, merge_eisenach=False)

        test_df = test_df.drop(
            test_df[test_df[dd.EngEng['population']] == 0].index)
        pd.testing.assert_frame_equal(
            test_df.astype('int64'), self.test_current_population_gender_result)

    @ patch('pandas.read_excel', return_value=test_counties)
    @ patch('pandas.read_excel', return_value=test_reg_key)
    @ patch('memilio.epidata.getDataIntoPandasDataFrame.loadCsv',
            return_value=test_zensus)
    def test_load_population_data(
            self, mock_read_excel1, mock_read_excel2, mock_read_csv):
=======
    def test_read_population_data(self):
>>>>>>> e13e4c705c931ad490aea943f701f918bdc8803b

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
