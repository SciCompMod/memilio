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
from pyfakefs import fake_filesystem_unittest

import os
import pandas as pd
import numpy as np

from memilio.epidata import getPopulationData as gpd
from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from unittest.mock import patch


class Test_getPopulationData(fake_filesystem_unittest.TestCase):

    path = '/home/Population_Data'

    test_old_counties = np.zeros((18, 2))
    test_old_counties[:, 0] = [3152, 3156, 13056, 13002, 13055, 13052, 13051,
                               13053, 13061, 13005, 13057, 13006, 13058, 13059, 13062, 13001, 13054, 13060]
    test_old_counties[:, 1] = np.arange(len(test_old_counties))

    test_new_counties = np.zeros((7, 2))
    test_new_counties[:, 0] = [3159, 13071, 13072, 13073, 13074, 13075, 13076]
    test_new_counties[:, 1] = [1, 14, 13, 27, 23, 42, 33]

    data = np.zeros((5, 30))
    data[:, 0] = np.arange(1, 6)
    for i in range(len(data)):
        data[i, 3] = 22*(i+1)
        data[i, 4] = 11*(i+1)
        data[i, 5:-2] = 1*(i+1)
        data[i, 16] = 11*(i+1)

    test_zensus = pd.DataFrame(
        data,
        columns=["FID", "DES", "Name", "EWZ", "Gesamt_Maennlich",
                 'M_Unter_3', 'M_3_bis_5', 'M_6_bis_14', 'M_15_bis_17',
                 'M_18_bis_24', 'M_25_bis_29', 'M_30_bis_39', 'M_40_bis_49',
                 'M_50_bis_64', 'M_65_bis_74', 'M_75_und_aelter',
                 "Gesamt_Weiblich", 'W_Unter_3', 'W_3_bis_5', 'W_6_bis_14',
                 'W_15_bis_17', 'W_18_bis_24', 'W_25_bis_29', 'W_30_bis_39',
                 'W_40_bis_49', 'W_50_bis_64', 'W_65_bis_74',
                 'W_75_und_aelter', 'SHAPE_Length', 'SHAPE_Area'])
    test_zensus["DES"] = "Kreis"
    test_zensus["Name"] = ["Hogwarts", "Narnia",
                           "MittelErde", "Westeros", "Wakanda"]

    data = np.zeros((5, 3))
    data[:, 0] = [1001, 1002, 1003, 1004, 1005]
    data[:, 2] = [(x+1)*22/1000 for x in range(len(data))]

    test_reg_key = pd.DataFrame(data, columns=['AGS', 'NAME', 'Zensus_EWZ'])
    test_reg_key['NAME'] = ["Hogwarts", "Narnia",
                            "MittelErde", "Westeros", "Wakanda"]

    data = np.zeros((5, 2))
    data[:, 0] = [1001, 1002, 1003, 1004, 1005]
    data[:, 1] = [(x+1)*44 for x in range(len(data))]

    test_counties = pd.DataFrame(
        data, columns=['Schlüssel-nummer', 'Bevölkerung2)'])

    columns = [
        'ID_County', 'Population', '<3 years', '3-5 years', '6-14 years',
        '15-17 years', '18-24 years', '25-29 years', '30-39 years',
        '40-49 years', '50-64 years', '65-74 years', '>74 years']

    data = np.zeros((5, len(columns)))
    for i in range(len(data)):
        data[i, 0] = 1001 + i
        data[i, 1] = 22*(i+1)
        data[i, 2:] = 2*(i+1)
    test_population_result = pd.DataFrame(data, columns=columns)
    test_population_result = test_population_result.astype('int64')

    data = np.zeros((5, len(columns)))
    for i in range(len(data)):
        data[i, 0] = 1001 + i
        data[i, 1] = 44 * (i+1)
        data[i, 2:] = 4 * (i+1)

    test_current_population_result = pd.DataFrame(data, columns=columns)
    test_current_population_result = test_current_population_result.astype(
        'int64')

    columns_gender = [
        'ID_County', 'Population', 'M <3 years', 'M 3-5 years', 'M 6-14 years',
        'M 15-17 years', 'M 18-24 years', 'M 25-29 years', 'M 30-39 years',
        'M 40-49 years', 'M 50-64 years', 'M 65-74 years', 'M >74 years',
        'F <3 years', 'F 3-5 years', 'F 6-14 years', 'F 15-17 years',
        'F 18-24 years', 'F 25-29 years', 'F 30-39 years', 'F 40-49 years',
        'F 50-64 years', 'F 65-74 years', 'F >74 years']

    data = np.zeros((5, len(columns_gender)))
    for i in range(len(data)):
        data[i, 0] = 1001 + i
        data[i, 1] = 44 * (i+1)
        data[i, 2:] = 2 * (i+1)
    test_current_population_gender_result = pd.DataFrame(
        data, columns=columns_gender)
    test_current_population_gender_result = test_current_population_gender_result.astype(
        'int64')

    def setUp(self):
        self.setUpPyfakefs()
    def test_get_new_counties(self):
        test = gpd.get_new_counties(self.test_old_counties)
        self.assertTrue(np.array_equal(test, self.test_new_counties))
    
    @patch('memilio.epidata.getPopulationData.load_population_data',
           return_value=(test_counties, test_zensus, test_reg_key))
    def test_get_population(self, mock_data):

        gpd.get_population_data(
            self.path, read_data=True, file_format='json', 
            no_raw=False, split_gender=False, merge_eisenach=False)

        test_df = pd.read_json(os.path.join(
            self.path, 'Germany/', 'county_current_population_dim401.json'))
        test_df = test_df.drop(
            test_df[test_df[dd.EngEng['population']] == 0].index)
        pd.testing.assert_frame_equal(
            test_df.astype('int64'), self.test_current_population_result)

    @patch('memilio.epidata.getPopulationData.load_population_data',
           return_value=(test_counties, test_zensus, test_reg_key))
    def test_popul_split_gender(self, mock_data):

        test_df = gpd.get_population_data(
            self.path, read_data=False, file_format='json',
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

        directory = os.path.join(self.path, 'Germany/')

        counties_write, zensus_write, reg_key_write = gpd.load_population_data(
            self.path, read_data=False)
        self.assertEqual(len(os.listdir(directory)), 3)

        counties_read, zensus_read, reg_key_read = gpd.load_population_data(
            self.path, read_data=True)

        pd.testing.assert_frame_equal(
            counties_read, counties_write, check_dtype=False)
        pd.testing.assert_frame_equal(
            reg_key_read, reg_key_write, check_dtype=False)
        pd.testing.assert_frame_equal(
            zensus_read, zensus_write, check_dtype=False)

        # TODO: How to test hdf5 export?

    @ patch('memilio.epidata.getPopulationData.gd.loadCsv')
    @ patch('memilio.epidata.getPopulationData.gd.loadExcel')
    def test_errors(self, mocklexcel, mocklcsv):
        mocklexcel.side_effect = ValueError
        mocklcsv.side_effect = ValueError

        with self.assertRaises(FileNotFoundError) as error:
            gpd.load_population_data(self.path)
        error_message = "Error: The counties file does not exist."
        self.assertEqual(str(error.exception), error_message)

        mocklexcel.side_effect = None
        mocklexcel.return_value = self.test_counties.copy()
        with self.assertRaises(FileNotFoundError) as error:
            gpd.load_population_data(self.path)
        error_message = "Error: The zensus file does not exist."
        self.assertEqual(str(error.exception), error_message)


if __name__ == '__main__':
    unittest.main()
