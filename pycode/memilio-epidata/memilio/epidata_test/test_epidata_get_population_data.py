#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Patrick Lenz
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
import configparser
import json
import pandas as pd

from unittest.mock import patch
from pyfakefs import fake_filesystem_unittest

from memilio.epidata import getPopulationData as gpd
from memilio.epidata import progress_indicator


class Test_getPopulationData(fake_filesystem_unittest.TestCase):

    path = '/home/Population_Data'

    config_file_name = 'CredentialsRegio.ini'
    test_username = 'username_test'
    test_password = 'password_test'

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

    @patch('builtins.input', return_value=test_username)
    @patch('getpass.getpass', return_value=test_password)
    @patch('memilio.epidata.getDataIntoPandasDataFrame.user_choice', return_value=True)
    @patch('os.path.abspath', return_value='')
    @patch('memilio.epidata.getPopulationData.read_population_data', return_value=df_pop_raw)
    @patch('memilio.epidata.getPopulationData.assign_population_data', return_value=df_pop)
    @patch('memilio.epidata.getPopulationData.test_total_population')
    def test_config_write(self, mock_test, mock_export, mock_raw, mock_path, mock_choice, mock_pw, mock_un):
        # username and password should be written into the config file.
        # The download and assigning to counties of the population data is mocked.
        gpd.get_population_data(username=None, password=None)
        # Check if the file is written.
        self.assertTrue(self.config_file_name in os.listdir(os.getcwd()))
        # Check content of the file.
        # Read file.
        parser = configparser.ConfigParser()
        parser.read(os.path.join(os.getcwd(), self.config_file_name))
        # Test content.
        self.assertEqual(parser['CREDENTIALS']['Username'], self.test_username)
        self.assertEqual(parser['CREDENTIALS']['Password'], self.test_password)

    @patch('os.path.abspath', return_value='')
    @patch('memilio.epidata.getPopulationData.read_population_data', return_value=df_pop_raw)
    @patch('memilio.epidata.getPopulationData.assign_population_data', return_value=df_pop)
    @patch('memilio.epidata.getPopulationData.test_total_population')
    def test_config_read(self, mock_test, mock_export, mock_read, mock_path):
        # File should not exist yet.
        self.assertFalse(self.config_file_name in os.listdir(os.getcwd()))
        # Create config file.
        string = '[CREDENTIALS]\nUsername = ' + \
            self.test_username+'\nPassword = '+self.test_password
        path = os.path.join(os.getcwd(), self.config_file_name)
        with open(path, 'w+') as file:
            file.write(string)
        # Check if the file is written.
        self.assertTrue(self.config_file_name in os.listdir(os.getcwd()))
        # The download and assigning to counties of the population data is mocked.
        gpd.get_population_data(
            username=None, password=None, read_data=False, out_folder=self.path)
        # The file exist in the directory (mocked) and the credentials should be read.
        mock_read.assert_called_with(
            self.test_username, self.test_password, False, os.path.join(self.path, 'Germany'))


if __name__ == '__main__':
    unittest.main()
