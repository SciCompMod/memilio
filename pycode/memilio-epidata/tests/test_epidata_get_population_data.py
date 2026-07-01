#############################################################################
# Copyright (C) 2020-2026 MEmilio
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
import io
import unittest
import json
import pandas as pd

from unittest.mock import patch, Mock
from pyfakefs import fake_filesystem_unittest

from memilio.epidata import getPopulationData as gpd


class Test_getPopulationData(fake_filesystem_unittest.TestCase):
    """ """

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
        """ """
        self.setUpPyfakefs()

    def test_export_population_data(self):
        """ """

        result_df = gpd.export_population_dataframe(
            self.df_pop, self.path, 'json', True, 'newest')
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

    @patch('memilio.epidata.getPopulationData.read_population_data',
           return_value=(df_pop_raw, None))
    @patch('memilio.epidata.getPopulationData.assign_population_data', return_value=df_pop)
    @patch('memilio.epidata.getPopulationData.test_total_population')
    def test_get_population_data_full(self, mock_test, mock_assign, mock_download):
        """

        :param mock_test: 
        :param mock_assign: 
        :param mock_download: 

        """
        # should not raise any errors
        gpd.get_population_data(out_folder=self.path)
        # test ref_year
        gpd.get_population_data(out_folder=self.path, ref_year=2013)

    @patch('io.StringIO')
    @patch('pandas.read_csv')
    @patch('requests.get')
    def test_read_population_data(self, mock_req, mock_pd, mock_io):
        """

        :param mock_req: 
        :param mock_pd: 
        :param mock_io: 

        """
        # Test a year that does not have population Data. Function should throw a
        # warning, and download the newest data (ref_year = None)
        test_year = 2000
        # Create a mock response object for requests.get()
        mock_response = Mock()
        mock_response.text = "mocked csv data"
        mock_req.return_value = mock_response
        # Mock pandas.read_csv to raise a ParserError on the first call and return a DataFrame on the second
        mock_pd.side_effect = [pd.errors.ParserError, pd.DataFrame()]
        # Mock io.StringIO to return the StringIO object for pandas.read_csv
        mock_io.return_value = io.StringIO("mocked csv data")
        df, year = gpd.read_population_data(test_year)
        # Test results, ref_year should now be None
        self.assertTrue(df.empty)  # from mock
        self.assertIsNone(year)


if __name__ == '__main__':
    unittest.main()
