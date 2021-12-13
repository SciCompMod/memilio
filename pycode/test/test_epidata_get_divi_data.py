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
from pandas.core.indexes.base import Index
from pyfakefs import fake_filesystem_unittest
from freezegun import freeze_time
from datetime import date, datetime, time, timedelta

import os
import pandas as pd

from epidemiology.epidata import getDIVIData as gdd
from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from unittest.mock import patch, call


class TestGetDiviData(fake_filesystem_unittest.TestCase):

    maxDiff = None

    path = '/home/DiviData'

    test_df = pd.DataFrame(
        {
            'Date':
            ['2021-09-08', '2021-09-08', '2021-09-08', '2021-09-08', '2021-09-08',
             '2021-09-08', '2021-09-08', '2021-09-08', '2021-09-08', '2021-09-08',
             '2021-09-08', '2021-09-08', '2021-09-08', '2021-09-08', '2021-09-08',
             '2021-09-08'],
            'ICU':
            [16, 52, 111, 7, 432, 126, 74, 175, 208, 33, 79, 16, 11, 27, 5, 15],
            'ICU_ventilated':
            [13, 34, 63, 5, 220, 53, 38, 79, 111, 15, 53, 8, 7, 13, 2, 9],
            'ID_State': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
            'State':
            ['Schleswig-Holstein', 'Hamburg', 'Niedersachsen', 'Bremen',
                'Nordrhein-Westfalen', 'Hessen', 'Rheinland-Pfalz',
                'Baden-Württemberg', 'Bayern', 'Saarland', 'Berlin', 'Brandenburg',
                'Mecklenburg-Vorpommern', 'Sachsen', 'Sachsen-Anhalt', 'Thüringen']})

    (df_raw, df_counties, df_states, df_ger) = gdd.get_divi_data(out_folder=path)

    def setUp(self):
        self.setUpPyfakefs()

    def gdd_calls(self,text=''):
        directory = os.path.join(self.path, 'Germany/')
        gdd_calls= [
        call('Information: Data has been written to',
            os.path.join(directory, 'FullData_DIVI.json')),
        call('Information: Data has been written to',
            os.path.join(directory, 'county_divi'+text+'.json')),
        call(
            'Information: Data has been written to',
            os.path.join(directory, 'state_divi'+text+'.json')),
        call(
            'Information: Data has been written to',
            os.path.join(directory, 'germany_divi'+text+'.json'))]
        return gdd_calls 

    def test_cut_of_dates(self):
        # test if only dates from 08-09-2021 are returned
        df_state_testdate = gdd.cut_of_dates(
            self.df_states, date(2021, 9, 8),
            date(2021, 9, 8))
        pd.testing.assert_frame_equal(self.test_df, df_state_testdate)

    @patch('epidemiology.epidata.getDIVIData.pd.read_json')
    @patch('epidemiology.epidata.getDataIntoPandasDataFrame.loadCsv')
    def test_exit_strings(self, mocklcsv, mockrjson):

        # test read_data Error call if json file is not found
        mockrjson.side_effect = ValueError
        with self.assertRaises(SystemExit) as cm:
            gdd.get_divi_data(read_data=True, out_folder=self.path)
        file_in = os.path.join(self.path, "Germany/FullData_DIVI.json")
        exit_string = "Error: The file: " + file_in + " does not exist. "\
            "Call program without -r flag to get it."
        self.assertEqual(cm.exception.code, exit_string)

        # test loadCsv Error if file can't be downloaded
        mocklcsv.side_effect = Exception
        with self.assertRaises(SystemExit) as cm:
            gdd.get_divi_data(read_data=False)
        exit_string = "Error: Download link for Divi data has changed."
        self.assertEqual(cm.exception.code, exit_string)

    @patch('epidemiology.epidata.getDataIntoPandasDataFrame.loadCsv')
    def test_df_empty(self, mocklcsv):

        # new test function because of the new mock value
        # test Error for empty returned dataframe
        mocklcsv.value = pd.DataFrame()
        with self.assertRaises(SystemExit) as cm:
            gdd.get_divi_data(read_data=False)
        exit_string = "Something went wrong, dataframe is empty."
        self.assertEqual(cm.exception.code, exit_string)

    @patch('builtins.print')
    @patch('epidemiology.epidata.getDataIntoPandasDataFrame.loadCsv', return_value = df_raw.copy())
    def test_det_divi_data_prints(self, mocklcsv, mock_print):

        # case with start_date before 2020-04-24
        gdd.get_divi_data(out_folder=self.path, start_date=date(2020, 1, 1))
        expected_call = [
            call(
                'Warning: First data available on 2020-04-24. You asked for 2020-01-01.')]
        gdd_calls = self.gdd_calls(text='')
        expected_calls = expected_call + gdd_calls
        mock_print.assert_has_calls(expected_calls)

    @patch('builtins.print')
    @patch('epidemiology.epidata.getDIVIData.pd.read_json', return_value= df_raw.copy())
    def test_gdd_ma(self, mockrjson2, mock_print):

        # test case with moving average
        gdd.get_divi_data(
            read_data=True, out_folder=self.path, moving_average=3)
        expected_calls = self.gdd_calls('_ma3')
        mock_print.assert_has_calls(expected_calls)

    @patch('builtins.print')
    @patch('epidemiology.epidata.getDIVIData.pd.read_json', return_value= df_raw.copy())
    def test_gdd_all_dates(self, mockrjson2, mock_print):
        # test case with impute dates is True
        gdd.get_divi_data(
            read_data=True, out_folder=self.path, impute_dates=True)
        expected_calls = self.gdd_calls(text='_all_dates')
        mock_print.assert_has_calls(expected_calls)

    @patch('epidemiology.epidata.getDIVIData.pd.read_json', return_value=test_df.copy())
    def test_gdd_sanity_checks(self, mockrjson3):

        # first test
        # get random dataframe
        df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
        with self.assertRaises(SystemExit) as cm:
            gdd.gdd_sanity_checks(df)
        exit_string = "Error: Number of data categories changed."
        self.assertEqual(cm.exception.code, exit_string)

        # second test
        # get dataframe with 11 columns but different names
        df = pd.DataFrame(
            {'date_fake': [1, 2, 3],
             '6': [6, 7, 8],
             '7': [7, 8, 9],
             '8': [8, 9, 0],
             'bundesland': [2, 3, 4],
             '9': [9, 0, 1],
             '10': [0, 1, 2],
             '11': [0, 0, 0],
             'gemeindeschluessel_fake': [3, 4, 5],
             'faelle_covid_aktuell_fake': [4, 5, 6],
             'faelle_covid_aktuell_invasiv_beatmet': [5, 6, 7]})
        with self.assertRaises(SystemExit) as cm:
            gdd.gdd_sanity_checks(df)
        exit_string = "Error: Data categories have changed."
        self.assertEqual(cm.exception.code, exit_string)

        # third test
        # get dataframe with 11 columns and same headers but only 3 rows
        df = pd.DataFrame(
            {'date': [1, 2, 3],
             '6': [6, 7, 8],
             '7': [7, 8, 9],
             '8': [8, 9, 0],
             'bundesland': [2, 3, 4],
             '9': [9, 0, 1],
             '10': [0, 1, 2],
             '11': [0, 0, 0],
             'gemeindeschluessel': [3, 4, 5],
             'faelle_covid_aktuell': [4, 5, 6],
             'faelle_covid_aktuell_invasiv_beatmet': [5, 6, 7]})
        with self.assertRaises(SystemExit) as cm:
            gdd.gdd_sanity_checks(df)
        exit_string = "Error: unexpected length of dataframe."
        self.assertEqual(cm.exception.code, exit_string)

        # test if it works in main
        with self.assertRaises(SystemExit) as cm:
            gdd.get_divi_data(read_data=True, out_folder=self.path)
        exit_string = "Error: Number of data categories changed."
        self.assertEqual(cm.exception.code, exit_string)


if __name__ == '__main__':
    unittest.main()
