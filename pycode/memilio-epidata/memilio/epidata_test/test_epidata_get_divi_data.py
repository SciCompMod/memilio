#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Patrick Lenz, Annette Lutz
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
from datetime import date

import os
import pandas as pd

from memilio.epidata import getDIVIData as gdd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import modifyDataframeSeries as mDfS
from unittest.mock import patch, call


class TestGetDiviData(fake_filesystem_unittest.TestCase):

    maxDiff = None

    path = '/home/DiviData'

    here = os.path.dirname(os.path.abspath(__file__))
    test_data_file = os.path.join(
        here, "test_data", "TestSetFullDIVIData.json")

    df_test = pd.read_json(test_data_file)

    df_test_error = pd.DataFrame(
        {'dates': ['2021-10-30', '2021-11-01'],
         'numbers': [542, 22],
         'strings': ['one', 'two']})

    def setUp(self):
        self.setUpPyfakefs()

    def gdd_calls(self, text=''):
        directory = os.path.join(self.path, 'Germany/')
        gdd_calls = [
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

    @patch('memilio.epidata.getDIVIData.pd.read_json')
    @patch('memilio.epidata.getDataIntoPandasDataFrame.loadCsv')
    def test_exit_strings(self, mocklcsv, mockrjson):

        # test read_data Error call if json file is not found
        mockrjson.side_effect = ValueError
        with self.assertRaises(FileNotFoundError) as error:
            gdd.get_divi_data(self.path, read_data=True)
        file_in = os.path.join(self.path, "Germany/FullData_DIVI.json")
        error_message = "Error: The file: " + file_in + " does not exist. "\
            "Call program without -r flag to get it."
        self.assertEqual(str(error.exception), error_message)

        # test loadCsv Error if file can't be downloaded
        mocklcsv.side_effect = Exception
        with self.assertRaises(FileNotFoundError) as error:
            gdd.get_divi_data(self.path, read_data=False)
        error_message = "Error: Download link for Divi data has changed."
        self.assertEqual(str(error.exception), error_message)

    @patch('memilio.epidata.getDataIntoPandasDataFrame.loadCsv')
    def test_df_empty(self, mocklcsv):

        # new test function because of the new mock value
        # test Error for empty returned dataframe
        mocklcsv.value = pd.DataFrame()
        with self.assertRaises(gd.DataError) as error:
            gdd.get_divi_data(self.path, read_data=False)
        error_message = "Something went wrong, dataframe is empty."
        self.assertEqual(str(error.exception), error_message)

    @patch('memilio.epidata.getDIVIData.divi_data_sanity_checks')
    @patch('memilio.epidata.getDIVIData.gd.loadCsv')
    @patch('builtins.print')
    def test_get_divi_data_prints(self, mock_print, mock_csv, mock_san):
        mock_csv.return_value = self.df_test
        # case with start_date before 2020-04-24
        gdd.get_divi_data(self.path, start_date=date(2020, 1, 1))
        expected_call = [
            call(
                'Warning: First data available on 2020-04-24. You asked for 2020-01-01.')]
        gdd_calls = self.gdd_calls(text='')
        expected_calls = expected_call + gdd_calls
        mock_print.assert_has_calls(expected_calls)
        mock_san.assert_has_calls([call(self.df_test)])

    @patch('memilio.epidata.getDIVIData.divi_data_sanity_checks')
    @patch('memilio.epidata.getDIVIData.gd.loadCsv')
    @patch('builtins.print')
    def test_get_divi_data(self, mock_print, mock_csv, mock_san):
        mock_csv.return_value = self.df_test
        # test case with standard parameters
        (df, df_county, df_states, df_ger) = gdd.get_divi_data(self.path)
        mock_san.assert_has_calls([call(self.df_test)])
        pd.testing.assert_frame_equal(df, self.df_test)
        self.assertEqual(
            df_ger[df_ger["Date"] == "2021-05-29"]["ICU"].item(), 8)
        self.assertEqual(
            df_ger[df_ger["Date"] == "2021-05-29"]
            ["ICU_ventilated"].item(),
            2)
        self.assertEqual(
            df_states
            [(df_states["Date"] == "2021-05-29") &
             (df_states["ID_State"] == 3)]["ICU"].item(),
            6)
        self.assertEqual(
            df_states
            [(df_states["Date"] == "2021-05-29") & (df_states["ID_State"] == 3)]
            ["ICU_ventilated"].item(),
            2)
        self.assertEqual(
            df_county
            [(df_county["Date"] == "2021-05-29") &
             (df_county["ID_County"] == 3151)]["ICU"].item(),
            5)
        self.assertEqual(
            df_county
            [(df_county["Date"] == "2021-05-29") &
             (df_county["ID_County"] == 3151)]["ICU_ventilated"].item(),
            2)

    @patch('memilio.epidata.getDIVIData.divi_data_sanity_checks')
    @patch('memilio.epidata.getDIVIData.gd.loadCsv')
    @patch('builtins.print')
    def test_gdd_ma(self, mock_print, mock_csv, mock_san):
        mock_csv.return_value = self.df_test
        # test case with moving average
        (df, df_county, df_states, df_ger) = gdd.get_divi_data(
            self.path, moving_average=3)
        mock_san.assert_has_calls([call(self.df_test)])
        pd.testing.assert_frame_equal(df, self.df_test)
        self.assertAlmostEqual(
            df_ger[df_ger["Date"] == "2021-05-27"]["ICU"].item(), 143+(1/3.0))
        self.assertAlmostEqual(
            df_ger[df_ger["Date"] == "2021-05-27"]
            ["ICU_ventilated"].item(),
            108+(2/3.0))
        self.assertAlmostEqual(
            df_states
            [(df_states["Date"] == "2021-05-27") &
             (df_states["ID_State"] == 3)]["ICU"].item(),
            16+(1/3.0))
        self.assertAlmostEqual(
            df_states
            [(df_states["Date"] == "2021-05-27") & (df_states["ID_State"] == 3)]
            ["ICU_ventilated"].item(),
            10+(2/3.0))
        self.assertAlmostEqual(
            df_county
            [(df_county["Date"] == "2021-05-27") &
             (df_county["ID_County"] == 1001)]["ICU"].item(),
            18+(1/3.0))
        self.assertAlmostEqual(
            df_county
            [(df_county["Date"] == "2021-05-27") &
             (df_county["ID_County"] == 1001)]["ICU_ventilated"].item(),
            13+(1/3.0))

    @patch('memilio.epidata.getDIVIData.divi_data_sanity_checks')
    @patch('memilio.epidata.getDIVIData.gd.loadCsv')
    @patch('builtins.print')
    def test_gdd_all_dates(self, mock_print, mock_csv, mock_san):
        mock_csv.return_value = self.df_test
        # test case with impute dates is True
        (df, df_county, df_states, df_ger) = gdd.get_divi_data(
            self.path, impute_dates=True)
        mock_san.assert_has_calls([call(self.df_test)])
        pd.testing.assert_frame_equal(df, self.df_test)
        self.assertEqual(
            df_ger[df_ger["Date"] == "2021-05-26"]["ICU"].item(), 119)
        self.assertEqual(
            df_ger[df_ger["Date"] == "2021-05-26"]
            ["ICU_ventilated"].item(),
            95)
        self.assertEqual(
            df_states
            [(df_states["Date"] == "2021-05-27") &
             (df_states["ID_State"] == 3)]["ICU"].item(),
            20)
        self.assertEqual(
            df_states
            [(df_states["Date"] == "2021-05-27") & (df_states["ID_State"] == 1)]
            ["ICU_ventilated"].item(),
            9)
        self.assertEqual(
            df_county
            [(df_county["Date"] == "2021-05-27") &
             (df_county["ID_County"] == 3151)]["ICU"].item(),
            0)
        self.assertEqual(
            df_county
            [(df_county["Date"] == "2021-05-27") &
             (df_county["ID_County"] == 1051)]["ICU_ventilated"].item(),
            6)

    @patch('memilio.epidata.getDIVIData.pd.read_json',
           return_value=df_test_error)
    def test_divi_data_sanity_checks(self, mockrjson3):

        # first test
        # get random dataframe
        df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
        with self.assertRaises(gd.DataError) as error:
            gdd.divi_data_sanity_checks(df)
        error_message = "Error: Number of data categories changed."
        self.assertEqual(str(error.exception), error_message)

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
        with self.assertRaises(gd.DataError) as error:
            gdd.divi_data_sanity_checks(df)
        error_message = "Error: Data categories have changed."
        self.assertEqual(str(error.exception), error_message)

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
        with self.assertRaises(gd.DataError) as error:
            gdd.divi_data_sanity_checks(df)
        error_message = "Error: unexpected length of dataframe."
        self.assertEqual(str(error.exception), error_message)

        # test if it works in main
        with self.assertRaises(gd.DataError) as error:
            gdd.get_divi_data(self.path, read_data=True)
        error_message = "Error: Number of data categories changed."
        self.assertEqual(str(error.exception), error_message)


if __name__ == '__main__':
    unittest.main()
