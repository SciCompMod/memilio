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
import random
import unittest
from unittest.mock import call, patch

import numpy as np
import pandas as pd
from pyfakefs import fake_filesystem_unittest

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import getHospitalizationData as ghd
from memilio.epidata import progress_indicator


class TestGetHospitalizationData(fake_filesystem_unittest.TestCase):

    maxDiff = None

    path = '/home/HospitalizationData'

    here = os.path.dirname(os.path.abspath(__file__))
    test_data_file = os.path.join(
        here, "test_data", "TestSetFullHospitalizationData.json")

    df_test = pd.read_json(test_data_file)

    six_day_offset = [0, 0, 0, 0, 0, 0]
    array_to_test_constant = np.array(six_day_offset + [1 for i in range(15)])
    output_constant = [1/7 for i in range(21)]

    # start with a seven day offset
    output_random = [random.randint(0, 5)for i in range(200)]
    array_to_test_random = np.array([0 for i in range(len(output_random))])
    for j in range(len(array_to_test_random)):
        for i in range(7):
            try:
                array_to_test_random[i+j] += output_random[j]
            except:
                pass
    array_to_test_random[:6] = 0

    output_standart = six_day_offset + [
        3, 2, 0, 5, 0, 1, 4, 5, 4, 1, 6, 5, 4, 0, 4, 8, 1, 0, 0, 1, 0, 1, 0, 5,
        2, 20, 8, 0, 8, 9, 6, 5, 7, 1, 8, 4, 0, 0, 1, 0, 1, 0, 2, 0, 5, 6, 4,
        0, 6, 0, 4, 1, 5, 0, 4, 0, 6, 0, 4, 9, 0, 4, 7, 4, 4, 12, 2, 8, 2, 4,
        8, 2, 6, 9, 2, 4, 2, 4, 7, 0, 5, 0, 2, 0, 4, 1]

    array_to_test_standart = np.array([0 for i in range(len(output_standart))])
    for j in range(len(array_to_test_standart)):
        for i in range(7):
            try:
                array_to_test_standart[i+j] += output_standart[j]
            except:
                pass

    def setUp(self):
        self.setUpPyfakefs()
        progress_indicator.ProgressIndicator.disable_indicators(True)

    @patch('builtins.print')
    def test_divi_data_hospit_sanity_checks(self, mock_print):
        # first test
        # test empty Dataframe
        df = pd.DataFrame()
        with self.assertRaises(gd.DataError) as error:
            ghd.hospit_sanity_checks(df)
        error_message = "Download of Hospitalization Data failed. File is empty."
        self.assertEqual(str(error.exception), error_message)

        # second test
        # get dataframe with wanted columns and one additional column
        df = pd.DataFrame(
            {'Datum': [1, 2, 3],
             'Bundesland': ['A', 'B', 'C'],
             'Altersgruppe': ['00', '2+', 'alle'],
             'Bundesland_Id': [4, 5, 6],
             'a': [1, 2, 3],
             '7T_Hospitalisierung_Inzidenz': [4, 5, 6],
             '7T_Hospitalisierung_Faelle': [100, 1001, 100]})
        ghd.hospit_sanity_checks(df)
        expected_print = [call("Warning: Number of data categories changed.")]
        mock_print.assert_has_calls(expected_print)

        # third test
        # get dataframe with 6 columns but different names
        df = pd.DataFrame(
            {'date_fake': [1, 2, 3],
             '6': [6, 7, 8],
             '7': [7, 8, 9],
             '8': [8, 9, 0],
             'bundesland': [2, 3, 4],
             '9': [9, 0, 1]})
        with self.assertRaises(gd.DataError) as error:
            ghd.hospit_sanity_checks(df)
        error_message = "Error: Data categories have changed."
        self.assertEqual(str(error.exception), error_message)

    @patch('builtins.input', return_value='Y')
    @patch('memilio.epidata.getHospitalizationData.pd.read_csv',
           return_value=df_test)
    def test_get_hospitalization_data(self, mock_file, mock_in):
        # this should not raise any errors
        ghd.get_hospitalization_data(out_folder=self.path)

        # check if all files are written
        self.assertEqual(
            len(os.listdir(os.path.join(self.path, 'Germany'))), 5)

    def test_compute_hospitailzations_per_day(self):
        # test constant input array
        inputarray_c = self.array_to_test_constant
        result_c = ghd.get_hospitailzations_per_day(inputarray_c)
        self.assertEqual(list(result_c), self.output_constant)

        # test standart input
        inputarray_s = self.array_to_test_standart
        result_s = ghd.get_hospitailzations_per_day(inputarray_s)
        self.assertEqual(list(result_s), self.output_standart)

        # test random input
        inputarray_r = self.array_to_test_random
        result_r = ghd.get_hospitailzations_per_day(inputarray_r)
        # for random inputs the values in the result doesn't necessarily have to be integers.
        # Therefore we check the 7 day sum:
        control_sum = np.array([0. for i in range(len(result_r))])
        for j in range(len(control_sum)):
            for i in range(7):
                try:
                    control_sum[i+j] += result_r[j]
                except:
                    pass
        # should be equal to the input array
        np.testing.assert_array_almost_equal(
            inputarray_r[6:], control_sum[6:], 10)


if __name__ == '__main__':
    unittest.main()
