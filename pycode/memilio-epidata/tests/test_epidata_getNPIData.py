######################################################################
# Copyright (C) 2020-2026 MEmilio
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
import unittest
from pyfakefs import fake_filesystem_unittest
from unittest.mock import patch

import os
import pandas as pd
import numpy as np

from datetime import date

from memilio.epidata import getNPIData as gnd
from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gd


class TestGetNPIData(fake_filesystem_unittest.TestCase):
    """ """
    maxDiff = None

    path = '/home/NPIData'

    # load test data from test_data folder
    here = os.path.dirname(os.path.abspath(__file__))
    df_npis_old = pd.read_json(os.path.join(
        here, 'test_data', 'TestSetNPIsUnterkategorien.json'))
    df_npis_desc = pd.read_json(os.path.join(
        here, 'test_data', 'TestSetNPIsDescription.json'))
    df_npis_combinations_pre = pd.read_json(os.path.join(
        here, 'test_data', 'TestSetNPIsCombinations.json'))
    df_cases = pd.read_json(os.path.join(
        here, 'test_data', 'TestSetNPIsCaseData.json'))
    df_pop = pd.read_json(os.path.join(
        here, 'test_data', 'TestSetNPIsPopulationData.json'))

    df_npis_old_renamed = df_npis_old.rename(dd.GerEng, axis=1, inplace=False)

    threshold = 1.5
    incid = pd.Series(
        [1, 2, 3, 4, 1, 1, 0, 5, 0, 1, 2, 3, 5, 1, 0, 2, 4, 5, 6])
    active_11 = pd.Series(
        [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1])
    active_32 = pd.Series(
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0])
    active_35 = pd.Series(
        [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

    incid_start_above_threshold = pd.Series(
        [2, 2, 3, 4, 1, 1, 0, 5, 0, 1, 2, 3, 5, 1, 0, 2, 4, 5, 6])
    active_start_above_threshold_11 = pd.Series(
        [0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1])
    active_start_above_threshold_32 = pd.Series(
        [0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0])
    active_start_above_threshold_35 = pd.Series(
        [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

    codes_to_correct = [    # should be renamed:
        'M04_2', 'M04_2', 'M04_2', 'M04_2', 'M04_2',
        'M04_3', 'M04_3', 'M04_3', 'M04_3', 'M04_3',
        'M04_4', 'M04_4', 'M04_4', 'M04_4', 'M04_4',
        'M04_5', 'M04_5', 'M04_5', 'M04_5', 'M04_5',
        'M05_1', 'M05_1', 'M05_1', 'M05_1', 'M05_1', 'M05_1', 'M05_1',
        'M05_2', 'M05_2', 'M05_2', 'M05_2', 'M05_2', 'M05_2', 'M05_2',
        'M05_3', 'M05_3', 'M05_3', 'M05_3', 'M05_3', 'M05_3', 'M05_3',
        'M05_4', 'M05_4', 'M05_4', 'M05_4', 'M05_4', 'M05_4', 'M05_4',
        'M05_5', 'M05_5', 'M05_5', 'M05_5', 'M05_5', 'M05_5', 'M05_5',
        'M16_200_2',
                            # should not be renamed/dropped:
        'M04_120', 'M04_110', 'M04_100', 'M04_130', 'M04_140',
        'M01_100', 'M01_110', 'M01_120', 'M01_130', 'M01_140',
        'M01_100_1', 'M01_110_1', 'M01_120_1', 'M01_130_1', 'M01_140_1',
        'M01_100_2', 'M01_110_2', 'M01_120_2', 'M01_130_2', 'M01_140_2',
        'M01_100_3', 'M01_110_3', 'M01_120_3', 'M01_130_3', 'M01_140_3',
        'M01_100_4', 'M01_110_4', 'M01_120_4', 'M01_130_4', 'M01_140_4',
        'M01_100_5', 'M01_110_5', 'M01_120_5', 'M01_130_5', 'M01_140_5',
    ]

    corrected_codes = [     # renamed codes:
        'M04_120_2', 'M04_110_2', 'M04_100_2', 'M04_130_2', 'M04_140_2',
        'M04_120_3', 'M04_110_3', 'M04_100_3', 'M04_130_3', 'M04_140_3',
        'M04_120_4', 'M04_110_4', 'M04_100_4', 'M04_130_4', 'M04_140_4',
        'M04_120_5', 'M04_110_5', 'M04_100_5', 'M04_130_5', 'M04_140_5',
        'M05_130_1', 'M05_150_1', 'M05_120_1', 'M05_140_1', 'M05_110_1',
        'M05_100_1', 'M05_160_1', 'M05_130_2', 'M05_150_2', 'M05_120_2',
        'M05_140_2', 'M05_110_2', 'M05_100_2', 'M05_160_2', 'M05_130_3',
        'M05_150_3', 'M05_120_3', 'M05_140_3', 'M05_110_3', 'M05_100_3',
        'M05_160_3', 'M05_130_4', 'M05_150_4', 'M05_120_4', 'M05_140_4',
        'M05_110_4', 'M05_100_4', 'M05_160_4', 'M05_130_5', 'M05_150_5',
        'M05_120_5', 'M05_140_5', 'M05_110_5', 'M05_100_5', 'M05_160_5',
        'M16_100_2',
                            # not changed codes:
        'M04_120', 'M04_110', 'M04_100', 'M04_130', 'M04_140',
        'M01_100', 'M01_110', 'M01_120', 'M01_130', 'M01_140',
        'M01_100_1', 'M01_110_1', 'M01_120_1', 'M01_130_1', 'M01_140_1',
        'M01_100_2', 'M01_110_2', 'M01_120_2', 'M01_130_2', 'M01_140_2',
        'M01_100_3', 'M01_110_3', 'M01_120_3', 'M01_130_3', 'M01_140_3',
        'M01_100_4', 'M01_110_4', 'M01_120_4', 'M01_130_4', 'M01_140_4',
        'M01_100_5', 'M01_110_5', 'M01_120_5', 'M01_130_5', 'M01_140_5', ]

    codes_to_drop = ['M02_120', 'M02_120_1', 'M02_120_2',
                     'M02_120_3', 'M02_120_4', 'M02_120_5',
                     'M07_100', 'M07_100_1', 'M07_100_2',
                     'M07_100_3', 'M07_100_4', 'M07_100_5']

    missing_codes = ['M02_120_1', 'M02_120_2',
                     'M02_120_3', 'M02_120_4', 'M02_120_5',
                     'M07_100_1', 'M07_100_2',
                     'M07_100_3', 'M07_100_4', 'M07_100_5']

    def setUp(self):
        """ """
        self.setUpPyfakefs()

    def test_activate_npis_based_on_threshold(self):
        """ """

        # test with delay = 1; should be active two days after incid > 1.5
        npi_activation_days = 1
        npi_lifting_days = 1
        int_active = gnd.activate_npis_based_on_incidence(
            self.incid, npi_lifting_days, npi_activation_days, self.threshold)
        self.assertEqual(
            int_active.to_list(),
            self.active_11.to_list())

        # test same data set with first value of incid_delay_1 = 2 > 1.5
        int_active_start_above_threshold = gnd.activate_npis_based_on_incidence(
            self.incid_start_above_threshold, npi_lifting_days,
            npi_activation_days, self.threshold)
        self.assertEqual(
            int_active_start_above_threshold.to_list(),
            self.active_start_above_threshold_11.to_list())

        # tests with day values larger 1
        npi_activation_days = 3
        npi_lifting_days = 2
        int_active = gnd.activate_npis_based_on_incidence(
            self.incid, npi_lifting_days, npi_activation_days, self.threshold)
        self.assertEqual(
            int_active.to_list(),
            self.active_32.to_list())
        int_active_start_above_threshold = gnd.activate_npis_based_on_incidence(
            self.incid_start_above_threshold, npi_lifting_days,
            npi_activation_days, self.threshold)
        self.assertEqual(
            int_active_start_above_threshold.to_list(),
            self.active_start_above_threshold_32.to_list())

        npi_activation_days = 3
        npi_lifting_days = 5
        int_active = gnd.activate_npis_based_on_incidence(
            self.incid, npi_lifting_days, npi_activation_days, self.threshold)
        self.assertEqual(
            int_active.to_list(),
            self.active_35.to_list())
        int_active_start_above_threshold = gnd.activate_npis_based_on_incidence(
            self.incid_start_above_threshold, npi_lifting_days,
            npi_activation_days, self.threshold)
        self.assertEqual(
            int_active_start_above_threshold.to_list(),
            self.active_start_above_threshold_35.to_list())

    def test_drop_codes_and_categories(self):
        """ """
        # test with no dropped codes or categories
        # just rename
        fine_resolution = 2
        npi_codes_prior_test = pd.Series(self.codes_to_correct)
        # create dataframe with corrected columns
        # just one column, no data
        df_npis_old_test = pd.DataFrame(
            columns=[dd.EngEng['npiCode']],
            data=self.corrected_codes)
        # for this test no categories should be removed
        # no 'Platzhalter' in the description
        npi_codes_prior_desc_test = pd.Series(
            ['-' for i in range(len(self.codes_to_correct))])
        # work with copies to not change the original data
        codes_dropped, npi_codes_prior, df_npis_old = gnd.drop_codes_and_categories(
            npi_codes_prior_test.copy(), npi_codes_prior_desc_test.copy(), df_npis_old_test.copy(), fine_resolution)
        # no codes should be dropped
        self.assertEqual(codes_dropped, [])
        # codes should now be corrected
        self.assertEqual(npi_codes_prior.to_list(), self.corrected_codes)
        # dataframe should not have changed
        pd.testing.assert_frame_equal(df_npis_old, df_npis_old_test)

        # now test with codes dropped
        # npi_codes_prior and npi_codes_prior_desc same as above
        # only some codes added to df_npis_old which should be removed.
        df_npis_old_test_cd = pd.DataFrame(
            columns=[dd.EngEng['npiCode']],
            data=self.corrected_codes+self.codes_to_drop)
        # work with copies to not change the original data
        codes_dropped, npi_codes_prior, df_npis_old = gnd.drop_codes_and_categories(
            npi_codes_prior_test.copy(), npi_codes_prior_desc_test.copy(), df_npis_old_test_cd.copy(), fine_resolution)
        # codes should be dropped
        self.assertEqual(codes_dropped, self.codes_to_drop)
        # codes should now be corrected as above
        self.assertEqual(npi_codes_prior.to_list(), self.corrected_codes)
        # dataframe should have changed (expect AssertionError)
        with self.assertRaises(AssertionError):
            pd.testing.assert_frame_equal(df_npis_old, df_npis_old_test_cd)

        # test handling of missing codes
        # more codes in npi_codes_prior than in df
        npi_codes_prior_test_mc = pd.Series(
            self.codes_to_correct + self.missing_codes)
        # create dataframe with corrected columns
        # just one column, no data
        df_npis_old_test = pd.DataFrame(
            columns=[dd.EngEng['npiCode']],
            data=self.corrected_codes)
        # for this test no categories should be removed
        # no 'Platzhalter' in the description
        npi_codes_prior_desc_test_mc = pd.Series(
            ['-'
             for i in range(
                 len(self.codes_to_correct + self.missing_codes))])

        # with fine_resolution = 2 every missing code should raise a DataError
        with self.assertRaises(gd.DataError):
            codes_dropped, npi_codes_prior, df_npis_old = gnd.drop_codes_and_categories(
                npi_codes_prior_test_mc.copy(), npi_codes_prior_desc_test_mc.copy(), df_npis_old_test.copy(), fine_resolution)

        # fine_resolution = 1 should handle missing codes
        fine_resolution = 1
        codes_dropped, npi_codes_prior, df_npis_old = gnd.drop_codes_and_categories(
            npi_codes_prior_test_mc.copy(), npi_codes_prior_desc_test_mc.copy(), df_npis_old_test.copy(), fine_resolution)
        # no codes should be dropped
        self.assertEqual(codes_dropped, [])
        # codes should now be corrected
        self.assertEqual(
            npi_codes_prior.to_list(),
            self.corrected_codes + self.missing_codes)
        self.assertNotEqual(npi_codes_prior.to_list(), self.corrected_codes)
        # dataframe should not have changed
        pd.testing.assert_frame_equal(df_npis_old, df_npis_old_test)

        # create dataframe with Platzhalter categories
        npi_codes_prior_test = pd.Series(self.corrected_codes[0:12])
        df_npis_old_test = pd.DataFrame(
            columns=[dd.EngEng['npiCode']],
            data=self.corrected_codes[0:12])
        # add 'Platzhalter' every second element
        npi_codes_prior_desc_test = pd.Series(
            ['-' for i in range(0, 12)])
        npi_codes_prior_desc_test[0:12:2] = 'Platzhalter'
        # work with copies to not change the original data
        codes_dropped, npi_codes_prior, df_npis_old = gnd.drop_codes_and_categories(
            npi_codes_prior_test.copy(), npi_codes_prior_desc_test.copy(), df_npis_old_test.copy(), fine_resolution)
        # check dropped codes
        self.assertEqual(codes_dropped, np.sort(
            self.corrected_codes[0:12:2]).tolist())
        # every second category dropped but naming should be the same
        self.assertEqual(
            df_npis_old[dd.EngEng['npiCode']].tolist(),
            self.corrected_codes[1: 12: 2])

    # Test full functionality only for 13 days and one county.
    # Use side_effect to return first case data, then population data with 1st
    # and 2nd call to pd.read_json().
    @patch('pandas.read_json', side_effect=[df_cases, df_pop])
    @patch('memilio.epidata.getNPIData.read_files',
           return_value=[df_npis_old, df_npis_desc, df_npis_combinations_pre])
    @patch('memilio.epidata.getNPIData.drop_codes_and_categories',
           return_value=[[],
                         df_npis_desc['Variablenname'],
                         df_npis_old_renamed])
    @patch('memilio.epidata.getNPIData.plot_interaction_matrix')
    def test_get_npi_data(self, mock_plot, mock_codes, mock_read, mock_data):
        """

        :param mock_plot: 
        :param mock_codes: 
        :param mock_read: 
        :param mock_data: 

        """
        # print 'Additional errors in consistent naming' is expected.
        # print 'WARNING: DataFrame starts with reported cases > 0 for more than 5 percent...' is expected.
        npis_test = gnd.get_npi_data(
            fine_resolution=2, out_folder=self.path,
            counties_considered=[1001],
            end_date=date(2024, 1, 1),
            npi_activation_days_threshold=3, npi_lifting_days_threshold=2)
        # test if mocks have been called
        self.assertEqual(mock_data.call_count, 2)
        self.assertEqual(mock_read.call_count, 1)
        self.assertEqual(mock_codes.call_count, 1)
        # some columns should be empty
        # either because they're not mentioned or because the incidence is not exceeded.
        self.assertEqual(
            npis_test.iloc[:, [4, 5, 6, 7, 9, 11, 13, 15, 17, 19]].values.sum(), 0)
        # incidence independent NPIs should not have changed
        self.assertEqual(
            npis_test.M1_1.to_list(),
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])
        self.assertEqual(
            npis_test.M1_2.to_list(),
            [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0])
        self.assertEqual(
            npis_test.M1_3.to_list(),
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        # test remaining_columns
        # incidence depending NPIs can first be activated on day 4 due to activation_days_threshold=3
        # incidence is constantly > 20
        # M1_2_1,M1_3_1,M1_1_2,M1_3_2 always 0
        # M1_2_2 is 0 if M1_2 is 1
        self.assertEqual(
            npis_test.M1_2_2.tolist(),
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

        self.assertEqual(
            npis_test.M1_3_2.tolist(),
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        # _4 -> Incidence > 50
        # M1_3_4 is always 0 - cant be simultaneously active with M1_3_2
        self.assertEqual(
            npis_test.M1_3_4.to_list(),
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(
            npis_test.M1_2_4.to_list(),
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0])
        # _5 is always 0 since incidence is never > 100 in this test

        # M1_1_1 should not be active when M2,3_2,3,4,5 is active
        self.assertEqual(
            npis_test.M1_1_1.to_list(),
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

        self.assertEqual(mock_plot.call_count, 4)


if __name__ == '__main__':
    unittest.main()
