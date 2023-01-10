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
import unittest
from pyfakefs import fake_filesystem_unittest

import pandas as pd

from memilio.epidata import getNPIData as gnd
from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gd

class TestGetNPIData(fake_filesystem_unittest.TestCase):
    maxDiff = None

    path = '/home/NPIData'

    incid = pd.Series(
        [1, 2, 3, 4, 1, 1, 0, 5, 0, 1, 2, 3, 5, 1, 0, 2, 4, 5, 6])
    active = pd.Series(
        [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1])

    incid_start_above_threshold = pd.Series(
        [2, 2, 3, 4, 1, 1, 0, 5, 0, 1, 2, 3, 5, 1, 0, 2, 4, 5, 6])
    active_start_above_threshold = pd.Series(
        [0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1])

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
        self.setUpPyfakefs()

    def test_activate_npis_based_on_threshold(self):
        threshold = 1.5
        # test with delay = 1 ; int should be active two days after incid > 1.5
        # should start with 0 since first value of incid_delay_1 = 1 < 1.5
        npi_activation_delay = 1
        npi_lifting_delay = 1
        int_active = gnd.activate_npis_based_on_incidence(
            self.incid, npi_lifting_delay, npi_activation_delay, threshold)
        self.assertEqual(
            int_active.to_list(),
            self.active.to_list())

        # test same data set with first value of incid_delay_1 = 2 > 1.5
        int_active_start_above_threshold = gnd.activate_npis_based_on_incidence(
            self.incid_start_above_threshold, npi_lifting_delay,
            npi_activation_delay, threshold)
        self.assertEqual(
            int_active_start_above_threshold.to_list(),
            self.active_start_above_threshold.to_list())
        
        # TODO maybe test bigger data sets and npi_activation_delay != npi_lifting_delay but both >= 3

    def test_drop_codes_and_categories(self):
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
        npi_codes_prior_test_mc = pd.Series(self.codes_to_correct + self.missing_codes)
        # create dataframe with corrected columns
        # just one column, no data
        df_npis_old_test = pd.DataFrame(
            columns=[dd.EngEng['npiCode']],
            data=self.corrected_codes)
        # for this test no categories should be removed
        # no 'Platzhalter' in the description
        npi_codes_prior_desc_test_mc = pd.Series(
            ['-' for i in range(len(self.codes_to_correct+ self.missing_codes))])
            
        # with fine_resolution = 2 every missing code should raise a DataError
        with self.assertRaises(gd.DataError):
            codes_dropped, npi_codes_prior, df_npis_old = gnd.drop_codes_and_categories(
            npi_codes_prior_test_mc.copy(), npi_codes_prior_desc_test_mc.copy(), df_npis_old_test.copy(), fine_resolution)

        # fine_resolution = 1 should handle missing codes and 
        fine_resolution=1
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

        # TODO test Platzhalter


if __name__ == '__main__':
    unittest.main()
