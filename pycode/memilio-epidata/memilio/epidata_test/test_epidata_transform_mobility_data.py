######################################################################
# Copyright (C) 2020-2025 MEmilio
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
import os
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from pyfakefs import fake_filesystem_unittest

from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import transformMobilityData as tfmd


class TestTransformMobilityData(fake_filesystem_unittest.TestCase):
    """ """
    maxDiff = None

    path = '/home/Mobility/'

    counties = geoger.get_county_ids(merge_eisenach=False)

    df_401 = pd.DataFrame(np.zeros((len(counties), len(counties))))

    # two neighboring counties (state 0) to one county (state 1)
    df_401.iloc[0, 15] = 2
    df_401.iloc[1, 15] = 3
    # two neighboring counties (state 2) to one county (state 5)
    df_401.iloc[44, 122] = 5
    df_401.iloc[45, 122] = 4
    # one county (state 2) to two neighboring counties (state 11)
    df_401.iloc[44, 333] = 1
    df_401.iloc[44, 335] = 7
    # Eisenach (state 15) to ... (state 2 and 15)
    df_401.iloc[383, 27] = 3
    df_401.iloc[383, 399] = 14
    # index between Eisenach and Wartburgkreis (state 15) to ... (state 7)
    df_401.iloc[384, 217] = 22
    # Wartburgkreis (state 15) to ... (state 2, 15, 2)
    df_401.iloc[386, 27] = 3
    df_401.iloc[399, 399] = 27
    df_401.iloc[399, 24] = 11
    # state 2 to state 15
    df_401.iloc[55, 387] = 11

    # take reduced dataframe to check functionality
    df_400 = df_401.iloc[0:400, 0:400].copy()

    def setUp(self):
        """ """
        self.setUpPyfakefs()

    # Test that correct data frame is returned by get method

    @patch('memilio.epidata.transformMobilityData.pd.read_csv',
           return_value=df_401.copy())
    def test_get_mobility_from_file(self, mock_load_file):
        """

        :param mock_load_file: 

        """
        df_return = tfmd.getMobilityFromFile(
            self.path, mobility_file='mobility')
        self.assertTrue(df_return.equals(self.df_401))

    # Test that data frame of size 401 is reduced
    @patch('memilio.epidata.transformMobilityData.getMobilityFromFile',
           return_value=df_401.copy())
    def test_update_mobility_reduction_401to400(self, mock_load_file):
        """

        :param mock_load_file: 

        """
        # read is mocked
        df_read = tfmd.getMobilityFromFile(self.path, 'mobility')
        # create folder where new file can be stored and run updateMobility
        gd.check_dir(self.path)
        df_return = tfmd.updateMobility2022(
            self.path, mobility_file='mobility')
        # different frame is returned
        self.assertFalse(df_return.equals(df_read))
        # reduced frame is returned
        self.assertEqual(len(df_return), 400)
        self.assertEqual(len(df_return.columns), 400)
        # new file is written / file is overwritten
        self.assertTrue(os.path.isfile(
            os.path.join(self.path, 'mobility.txt')))
        # Eisenach is index 383 and Wartburgkreis 386 in old county list
        self.assertTrue(abs(df_return.iloc[0:383, 0:383].values -
                        df_read.iloc[0:383, 0:383].values).max().max() < 1e-10)
        self.assertTrue(abs(df_return.iloc[383:385, :383].values -
                        df_read.iloc[384:386, :383].values).max().max() < 1e-10)
        self.assertTrue(abs(
            df_return.iloc[386:, 386:].values-df_read.iloc[387:, 387:].values).max().max() < 1e-10)
        self.assertEqual(df_return.iloc[398, 24], 11)
        self.assertEqual(df_return.iloc[55, 386], 11)
        self.assertEqual(df_return.iloc[385, 27], 6)

    # Test that data frame of size 400 is not reduced
    @patch('memilio.epidata.transformMobilityData.getMobilityFromFile',
           return_value=df_400.copy())
    def test_update_mobility_reduction_400to400(self, mock_load_file):
        """

        :param mock_load_file: 

        """
        # read is mocked
        df_read = tfmd.getMobilityFromFile(self.path, 'mobility')
        # create folder where new file can be stored and run updateMobility
        gd.check_dir(self.path)
        df_return = tfmd.updateMobility2022(
            self.path, mobility_file='mobility')
        # identical frame is returned
        self.assertTrue(df_return.equals(df_read))
        # no file is written
        self.assertFalse(os.path.isfile(
            os.path.join(self.path, 'mobility.txt')))

    # Test that empty data frame is returned if wrong sized county mobility information is given
    @patch('memilio.epidata.transformMobilityData.getMobilityFromFile',
           return_value=df_401.copy())
    def test_create_federal_states_mobility_wrong_input(self, mock_load_file):
        """

        :param mock_load_file: 

        """
        # run createFederalStatesMobility
        df_return = tfmd.createFederalStatesMobility(
            self.path, mobility_file='mobility')
        # empty frame is returned
        self.assertTrue(df_return.equals(pd.DataFrame()))
        # no file is written
        self.assertFalse(os.path.isfile(
            os.path.join(self.path, 'mobility_states.txt')))

    # Test that data frame of size 400 will be aggregated to federal state mobility
    @patch('memilio.epidata.transformMobilityData.getMobilityFromFile',
           return_value=df_400.copy())
    def test_create_federal_states_mobility(self, mock_load_file):
        """

        :param mock_load_file: 

        """
        # create folder where new file can be stored and run createFederalStatesMobility
        gd.check_dir(self.path)
        df_return = tfmd.createFederalStatesMobility(
            self.path, mobility_file='mobility')
        # data frame iof dimensions of the federal states is returned
        self.assertEqual(len(df_return), 16)
        self.assertEqual(len(df_return.columns), 16)
        # file is written
        self.assertTrue(os.path.isfile(
            os.path.join(self.path, 'mobility_states.txt')))
        # entries are aggregated, no other entries were added, no NaNs found
        self.assertEqual(df_return.iloc[0, 1], 2+3)
        self.assertEqual(df_return.iloc[2, 5], 5+4)
        self.assertEqual(df_return.iloc[2, 11], 1+7)
        self.assertEqual(df_return.iloc[15, 2], 3+3+11)
        self.assertEqual(df_return.iloc[14, 15], 0)
        self.assertEqual(df_return.iloc[15, 7], 22)
        self.assertEqual(df_return.iloc[2, 15], 11)
        self.assertEqual(df_return.sum().sum(), 72)
        self.assertEqual(len(np.where(np.isnan(df_return) == True)[0]), 0)


if __name__ == '__main__':
    unittest.main()
