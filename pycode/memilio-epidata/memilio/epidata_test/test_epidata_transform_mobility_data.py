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
import os
import unittest
from unittest.mock import patch
from pyfakefs import fake_filesystem_unittest

import pandas as pd
import numpy as np

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import transformMobilityData as tfmd


class TestTransformMobilityData(fake_filesystem_unittest.TestCase):
    maxDiff = None

    path = '/home/Mobility/'

    counties = geoger.get_county_ids(merge_eisenach=False)
    idx_cols = [i for i in range(len(counties))]

    df_401 = pd.DataFrame(columns=idx_cols)
    df_401[idx_cols] = np.zeros((len(counties), len(counties)))

    # two neighboring counties to one county (and federal state, i.e., Hamburg)
    df_401.iloc[0, 15] = 2
    df_401.iloc[1, 15] = 3
    # two neighboring counties to one county
    df_401.iloc[44, 122] = 5
    df_401.iloc[45, 122] = 4
    # one county to two neighboring counties
    df_401.iloc[44, 333] = 2
    df_401.iloc[44, 335] = 7  
    # Eisenach to ...
    df_401.iloc[383, 27] = 3
    df_401.iloc[383, 399] = 14
    # index between Eisenach and Wartburgkreis to ...
    df_401.iloc[384, 217] = 22
    # Wartburgkreis to ...
    df_401.iloc[386, 27] = 3
    df_401.iloc[399, 399] = 27
    df_401.iloc[399, 24] = 11
    df_401.iloc[55, 387] = 11    

    def setUp(self):
        self.setUpPyfakefs()
    
    # Test that data frame of size 401 is reduced
    @patch('memilio.epidata.transformMobilityData.getMobilityFromFile',
           return_value=df_401.copy())
    def test_update_mobility_reduction_401to400(self, mock_load_file):
        # read is mocked
        df_read = tfmd.getMobilityFromFile(self.path, 'mobility')        
        # create folder where new file can be stored and run updateMobility
        gd.check_dir(self.path)
        df_return = tfmd.updateMobility2022(self.path, mobility_file='mobility')
        # different frame is returned
        self.assertFalse(df_return.equals(df_read))
        # reduced frame is returned
        self.assertEqual(len(df_return), 400)
        self.assertEqual(len(df_return.columns), 400)
        # new file is written / file is overwritten
        self.assertTrue(os.path.isfile(os.path.join(self.path, 'mobility.txt')))
        # Eisenach is index 383 and Wartburgkreis 386 in old county list
        self.assertTrue(abs(df_return.iloc[0:383, 0:383].values-df_read.iloc[0:383, 0:383].values).max().max() < 1e-10)
        self.assertTrue(abs(df_return.iloc[383:385, :383].values-df_read.iloc[384:386, :383].values).max().max() < 1e-10)
        self.assertTrue(abs(df_return.iloc[386:, 386:].values-df_read.iloc[387:, 387:].values).max().max() < 1e-10)
        self.assertTrue(df_return.iloc[398, 24] == 11)
        self.assertTrue(df_return.iloc[55, 386] == 11)
        self.assertTrue(df_return.iloc[385, 27] == 6)                     

    # Test that data frame of size 400 is not reduced
    @patch('memilio.epidata.transformMobilityData.getMobilityFromFile',
           return_value=df_401.iloc[0:400,0:400].copy())
    def test_update_mobility_reduction_400to400(self, mock_load_file):
        # read is mocked
        df_read = tfmd.getMobilityFromFile(self.path, 'mobility')
        # create folder where new file can be stored and run updateMobility
        gd.check_dir(self.path)
        df_return = tfmd.updateMobility2022(self.path, mobility_file='mobility')
        # identical frame is returned
        self.assertTrue(df_return.equals(df_read))
        # no file is written
        self.assertFalse(os.path.isfile(os.path.join(self.path, 'mobility.txt')))


if __name__ == '__main__':
    unittest.main()
