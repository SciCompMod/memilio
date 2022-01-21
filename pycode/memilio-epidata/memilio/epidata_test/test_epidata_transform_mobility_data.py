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
from memilio.epidata import transformMoblityData as tfmd


class TestTransformMobilityData(fake_filesystem_unittest.TestCase):
    maxDiff = None

    path = '/home/Mobility/'

    counties = geoger.get_county_ids(merge_eisenach=False)
    idx_cols = [i for i in range(len(counties))]

    df = pd.DataFrame(columns=idx_cols)
    df[idx_cols] = np.zeros((len(counties), len(counties)))

    df.iloc[0, 15] = 2
    df.iloc[1, 15] = 3

    def setUp(self):
        self.setUpPyfakefs()

    def write_mobility_data(self, out_folder):
        file_mobility = 'mobility.txt'
        file_mobility = os.path.join(out_folder, file_mobility)
        self.df.to_csv(file_mobility, header=None, index=False)
    
    def test_update_mobility_reduction(self):
        tfmd.updateMobility2022(self.path, mobility_file='mobility')

    

if __name__ == '__main__':
    unittest.main()
