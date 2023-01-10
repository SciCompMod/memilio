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
from unittest.mock import patch, call
from pyfakefs import fake_filesystem_unittest

import pandas as pd
import numpy as np

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import getNPIData as gnd
from memilio.epidata import defaultDict as dd


class TestGetNPIData(fake_filesystem_unittest.TestCase):
    maxDiff = None

    path = '/home/NPIData'

    incid_no_delay =        pd.Series([1,2,3,4,1,1,0,5,0,1,2,3,5,1,0,2,4,5,6])
    test_series_no_delay =  pd.Series([0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,1,1])

    def setUp(self):
        self.setUpPyfakefs()

    # first test of everything
    # should not raise any errors
    # to be deleted if enough tests are written
    def test_get_npi_data(self, mockv):
        gnd.get_npi_data(out_folder=self.path)
    
    def test_activate_npis_based_on_threshold(self):
        threshold=1.5
        # test with delay = 1 ; int should be active two days after incid > 1.5
        # should start with 0 since first value of incid no delay = 1 < 1.5
        npi_activation_delay = 1
        npi_lifting_deay = 1
        int_active_no_delay = gnd.activate_npis_based_on_incidence(
            self.incid_no_delay, npi_lifting_deay, npi_activation_delay, threshold)
        self.assertEqual(int_active_no_delay.to_list(), self.test_series_no_delay.to_list())


if __name__ == '__main__':
    unittest.main()
