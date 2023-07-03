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
import os
import unittest
from unittest.mock import patch

import twill
import numpy as np
import pandas as pd
from pyfakefs import fake_filesystem_unittest

from memilio.epidata import defaultDict as dd
from memilio.epidata import getPopulationData as gpd
from memilio.epidata import progress_indicator

progress_indicator.ProgressIndicator.disable_indicators(True)


class Test_getPopulationData(fake_filesystem_unittest.TestCase):

    def setUp(self):
        self.setUpPyfakefs()

    def test_read_population_data(self):

        directory = os.path.join(self.path, 'Germany/')

        # test file not found
        with self.assertRaises(FileNotFoundError) as error:
            df = gpd.read_population_data(
                username='', password='', read_data=True, directory=directory)



if __name__ == '__main__':
    unittest.main()
