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

import numpy as np
from numpy.testing import assert_array_equal

import memilio.simulation as mio


class Test_TimeSeries(unittest.TestCase):
    def test_add_time_point(self):
        ts = mio.TimeSeries(1)
        ts.add_time_point(2, np.r_[1])
        self.assertEqual(ts.get_num_time_points(), 1)
        self.assertEqual(ts.get_time(0), 2)
        self.assertEqual(ts.get_value(0), np.r_[1])

        ts.add_time_point(3.5)
        self.assertEqual(ts.get_num_time_points(), 2)
        self.assertEqual(ts.get_last_time(), 3.5)

    def test_set_value(self):
        ts = mio.TimeSeries(1)
        ts.add_time_point(0.0, np.r_[1.0])
        ts.get_value(0)[:] = np.r_[2.0]
        self.assertEqual(ts.get_value(0), np.r_[2.0])

    def test_ndarray(self):
        ts = mio.TimeSeries(2)
        ts.add_time_point()
        ts.add_time_point()
        arr = ts.as_ndarray()
        self.assertEqual(arr.shape, (3, 2))
        arr[:, 1] = np.r_[1.0, 1.1, 1.2]
        assert_array_equal(ts.get_last_value(), np.r_[1.1, 1.2])
        assert_array_equal(ts.get_last_time(), 1.0)


if __name__ == '__main__':
    unittest.main()
