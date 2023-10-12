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

import memilio.simulation as mio


class Test_UncertainMatrix(unittest.TestCase):
    def test_dampings(self):
        m = mio.UncertainContactMatrix(
            mio.ContactMatrixGroup(num_matrices=2, size=2))
        d = mio.DampingSampling(
            value=mio.UncertainValue(3.0),
            level=1,
            type=2,
            time=10.0,
            matrix_indices=[0, 1],
            group_weights=np.r_[2.0, 1.0])
        m.dampings = [d]
        self.assertEqual(m.dampings[0].value.value, 3.0)
        self.assertEqual(m.dampings[0].level, 1)
        self.assertEqual(m.dampings[0].type, 2)
        self.assertEqual(m.dampings[0].time, 10.0)
        self.assertEqual(m.dampings[0].matrix_indices, [0, 1])
        self.assertTrue((m.dampings[0].group_weights == np.r_[2.0, 1.0]).all())


if __name__ == '__main__':
    unittest.main()
