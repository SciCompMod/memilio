#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Henrik Zunker, Maximilian Betz
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
from memilio.simulation.osecir import Model


class Test_DynamicNPIs(unittest.TestCase):
    """ """

    def test_dynamic_npis(self):
        """ """
        model = Model(0)
        dynamic_npis = model.parameters.DynamicNPIsInfectedSymptoms
        dynamic_npis.interval = 3.0
        dynamic_npis.duration = 14.0
        dynamic_npis.base_value = 100000
        self.assertEqual(dynamic_npis.interval, 3.0)
        self.assertEqual(dynamic_npis.duration, 14.0)
        self.assertEqual(dynamic_npis.base_value, 100000)
        self.assertEqual(len(dynamic_npis.threshold), 0)

        dynamic_npis.set_threshold(
            0.5, [mio.DampingSampling(
                mio.UncertainValue(123),
                0,
                0,
                0.0,
                [],
                np.array([1]))])
        dynamic_npis.set_threshold(
            1.0, [mio.DampingSampling(
                mio.UncertainValue(543),
                0,
                0,
                0.0,
                [],
                np.array([1]))])
        self.assertEqual(len(dynamic_npis.threshold), 2)
        self.assertEqual(dynamic_npis.threshold[0][0], 1.0)
        self.assertEqual(dynamic_npis.threshold[1][0], 0.5)
        self.assertEqual(dynamic_npis.threshold[1][1][0].value.value, 123)
        self.assertEqual(dynamic_npis.threshold[0][1][0].value.value, 543)


if __name__ == '__main__':
    unittest.main()
