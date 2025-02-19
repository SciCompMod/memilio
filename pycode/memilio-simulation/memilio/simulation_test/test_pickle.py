#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Maximilian Betz
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
import pickle
import unittest

import numpy as np

import memilio.simulation as msim


class Test_Pickle(unittest.TestCase):
    def test_date(self):
        test = msim.Date(1, 2, 3)

        data = pickle.dumps(test)
        pickle_test = pickle.loads(data)

        self.assertEqual(pickle_test.year, 1)
        self.assertEqual(pickle_test.month, 2)
        self.assertEqual(pickle_test.day, 3)

    def test_uncertain_value(self):
        test = msim.UncertainValue(2.2)

        data = pickle.dumps(test)
        pickle_test = pickle.loads(data)

        self.assertEqual(pickle_test.value, 2.2)

    def test_distribution(self):

        test = msim.ParameterDistributionNormal(0, 1, 0.4, 0.1)

        data = pickle.dumps(test)
        pickle_test = pickle.loads(data)

        self.assertEqual(pickle_test.mean, 0.4)
        self.assertEqual(pickle_test.standard_dev, 0.1)

    def test_damping_sampling(self):
        test = msim.UncertainValue(2.2)
        test.set_distribution(msim.ParameterDistributionNormal(0, 1, 0.4, 0.1))
        test = msim.DampingSampling(test, 1, 2, 3, [0, 1], np.arange(2))

        data = pickle.dumps(test)
        pickle_test = pickle.loads(data)

        self.assertEqual(pickle_test.value.value, 2.2)

        self.assertEqual(pickle_test.value.get_distribution().mean, 0.4)
        self.assertEqual(
            pickle_test.value.get_distribution().standard_dev, 0.1)

        self.assertEqual(pickle_test.level, 1)
        self.assertEqual(pickle_test.type, 2)
        self.assertEqual(pickle_test.time, 3)

        self.assertEqual(pickle_test.matrix_indices, [0, 1])
        self.assertEqual(pickle_test.group_weights[0], 0)
        self.assertEqual(pickle_test.group_weights[1], 1)


if __name__ == '__main__':
    unittest.main()
