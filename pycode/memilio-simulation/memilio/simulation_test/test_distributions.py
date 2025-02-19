#############################################################################
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
#############################################################################
import unittest

import memilio.simulation as mio


class Test_ParameterDistribution(unittest.TestCase):
    def test_uniform(self):
        U = mio.ParameterDistributionUniform(1.0, 2.0)
        # properties
        self.assertEqual(U.lower_bound, 1.0)
        self.assertEqual(U.upper_bound, 2.0)
        # sample
        u = U.get_sample()
        self.assertGreaterEqual(u, 1.0)
        self.assertLessEqual(u, 2.0)

    def test_normal(self):
        N = mio.ParameterDistributionNormal(-1.0, 1.0, 0.0, 1.0)
        # properties
        self.assertEqual(N.mean, 0.0)  # std_dev automatically adapted

        # sample
        n = N.get_sample()
        self.assertGreaterEqual(n, -1)
        self.assertLessEqual(n, 1)

    def test_polymorphic(self):
        uv = mio.UncertainValue()
        N = mio.ParameterDistributionNormal(0, 2, 1.0)
        uv.set_distribution(N)
        self.assertIsNot(uv.get_distribution(), N)
        self.assertIsInstance(uv.get_distribution(),
                              mio.ParameterDistributionNormal)
        self.assertIsInstance(uv.get_distribution(), mio.ParameterDistribution)
        self.assertEqual(uv.get_distribution().mean, 1.0)


if __name__ == '__main__':
    unittest.main()
