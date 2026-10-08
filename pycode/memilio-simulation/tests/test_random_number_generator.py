#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Carlotta Gerstein
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


class Test_RandomNumberGenerator(unittest.TestCase):
    """ """

    def test_advance_counter(self):
        """ """
        RNG = mio.RandomNumberGenerator()

        RNG.seed([42])
        counter = RNG.counter

        mio.DiscreteDistribution.get_instance()(RNG, [0.2, 0.4, 0.3, 0.1])
        self.assertEqual(counter + 1, RNG.counter)

    def test_discrete_distribution(self):
        RNG = mio.RandomNumberGenerator()

        sample = mio.DiscreteDistribution.get_instance()(RNG, [1.])
        self.assertEqual(sample, 0)

        samples = [mio.DiscreteDistribution.get_instance()(
            RNG, [0.2, 0.4, 0.3, 0.1]) for _ in range(100)]

        for s in samples:
            self.assertGreaterEqual(s, 0)
            self.assertLessEqual(s, 3)


if __name__ == '__main__':
    unittest.main()
