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


class Test_TimeSeriesFunctor(unittest.TestCase):
    """ """

    def test_timeseries_functor(self):
        """ """
        functor = mio.TimeSeriesFunctor([[0, 0.0], [2, 0.67], [4, 0.4]])

        self.assertEqual(functor(0), 0.0)
        self.assertGreaterEqual(functor(1), 0.0)
        self.assertLessEqual(functor(1), 0.67)
        self.assertEqual(functor(2), 0.67)
        self.assertGreaterEqual(functor(3), 0.4)
        self.assertLessEqual(functor(3), 0.67)
        self.assertEqual(functor(4), 0.4)


if __name__ == '__main__':
    unittest.main()
