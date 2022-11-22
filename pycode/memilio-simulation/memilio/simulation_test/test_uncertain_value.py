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

import memilio.simulation as mio


class Test_UncertainValue(unittest.TestCase):
    def test_value(self):
        uv = mio.UncertainValue(0)
        self.assertEqual(uv.value, 0.0)
        uv.value = 1.0
        self.assertEqual(uv.value, 1.0)

    def test_distribution(self):
        uv = mio.UncertainValue(0)
        uv.set_distribution(mio.ParameterDistributionUniform(2, 2))
        uv.draw_sample()
        self.assertEqual(uv.value, 2.0)


if __name__ == '__main__':
    unittest.main()
