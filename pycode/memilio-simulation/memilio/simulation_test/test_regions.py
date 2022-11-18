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


class Test_Regions(unittest.TestCase):
    def test_get_holidays(self):
        holidays = mio.get_holidays_de(9, start_date=mio.Date(
            2020, 10, 15), end_date=mio.Date(2020, 11, 15))
        self.assertEqual(len(holidays), 1)
        self.assertEqual(holidays[0][0], mio.Date(2020, 10, 31))
        self.assertEqual(holidays[0][1], mio.Date(2020, 11, 7))

    def test_state_id(self):
        self.assertEqual(mio.get_state_id_de(1001), 1)
        self.assertEqual(mio.get_state_id_de(2000), 2)
        self.assertEqual(mio.get_state_id_de(9161), 9)


if __name__ == '__main__':
    unittest.main()
