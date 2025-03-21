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

import unittest

from memilio.simulation import Date


class TestCustomIndexArray(unittest.TestCase):
    """ """

    def test_init(self):
        """ """
        year = 2020
        month = 3
        day = 15

        date = Date(year, month, day)
        self.assertEqual(date.year, year)
        self.assertEqual(date.month, month)
        self.assertEqual(date.day, day)
        self.assertEqual(date.day_in_year, 75)

    def test_comparison(self):
        """ """
        self.assertEqual(Date(2021, 3, 12), Date(2021, 3, 12))
        self.assertNotEqual(Date(2021, 5, 11), Date(2021, 5, 12))
        self.assertNotEqual(Date(2021, 5, 11), Date(2021, 6, 11))
        self.assertNotEqual(Date(2021, 5, 11), Date(2022, 5, 11))

        self.assertTrue(Date(2020, 5, 10) < Date(2021, 1, 1))
        self.assertTrue(Date(2020, 5, 10) < Date(2020, 6, 1))
        self.assertTrue(Date(2020, 5, 10) < Date(2020, 5, 11))
        self.assertFalse(Date(2021, 5, 10) < Date(2020, 5, 11))
        self.assertFalse(Date(2020, 5, 10) < Date(2020, 5, 10))
        self.assertTrue(Date(2020, 5, 10) <= Date(2020, 5, 10))
        self.assertTrue(Date(2020, 5, 10) <= Date(2020, 5, 11))
        self.assertFalse(Date(2020, 5, 10) <= Date(2020, 5, 9))

        self.assertFalse(Date(2020, 5, 10) > Date(2021, 1, 1))
        self.assertFalse(Date(2020, 5, 10) > Date(2020, 6, 1))
        self.assertFalse(Date(2020, 5, 10) > Date(2020, 5, 11))
        self.assertTrue(Date(2021, 5, 10) > Date(2020, 5, 11))
        self.assertFalse(Date(2020, 5, 10) > Date(2020, 5, 10))
        self.assertTrue(Date(2020, 5, 10) >= Date(2020, 5, 10))
        self.assertFalse(Date(2020, 5, 10) >= Date(2020, 5, 11))
        self.assertTrue(Date(2020, 5, 10) >= Date(2020, 5, 9))

    def test_calculation(self):
        """ """
        self.assertEqual(Date(2020, 8, 30) - Date(2020, 8, 15), 15)
        self.assertEqual(Date(2020, 8, 30) - Date(2020, 8, 31), -1)
        self.assertEqual(Date(2020, 7, 25) - Date(2020, 5, 25), 61)

        date1 = Date(2020, 8, 15)
        date2 = date1 + 15
        self.assertEqual(date1, Date(2020, 8, 15))
        self.assertEqual(date2, Date(2020, 8, 30))

        date1 += 15
        self.assertEqual(date1, Date(2020, 8, 30))

        date1 = Date(2020, 7, 25)
        date2 = date1 - 61
        self.assertEqual(date1, Date(2020, 7, 25))
        self.assertEqual(date2, Date(2020, 5, 25))

        date1 -= 61
        self.assertEqual(date1, Date(2020, 5, 25))


if __name__ == '__main__':
    unittest.main()
