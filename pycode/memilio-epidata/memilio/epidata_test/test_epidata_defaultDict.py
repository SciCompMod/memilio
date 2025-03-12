#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Patrick Lenz
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
from datetime import date

from memilio.epidata import defaultDict as dd


class Test_defaultDict(unittest.TestCase):
    """ """
    test_dict_county = {
        'Flensburg, Stadt': 1001,
        'Berlin': 11000,
        'Eisenach, Stadt': 16056,
        'Berlin Neukölln': 11008
    }
    test_dict_state = {
        'Schleswig-Holstein': 1,
        'Hamburg': 2,
        'Niedersachsen': 3,
        'Bremen': 4,
        'Nordrhein-Westfalen': 5,
        'Hessen': 6,
        'Rheinland-Pfalz': 7,
        'Baden-Württemberg': 8,
        'Bayern': 9,
        'Saarland': 10,
        'Berlin': 11,
        'Brandenburg': 12,
        'Mecklenburg-Vorpommern': 13,
        'Sachsen': 14,
        'Sachsen-Anhalt': 15,
        'Thüringen': 16
    }
    test_dict_dd = {
        date.today(): 'end_date',
        date(2020, 1, 1): 'start_date',
        'json_timeasstring': 'file_format'
    }
    test_values_notin_dd = ('no_raw', 'impute_dates', 'split_berlin',
                            'make_plot')

    def test_invert_dict(self):
        """ """
        inv_county = dd.invert_dict(dd.County)
        # check that test_dict_county is a subset of inv_county
        for county in self.test_dict_county.keys():
            self.assertIn(county, inv_county.keys())
            self.assertEqual(self.test_dict_county[county],
                             inv_county[county])
        # check inverted state dictionary complete
        inv_state = dd.invert_dict(dd.State)
        self.assertDictEqual(self.test_dict_state, inv_state)
        # check that test_dict_dd is a subset of inv_dd
        inv_dd = dd.invert_dict(dd.defaultDict)
        for defaultvalue in self.test_dict_dd.keys():
            self.assertIn(defaultvalue, inv_dd.keys())
            self.assertEqual(self.test_dict_dd[defaultvalue],
                             inv_dd[defaultvalue])
        for value in self.test_values_notin_dd:
            self.assertNotIn(value, inv_dd.values())

    def test_uniqueness(self):
        """ """
        # tests whether the names of counties and states are unique
        self.assertEqual(dd.invert_dict(dd.invert_dict(dd.County)), dd.County)
        self.assertEqual(dd.invert_dict(dd.invert_dict(dd.State)), dd.State)


if __name__ == '__main__':
    unittest.main()
