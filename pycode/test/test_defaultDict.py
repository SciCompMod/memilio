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
from collections import defaultdict
from typing import DefaultDict
import unittest
import os
from datetime import date
from epidemiology.epidata import defaultDict as dd
from pyfakefs import fake_filesystem_unittest


class Test_defaultDict(fake_filesystem_unittest.TestCase):
    test_dict_county = {
        'Flensburg, Stadt': 1001}, {
        'Berlin': 11000}, {
        'Eisenach, Stadt': 16056}, {
        'Berlin Neukölln': 11008}
    test_dict_state = {
        'Schleswig-Holstein': 1}, {
        'Hamburg': 2}, {
        'Niedersachsen': 3}, {
        'Bremen': 4}, {
        'Nordrhein-Westfalen': 5}, {
        'Hessen': 6}, {
        'Rheinland-Pfalz': 7}, {
        'Baden-Württemberg': 8}, {
        'Bayern': 9}, {
        'Saarland': 10}, {
        'Berlin': 11}, {
        'Brandenburg': 12}, {
        'Mecklenburg-Vorpommern': 13}, {
        'Sachsen': 14}, {
        'Sachsen-Anhalt': 15}, {
        'Thüringen': 16}
    test_dict_dd = {
        date.today(): 'end_date'}, {
        date(2020, 4, 24): 'start_date'}, {
        'json_timeasstring': 'file_format'}
    test_values_notin_dd = ('no_raw', 'impute_dates', 'split_berlin', 'update_data', 'make_plot')

    def setUp(self):
        self.setUpPyfakefs()

    def test_invert_dict(self):
        inv_county = dd.invert_dict(dd.County)
        for item in self.test_dict_county:
            self.assertDictContainsSubset(item, inv_county)
        inv_state = dd.invert_dict(dd.State)
        for item in self.test_dict_state:
            self.assertDictContainsSubset(item, inv_state)
        inv_dd = dd.invert_dict(dd.defaultDict)
        for item in self.test_dict_dd:
            self.assertDictContainsSubset(item,inv_dd)
        for value in self.test_values_notin_dd:
            self.assertNotIn(value, inv_dd.values())

if __name__ == '__main__':
    unittest.main()
