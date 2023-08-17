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
import os
import json
import collections
import unittest
import pandas as pd
from unittest.mock import patch

from pyfakefs import fake_filesystem_unittest

from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getCommuterMobility as gcm
from memilio.epidata import progress_indicator


class TestCommuterMigration(fake_filesystem_unittest.TestCase):

    path = '/home/CMData/'

    here = os.path.dirname(os.path.abspath(__file__))

    test_govkey_list = [
        '01', '02', '031', '032', '033', '034', '04', '051', '053', '055',
        '057', '059']

    countykey_list = geoger.get_county_ids(merge_eisenach=False, zfill=True)
    govkey_list = geoger.get_governing_regions()

    test_countykey_list = [1001, 1002, 1003, 1004,
                           1051, 1053, 1054, 1055, 1056, 1057, 1058, 3159]

    test_countykey_list2 = [1001, 1002, 1003, 1051,
                            1053, 1004, 1055, 1054, 1056, 1057, 1058]

    test_countykey2govkey = {'01055': 0, '02000': 1, '03101': 2}

    test_countykey2localnumlist = {
        '01001': 0, '01004': 3, '01060': 12, '10045': 4, '11000': 0}

    test_state_gov_table = [['01'], ['02'], ['031', '032', '033', '034'], [
                            '04'], ['051', '053', '055', '057', '059']]

    test_gov_county_table = (
        ['02000'],
        ['05512', '05513', '05515', '05554', '05558', '05562', '05566',
         '05570'],
        ['05711', '05754', '05758', '05762', '05766', '05770', '05774'])

    filename = os.path.join(
        here, 'test_data', 'TestSetPopulationFinal.json')
    with open(filename) as file_object:
        df_pop = pd.DataFrame(json.load(file_object))

    def setUp(self):
        self.setUpPyfakefs()
        progress_indicator.ProgressIndicator.disable_indicators(True)

    @patch('builtins.print')
    def test_verify_sorted(self, mock_print):
        self.assertEqual(True, gcm.verify_sorted(self.test_countykey_list))
        self.assertEqual(False, gcm.verify_sorted(self.test_countykey_list2))
        Errorcall = ('Error. Input list not sorted.')
        mock_print.assert_called_with(Errorcall)

    @patch('builtins.print')
    def test_assign_geographical_entities(self, mock_print):
        (
            countykey2govkey, countykey2localnumlist, gov_county_table,
            state_gov_table) = gcm.assign_geographical_entities(
            self.countykey_list, self.govkey_list)
        for item in self.test_countykey2govkey.keys():
            self.assertEqual(
                self.test_countykey2govkey.get(item),
                countykey2govkey.get(item))

        # check if all countyIDs are in countykey2govkey
        for key in geoger.get_county_ids(True, False, True):
            self.assertIn(key, countykey2govkey.keys())
            self.assertIn(key, countykey2localnumlist.keys())
        for item in self.test_countykey2localnumlist.keys():
            self.assertEqual(self.test_countykey2localnumlist.get(
                item), countykey2localnumlist.get(item))
        for item in self.test_state_gov_table:
            self.assertIn(item, state_gov_table)
        for item in self.test_gov_county_table:
            self.assertIn(item, gov_county_table)

        # test case with not matching countykey and govkey lists
        (countykey2govkey, countykey2localnumlist, gov_county_table,
         state_gov_table) = gcm.assign_geographical_entities(
            self.test_countykey_list, self.test_govkey_list)
        self.assertEqual(countykey2govkey, collections.OrderedDict())
        self.assertEqual(countykey2localnumlist, collections.OrderedDict())
        self.assertEqual(gov_county_table, [
                         [], [], [], [], [], [], [], [], [], [], [], []])
        self.assertEqual(state_gov_table, self.test_state_gov_table)

        # test case with different number of data
        gcm.assign_geographical_entities(
            self.test_countykey_list, self.govkey_list)
        Errorcall = ('Error. Number of government regions wrong.')
        mock_print.assert_called_with(Errorcall)

    @patch('memilio.epidata.getPopulationData.get_population_data', return_value=df_pop)
    @patch('builtins.input', return_value='y')
    def test_commuter_data(self, mock_input, mock_popul):
        """! Tests migration data by some randomly chosen tests.
        """

        df_commuter_migration = gcm.get_commuter_data(
            out_folder=self.path, ref_year=2022)

        # just do some tests on randomly chosen migrations

        # check migration from Leverkusen (averaged from NRW, 05) to Hildburghausen
        city_from = 5316
        city_to = 16069
        self.assertAlmostEqual(
            df_commuter_migration.loc[city_from, city_to], 0.451, 2)

        # check migration from Duisburg to Oberspreewald-Lausitz
        city_from = 5112
        city_to = 12066
        self.assertAlmostEqual(
            df_commuter_migration.loc[city_from, city_to], 3.143, 2)

        # check migration from Lahn-Dill-Kreis to Hamburg
        city_from = 6532
        city_to = 2000
        self.assertEqual(df_commuter_migration.loc[city_from, city_to], 100)

        # check migration from Herzogtum Lauenburg to Flensburg, Stadt
        city_from = 1001
        city_to = 1053
        self.assertEqual(df_commuter_migration.loc[city_from, city_to], 29)

    @patch('memilio.epidata.getPopulationData.get_population_data', return_value=df_pop)
    @patch('builtins.input', return_value='y')
    @patch('builtins.print')
    def test_get_neighbors_mobility(self, mock_print, mock_input, mock_popul):

        testcountyid = 1051
        # direction = both
        (countykey_list, commuter_all) = gcm.get_neighbors_mobility(
            testcountyid, direction='both', abs_tol=0, rel_tol=0,
            tol_comb='or', out_folder=self.path)
        self.assertEqual(len(countykey_list), 398)
        self.assertEqual(271, commuter_all[0])
        self.assertEqual(2234, commuter_all[9])
        self.assertEqual(366, commuter_all[11])
        self.assertEqual(34, commuter_all[347])

        # direction = in
        (countykey_list, commuter_all) = gcm.get_neighbors_mobility(
            testcountyid, direction='in', abs_tol=0, rel_tol=0,
            tol_comb='or', out_folder=self.path)
        self.assertEqual(len(countykey_list), 393)
        self.assertEqual(70, commuter_all[0])
        self.assertEqual(892, commuter_all[9])
        self.assertEqual(112, commuter_all[11])

        # direction = out
        (countykey_list, commuter_all) = gcm.get_neighbors_mobility(
            testcountyid, direction='out', abs_tol=0, rel_tol=0,
            tol_comb='or', out_folder=self.path)
        self.assertEqual(len(countykey_list), 378)
        self.assertEqual(201, commuter_all[0])
        self.assertEqual(1342, commuter_all[9])
        self.assertEqual(254, commuter_all[11])


if __name__ == '__main__':
    unittest.main()
