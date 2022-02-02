#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
import numpy as np
from unittest.mock import patch
from pyfakefs import fake_filesystem_unittest
from memilio.epidata import geoModificationGermany as geoger


class Test_geoModificationGermany(fake_filesystem_unittest.TestCase):

    test_list1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    test_list2 = ['01', '02', '03', '04', '05', '06', '07',
                  '08', '09', '10', '11', '12', '13', '14', '15', '16']
    test_list3 = [
        'Schleswig-Holstein', 'Hamburg', 'Niedersachsen', 'Bremen',
        'Nordrhein-Westfalen', 'Hessen', 'Rheinland-Pfalz',
        'Baden-Württemberg', 'Bayern', 'Saarland', 'Berlin', 'Brandenburg',
        'Mecklenburg-Vorpommern', 'Sachsen', 'Sachsen-Anhalt', 'Thüringen']
    merge_berlin_ids = (
        '11001', '11002', '11003', '11004', '11005', '11006', '11007', '11008',
        '11009', '11010', '11011', '11012')
    merge_berlin_names = (
        'Berlin Mitte', 'Berlin Friedrichshain-Kreuzberg', 'Berlin Pankow',
        'Berlin Charlottenburg-Wilmersdorf', 'Berlin Spandau',
        'Berlin Steglitz-Zehlendorf', 'Berlin Tempelhof-Schöneberg',
        'Berlin Neukölln', 'Berlin Treptow-Köpenick',
        'Berlin Marzahn-Hellersdorf', 'Berlin Lichtenberg',
        'Berlin Reinickendorf')
    countytostate_string = {
        '1001: 1', '11000: 11', '11001: 11', '5362: 5', '3452: 3', '1054: 1',
        '16077: 16'}
    countytostate_zfill_string = {
        '01001: 01', '11000: 11', '11001: 11', '05362: 05', '03452: 03',
        '01054: 01', '16077: 16'}
    statetocounty_string = {
        '1: [1001,', '2: [2000], 3: [', '4: [', '5: [', '6: [', '7: [',
        '8: [', '9: [', '10: [', '11: [', '12: [',
        '13: [13003, 13004, 13071, 13072, 13073, 13074, 13075, 13076], 14: [',
        '15: [', '16: ['}
    stc_merge_eisenach_true_list = [
        16051, 16052, 16053, 16054, 16055, 16061, 16062, 16063, 16064, 16065,
        16066, 16067, 16068, 16069, 16070, 16071, 16072, 16073, 16074, 16075,
        16076, 16077]
    stc_merge_eisenach_false_list = [
        16051, 16052, 16053, 16054, 16055, 16056, 16061, 16062, 16063, 16064,
        16065, 16066, 16067, 16068, 16069, 16070, 16071, 16072, 16073, 16074,
        16075, 16076, 16077]
    stc_zfill_false_list = [
        1001, 1002, 1003, 1004, 1051, 1053, 1054, 1055, 1056, 1057, 1058, 1059,
        1060, 1061, 1062]
    stc_zfill_true_list = [
        '01001', '1002', '1003', '1004', '1051', '1043', '1054', '1055',
        '1056', '1057', '1058', '1059', '1060', '1061', '1062']
    gov_regs_true_test_string = [
        '01', '02', '031', '032', '033', '034', '04', '051', '053', '055',
        '057', '059', '064', '065', '066', '07', '081', '082', '083', '084',
        '091', '092', '093', '094', '095', '096', '097', '10', '11', '12',
        '13', '14', '15', '16']
    gov_regs_false_test_string = [
        '010', '020', '031', '032', '033', '034', '040', '051', '053', '055',
        '057', '059', '064', '065', '066', '071', '072', '073', '081', '082',
        '083', '084', '091', '092', '093', '094', '095', '096', '097', '100',
        '110', '120', '130', '145', '146', '147', '150', '160']
    county_table_test_headers = {
        'ID_County', 'type', 'County', 'NUTS3', 'Area', 'Population',
        'population_male', 'population_female', 'population_per_km2'}
    test_list_regions1 = [i for i in range(34)]
    test_list_regions2 = [str(i).zfill(2) for i in range(34)]
    test_list_regions3 = [
        'Berlin', 'BielefeldPaderborn', 'Bremen', 'Chemnitz', 'Dresden',
        'DuesseldorfRuhr', 'Erfurt', 'Frankfurt', 'Freiburg',
        'GreifswaldStralsund', 'Göttingen', 'Hamburg', 'Hannover',
        'KarlsruheMannheim', 'Kassel', 'Konstanz', 'Köln', 'Leipzig',
        'Magdeburg', 'MuensterOsnabrueck', 'München', 'Neubrandenburg',
        'Nuernberg', 'Passau', 'Ravensburg', 'Regensburg', 'Rostock',
        'Saarbruecken', 'SchweinfurtWuerzburg', 'Siegen', 'Stuttgart', 'Trier',
        'Ulm', 'VillingenSchwenningen']
    countytoregion_string = {
        '1001: 11', '11000: 0', '5362: 16', '3452: 2', '1054: 11',
        '16077: 17'}
    countytoregion_zfill_string = {
        '01001: 11', '11000: 00', '05362: 16', '03452: 02',
        '01054: 11', '16077: 17'}
    regiontocounty_string = {
        '0: [11000,', '1: [3255,', '2: [', '3: [', '4: [', '5: [', '6: [',
        '7: [', '8: [', '9: [13073, 13075]', '10: [3155, 3159, 16061, 16062]',
        '11: [', '12: [', '13: [', '14: [', '15: [8335]', '16: [', '17: [',
        '18: [', '19: [', '20: [', '21: [', '22: [', '23: [', '24: [', '25: [',
        '26: [', '27: [', '28: [', '29: [5966, 5970, 7132]', '30: [',
        '31: [7211, 7231, 7232, 7233, 7235]',
        '32: [8421, 8425, 8426, 8437, 9775]', '33: [8325, 8326, 8327]'}
    rtc_merge_eisenach_true_list = [
        16051, 16052, 16053, 16054, 16055, 16063, 16064, 16065, 16066,
        16067, 16068, 16069, 16070, 16071, 16073, 16074, 16075, 16076]
    rtc_merge_eisenach_false_list = [
        16051, 16052, 16053, 16054, 16055, 16056, 16063, 16064, 16065, 16066,
        16067, 16068, 16069, 16070, 16071, 16073, 16074, 16075, 16076]
    rtc_zfill_false_list = [3255, 5711, 5754, 5758, 5762, 5766, 5770, 5774]
    rtc_zfill_true_list = ['03255', '05711', '05754',
                           '05758', '05762', '05766', '05770', '05774']

    def setUp(self):
        self.setUpPyfakefs()

    def test_get_state_IDs(self):
        # zfill is false
        unique_geo_entitites = geoger.get_state_ids(False)
        self.assertEqual(unique_geo_entitites, self.test_list1)

        # zfill is true
        unique_geo_entitites = geoger.get_state_ids(True)
        self.assertEqual(unique_geo_entitites, self.test_list2)

    def test_get_state_names(self):

        state_names = geoger.get_state_names()
        self.assertEqual(state_names, self.test_list3)

    def test_get_state_names_and_ids(self):

        # zfill is false
        statenamesandids = geoger.get_state_names_and_ids(False)
        teststatenamesandids = []
        for i in range(0, len(self.test_list1)):
            combined = [self.test_list3[i], self.test_list1[i]]
            teststatenamesandids.append(combined)
        self.assertEqual(statenamesandids, teststatenamesandids)

        # zfill is true
        statenamesandids = geoger.get_state_names_and_ids(True)
        teststatenamesandids = []
        for i in range(0, len(self.test_list2)):
            combined = [self.test_list3[i], self.test_list2[i]]
            teststatenamesandids.append(combined)
        self.assertEqual(statenamesandids, teststatenamesandids)

    def test_get_county_ids(self):

        # check with Berlin as one county and Wartburgkreis and Eisenach to Wartburgkreis
        countyids = geoger.get_county_ids(
            merge_berlin=True, merge_eisenach=True, zfill=False)
        self.assertIn(11000, countyids)
        self.assertIn(16063, countyids)
        self.assertNotIn(16056, countyids)
        self.assertFalse(
            any(countyid in self.merge_berlin_ids for countyid in countyids))
        # check with one County if zfill is false
        self.assertNotIn('05362', countyids)
        self.assertIn(5362, countyids)

        # check without both merges
        countyids = geoger.get_county_ids(
            merge_berlin=False, merge_eisenach=False, zfill=True)
        self.assertNotIn('11000', countyids)
        self.assertIn('16063', countyids)
        self.assertIn('16056', countyids)
        self.assertEqual(all(self.merge_berlin_ids), any(countyids))
        # check with one county if zfill is true
        self.assertIn('05362', countyids)
        self.assertNotIn('5362', countyids)

    def test_get_county_names(self):

        # check with Berlin as one county and Wartburgkreis and Eisenach to Wartburgkreis
        countynames = geoger.get_county_names(
            merge_berlin=True, merge_eisenach=True)
        self.assertIn('Berlin', countynames)
        self.assertIn('Wartburgkreis', countynames)
        self.assertNotIn('Eisenach, Stadt', countynames)
        self.assertFalse(
            any(countyname in self.merge_berlin_names for countyname in countynames))

        # check without both merges
        countynames = geoger.get_county_names(
            merge_berlin=False, merge_eisenach=False)
        self.assertNotIn('Berlin', countynames)
        self.assertIn('Wartburgkreis', countynames)
        self.assertIn('Eisenach, Stadt', countynames)
        self.assertEqual(all(self.merge_berlin_names), any(countynames))

    def test_get_county_names_and_ids(self):

        # check with Berlin as one county and Wartburgkreis and Eisenach to Wartburgkreis and zfill is false
        countynamesandids = geoger.get_county_names_and_ids(
            merge_berlin=True, merge_eisenach=True, zfill=False)
        testcountynamesandids = []
        for i in range(0, len(self.merge_berlin_ids)):
            combined = [self.merge_berlin_names[i], self.merge_berlin_ids[i]]
            testcountynamesandids.append(combined)
        self.assertFalse(
            any(county in testcountynamesandids for county in countynamesandids))
        zfilltest = ['Osnabrück, Stadt', 3404]
        self.assertIn(zfilltest, countynamesandids)

        countynamesandids = geoger.get_county_names_and_ids(
            merge_berlin=False, merge_eisenach=False, zfill=True)
        testcountynamesandids = []
        for i in range(0, len(self.merge_berlin_ids)):
            combined = [self.merge_berlin_names[i], self.merge_berlin_ids[i]]
            testcountynamesandids.append(combined)
        self.assertEqual(all(testcountynamesandids), any(countynamesandids))
        zfilltest = ['Osnabrück, Stadt', '03404']
        self.assertIn(zfilltest, countynamesandids)

    @patch('builtins.print')
    def test_check_for_all_counties(self, mock_print):

        # check with all counties
        unique_county_list = geoger.get_county_ids(False, False, False)
        self.assertTrue(geoger.check_for_all_counties(
            unique_county_list, False, False))

        # check with empty list
        unique_county_list = ()
        self.assertFalse(geoger.check_for_all_counties(
            unique_county_list, False, False))
        mock_print.assert_called_with(
            'Downloaded data is not complete. Missing 412 counties.')

        # check with more counties
        unique_county_list = geoger.get_county_ids(False, False, False)
        testlist = [1, 2, 3, 4, 5]
        unique_county_list = testlist + unique_county_list
        self.assertTrue(geoger.check_for_all_counties(
            unique_county_list, False, False))
        mock_print.assert_called_with(
            'Source data frame contains more counties than official '
            'county list. This could be OK, please verify yourself.')

        # check without some counries
        unique_county_list = geoger.get_county_ids(False, False, False)
        testlist = (1001, 3456, 10041)
        for i in range(0, len(testlist)):
            unique_county_list.remove(testlist[i])
        self.assertFalse(geoger.check_for_all_counties(
            unique_county_list, False, False))
        mock_print.assert_called_with('Missing counties: [3456, 10041, 1001]')

        # check without merged counties
        unique_county_list = geoger.get_county_ids(True, True, False)
        self.assertFalse(geoger.check_for_all_counties(
            unique_county_list, False, False))
        mock_print.assert_called_with(
            'Downloaded data is not complete. Missing 12 counties.')

    def test_get_countyid_to_stateid_map(self):
        countytostate = geoger.get_countyid_to_stateid_map(
            merge_eisenach=True, zfill=False)
        self.assertEqual(all(self.countytostate_string), any(countytostate))
        self.assertEqual(any('16063: 16'), any(countytostate))
        self.assertNotEqual(('16056: 16'), any(countytostate))

        countytostate = geoger.get_countyid_to_stateid_map(
            merge_eisenach=False, zfill=True)
        self.assertEqual(
            all(self.countytostate_zfill_string),
            any(countytostate))
        self.assertEqual(any('16063: 16'), any(countytostate))
        self.assertEqual(any('16056: 16'), any(countytostate))

    def test_get_stateid_to_countyids_map(self):

        # test merge_eisenach = true and zfill = false
        statetocounty = geoger.get_stateid_to_countyids_map(
            merge_eisenach=True, zfill=False)
        self.assertEqual(all(self.statetocounty_string), any(statetocounty))
        self.assertEqual(self.stc_merge_eisenach_true_list, statetocounty[16])
        self.assertNotEqual(
            self.stc_merge_eisenach_false_list, statetocounty[16])
        self.assertEqual(self.stc_zfill_false_list, statetocounty[1])
        self.assertNotEqual(self.stc_zfill_true_list, statetocounty[1])

        # test without merge_eisenach and zfill = true
        statetocounty = geoger.get_countyid_to_stateid_map(
            merge_eisenach=False, zfill=True)
        self.assertEqual(all(self.statetocounty_string), any(statetocounty))
        self.assertEqual(
            all(self.countytostate_zfill_string),
            any(statetocounty))
        self.assertEqual(any('16063: 16'), any(statetocounty))
        self.assertEqual(any('16056: 16'), any(statetocounty))

    def test_get_governing_regions(self):

        # test currently governing regions
        gov_regs = geoger.get_governing_regions(strict=True)
        self.assertEqual(gov_regs, self.gov_regs_true_test_string)

        # test non strict governing regions
        gov_regs = geoger.get_governing_regions(strict=False)
        self.assertEqual(gov_regs, self.gov_regs_false_test_string)

    def test_get_official_county_table(self):

        county_table = geoger.get_official_county_table()
        # test headers of df
        for name in self.county_table_test_headers:
            if(name not in county_table.columns.tolist()):
                self.assertFalse("headers have changed.")

    @patch('builtins.print')
    def test_get_nuts3_county_id_map(self, mock_print):
        # merge_berlin = True, merge_eisenach = False
        nuts_key_dict = geoger.get_nuts3_county_id_map()
        assert 16056 in nuts_key_dict.values()
        assert 11000 in nuts_key_dict.values()
        for id in self.merge_berlin_ids:
            assert int(id) not in nuts_key_dict.values()
        mock_print.assert_called_with(
            'Source data frame contains more counties than official'
            ' county list. This could be OK, please verify yourself.')

    def test_get_intermediateregion_IDs(self):
        # zfill is false
        unique_geo_entitites = geoger.get_intermediateregion_ids(False)
        self.assertEqual(unique_geo_entitites, self.test_list_regions1)

        # zfill is true
        unique_geo_entitites = geoger.get_intermediateregion_ids(True)
        self.assertEqual(unique_geo_entitites, self.test_list_regions2)

    def test_get_intermediateregion_names(self):

        region_names = geoger.get_intermediateregion_names()
        self.assertEqual(region_names, self.test_list_regions3)

    def get_intermediateregion_to_name(self):

        region_to_name = geoger.get_intermediateregion_to_name()
        self.assertEqual(region_to_name, dict(zip(self.test_list_regions1,
                                                  self.test_list_regions3)))

    def test_get_intermediateregion_names_and_ids(self):

        # zfill is false
        regionnamesandids = geoger.get_intermediateregion_names_and_ids(False)
        testregionnamesandids = []
        for i in range(0, len(self.test_list_regions1)):
            combined = [self.test_list_regions3[i], self.test_list_regions1[i]]
            testregionnamesandids.append(combined)
        self.assertEqual(regionnamesandids, testregionnamesandids)

        # zfill is true
        regionnamesandids = geoger.get_intermediateregion_names_and_ids(True)
        testregionnamesandids = []
        for i in range(0, len(self.test_list_regions2)):
            combined = [self.test_list_regions3[i], self.test_list_regions2[i]]
            testregionnamesandids.append(combined)
        self.assertEqual(regionnamesandids, testregionnamesandids)

    def test_get_countyid_to_intermediateregionid_map(self):

        countytoregion = geoger.get_countyid_to_intermediateregionid_map(
            merge_eisenach=True, zfill=False)
        self.assertEqual(all(self.countytoregion_string), any(countytoregion))
        self.assertEqual(any('16063: 6'), any(countytoregion))
        self.assertNotEqual(('16056: 6'), any(countytoregion))

        countytoregion = geoger.get_countyid_to_stateid_map(
            merge_eisenach=False, zfill=True)
        self.assertEqual(
            all(self.countytoregion_zfill_string),
            any(countytoregion))
        self.assertEqual(any('16063: 06'), any(countytoregion))
        self.assertEqual(any('16056: 06'), any(countytoregion))

    def test_get_intermediateregionid_to_countyids_map(self):

        # test merge_eisenach = true and zfill = false
        regiontocounty = geoger.get_intermediateregionid_to_countyids_map(
            merge_eisenach=True, zfill=False)
        self.assertEqual(all(self.regiontocounty_string), any(regiontocounty))
        self.assertEqual(self.rtc_merge_eisenach_true_list, regiontocounty[6])
        self.assertNotEqual(
            self.rtc_merge_eisenach_false_list, regiontocounty[6])
        self.assertEqual(self.rtc_zfill_false_list, regiontocounty[1])
        self.assertNotEqual(self.rtc_zfill_true_list, regiontocounty[1])

        # test without merge_eisenach and zfill = true
        regiontocounty = geoger.get_countyid_to_stateid_map(
            merge_eisenach=False, zfill=True)
        self.assertEqual(all(self.regiontocounty_string), any(regiontocounty))
        self.assertEqual(
            all(self.countytostate_zfill_string),
            any(regiontocounty))
        self.assertEqual(any('16063: 16'), any(regiontocounty))
        self.assertEqual(any('16056: 16'), any(regiontocounty))


if __name__ == '__main__':
    unittest.main()
