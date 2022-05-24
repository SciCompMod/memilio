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

from pyfakefs import fake_filesystem_unittest
from unittest.mock import patch
import unittest
import pandas as pd
import os
import collections
from os.path import dirname as up
from memilio.epidata import getCommuterMobility as gcm
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gD
from memilio.epidata import defaultDict as dd
from memilio.epidata import getPopulationData as gpd


class TestCommuterMigration(fake_filesystem_unittest.TestCase):
    path = '/home/CMData/'

    setup_dict = {'num_counties': 2, 'abs_tol': 100, 'rel_tol': 0.01,
                  'num_govregions': 2, 'counties': 'empty', 'path': path}

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
    
    memilio_path = up(up(up(up(up(__file__)))))
    center_coordinates = pd.read_json(os.path.join(
        memilio_path, 'data', 'mobility' , 'county_centers_dim400.json'))

    def setUp(self):
        self.setUpPyfakefs()

    def write_kreise_deu_data(self, out_folder):
        # sheet 0 is unused in commuter_migration_bfa, but other one has to have index 1
        sheet0 = pd.DataFrame(
            {'0': ['0', '0', '0', '0'], '1': ['1', '2', '3', '4']})
        # 'nothing' strings used to interpret keys as Strings instead of numbers
        sheet1 = pd.DataFrame({'Schlüssel-nummer': ['nothing', '01', '01001', '01053', '02', '02000', '03', '03101',
                                                    '04', '04012', '05', '051', '05112', '053', '05316', '06', '06532',
                                                    '06632', '07', '071', '07141',
                                                    '07235', '08', '081', '0811', '08111', '09', '091', '09181', '10',
                                                    '10044', '11',
                                                    '11000', '12', '12051', '12066', '13', '13003', '14', '145',
                                                    '14511', '15',
                                                    '15082',
                                                    '16', '16069'],

                               'Regionale Bezeichnung': ['', 'Schleswig-Holstein', 'Kreisfreie Stadt', 'Kreis',
                                                         'Hamburg',
                                                         'Kreisfreie Stadt', 'Niedersachsen', 'Kreisfreie Stadt',
                                                         'Bremen',
                                                         'Kreisfreie Stadt', 'Nordrhein-Westfalen',
                                                         'Reg.-Bez. Düsseldorf',
                                                         'Kreisfreie Stadt', 'Reg.-Bez. Köln', 'Kreisfreie Stadt',
                                                         'Hessen', 'Landkreis',
                                                         'Landkreis', 'Rheinland-Pfalz', 'früher: Reg.-Bez. Koblenz',
                                                         'Landkreis',
                                                         'Landkreis', 'Baden-Württemberg', 'Reg.-Bez. Stuttgart',
                                                         'Region Stuttgart',
                                                         'Stadtkreis', 'Bayern', 'Reg.-Bez. Oberbayern', 'Landkreis',
                                                         'Saarland',
                                                         'Landkreis', 'Berlin', 'Kreisfreie Stadt', 'Brandenburg',
                                                         'Kreisfreie Stadt',
                                                         'Landkreis', 'Mecklenburg-Vorpommern', 'Kreisfreie Stadt',
                                                         'Sachsen', 'Direktionsbezirk Chemnitz',
                                                         'Kreisfreie Stadt', 'Sachsen-Anhalt', 'Landkreis', 'Thüringen',
                                                         'Landkreis'],

                               'Kreis / Landkreis': ['', '', 'Flensburg, Stadt', 'Herzogtum Lauenburg', '',
                                                     'Hamburg, Freie und Hansestadt', '', 'Braunschweig, Stadt', '',
                                                     'Bremerhaven, Stadt', '', '', 'Duisburg, Stadt', '',
                                                     'Leverkusen, Stadt',
                                                     '', 'Lahn-Dill-Kreis', 'Hersfeld-Rotenburg', '', '',
                                                     'Rhein-Lahn-Kreis',
                                                     'Trier-Saarburg', '', '', '', 'Stuttgart, Stadtkreis', '',
                                                     '', 'Landsberg am Lech', '', 'Saarlouis', '', 'Berlin, Stadt', '',
                                                     'Brandenburg an der Havel, Stadt', 'Oberspreewald-Lausitz', '',
                                                     'Rostock', '', '', 'Chemnitz, Stadt', '', 'Anhalt-Bitterfeld', '',
                                                     'Hildburghausen'],

                               'NUTS3': ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                                         '0', '0', '0', '0',
                                         '0', '0',
                                         '0',
                                         '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                                         '0', '0', '0',
                                         '0', '0',
                                         '0', '0', '0'],
                               'Fläche in km2': ['', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                                                 '0', '0', '0',
                                                 '0', '0', '0',
                                                 '0', '0', '0',
                                                 '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                                                 '0', '0', '0', '0',
                                                 '0', '0',
                                                 '0', '0', '0'],

                               'insgesamt': ['', '', 90164, 198019, '', 1847253, '', 249406, '', 113643, '', '', 498686,
                                             '', 163729, '',
                                             253319, 120719, '', '', 122297, 149398, '', '', '', 635911, '', '', 120302,
                                             '', 194319, '', 3669491, '', 72184, 109371, '', 209191, '', '', 246334, '',
                                             158486,
                                             '', 63197]})

        sheets = {'Deckblatt': sheet0,
                  'Kreisfreie Städte u. Landkreise': sheet1}
        data_path = os.path.join(out_folder, 'kreise_deu.xlsx')
        dummy = pd.ExcelWriter(data_path)
        for sheet_name in sheets.keys():
            sheets[sheet_name].to_excel(
                dummy, sheet_name=sheet_name, index=False)
        dummy.save()
        dummy.close()

    def write_commuter_all_federal_states(self, out_folder):
        # just 3rd sheet is interesting
        sheet0 = pd.DataFrame(
            {'0': ['0', '0', '0', '0'], '1': ['1', '2', '3', '4']})

        sheet1 = pd.DataFrame({'Arbeitsort': ['nothing', '01001', '', '', '', '', '', '', '', '', '', '', '01053', ''],
                               'Arbeitsort2': ['', 'Flensburg, Stadt', '', '', '', '', '', '', '', '', '', '',
                                               'Herzogtum Lauenburg', ''],
                               'Wohnort': ['nothing', '', '01053', '02', '02000', '03', '04', '05',
                                           '051', '05112', '053', '05316', '', '01001'],
                               'Wohnort2': ['', '', 'Herzogtum Lauenburg', 'Hamburg',
                                            'Hamburg, Freie und Hansestadt', 'Niedersachsen', 'Bremen',
                                            'Nordrhein-Westfalen', 'Reg.-Bez. Düsseldorf', 'Duisburg, Stadt',
                                            'Reg.-Bez. Köln', 'Leverkusen, Stadt', '', 'Flensburg, Stadt'],
                               'Insgesamt': ['', '', 30, 272, 272, 543, 48, 251, 76, 13, 52, 17, '', 17]})
        sheet2 = pd.DataFrame({'Arbeitsort': ['nothing', '02000', '', '', ''],
                               'Arbeitsort2': ['', 'Hamburg, Freie und Hansestadt', '', '', ''],
                               'Wohnort': ['nothing', '', '01053', '04', '06532'],
                               'Wohnort2': ['', '', 'Herzogtum Lauenburg', 'Bremen',
                                            'Lahn-Dill-Kreis'],
                               'Insgesamt': ['', '', 23230, 4628, 92]})
        sheet3 = pd.DataFrame({'Arbeitsort': ['nothing', '03101', '', '', ''],
                               'Arbeitsort2': ['', 'Braunschweig, Stadt', '', '', ''],
                               'Wohnort': ['nothing', '', '01053', '07', '14511'],
                               'Wohnort2': ['', '', 'Herzogtum Lauenburg', 'Rheinland-Pfalz',
                                            'Chemnitz, Stadt'],
                               'Insgesamt': ['', '', 17, 123, 15]})
        sheet4 = pd.DataFrame({'Arbeitsort': ['nothing', '04012', '', '', ''],
                               'Arbeitsort2': ['', 'Bremerhaven, Stadt', '', '', ''],
                               'Wohnort': ['nothing', '', '06532', '05', '12051'],
                               'Wohnort2': ['', '', 'Lahn-Dill-Kreis', 'Nordrhein-Westfalen',
                                            'Brandenburg an der Havel, St.'],
                               'Insgesamt': ['', '', 22, 299, 10]})
        sheet5 = pd.DataFrame({'Arbeitsort': ['nothing', '05112', '', '', '', '05316', '', ''],
                               'Arbeitsort2': ['', 'Duisburg, Stadt', '', '', '', 'Leverkusen, Stadt', '', ''],
                               'Wohnort': ['nothing', '', '12066', '13003', 'ZZ', '', '14511', 'ZZ'],
                               'Wohnort2': ['', '', 'Oberspreewald-Lausitz', 'Rostock, Hansestadt',
                                            'Einpendler aus dem Bundesgebiet',
                                            '', 'Chemnitz, Stadt', 'Einpendler aus dem Bundesgebiet'],
                               'Insgesamt': ['', '', 34, 299, 305, '', 12, 12]})
        sheet6 = pd.DataFrame({'Arbeitsort': ['nothing', '06532', '', '', '06632', '', ''],
                               'Arbeitsort2': ['', 'Lahn-Dill-Kreis', '', '', 'Hersfeld-Rotenburg', '', ''],
                               'Wohnort': ['nothing', '', '02000', '01053', '', '02000', '091'],
                               'Wohnort2': ['', '', 'Hamburg, Freie und Hansestadt', 'Herzogtum Lauenburg',
                                            '', 'Hamburg, Freie und Hansestadt', 'Übrige Kreise (Regierungsbezirk)'],
                               'Insgesamt': ['nothing', '', 112, 13, '', 28, 47]})
        sheet7 = pd.DataFrame({'Arbeitsort': ['nothing', '07141', '', '', '07235', ''],
                               'Arbeitsort2': ['', 'Rhein-Lahn-Kreis', '', '', 'Trier-Saarburg', ''],
                               'Wohnort': ['nothing', '', '02000', '07235', '', '02000'],
                               'Wohnort2': ['', '', 'Hamburg, Freie und Hansestadt', 'Trier-Saarburg',
                                            '', 'Hamburg, Freie und Hansestadt'],
                               'Insgesamt': ['', '', 18, 19, '', 10]})
        sheet8 = pd.DataFrame({'Arbeitsort': ['nothing', '08111', '', ''],
                               'Arbeitsort2': ['', 'Stuttgart, Landeshauptstadt', '', ''],
                               'Wohnort': ['nothing', '', '01053', '13003'],
                               'Wohnort2': ['', '', 'Herzogtum Lauenburg', 'Rostock, Hansestadt'],
                               'Insgesamt': ['', '', 32, 38]})
        sheet9 = pd.DataFrame({'Arbeitsort': ['nothing', '09181', '', ''],
                               'Arbeitsort2': ['', 'Landsberg am Lech', '', ''],
                               'Wohnort': ['nothing', '', '08111', '11000'],
                               'Wohnort2': ['', '', 'Stuttgart, Landeshauptstadt', 'Berlin, Stadt'],
                               'Insgesamt': ['', '', 57, 43]})
        sheet10 = pd.DataFrame({'Arbeitsort': ['nothing', '10044', ''],
                                'Arbeitsort2': ['', 'Saarlouis', ''],
                                'Wohnort': ['nothing', '', '07235'],
                                'Wohnort2': ['', '', 'Trier-Saarburg'],
                                'Insgesamt': ['', '', 220]})
        sheet11 = pd.DataFrame({'Arbeitsort': ['nothing', '11000', '', ''],
                                'Arbeitsort2': ['', 'Berlin, Stadt', '', ''],
                                'Wohnort': ['nothing', '', '01053', '08111'],
                                'Wohnort2': ['', '', 'Herzogtum Lauenburg', 'Stuttgart, Landeshauptstadt'],
                                'Insgesamt': ['', '', 297, 1186]})
        sheet12 = pd.DataFrame({'Arbeitsort': ['nothing', '12051', '', '', '12066', ''],
                                'Arbeitsort2': ['', 'Brandenburg an der Havel, St.', '', '', 'Oberspreewald-Lausitz',
                                                ''],
                                'Wohnort': ['nothing', '', '11000', '12066', '', '05112'],
                                'Wohnort2': ['', '', 'Berlin, Stadt', 'Oberspreewald-Lausitz', '', 'Duisburg, Stadt'],
                                'Insgesamt': ['', '', 571, 24, '', 10]})
        sheet13 = pd.DataFrame({'Arbeitsort': ['nothing', '13003', '', ''],
                                'Arbeitsort2': ['', 'Rostock, Hansestadt', '', ''],
                                'Wohnort': ['nothing', '', '01053', '12066'],
                                'Wohnort2': ['', '', 'Herzogtum Lauenburg', 'Oberspreewald-Lausitz'],
                                'Insgesamt': ['', '', 42, 12]})
        sheet14 = pd.DataFrame({'Arbeitsort': ['nothing', '14511', '', ''],
                                'Arbeitsort2': ['', 'Chemnitz, Stadt', '', ''],
                                'Wohnort': ['nothing', '', '02000', '11000'],
                                'Wohnort2': ['', '', 'Hamburg, Freie und Hansestadt', 'Berlin, Stadt'],
                                'Insgesamt': ['', '', 74, 524]})
        sheet15 = pd.DataFrame({'Arbeitsort': ['nothing', '15082', '', ''],
                                'Arbeitsort2': ['', 'Anhalt-Bitterfeld', '', ''],
                                'Wohnort': ['nothing', '', '12066', '14511'],
                                'Wohnort2': ['', '', 'Oberspreewald-Lausitz', 'Chemnitz, Stadt'],
                                'Insgesamt': ['', '', 27, 26]})
        sheet16 = pd.DataFrame({'Arbeitsort': ['nothing', '16069', '', ''],
                                'Arbeitsort2': ['', 'Hildburghausen', '', ''],
                                'Wohnort': ['nothing', '', '03', '05'],
                                'Wohnort2': ['', '', 'Niedersachsen', 'Nordrhein-Westfalen'],
                                'Insgesamt': ['', '', 29, 34]})
        states = (
            sheet1, sheet2, sheet3, sheet4, sheet5, sheet6, sheet7, sheet8, sheet9, sheet10, sheet11, sheet12, sheet13,
            sheet14, sheet15, sheet16)
        for i in range(16):
            sheets = {'Deckblatt': sheet0, 'Impressum': sheet0, 'Auspendler Kreise': sheet0,
                      'Einpendler Kreise': states[i]}
            if i < 9:
                name = 'krpend_0' + str(i + 1) + "_0.xlsx"
            else:
                name = 'krpend_' + str(i + 1) + "_0.xlsx"
            data_path = os.path.join(out_folder, name)
            dummy = pd.ExcelWriter(data_path)
            for sheet_name in sheets.keys():
                sheets[sheet_name].to_excel(
                    dummy, sheet_name=sheet_name, index=False)
            dummy.save()
            dummy.close()

    @patch('builtins.print')
    def test_verify_sorted(self, mock_print):
        self.assertEqual(True, gcm.verify_sorted(self.test_countykey_list))
        self.assertEqual(False, gcm.verify_sorted(self.test_countykey_list2))
        Errorcall = ('Error. Input list not sorted.')
        mock_print.assert_called_with(Errorcall)

    @patch('builtins.print')
    def test_assign_geographical_entities(self, mock_print):
        (countykey2govkey, countykey2localnumlist, gov_county_table,
         state_gov_table) = gcm.assign_geographical_entities(self.countykey_list, self.govkey_list)
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

    @patch('builtins.print')
    def test_some_errors(self, mock_print):
        gD.check_dir(self.path)
        self.write_kreise_deu_data(self.path)
        self.write_commuter_all_federal_states(self.path)
        self.assertEqual(len(os.listdir(self.path)), 17)

    def test_commuter_data(self):
        """! Tests migration data by some randomly chosen tests.
        """

        df_commuter_migration = gcm.get_commuter_data(
            out_folder=self.path, center_coordinates=self.center_coordinates)
        # the first column and first row are just the county IDs
        # mat_commuter_migration is the real Data that should be tested
        mat_commuter_migration = df_commuter_migration.iloc[:, 0:]
        mat_commuter_migration = mat_commuter_migration.iloc[0:, :]

        countykey2numlist = collections.OrderedDict(
            zip(self.countykey_list, list(range(0, len(self.countykey_list)))))

        gD.check_dir(self.path)
        self.write_kreise_deu_data(self.path)
        self.write_commuter_all_federal_states(self.path)
        self.assertEqual(len(os.listdir(self.path)), 18)

        # just do some tests on randomly chosen migrations

        # check migration from Leverkusen (averaged from NRW, 05) to Hildburghausen
        city_from = countykey2numlist['05316']
        city_to = countykey2numlist['16069']
        population = gpd.get_population_data(
            out_folder=self.path, merge_eisenach=False)
        countypop_list = list(population[dd.EngEng['population']])
        self.assertEqual(countypop_list[city_from], 163905)
        self.assertAlmostEqual(
            mat_commuter_migration.iat[city_from, city_to], 0.305, 2)

        # check migration from Duisburg to Oberspreewald-Lausitz
        city_from = countykey2numlist['05112']
        city_to = countykey2numlist['12066']
        self.assertEqual(mat_commuter_migration.iat[city_from, city_to], 2.2)

        # check migration from Lahn-Dill-Kreis to Hamburg
        city_from = countykey2numlist['06532']
        city_to = countykey2numlist['02000']
        self.assertAlmostEqual(mat_commuter_migration.iat[city_from, city_to], 19.2, 2)

        # check migration from Herzogtum Lauenburg to Flensburg, Stadt
        city_from = countykey2numlist['01001']
        city_to = countykey2numlist['01053']
        self.assertEqual(mat_commuter_migration.iat[city_from, city_to], 7)

    @patch('builtins.print')
    def test_get_neighbors_mobility(self, mock_print):

        testcountyid = 1051
        tci = testcountyid
        #direction = both
        (countykey_list, commuter_all) = gcm.get_neighbors_mobility(
            tci, direction='both', abs_tol=0, rel_tol=0,
            tol_comb='or', merge_eisenach=True, out_folder=self.path, 
            center_coordinates=self.center_coordinates)
        self.assertEqual(len(countykey_list), 398)
        self.assertAlmostEqual(228, commuter_all[0], 2)
        self.assertAlmostEqual(2146, commuter_all[9], 2)
        self.assertAlmostEqual(146.5, commuter_all[11], 2)
        self.assertAlmostEqual(.2, commuter_all[397], 2)

        # direction = in
        (countykey_list, commuter_all) = gcm.get_neighbors_mobility(
            tci, direction='in', abs_tol=0, rel_tol=0,
            tol_comb='or', merge_eisenach=True, out_folder=self.path,
            center_coordinates=self.center_coordinates)
        self.assertEqual(len(countykey_list), 393)
        self.assertAlmostEqual(48, commuter_all[0], 2)
        self.assertAlmostEqual(842, commuter_all[9], 2)
        self.assertAlmostEqual(92, commuter_all[11], 2)

        # direction = out
        (countykey_list, commuter_all) = gcm.get_neighbors_mobility(
            tci, direction='out', abs_tol=0, rel_tol=0,
            tol_comb='or', merge_eisenach=True, out_folder=self.path,
            center_coordinates=self.center_coordinates)
        self.assertEqual(len(countykey_list), 375)
        self.assertAlmostEqual(180, commuter_all[0], 2)
        self.assertAlmostEqual(1304, commuter_all[9], 2)
        self.assertAlmostEqual(201, commuter_all[11], 2)


if __name__ == '__main__':
    unittest.main()
