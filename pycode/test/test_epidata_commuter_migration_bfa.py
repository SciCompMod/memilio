import unittest
from pyfakefs import fake_filesystem_unittest
from mock import patch
import pandas as pd
import os
from epidemiology.epidata import commuter_migration_bfa as cm
from epidemiology.epidata import getDataIntoPandasDataFrame as gD


class TestCommuterMigration(fake_filesystem_unittest.TestCase):
    path = '/home/CMData/'

    def setUp(self):
        self.setUpPyfakefs()

    def write_kreise_deu_daten(self, out_folder):
        # sheet 0 is unused in commuter_migration_bfa, but other one has to have index 1
        sheet0 = pd.DataFrame({'0': ['0', '0', '0', '0'], '1': ['1', '2', '3', '4']})
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

        sheets = {'Deckblatt': sheet0, 'Kreisfreie Städte u. Landkreise': sheet1}
        data_path = os.path.join(out_folder, 'kreise_deu.xlsx')
        dummy = pd.ExcelWriter(data_path)
        for sheet_name in sheets.keys():
            sheets[sheet_name].to_excel(dummy, sheet_name=sheet_name, index=False)
        dummy.save()

    def write_commuter_all_federal_states(self, out_folder):
        # just 3rd sheet is interesting
        sheet0 = pd.DataFrame({'0': ['0', '0', '0', '0'], '1': ['1', '2', '3', '4']})

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
                sheets[sheet_name].to_excel(dummy, sheet_name=sheet_name, index=False)
            dummy.save()

    @patch('builtins.print')
    def test_some_errors(self, mock_print):
        gD.check_dir(self.path)
        self.write_kreise_deu_daten(self.path)
        self.write_commuter_all_federal_states(self.path)
        self.assertEqual(len(os.listdir(self.path)), 17)
        counties = gD.loadExcel('kreise_deu', apiUrl=self.path, extension='.xlsx', sheet_name=1)
        setup_dict = {'num_counties': 2,
                      'abs_tol': 100,
                      'rel_tol': 0.01,
                      'num_govregions': 2,
                      'counties': counties,
                      'path': self.path}
        cm.map_keys_to_numlists(setup_dict)
        mock_print.assert_any_call('Error. Number of government regions wrong. Having', 17, 'instead of', 2)
        mock_print.assert_any_call("Error. Number of counties wrong.")
        mock_print.assert_any_call('Error. Number of governing regions wrong.')

        mock_print.call_args_list = []
        cm.assign_geographical_entities(setup_dict)
        mock_print.assert_any_call('Error. Number of government regions wrong.')

        mock_print.call_args_list = []
        cm.verify_sorted(countykey_list=pd.Series(['01001', '01053', '02000', '03101', '06632', '04012', '05112',
                                                   '05316', '06532', '07141']))
        mock_print.assert_any_call('Error. Input list not sorted, population per county list had to be sorted '
                                   'accordingly.')

        mock_print.call_args_list = []
        setup_dict['num_counties'] = 21
        setup_dict['num_govregions'] = 17
        cm.get_matrix_commuter_migration_patterns(setup_dict)
        mock_print.assert_any_call('Error in calculations for county ', 'Duisburg, Stadt', '\nAccumulated values:',
                                   333.0, ', correct sum:', 305.0)

    def test_migration_data(self):
        """! Tests migration data by some randomly chosen tests.
        """
        gD.check_dir(self.path)
        self.write_kreise_deu_daten(self.path)
        self.write_commuter_all_federal_states(self.path)
        self.assertEqual(len(os.listdir(self.path)), 17)
        counties = gD.loadExcel('kreise_deu', apiUrl=self.path, extension='.xlsx', sheet_name=1)
        setup_dict = {'num_counties': 21,
                      'abs_tol': 100,
                      'rel_tol': 0.01,
                      'num_govregions': 17,
                      'counties': counties,
                      'path': self.path}
        (countykey_list, countypop_list, govkey_list, countykey2numlist, govkey2numlist, gov_county_table,
         countykey2govkey,
         countykey2localnumlist, state_gov_table, mat_commuter_migration) = cm.get_data(setup_dict)
        # just do some tests on randomly chosen migrations

        # check migration from Leverkusen (averaged from NRW, 05) to Hildburghausen
        city_from = countykey2numlist['05316']
        city_to = countykey2numlist['16069']
        self.assertEqual(countypop_list[city_from], 163729)
        self.assertEqual(mat_commuter_migration[city_from][city_to], 34 * countypop_list[city_from] / (498686 + 163729))

        # check migration from Duisburg to Oberspreewald-Lausitz
        city_from = countykey2numlist['05112']
        city_to = countykey2numlist['12066']
        self.assertEqual(mat_commuter_migration[city_from][city_to], 10)

        # check migration from Lahn-Dill-Kreis to Hamburg
        city_from = countykey2numlist['06532']
        city_to = countykey2numlist['02000']
        self.assertEqual(mat_commuter_migration[city_from][city_to], 92)

        # check migration from Landsberg am Lech (averaged from 091) to Hersfeld-Rotenburg
        city_from = countykey2numlist['09181']
        city_to = countykey2numlist['06632']
        self.assertEqual(mat_commuter_migration[city_from][city_to], 47)

        # check migration from Herzogtum Lauenburg to Flensburg, Stadt
        city_from = countykey2numlist['01001']
        city_to = countykey2numlist['01053']
        self.assertEqual(mat_commuter_migration[city_from][city_to], 17)


if __name__ == '__main__':
    unittest.main()
