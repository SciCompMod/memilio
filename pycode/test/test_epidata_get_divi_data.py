import unittest
from pyfakefs import fake_filesystem_unittest
from datetime import date, timedelta


import os
import pandas as pd

from epidemiology.epidata import getDIVIData as gdd
from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from unittest.mock import patch, mock_open


class Test_getDiviData(fake_filesystem_unittest.TestCase):

    # construct fake directory for testing
    path = '/home/DiviData'

    # string for read and download data
    test_string = ("""[\
{"ID_State":1,"ID_County":1001,"anzahl_meldebereiche":2,"ICU":0,"ICU_ventilated":0,\
"reporting_hospitals":2,"ICU_free":48,"ICU_occupied":34,"Date":"2020-07-07 12:15:00"},\
{"ID_State":2,"ID_County":2000,"anzahl_meldebereiche":28,"ICU":7,"ICU_ventilated":6,\
"reporting_hospitals":24,"ICU_free":396,"ICU_occupied":574,"Date":"2020-07-07 12:15:00"},\
{"ID_State":3,"ID_County":3101,"anzahl_meldebereiche":5,"ICU":1,"ICU_ventilated":1,\
"reporting_hospitals":5,"ICU_free":60,"ICU_occupied":96,"Date":"2020-07-07 12:15:00"},\
{"ID_State":3,"ID_County":3103,"anzahl_meldebereiche":1,"ICU":4,"ICU_ventilated":1,\
"reporting_hospitals":1,"ICU_free":11,"ICU_occupied":23,"Date":"2020-07-07 12:15:00"},\
{"ID_State":2,"ID_County":2000,"anzahl_meldebereiche":28,"ICU":7,"ICU_ventilated":6,\
"reporting_hospitals":24,"ICU_free":397,"ICU_occupied":597,"Date":"2020-07-08 12:15:00"},\
{"ID_State":3,"ID_County":3101,"anzahl_meldebereiche":5,"ICU":1,"ICU_ventilated":1,\
"reporting_hospitals":5,"ICU_free":65,"ICU_occupied":91,"Date":"2020-07-08 12:15:00"}]""")

    # result string for counties
    test_stringr1 = ("""[\
{"County":"SK Flensburg","ID_County":1001,"ICU":0,"ICU_ventilated":0,"Date":1594124100000},\
{"County":"SK Hamburg","ID_County":2000,"ICU":7,"ICU_ventilated":6,"Date":1594124100000},\
{"County":"SK Braunschweig","ID_County":3101,"ICU":1,"ICU_ventilated":1,"Date":1594124100000},\
{"County":"SK Wolfsburg","ID_County":3103,"ICU":4,"ICU_ventilated":1,"Date":1594124100000},\
{"County":"SK Hamburg","ID_County":2000,"ICU":7,"ICU_ventilated":6,"Date":1594210500000},\
{"County":"SK Braunschweig","ID_County":3101,"ICU":1,"ICU_ventilated":1,"Date":1594210500000}]""")

    # result string for states
    test_stringr2 = ("""[\
{"Date":1594124100000,"ICU":0,"ICU_ventilated":0,"ID_State":1,"State":"Schleswig-Holstein"},\
{"Date":1594124100000,"ICU":7,"ICU_ventilated":6,"ID_State":2,"State":"Hamburg"},\
{"Date":1594210500000,"ICU":7,"ICU_ventilated":6,"ID_State":2,"State":"Hamburg"},\
{"Date":1594124100000,"ICU":5,"ICU_ventilated":2,"ID_State":3,"State":"Niedersachsen"},\
{"Date":1594210500000,"ICU":1,"ICU_ventilated":1,"ID_State":3,"State":"Niedersachsen"}]""")

    # result string for germany
    test_stringr3 = ("""[\
{"Date":1594124100000,"ICU":12,"ICU_ventilated":8},\
{"Date":1594210500000,"ICU":8,"ICU_ventilated":7}]""")

    # data for test dataframe to test adjust data
    d24 = {'bundesland': [1, 2], 'kreis': [1001, 2000], 'ICU': [0, 7]}
    d25 = {'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7]}
    d262728 = {'bundesland': [1, 2],'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7]}
    d29 = {'Unnamed: 0': [0,0],'bundesland': [1, 2],'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7]}
    liste_input = [d24,d25,d262728,d262728,d262728,d29]

    # data for result dataframe to test adjust data
    dr24 = {'bundesland': [1, 2],'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7], 'daten_stand':["2020-04-24 09:15:00","2020-04-24 09:15:00"]}
    dr25 = {'bundesland': [1, 2],'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7], 'daten_stand':["2020-04-25 09:15:00","2020-04-25 09:15:00"]}
    dr26 = {'bundesland': [1, 2],'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7], 'daten_stand':["2020-04-26 09:15:00","2020-04-26 09:15:00"]}
    dr27 = {'bundesland': [1, 2],'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7], 'daten_stand':["2020-04-27 09:15:00","2020-04-27 09:15:00"]}
    dr2829 = {'bundesland': [1, 2],'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7]}
    liste_result = [dr24,dr25,dr26,dr27,dr2829,dr2829]

    def setUp(self):
        self.setUpPyfakefs()

# TODO: test download data for one day

    def test_gdd_adjust_data(self):

        start_date = date(2020,4,24)
        for i in range(6):
            # construct test dataframe
            df = pd.DataFrame(self.liste_input[i])

            # adjust dataframe
            df = gdd.adjust_data(df, start_date)
            print(df)
            df.sort_index(axis=1, inplace=True)

            # construct result
            df_res = pd.DataFrame(self.liste_result[i])
            print(df_res)
            df_res.sort_index(axis=1, inplace=True)

            # check equality
            self.assertTrue((df == df_res).all().all())

            start_date += timedelta(days=1)


    @patch('getDIVIData.pandas.read_csv')
    def test_gdd_download_data(self, mock_read_csv):

        mock_read_csv.return_value = pd.read_json(self.test_string)

        [read_data, make_plot, out_form, out_folder] = [False, False, "json", self.path]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        file = "FullData_DIVI.json"

        file_out1 = "county_divi.json"
        file_out2 = "state_divi.json"
        file_out3 = "germany_divi.json"

        gdd.get_divi_data(read_data, make_plot, out_form, out_folder, start_date=date(2020,7,7), end_date=date(2020,7,7))

        self.assertEqual(len(os.listdir(directory)), 4)
        self.assertEqual(os.listdir(directory).sort(), [file, file_out1, file_out2, file_out3].sort())

        f_path = os.path.join(directory, file_out1)
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_stringr1)

        f_path = os.path.join(directory, file_out2)
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_stringr2)

        f_path = os.path.join(directory, file_out3)
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_stringr3)

    def test_gdd_read_data(self):

        [read_data, make_plot, out_form, out_folder] = [True, False, "json", self.path]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        file = "FullData_DIVI.json"

        file_out1 = "county_divi.json"
        file_out2 = "state_divi.json"
        file_out3 = "germany_divi.json"

        with open(os.path.join(directory, file), 'w') as f:
            f.write(self.test_string)

        gdd.get_divi_data(read_data, make_plot, out_form, out_folder)

        self.assertEqual(len(os.listdir(directory)), 4)
        self.assertEqual(os.listdir(directory).sort(), [file, file_out1, file_out2, file_out3].sort())

        f_path = os.path.join(directory, file_out1)
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_stringr1)

        f_path = os.path.join(directory, file_out2)
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_stringr2)
        
        f_path = os.path.join(directory, file_out3)
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_stringr3)

if __name__ == '__main__':
    unittest.main()
