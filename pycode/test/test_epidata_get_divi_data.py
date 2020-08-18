import unittest
from pyfakefs import fake_filesystem_unittest
from datetime import date, timedelta

import os
import pandas as pd

from epidemiology.epidata import getDIVIData as gdd
from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from unittest.mock import patch, call


class Test_getDiviData(fake_filesystem_unittest.TestCase):
    # construct fake directory for testing
    path = '/home/DiviData'

    # strings for read, download and update data
    test_string1 = ("""[\
{"ID_State":1,"ID_County":1001,"anzahl_meldebereiche":2,"ICU":0,"ICU_ventilated":0,\
"reporting_hospitals":2,"ICU_free":48,"ICU_occupied":34,"Date":"2020-07-07 12:15:00"},\
{"ID_State":2,"ID_County":2000,"anzahl_meldebereiche":28,"ICU":7,"ICU_ventilated":6,\
"reporting_hospitals":24,"ICU_free":396,"ICU_occupied":574,"Date":"2020-07-07 12:15:00"},\
{"ID_State":3,"ID_County":3101,"anzahl_meldebereiche":5,"ICU":1,"ICU_ventilated":1,\
"reporting_hospitals":5,"ICU_free":60,"ICU_occupied":96,"Date":"2020-07-07 12:15:00"},\
{"ID_State":3,"ID_County":3103,"anzahl_meldebereiche":1,"ICU":4,"ICU_ventilated":1,\
"reporting_hospitals":1,"ICU_free":11,"ICU_occupied":23,"Date":"2020-07-07 12:15:00"}""")

    test_string2 = ("""
{"ID_State":2,"ID_County":2000,"anzahl_meldebereiche":28,"ICU":7,"ICU_ventilated":6,\
"reporting_hospitals":24,"ICU_free":397,"ICU_occupied":597,"Date":"2020-07-08 12:15:00"},\
{"ID_State":3,"ID_County":3101,"anzahl_meldebereiche":5,"ICU":1,"ICU_ventilated":1,\
"reporting_hospitals":5,"ICU_free":65,"ICU_occupied":91,"Date":"2020-07-08 12:15:00"}]""")

    test_string = test_string1 + "," + test_string2

    test_string_update1 = test_string1 + "]"

    test_string_update2 = "[" + test_string2

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
    d262728 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7]}
    d29 = {'Unnamed: 0': [0, 0], 'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7]}
    liste_input = [d24, d25, d262728, d262728, d262728, d29]

    # data for result dataframe to test adjust data
    dr24 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7],
            'daten_stand': ["2020-04-24 09:15:00", "2020-04-24 09:15:00"]}
    dr25 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7],
            'daten_stand': ["2020-04-25 09:15:00", "2020-04-25 09:15:00"]}
    dr26 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7],
            'daten_stand': ["2020-04-26 09:15:00", "2020-04-26 09:15:00"]}
    dr27 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7],
            'daten_stand': ["2020-04-27 09:15:00", "2020-04-27 09:15:00"]}
    dr2829 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7]}
    liste_result = [dr24, dr25, dr26, dr27, dr2829, dr2829]

    def setUp(self):
        self.setUpPyfakefs()

    def test_gdd_adjust_data(self):

        start_date = date(2020, 4, 24)
        for i in range(6):
            # construct test dataframe
            df = pd.DataFrame(self.liste_input[i])

            # adjust dataframe
            df = gdd.adjust_data(df, start_date)
            df.sort_index(axis=1, inplace=True)

            # construct correct result
            df_res = pd.DataFrame(self.liste_result[i])
            df_res.sort_index(axis=1, inplace=True)

            # check equality
            self.assertTrue((df == df_res).all().all())

            start_date += timedelta(days=1)

    # note: patches have to have the right order!
    @patch('builtins.print')
    @patch('epidemiology.epidata.getDIVIData.pandas.read_csv')
    def test_gdd_call_call_url(self, mock_read_csv, mock_print):

        mock_read_csv.return_value = pd.read_json(self.test_string_update1)

        # test_string1 has data from 7.7.2020
        download_date = date(2020,7,7)
        call_date = download_date.strftime("%Y-%m-%d")
        call_time = "-12-15"
        ext = ""

        url_prefix = "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" \
                 + call_date + call_time + ext

        call_number = 3961
        df =  gdd.call_call_url(url_prefix, call_number)

        call_url = "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" +\
                   "2020-07-07-12-15/viewdocument/3961"

        mock_read_csv.assert_called_with(call_url)

        # check wether the read data is as expected
        self.assertTrue((pd.read_json(self.test_string_update1) == df).all().all())

        # add an additional print out, which should not change the data
        input_string = "An additional output."
        df = gdd.call_call_url(url_prefix, call_number, input_string)
        # check if the prints are as expected.
        mock_print.assert_called_with("New cal number found. " +
                                      "Please copy the following to call_number_dict to increase runtime: " +
                                      input_string)

        mock_read_csv.assert_called_with(
            "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" +
            "2020-07-07-12-15/viewdocument/3961")

        # check wether the read data is as expected
        self.assertTrue((pd.read_json(self.test_string_update1) == df).all().all())

        # TODO: check exceptions
        # check exception when OSError is returned

        mock_read_csv.side_effect = OSError
        #exit_string = "ERROR: URL " + call_url + " could not be opened. " \
        #             + "Hint: check your internet connection."



    @patch('epidemiology.epidata.getDIVIData.call_call_url')
    def test_gdd_download_data_for_one_day(self, mock_ccu):

        mock_ccu.return_value = pd.read_json(self.test_string)

        test_url_in_call_number_dict = {
            date(2020, 4, 24): ["2020-04-24-09-15", 3974],
            date(2020, 7, 7): ["2020-07-07-12-15",3961],
            date(2020, 7, 15): ["2020-07-15-12-15", 4108],
            }

        test_url_ending_else = {
            date(2020, 5, 7): ["2020-05-07-09-15", 3692],
            date(2020, 6, 4): ["2020-06-04-09-15", 3720],
            date(2020, 6, 7): ["2020-06-07-12-15", 3844],
            date(2020, 6, 15): ["2020-06-15-12-15-2", 3909],
        }

        last_number = 0

        # test cases, where date is part of call_number_dict
        for test_date in test_url_in_call_number_dict.keys():

            df = gdd.download_data_for_one_day(last_number, test_date)

            # cases where date is in call_number_dict

            mock_ccu.assert_called_with(
                "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" +
                test_url_in_call_number_dict[test_date][0], test_url_in_call_number_dict[test_date][1])


        # test cases, where date has difference 1 to given call_number
        for test_date in test_url_ending_else.keys():

            df = gdd.download_data_for_one_day(test_url_ending_else[test_date][1]-1, test_date)

            mock_ccu.assert_called_with(
                "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" +
                 test_url_ending_else[test_date][0], test_url_ending_else[test_date][1], "")

        # test cases, where date has difference 2 to given call_number
        for test_date in test_url_ending_else.keys():
            df = gdd.download_data_for_one_day(test_url_ending_else[test_date][1] - 2, test_date)

            expected_calls = [call.call_call_url("https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" +
                test_url_ending_else[test_date][0], test_url_ending_else[test_date][1]-1, "") ,call.call_call_url("https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" +
                test_url_ending_else[test_date][0], test_url_ending_else[test_date][1], "")]

            # TODO: assert the loop calls when difference is == 2 AND afterwards if it is larger
            #mock_ccu.assert_has_calls(expected_calls)

        #self.assertTrue((pd.read_json(self.test_string) == df).all().all())

    @patch('epidemiology.epidata.getDIVIData.pandas.read_csv')
    def test_gdd_download_data(self, mock_read_csv):

        mock_read_csv.return_value = pd.read_json(self.test_string)

        [read_data, update_data, make_plot, out_form, out_folder] = [False, False, False, "json", self.path]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        file = "FullData_DIVI.json"

        file_out1 = "county_divi.json"
        file_out2 = "state_divi.json"
        file_out3 = "germany_divi.json"

        gdd.get_divi_data(read_data, update_data, make_plot, out_form, out_folder,
                          start_date=date(2020, 7, 7), end_date=date(2020, 7, 7))

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

        [read_data, update_data, make_plot, out_form, out_folder] = [True, False, False, "json", self.path]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        file = "FullData_DIVI.json"

        file_out1 = "county_divi.json"
        file_out2 = "state_divi.json"
        file_out3 = "germany_divi.json"

        with open(os.path.join(directory, file), 'w') as f:
            f.write(self.test_string)

        gdd.get_divi_data(read_data, update_data, make_plot, out_form, out_folder)

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

    #@patch('epidemiology.epidata.getPopulationData.gd.loadCsv')
    #def test_gdd_update_data(self, mock_loadCSV):

     #   mock_loadCSV.return_value = pd.read_json(self.test_string_update2)

      #  [read_data, update_data, make_plot, out_form, out_folder] = [False, True, False, "json", self.path]

       # directory = os.path.join(out_folder, 'Germany/')
        #gd.check_dir(directory)

        #file = "FullData_DIVI.json"

        #file_out1 = "county_divi.json"
        #file_out2 = "state_divi.json"
        #file_out3 = "germany_divi.json"

        #with open(os.path.join(directory, file), 'w') as f:
        #    f.write(self.test_string_update1)

        #gdd.get_divi_data(read_data, update_data, make_plot, out_form, out_folder)

        #self.assertEqual(len(os.listdir(directory)), 4)
        #self.assertEqual(os.listdir(directory).sort(), [file, file_out1, file_out2, file_out3].sort())

        #f_path = os.path.join(directory, file_out1)
        #f = open(f_path, "r")
        #self.assertEqual(f.read(), self.test_stringr1)

        #f_path = os.path.join(directory, file_out2)
        #f = open(f_path, "r")
        #self.assertEqual(f.read(), self.test_stringr2)

        #f_path = os.path.join(directory, file_out3)
        #f = open(f_path, "r")
        #self.assertEqual(f.read(), self.test_stringr3)

if __name__ == '__main__':
    unittest.main()
