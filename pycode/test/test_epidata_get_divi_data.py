import unittest
from pyfakefs import fake_filesystem_unittest
from freezegun import freeze_time
from datetime import date, timedelta

import os
import pandas as pd

from epidemiology.epidata import getDIVIData as gdd
from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from unittest.mock import patch, call


# The following lines are commented to remember a solution to write an output without using the function print()
# This is important, because the usage of print would alter the test results
# import sys
# sys.stdout.write(str())


class TestGetDiviData(fake_filesystem_unittest.TestCase):
    # construct fake directory for testing
    maxDiff = None

    path = '/home/DiviData'

    # strings for read, download and update data
    test_string1 = ("""{"bundesland":1,"gemeindeschluessel":1001,"anzahl_meldebereiche":2,"faelle_covid_aktuell":0,\
"faelle_covid_aktuell_beatmet":0,"anzahl_standorte":2,"betten_frei":48,"betten_belegt":34,\
"daten_stand":"2020-07-07 12:15:00"},\
{"bundesland":2,"gemeindeschluessel":2000,"anzahl_meldebereiche":28,"faelle_covid_aktuell":7,\
"faelle_covid_aktuell_beatmet":6,"anzahl_standorte":24,"betten_frei":396,"betten_belegt":574,\
"daten_stand":"2020-07-07 12:15:00"},\
{"bundesland":3,"gemeindeschluessel":3101,"anzahl_meldebereiche":5,"faelle_covid_aktuell":1,\
"faelle_covid_aktuell_beatmet":1,"anzahl_standorte":5,"betten_frei":60,"betten_belegt":96,\
"daten_stand":"2020-07-07 12:15:00"},\
{"bundesland":3,"gemeindeschluessel":3103,"anzahl_meldebereiche":1,"faelle_covid_aktuell":4,\
"faelle_covid_aktuell_beatmet":1,"anzahl_standorte":1,"betten_frei":11,"betten_belegt":23,\
"daten_stand":"2020-07-07 12:15:00"}""")

    test_string2 = ("""{"bundesland":2,"gemeindeschluessel":2000,"anzahl_meldebereiche":28,"faelle_covid_aktuell":7,\
"faelle_covid_aktuell_beatmet":6,"anzahl_standorte":24,"betten_frei":397,"betten_belegt":579,\
"daten_stand":"2020-07-08 12:15:00"},\
{"bundesland":3,"gemeindeschluessel":3101,"anzahl_meldebereiche":5,"faelle_covid_aktuell":1,\
"faelle_covid_aktuell_beatmet":1,"anzahl_standorte":5,"betten_frei":65,"betten_belegt":91,\
"daten_stand":"2020-07-08 12:15:00"}""")

    test_string3 = ("""{"bundesland":1,"gemeindeschluessel":1001,"anzahl_meldebereiche":2,"faelle_covid_aktuell":0,\
"faelle_covid_aktuell_beatmet":0,"anzahl_standorte":2,"betten_frei":48,"betten_belegt":34,\
"daten_stand":"2020-07-06 12:15:00"}""")

# The following are the same as test_string3 but with changed dates.
    test_string4 = ("""{"bundesland":1,"gemeindeschluessel":1001,"anzahl_meldebereiche":2,"faelle_covid_aktuell":0,\
"faelle_covid_aktuell_beatmet":0,"anzahl_standorte":2,"betten_frei":48,"betten_belegt":34,\
"daten_stand":"2020-04-27 09:15:00"}""")

    test_string5 = ("""{"bundesland":1,"gemeindeschluessel":1001,"anzahl_meldebereiche":2,"faelle_covid_aktuell":0,\
"faelle_covid_aktuell_beatmet":0,"anzahl_standorte":2,"betten_frei":48,"betten_belegt":34,\
"daten_stand":"2020-04-28 09:15:00"}""")

    test_string_read_fulldata = "[" + test_string1 + "," + test_string2 + "]"
    test_string_read_fulldata_update = "[" + test_string3 + "," + test_string1 + "," + test_string2 + "]"
    test_string_read_fulldata_update_nondic = "[" + test_string4 + "," + test_string5 + "]"

    test_string1 = "[" + test_string1 + "]"
    test_string2 = "[" + test_string2 + "]"
    test_string3 = "[" + test_string3 + "]"
    test_string4 = "[" + test_string4 + "]"
    test_string5 = "[" + test_string5 + "]"

    # result string for counties
    test_stringr1_county = ("""\
{"County":"SK Flensburg","ID_County":1001,"ICU":0,"ICU_ventilated":0,"Date":1594124100000},\
{"County":"SK Hamburg","ID_County":2000,"ICU":7,"ICU_ventilated":6,"Date":1594124100000},\
{"County":"SK Braunschweig","ID_County":3101,"ICU":1,"ICU_ventilated":1,"Date":1594124100000},\
{"County":"SK Wolfsburg","ID_County":3103,"ICU":4,"ICU_ventilated":1,"Date":1594124100000}""")

    test_stringr2_county = ("""\
{"County":"SK Hamburg","ID_County":2000,"ICU":7,"ICU_ventilated":6,"Date":1594210500000},\
{"County":"SK Braunschweig","ID_County":3101,"ICU":1,"ICU_ventilated":1,"Date":1594210500000}""")

    test_stringr3_county = ("""\
{"County":"SK Flensburg","ID_County":1001,"ICU":0,"ICU_ventilated":0,"Date":1594037700000}""")

    test_stringr4_county = ("""\
{"County":"SK Flensburg","ID_County":1001,"ICU":0,"ICU_ventilated":0,"Date":1587978900000}""")

    test_stringr5_county = ("""\
{"County":"SK Flensburg","ID_County":1001,"ICU":0,"ICU_ventilated":0,"Date":1588065300000}""")

    # result string for states
    test_stringr1_state = ("""\
{"Date":1594124100000,"ICU":0,"ICU_ventilated":0,"ID_State":1,"State":"Schleswig-Holstein"},\
{"Date":1594124100000,"ICU":7,"ICU_ventilated":6,"ID_State":2,"State":"Hamburg"},\
{"Date":1594124100000,"ICU":5,"ICU_ventilated":2,"ID_State":3,"State":"Niedersachsen"}""")

    test_stringr2_state = ("""\
{"Date":1594124100000,"ICU":0,"ICU_ventilated":0,"ID_State":1,"State":"Schleswig-Holstein"},\
{"Date":1594124100000,"ICU":7,"ICU_ventilated":6,"ID_State":2,"State":"Hamburg"},\
{"Date":1594210500000,"ICU":7,"ICU_ventilated":6,"ID_State":2,"State":"Hamburg"},\
{"Date":1594124100000,"ICU":5,"ICU_ventilated":2,"ID_State":3,"State":"Niedersachsen"},\
{"Date":1594210500000,"ICU":1,"ICU_ventilated":1,"ID_State":3,"State":"Niedersachsen"}""")

    test_stringr3_state = ("""\
{"Date":1594037700000,"ICU":0,"ICU_ventilated":0,"ID_State":1,"State":"Schleswig-Holstein"},\
{"Date":1594124100000,"ICU":0,"ICU_ventilated":0,"ID_State":1,"State":"Schleswig-Holstein"},\
{"Date":1594124100000,"ICU":7,"ICU_ventilated":6,"ID_State":2,"State":"Hamburg"},\
{"Date":1594210500000,"ICU":7,"ICU_ventilated":6,"ID_State":2,"State":"Hamburg"},\
{"Date":1594124100000,"ICU":5,"ICU_ventilated":2,"ID_State":3,"State":"Niedersachsen"},\
{"Date":1594210500000,"ICU":1,"ICU_ventilated":1,"ID_State":3,"State":"Niedersachsen"}""")

    test_stringr4_state = ("""\
{"Date":1587978900000,"ICU":0,"ICU_ventilated":0,"ID_State":1,"State":"Schleswig-Holstein"}""")

    test_stringr5_state = ("""\
{"Date":1588065300000,"ICU":0,"ICU_ventilated":0,"ID_State":1,"State":"Schleswig-Holstein"}""")

    # result string for germany
    test_stringr1_country = ("""\
{"Date":1594124100000,"ICU":12,"ICU_ventilated":8}""")

    test_stringr2_country = ("""\
{"Date":1594210500000,"ICU":8,"ICU_ventilated":7}""")

    test_stringr3_country = ("""\
{"Date":1594037700000,"ICU":0,"ICU_ventilated":0}""")

    test_stringr4_country = ("""\
{"Date":1587978900000,"ICU":0,"ICU_ventilated":0}""")

    test_stringr5_country = ("""\
{"Date":1588065300000,"ICU":0,"ICU_ventilated":0}""")

    # data for test dataframe to test the function adjust_data
    d24 = {'bundesland': [1, 2], 'kreis': [1001, 2000], 'ICU': [0, 7]}
    d25 = {'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7]}
    d262728 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7]}
    d29 = {'Unnamed: 0': [0, 0], 'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7]}
    list_input = [d24, d25, d262728, d262728, d262728, d29]

    # data for result dataframe to test the function adjust_data
    dr24 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7],
            'daten_stand': ["2020-04-24 09:15:00", "2020-04-24 09:15:00"],
            'faelle_covid_aktuell_beatmet': ['', '']}
    dr25 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7],
            'daten_stand': ["2020-04-25 09:15:00", "2020-04-25 09:15:00"],
            'faelle_covid_aktuell_beatmet': ['', '']}
    dr26 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7],
            'daten_stand': ["2020-04-26 09:15:00", "2020-04-26 09:15:00"]}
    dr27 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7],
            'daten_stand': ["2020-04-27 09:15:00", "2020-04-27 09:15:00"]}
    dr2829 = {'bundesland': [1, 2], 'gemeindeschluessel': [1001, 2000], 'ICU': [0, 7]}
    list_result = [dr24, dr25, dr26, dr27, dr2829, dr2829]

    def setUp(self):
        self.setUpPyfakefs()

    def test_gdd_adjust_data(self):

        start_date = date(2020, 4, 24)
        for i in range(6):
            # construct test dataframe
            df = pd.DataFrame(self.list_input[i])

            # adjust dataframe
            df = gdd.adjust_data(df, start_date)
            df.sort_index(axis=1, inplace=True)

            # construct correct result
            df_res = pd.DataFrame(self.list_result[i])
            df_res.sort_index(axis=1, inplace=True)

            # check equality
            self.assertTrue((df == df_res).all().all())

            # increase date by one day
            start_date += timedelta(days=1)

    # note: patches have to have the right order!
    @patch('epidemiology.epidata.getDIVIData.sys.exit')
    @patch('epidemiology.epidata.getDIVIData.pandas.read_csv')
    def test_gdd_call_call_url(self, mock_read_csv, mock_sys_error):

        mock_sys_error.return_value = None

        mock_read_csv.return_value = pd.read_json(self.test_string1)

        # test_string1 has data from 7.7.2020
        download_date = date(2020, 7, 7)
        call_date = download_date.strftime("%Y-%m-%d")
        call_time = "-12-15"
        ext = ""

        url_prefix = call_date + call_time + ext

        call_number = 3961
        df = gdd.call_call_url(url_prefix, call_number)

        call_url = "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" + \
                   "2020-07-07-12-15/viewdocument/3961"

        mock_read_csv.assert_called_with(call_url)

        # check whether the read data is as expected
        self.assertTrue((pd.read_json(self.test_string1) == df).all().all())

        # add an additional print out, which should not change the data
        input_string = "An additional output."
        df = gdd.call_call_url(url_prefix, call_number)

        mock_read_csv.assert_called_with(
            "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" +
            "2020-07-07-12-15/viewdocument/3961")

        # check whether the read data is as expected
        self.assertTrue((pd.read_json(self.test_string1) == df).all().all())

        # check exception when OSError is returned

        mock_read_csv.side_effect = OSError

        gdd.call_call_url(url_prefix, call_number)

        exit_string = "ERROR: URL " + call_url + " could not be opened. " \
                      + "Hint: check your internet connection."

        mock_sys_error.assert_called_with(exit_string)

        # check other exception different to OSError

        mock_read_csv.side_effect = BaseException

        df = gdd.call_call_url(url_prefix, call_number)
        # df should be empty
        self.assertTrue(df.empty)

    test_url_in_call_number_dict = {
        date(2020, 4, 24): ["2020-04-24-09-15", 3974],
        date(2020, 7, 7): ["2020-07-07-12-15", 3961],
        date(2020, 7, 15): ["2020-07-15-12-15", 4108],
    }

    # dates which are not in dict
    test_url_ending_else = {
        date(2020, 5, 7): ["2020-05-07-09-15", 3692],
        date(2020, 6, 4): ["2020-06-04-09-15", 3720],
        date(2020, 6, 7): ["2020-06-07-12-15", 3844],
        date(2020, 6, 15): ["2020-06-15-12-15-2", 3909],
    }

    def fake_call_call_url(self, url_prefix, call_number):

        date_string = url_prefix

        if [date_string, call_number] in self.test_url_in_call_number_dict.values():
            # at this point we do not care about the content, we just need a filled dataframe
            return pd.read_json(self.test_string1)

        elif [date_string, call_number] in self.test_url_ending_else.values():

            # at this point we do not care about the content, we just need a filled dataframe
            return pd.read_json(self.test_string1)

        # otherwise return an empty dataframe
        # this is comparable to call a wrong url
        return pd.DataFrame()

    def test_fake_call_call_url(self):

        df_content = pd.read_json(self.test_string1)
        df_empty = pd.DataFrame()

        # test for test_url_in_call_number_dict
        url_prefix = "2020-04-24-09-15"
        call_number = 3974

        df = self.fake_call_call_url(url_prefix, call_number)

        self.assertTrue((df_content == df).all().all())

        # wrong call number
        call_number = 3973
        df = self.fake_call_call_url(url_prefix, call_number)
        self.assertTrue((df_empty == df).all().all())

        # test for test_url_ending_else
        url_prefix = "2020-06-04-09-15"
        call_number = 3720

        df = self.fake_call_call_url(url_prefix, call_number)

        self.assertTrue((df_content == df).all().all())

        # wrong call number
        call_number = 3973
        df = self.fake_call_call_url(url_prefix, call_number)
        self.assertTrue((df_empty == df).all().all())

    @patch('epidemiology.epidata.getDIVIData.call_call_url')
    def test_gdd_download_data_for_one_day(self, mock_ccu):

        # real call of call_call_url is replaced by fake_call_call_url
        mock_ccu.side_effect = self.fake_call_call_url

        last_number = 0

        # test cases, where date is part of call_number_dict
        for test_date in self.test_url_in_call_number_dict.keys():
            [call_number, df, call_string] = gdd.download_data_for_one_day(last_number, test_date)

            self.assertTrue(call_number == self.test_url_in_call_number_dict[test_date][1])
            self.assertTrue(call_string == "")
            # dataframe should be filled
            self.assertFalse(df.empty)

            # Check of call_call_url would be called with correct parameters
            mock_ccu.assert_called_with(
                self.test_url_in_call_number_dict[test_date][0], self.test_url_in_call_number_dict[test_date][1])

        # test cases, where date has difference 1 to given call_number
        for test_date in self.test_url_ending_else.keys():

            [call_number, df, call_string] = gdd.download_data_for_one_day(self.test_url_ending_else[test_date][1] - 1,
                                                                           test_date)
            self.assertTrue(call_number == self.test_url_ending_else[test_date][1])
            self.assertTrue(call_string == "")
            self.assertFalse(df.empty)

            # Check of call_call_url would be called with correct parameters
            mock_ccu.assert_called_with(
                self.test_url_ending_else[test_date][0], self.test_url_ending_else[test_date][1])

        # test cases, where given call_number has a difference 2 to correct one
        test_date = date(2020, 5, 7)
        test_date_str = self.test_url_ending_else[test_date][0]
        test_call_number = self.test_url_ending_else[test_date][1]

        [call_number, df, call_string] = gdd.download_data_for_one_day(test_call_number - 2, test_date)

        self.assertTrue(call_number == test_call_number)
        self.assertTrue(call_string == "")
        self.assertFalse(df.empty)

        expected_calls = [call.call_call_url(test_date_str, test_call_number - 1),
                          call.call_call_url(test_date_str, test_call_number)]

        # check if expected_calls is part of calls
        mock_ccu.assert_has_calls(expected_calls)
        # check if expected calls are the last two calls
        self.assertTrue(mock_ccu.mock_calls[-2:] == expected_calls)

        # test case where given number has a difference larger than 2 to call_number, namely 4
        call_string_correct = "date(" + test_date.strftime("%Y, %-m, %-d") + "): " + str(test_call_number) + "," + "\n"

        [call_number, df, call_string] = gdd.download_data_for_one_day(self.test_url_ending_else[test_date][1] - 4, test_date)

        self.assertEqual(call_number, test_call_number)
        self.assertEqual(call_string, call_string_correct)
        self.assertFalse(df.empty)

        expected_calls = [call.call_call_url(test_date_str, test_call_number - 3),
                          call.call_call_url(test_date_str, test_call_number - 2),
                          call.call_call_url(test_date_str, test_call_number - 1),
                          call.call_call_url(test_date_str, test_call_number)]

        # check if expected_calls is part of calls
        mock_ccu.assert_has_calls(expected_calls)
        # check if expected calls are the last 4 calls
        self.assertTrue(mock_ccu.mock_calls[-4:] == expected_calls)

    def fake_download_data_for_one_day(self, last_number, start_date):
        data_dict = {
            date(2020, 4, 24): [0, pd.DataFrame(), ""],
            date(2020, 4, 28): [0, pd.read_json(self.test_string5), ""],
            date(2020, 7, 7): [3961, pd.read_json(self.test_string1), ""],
            date(2020, 7, 8): [3964, pd.read_json(self.test_string2), ""],
            date(2020, 7, 9): [0, pd.DataFrame(), ""],
            date(2020, 7, 10): [3987, pd.read_json(self.test_string2), "date(2020, 7, 9): 3987," + "\n"  ]
        }

        return data_dict[start_date]

    # The following tests all belong to function get_divi_data

    # Test download of the full data set
    @patch('builtins.print')
    @patch('epidemiology.epidata.getDIVIData.download_data_for_one_day')
    def test_gdd_download_data(self, mock_ddfod, mock_print):

        mock_ddfod.return_value = [3961, pd.read_json(self.test_string1), ""]

        # case simple download
        [read_data, update_data, make_plot, out_form, out_folder, start_date, end_date]\
            = [False, False, False, "json", self.path, date(2020, 7, 7), date(2020, 7, 7)]

        directory = os.path.join(out_folder, 'Germany/')

        file = "FullData_DIVI.json"

        file_out1 = "county_divi.json"
        file_out2 = "state_divi.json"
        file_out3 = "germany_divi.json"

        gdd.get_divi_data(read_data, out_form, out_folder,
                          end_date, start_date, update_data)

        # check if folder Germany exists now
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ['Germany'])

        # check that four files are generated
        self.assertEqual(len(os.listdir(directory)), 4)
        self.assertEqual(os.listdir(directory).sort(), [file, file_out1, file_out2, file_out3].sort())

        # check content of files
        f_path = os.path.join(directory, file)
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_string1)

        f_path = os.path.join(directory, file_out1)
        f = open(f_path, "r")
        fread = f.read()
        self.assertEqual(fread, "[" + self.test_stringr1_county + "]")

        f_path = os.path.join(directory, file_out2)
        f = open(f_path, "r")
        self.assertEqual(f.read(), "[" + self.test_stringr1_state + "]")

        f_path = os.path.join(directory, file_out3)
        f = open(f_path, "r")
        self.assertEqual(f.read(), "[" + self.test_stringr1_country + "]")

        # test for more than one day
        # TODO: It would be better to split the test here

        mock_ddfod.side_effect = self.fake_download_data_for_one_day

        start_date = date(2020, 7, 7)
        end_date = date(2020, 7, 9)

        gdd.get_divi_data(read_data, out_form, out_folder,
                          end_date, start_date, update_data)

        # check if folder Germany exists now
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ['Germany'])

        # check that four files are generated
        self.assertEqual(len(os.listdir(directory)), 4)
        self.assertEqual(os.listdir(directory).sort(), [file, file_out1, file_out2, file_out3].sort())

        # check content of files
        f_path = os.path.join(directory, file)
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_string_read_fulldata)

        f_path = os.path.join(directory, file_out1)
        f = open(f_path, "r")
        fread = f.read()
        self.assertEqual(fread, "[" + self.test_stringr1_county + "," + self.test_stringr2_county + "]")

        f_path = os.path.join(directory, file_out2)
        f = open(f_path, "r")
        self.assertEqual(f.read(), "[" + self.test_stringr2_state + "]")

        f_path = os.path.join(directory, file_out3)
        f = open(f_path, "r")
        self.assertEqual(f.read(), "[" + self.test_stringr1_country + "," + self.test_stringr2_country + "]")

        start_string = "Success: Data of date "
        end_string = " has been added to dataframe"
        expected_calls = [call(start_string + date(2020, 7, 7).strftime("%Y-%m-%d") + end_string),
                          call(start_string + date(2020, 7, 8).strftime("%Y-%m-%d") + end_string),
                          call("Warning: Data of date " + date(2020, 7, 9).strftime("%Y-%m-%d")
                               + " is not added to dataframe"),
                          call("Information: DIVI data has been written to", directory)]

        # check if expected_calls is part of calls
        mock_print.assert_has_calls(expected_calls)
        # check if expected calls are the last three calls
        self.assertTrue(mock_print.mock_calls[-4:] == expected_calls)

        # Check output if whole dataframe is empty
        with self.assertRaises(SystemExit) as cm:
            start_date = date(2020, 7, 9)
            end_date = date(2020, 7, 9)
            gdd.get_divi_data(read_data, out_form, out_folder,
                              end_date, start_date, update_data)
        self.assertEqual(cm.exception.code, "Something went wrong, dataframe is empty.")

        # Check warning for wrong start_date
        with self.assertRaises(SystemExit) as cm:
            start_date = date(2020, 4, 23)
            end_date = date(2020, 4, 24)
            gdd.get_divi_data(read_data, out_form, out_folder,
                              end_date, start_date, update_data)
        self.assertEqual(cm.exception.code, "Something went wrong, dataframe is empty.")

        expected_call = call("Warning: First data available on 2020-04-24. You asked for 2020-04-23.")

        # check if second last call is the expected warning
        self.assertTrue(mock_print.mock_calls[-2] == expected_call)

        # check case where call_number has to be searched for and should be saved to call_dict.
        start_date = date(2020, 7, 10)
        end_date = date(2020, 7, 10)
        gdd.get_divi_data(read_data, out_form, out_folder,
                          end_date, start_date, update_data)

        expected_calls = [call(start_string + date(2020, 7, 10).strftime("%Y-%m-%d") + end_string),
                          call("New drifting number in link found. "
                               "To decrease runtime, please copy the following "
                               "to the dcitionary \"call_number_dict\" in the function \"download_data_for_one_day\": "),
                          call("date(2020, 7, 9): 3987," + "\n"),
                          call("Information: DIVI data has been written to", directory)]

        # check if expected_calls is part of calls
        mock_print.assert_has_calls(expected_calls)

        # check if expected calls are the last three calls
        self.assertTrue(mock_print.mock_calls[-4:] == expected_calls)

    def test_gdd_read_data(self):

        [read_data, update_data, make_plot, out_form, out_folder, start_date, end_date] \
            = [True, False, False, "json", self.path, dd.defaultDict['start_date'], dd.defaultDict['end_date']]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        file = "FullData_DIVI.json"
        file_with_path = os.path.join(directory, file)

        file_out1 = "county_divi.json"
        file_out2 = "state_divi.json"
        file_out3 = "germany_divi.json"

        # Test case where file does not exist
        with self.assertRaises(SystemExit) as cm:
            gdd.get_divi_data(read_data, out_form, out_folder,
                              end_date, start_date, update_data)
        self.assertEqual(cm.exception.code, "Error: The file: " + file_with_path + " does not exist. "
                                                                                   "Call program without -r or -u flag to get it.")

        # Test case with empty file
        with open(file_with_path, 'w') as f:
            f.write("[]")

        with self.assertRaises(SystemExit) as cm:
            gdd.get_divi_data(read_data, out_form, out_folder,
                              end_date, start_date, update_data)
        self.assertEqual(cm.exception.code, "Something went wrong, dataframe is empty.")

        # Generate data to read in
        with open(file_with_path, 'w') as f:
            f.write(self.test_string_read_fulldata)

        gdd.get_divi_data(read_data, out_form, out_folder,
                          end_date, start_date, update_data)

        self.assertEqual(len(os.listdir(directory)), 4)
        self.assertEqual(os.listdir(directory).sort(), [file, file_out1, file_out2, file_out3].sort())

        f_path = os.path.join(directory, file_out1)
        f = open(f_path, "r")
        fread = f.read()
        self.assertEqual(fread, "[" + self.test_stringr1_county + "," + self.test_stringr2_county + "]")

        f_path = os.path.join(directory, file_out2)
        f = open(f_path, "r")
        self.assertEqual(f.read(), "[" + self.test_stringr2_state + "]")

        f_path = os.path.join(directory, file_out3)
        f = open(f_path, "r")
        self.assertEqual(f.read(), "[" + self.test_stringr1_country + "," + self.test_stringr2_country + "]")

    # freeze_time changes the output of date.today() from today to the given date
    @freeze_time("2020-7-8")
    @patch('epidemiology.epidata.getDIVIData.download_data_for_one_day')
    @patch('builtins.print')
    @patch('epidemiology.epidata.getDIVIData.gd.loadCsv')
    def test_gdd_update_data(self, mock_loadcsv, mock_print, mock_ddfod):

        [read_data, out_form, out_folder] = [False, "json", self.path]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        file = "FullData_DIVI.json"
        file_with_path = os.path.join(directory, file)

        file_out1 = "county_divi.json"
        file_out2 = "state_divi.json"
        file_out3 = "germany_divi.json"

        # Test case where file does not exist
        with self.assertRaises(SystemExit) as cm:
            gdd.get_divi_data(read_data, out_form, out_folder,
                              end_date=dd.defaultDict['end_date'], start_date=dd.defaultDict['start_date'],
                              update_data=True)
        self.assertEqual(cm.exception.code, "Error: The file: " + file_with_path + 
                         " does not exist. Call program without -r or -u flag to get it.")

        # Test case with empty file
        with open(file_with_path, 'w') as f:
            f.write("[]")

        with self.assertRaises(SystemExit) as cm:
            gdd.get_divi_data(read_data, out_form, out_folder,
                              end_date=dd.defaultDict['end_date'], start_date=dd.defaultDict['start_date'],
                              update_data=True)
        self.assertEqual(cm.exception.code, "Something went wrong, dataframe is empty.")

        # write file with data for one day
        with open(file_with_path, 'w') as f:
            f.write(self.test_string1)

        # test where "data of today" is missing, but data is not yet uploaded:
        mock_loadcsv.return_value = pd.read_json(self.test_string1)

        with self.assertRaises(SystemExit) as cm:
            gdd.get_divi_data(read_data, out_form, out_folder,
                              end_date=dd.defaultDict['end_date'], start_date=dd.defaultDict['start_date'],
                              update_data=True)
        self.assertEqual(cm.exception.code, "Data of today = " + date(2020, 7, 8).strftime("%Y-%m-%d")
                         + " has not yet uploaded. Please, try again later.")

        # case where data is online and just data of "today" is missing
        mock_loadcsv.return_value = pd.read_json(self.test_string2)

        gdd.get_divi_data(read_data, out_form, out_folder,
                          end_date=dd.defaultDict['end_date'], start_date=dd.defaultDict['start_date'],
                          update_data=True)

        mock_loadcsv.assert_called_with('DIVI-Intensivregister-Tagesreport', apiUrl='https://www.divi.de/')

        # check content of files
        f = open(file_with_path, "r")
        self.assertEqual(f.read(), self.test_string_read_fulldata)

        f_path = os.path.join(directory, file_out1)
        f = open(f_path, "r")
        fread = f.read()
        self.assertEqual(fread, "[" + self.test_stringr1_county + "," + self.test_stringr2_county + "]")

        f_path = os.path.join(directory, file_out2)
        f = open(f_path, "r")
        self.assertEqual(f.read(), "[" + self.test_stringr2_state + "]")

        f_path = os.path.join(directory, file_out3)
        f = open(f_path, "r")
        self.assertEqual(f.read(), "[" + self.test_stringr1_country + "," + self.test_stringr2_country + "]")

        start_string = "Success: Data of date "
        end_string = " has been added to dataframe"
        expected_calls = [call(start_string + date(2020, 7, 8).strftime("%Y-%m-%d") + end_string),
                          call("Information: DIVI data has been written to", directory)]

        # check if expected_calls is part of calls
        mock_print.assert_has_calls(expected_calls)
        self.assertTrue(mock_print.mock_calls == expected_calls)

        # check case where more than today is missing

        # write data of 6.7.2020 to FullData_DIVI.json
        with open(file_with_path, 'w') as f:
            f.write(self.test_string3)

        mock_ddfod.side_effect = self.fake_download_data_for_one_day

        gdd.get_divi_data(read_data, out_form, out_folder,
                          end_date=date.today(), start_date=dd.defaultDict['start_date'],
                          update_data=True)

        # check that four files are generated
        self.assertEqual(len(os.listdir(directory)), 4)
        self.assertEqual(os.listdir(directory).sort(), [file, file_out1, file_out2, file_out3].sort())

        # check content of files
        f = open(file_with_path, "r")
        fread = f.read()
        self.assertEqual(fread, self.test_string_read_fulldata_update)

        f_path = os.path.join(directory, file_out1)
        f = open(f_path, "r")
        fread = f.read()
        self.assertEqual(fread, "[" + self.test_stringr3_county + "," + self.test_stringr1_county + ","
                         + self.test_stringr2_county + "]")

        f_path = os.path.join(directory, file_out2)
        f = open(f_path, "r")
        self.assertEqual(f.read(), "[" + self.test_stringr3_state + "]")

        f_path = os.path.join(directory, file_out3)
        f = open(f_path, "r")
        self.assertEqual(f.read(),
                         "[" + self.test_stringr3_country + "," + self.test_stringr1_country + 
                         "," + self.test_stringr2_country + "]")

        start_string = "Success: Data of date "
        end_string = " has been added to dataframe"
        expected_calls = [call(start_string + date(2020, 7, 7).strftime("%Y-%m-%d") + end_string),
                          call(start_string + date(2020, 7, 8).strftime("%Y-%m-%d") + end_string),
                          call("Information: DIVI data has been written to", directory)]

        # check if expected_calls is part of calls
        mock_print.assert_has_calls(expected_calls)

        # check if expected calls are the last three calls
        self.assertTrue(mock_print.mock_calls[-3:] == expected_calls)

        # Check what happens if end_date (in real today) is not in dict
        # and the call_number has a difference > 2 to check that no message is printed

        # write data which will be read in
        with open(file_with_path, 'w') as f:
            f.write(self.test_string4)

        mock_ddfod.side_effect = self.fake_download_data_for_one_day

        gdd.get_divi_data(read_data, out_form, out_folder,
                          end_date=date(2020,4,28), start_date=date(2020,4,24),
                          update_data=True)

        # check that four files are generated
        self.assertEqual(len(os.listdir(directory)), 4)
        self.assertEqual(os.listdir(directory).sort(), [file, file_out1, file_out2, file_out3].sort())

        # check content of files
        f = open(file_with_path, "r")
        fread = f.read()
        self.assertEqual(fread, self.test_string_read_fulldata_update_nondic)

        # compare county data
        f_path = os.path.join(directory, file_out1)
        f = open(f_path, "r")
        fread = f.read()
        self.assertEqual(fread, "[" + self.test_stringr4_county + "," + self.test_stringr5_county + "]")

        # compare state data
        f_path = os.path.join(directory, file_out2)
        f = open(f_path, "r")
        self.assertEqual(f.read(), "[" + self.test_stringr4_state + "," + self.test_stringr5_state + "]")

        # compare country data
        f_path = os.path.join(directory, file_out3)
        f = open(f_path, "r")
        self.assertEqual(f.read(),
                         "[" + self.test_stringr4_country + "," + self.test_stringr5_country + "]")

        start_string = "Success: Data of date "
        end_string = " has been added to dataframe"
        expected_calls = [call(start_string + date(2020, 4, 28).strftime("%Y-%m-%d") + end_string),
                          call("Information: DIVI data has been written to", directory)]

        # check if expected_calls is part of calls
        mock_print.assert_has_calls(expected_calls)

        # check if expected calls are the last two calls
        self.assertTrue(mock_print.mock_calls[-2:] == expected_calls)


def test_nearest_earlier_date(self):

    list_of_dates = [date(2020, 4, 24),
                     date(2020, 5, 6),
                     date(2020, 6, 5),
                     date(2020, 6, 12),
                     date(2020, 6, 26),
                     date(2020, 6, 28),
                     date(2020, 6, 29),
                     date(2020, 6, 30)]

    date_to_search = date(2020, 6, 24)

    date_ce = gdd.nearest_earlier_date(list_of_dates, date_to_search)

    self.assertEqual(date_ce, date(2020, 6, 12))


if __name__ == '__main__':
    unittest.main()
