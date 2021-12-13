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
import unittest
from pyfakefs import fake_filesystem_unittest
import os
import sys
from io import StringIO
from unittest.mock import patch, call, mock_open
from datetime import date, datetime
import pandas as pd

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd

from epidemiology.epidata import getVaccinationData
from epidemiology.epidata import getPopulationData
from epidemiology.epidata import getRKIData
from epidemiology.epidata import getDIVIData
from epidemiology.epidata import getRKIDatawithEstimations
from epidemiology.epidata import getJHData


class Test_getDataIntoPandasDataFrame(fake_filesystem_unittest.TestCase):

    path = '/home/x/'

    data = ("""{"type": "FeatureCollection",\
"name": "RKI_COVID19",\
"features": [\
{ "type": "Feature", "properties": { "ObjectId": 1, "IdBundesland": 1, "Bundesland": "Schleswig-Holstein",\
"Landkreis": "SK Flensburg", "Altersgruppe": "A15-A34", "Geschlecht": "M", "AnzahlFall": 1, "AnzahlTodesfall": 0,\
"Meldedatum": "2021-03-26T00:00:00Z", "IdLandkreis": "01001", "Datenstand": "20.04.2021, 00:00 Uhr", "NeuerFall": 0, "NeuerTodesfall": -9, "Refdatum": "2021-03-22T00:00:00Z", "NeuGenesen": 0, "AnzahlGenesen": 1, "IstErkrankungsbeginn": 1, "Altersgruppe2": "Nicht übermittelt" }, "geometry": null },\
{ "type": "Feature", "properties": { "ObjectId": 2, "IdBundesland": 1, "Bundesland": "Schleswig-Holstein",
"Landkreis": "SK Flensburg", "Altersgruppe": "A15-A34", "Geschlecht": "M", "AnzahlFall": 7, "AnzahlTodesfall": 0,\
"Meldedatum": "2021-03-26T00:00:00Z", "IdLandkreis": "01001", "Datenstand": "20.04.2021, 00:00 Uhr", "NeuerFall": 0, "NeuerTodesfall": -9, "Refdatum": "2021-03-26T00:00:00Z", "NeuGenesen": 0, "AnzahlGenesen": 7, "IstErkrankungsbeginn": 0, "Altersgruppe2": "Nicht übermittelt" }, "geometry": null },\
{ "type": "Feature", "properties": { "ObjectId": 3, "IdBundesland": 1, "Bundesland": "Schleswig-Holstein",\
"Landkreis": "SK Flensburg", "Altersgruppe": "A15-A34", "Geschlecht": "M", "AnzahlFall": 1, "AnzahlTodesfall": 0,\
"Meldedatum": "2021-03-26T00:00:00Z", "IdLandkreis": "01001", "Datenstand": "20.04.2021, 00:00 Uhr", "NeuerFall": 0, "NeuerTodesfall": -9, "Refdatum": "2021-03-26T00:00:00Z", "NeuGenesen": -9, "AnzahlGenesen": 0, "IstErkrankungsbeginn": 0, "Altersgruppe2": "Nicht übermittelt" }, "geometry": null }\
]}""")

    def setUp(self):
        self.setUpPyfakefs()
        del sys.argv[1:] # In this unit tests parse_args is called when it is called through unittests. This has a lot of command lines which lead to errors in getDataIntoPandasDataframe.
        # TODO:Is this is a good way to solve this?

    @patch('epidemiology.epidata.getDataIntoPandasDataFrame.urlopen')
    def test_load_geojson_error(self, mock_urlopen):

        mock_urlopen.side_effect = OSError

        with self.assertRaises(SystemExit) as cm:
            gd.loadGeojson("targetFileName")

        exit_string = "ERROR: URL " + 'https://opendata.arcgis.com/datasets/' + \
            "targetFileName" + '.geojson' + " could not be opened."

        self.assertEqual(cm.exception.code, exit_string)

    @patch('epidemiology.epidata.getDataIntoPandasDataFrame.pd.read_excel')
    def test_load_Excel_error(self, mock_urlopen):

        mock_urlopen.side_effect = OSError

        with self.assertRaises(SystemExit) as cm:
            gd.loadExcel("targetFileName")

        exit_string = "ERROR: URL " + 'https://opendata.arcgis.com/datasets/' + \
            "targetFileName" + '.xls' + " could not be opened."

        self.assertEqual(cm.exception.code, exit_string)

    @patch('epidemiology.epidata.getDataIntoPandasDataFrame.pd.read_csv')
    def test_load_csv_error(self, mock_csv):
        # return an empty dataframe
        mock_csv.return_value = pd.DataFrame()

        df_test = gd.loadCsv("targetFileName")
        expected_call = [
            call(
                'https://opendata.arcgis.com/datasets/' +
                "targetFileName" + '.csv')]
        mock_csv.assert_has_calls(expected_call)

        assert df_test.empty

        df_test = gd.loadCsv(
            "targetFileName",
            apiUrl='https://opendata.arcgis.com/datasets/different/',
            extension='.notcsv')
        expected_call = [
            call(
                'https://opendata.arcgis.com/datasets/different/' +
                "targetFileName" + '.notcsv')]
        mock_csv.assert_has_calls(expected_call)

        assert df_test.empty

    @patch('epidemiology.epidata.getDataIntoPandasDataFrame.pd.read_csv')
    def test_load_csv_working(self, mock_csv):
        mock_csv.side_effect = OSError

        with self.assertRaises(SystemExit) as cm:
            gd.loadCsv("targetFileName")

        exit_string = "ERROR: URL " + 'https://opendata.arcgis.com/datasets/' + \
            "targetFileName" + '.csv' + " could not be opened."

        self.assertEqual(cm.exception.code, exit_string)

    def test_cli_correct_default(self):
        
        out_path_default = dd.defaultDict['out_folder']
        
        arg_dict = gd.cli("population")
        [read_data, file_format, out_folder, no_raw] = [
            arg_dict["read_data"],
            arg_dict["file_format"],
            arg_dict["out_folder"],
            arg_dict["no_raw"]]

        assert read_data == dd.defaultDict['read_data']
        assert file_format == dd.defaultDict['file_format']
        assert out_folder == out_path_default
        assert no_raw == dd.defaultDict['no_raw']
        
        arg_dict = gd.cli("jh")
        [read_data, file_format, out_folder, no_raw] = [arg_dict["read_data"],
                                                        arg_dict["file_format"], arg_dict["out_folder"], arg_dict["no_raw"]]

        assert read_data == dd.defaultDict['read_data']
        assert file_format == dd.defaultDict['file_format']
        assert out_folder == out_path_default
        assert no_raw == dd.defaultDict['no_raw']

        arg_dict = gd.cli("rki")
        [read_data, file_format, out_folder, impute_dates, make_plot,
         moving_average, split_berlin, no_raw, rep_date] = [
            arg_dict["read_data"],
            arg_dict["file_format"],
            arg_dict["out_folder"],
            arg_dict["impute_dates"],
            arg_dict["make_plot"],
            arg_dict["moving_average"],
            arg_dict["split_berlin"],
            arg_dict["no_raw"],
            arg_dict["rep_date"]]

        assert read_data == dd.defaultDict['read_data']
        assert file_format == dd.defaultDict['file_format']
        assert out_folder == out_path_default
        assert make_plot == dd.defaultDict['make_plot']
        assert split_berlin == dd.defaultDict['split_berlin']
        assert moving_average == dd.defaultDict['moving_average']
        assert impute_dates == dd.defaultDict['impute_dates']
        assert no_raw == dd.defaultDict['no_raw']
        assert rep_date == dd.defaultDict['rep_date']

        arg_dict = gd.cli("rkiest")
        [read_data, file_format, make_plot, out_folder, no_raw] = [
            arg_dict["read_data"],
            arg_dict["file_format"],
            arg_dict["make_plot"],
            arg_dict["out_folder"],
            arg_dict["no_raw"]]

        assert read_data == dd.defaultDict['read_data']
        assert file_format == dd.defaultDict['file_format']
        assert out_folder == out_path_default
        assert no_raw == dd.defaultDict['no_raw']
        assert make_plot == dd.defaultDict['make_plot']

        arg_dict = gd.cli("commuter_official")
        [read_data, file_format, out_folder, no_raw] = [
            arg_dict["read_data"],
            arg_dict["file_format"],
            arg_dict["out_folder"],
            arg_dict["no_raw"]]

        assert read_data == dd.defaultDict['read_data']
        assert file_format == dd.defaultDict['file_format']
        assert out_folder == out_path_default
        assert no_raw == dd.defaultDict['no_raw']

        arg_dict = gd.cli("divi")
        [read_data, file_format, out_folder, end_date, start_date,
         impute_dates, moving_average, no_raw] = [arg_dict["read_data"],
                                                  arg_dict["file_format"],
                                                  arg_dict["out_folder"],
                                                  arg_dict["end_date"],
                                                  arg_dict["start_date"],
                                                  arg_dict["impute_dates"],
                                                  arg_dict["moving_average"],
                                                  arg_dict["no_raw"]]

        assert read_data == dd.defaultDict['read_data']
        assert file_format == dd.defaultDict['file_format']
        assert out_folder == out_path_default
        assert end_date == dd.defaultDict['end_date']
        assert start_date == dd.defaultDict['start_date']
        assert impute_dates == dd.defaultDict['impute_dates']
        assert moving_average == dd.defaultDict['moving_average']
        assert no_raw == dd.defaultDict['no_raw']

        arg_dict = gd.cli("sim")
        [read_data, file_format, out_folder, end_date, make_plot, impute_dates,
         moving_average, split_berlin, start_date, no_raw] = [
            arg_dict["read_data"],
            arg_dict["file_format"],
            arg_dict["out_folder"],
            arg_dict["end_date"],
            arg_dict["impute_dates"],
            arg_dict["make_plot"],
            arg_dict["moving_average"],
            arg_dict["split_berlin"],
            arg_dict["start_date"],
            arg_dict["no_raw"]]

        assert read_data == dd.defaultDict['read_data']
        assert file_format == dd.defaultDict['file_format']
        assert out_folder == out_path_default
        assert no_raw == dd.defaultDict['no_raw']
        assert end_date == dd.defaultDict['end_date']
        assert impute_dates == dd.defaultDict['impute_dates']
        assert make_plot == dd.defaultDict['make_plot']
        assert moving_average == dd.defaultDict['moving_average']
        assert split_berlin == dd.defaultDict['split_berlin']
        assert start_date == dd.defaultDict['start_date']

    @patch('sys.stderr', new_callable=StringIO)
    def test_cli_correct_raise_exit(self, mock_stderr):

        with self.assertRaises(SystemExit) as cm:
            gd.cli("wrong_key")

        the_exception = cm.exception
        self.assertEqual(the_exception.code, "Wrong key or cli_dict.")

        test_args = ["prog", '-f', 'wrong_format']
        with patch.object(sys, 'argv', test_args):

            with self.assertRaises(SystemExit) as cm:
                gd.cli("rki")
            self.assertRegexpMatches(mock_stderr.getvalue(), r"invalid choice")

        test_args = ["prog", '--update-data', ]
        with patch.object(sys, 'argv', test_args):
            with self.assertRaises(SystemExit) as cm:
                gd.cli("rki")
            self.assertRegexpMatches(
                mock_stderr.getvalue(),
                r"unrecognized arguments")

        test_args = ["prog", '--start_date']
        with patch.object(sys, 'argv', test_args):
            with self.assertRaises(SystemExit) as cm:
                gd.cli("rki")
            self.assertRegexpMatches(
                mock_stderr.getvalue(),
                r"unrecognized arguments")

        test_args = ["prog", '--end_date']
        with patch.object(sys, 'argv', test_args):
            with self.assertRaises(SystemExit) as cm:
                gd.cli("rki")
            self.assertRegexpMatches(
                mock_stderr.getvalue(),
                r"unrecognized arguments")

        test_args = ["prog", '--make-plot']
        with patch.object(sys, 'argv', test_args):
            with self.assertRaises(SystemExit) as cm:
                gd.cli("divi")
            self.assertRegexpMatches(
                mock_stderr.getvalue(),
                r"unrecognized arguments")

        test_args = ["prog", '--start-date', '2020,3,24']
        with patch.object(sys, 'argv', test_args):
            with self.assertRaises(SystemExit) as cm:
                gd.cli("divi")

        test_args = ["prog", '--end-date', '24-3-2020']
        with patch.object(sys, 'argv', test_args):
            with self.assertRaises(SystemExit) as cm:
                gd.cli("divi")

    def test_cli_set_different_values(self):

        folder = "some_folder"

        test_args = ["prog", '--read-data', '--out-folder',
                     folder, '--file-format', 'hdf5', '--no-raw']

        with patch.object(sys, 'argv', test_args):

            arg_dict = gd.cli("population")
            [read_data, file_format, out_folder, no_raw] = [
                arg_dict["read_data"],
                arg_dict["file_format"],
                arg_dict["out_folder"],
                arg_dict["no_raw"]]

            assert read_data == True
            assert file_format == 'hdf5'
            assert out_folder == "some_folder"
            assert no_raw == True

            arg_dict = gd.cli("jh")
            [read_data, file_format, out_folder, no_raw] = [
                arg_dict["read_data"],
                arg_dict["file_format"],
                arg_dict["out_folder"],
                arg_dict["no_raw"]]

            assert read_data == True
            assert file_format == 'hdf5'
            assert out_folder == "some_folder"
            assert no_raw == True

        test_args = ["prog", '--read-data', '--out-folder', folder,
                     '--file-format', 'hdf5', '--make-plot', '--split-berlin',
                     '--moving-average', 0, '--no-raw', '--impute-dates']

        with patch.object(sys, 'argv', test_args):
            arg_dict = gd.cli("rki")
            [read_data, file_format, out_folder, impute_dates, make_plot,
             moving_average, split_berlin, no_raw, rep_date] = [
                arg_dict["read_data"],
                arg_dict["file_format"],
                arg_dict["out_folder"],
                arg_dict["impute_dates"],
                arg_dict["make_plot"],
                arg_dict["moving_average"],
                arg_dict["split_berlin"],
                arg_dict["no_raw"],
                arg_dict["rep_date"]]

            assert read_data == True
            assert file_format == 'hdf5'
            assert out_folder == "some_folder"
            assert impute_dates == True
            assert split_berlin == True
            assert moving_average == dd.defaultDict['moving_average']
            assert make_plot == True
            assert no_raw == True
            assert rep_date == False

        test_args = ["prog", '--read-data', '--out-folder',
                     folder, '--file-format', 'json', '--make-plot']

        with patch.object(sys, 'argv', test_args):
            arg_dict = gd.cli("rkiest")
            [read_data, file_format, out_folder, no_raw, make_plot] = [
                arg_dict["read_data"],
                arg_dict["file_format"],
                arg_dict["out_folder"],
                arg_dict["no_raw"],
                arg_dict["make_plot"]]

            assert read_data == True
            assert file_format == 'json'
            assert out_folder == "some_folder"
            assert make_plot == True
            assert no_raw == False

        test_args = [
            "prog", '--out-folder', folder, '--file-format', 'json',
            '--start-date', '2020-11-24', '--end-date', '2020-11-26', '-n']

        with patch.object(sys, 'argv', test_args):
            arg_dict = gd.cli("divi")
            [read_data, file_format, out_folder, end_date, start_date,
             no_raw] = [arg_dict["read_data"],
                        arg_dict["file_format"],
                        arg_dict["out_folder"],
                        arg_dict["end_date"],
                        arg_dict["start_date"],
                        arg_dict["no_raw"]]

            assert read_data == dd.defaultDict['read_data']
            assert file_format == 'json'
            assert out_folder == "some_folder"
            assert end_date == date(2020, 11, 26)
            assert start_date == date(2020, 11, 24)
            assert no_raw == True

        test_args = [
            "prog", '--out-folder', folder, '--file-format', 'json',
            '--make-plot', '--start-date', '2020-11-24', '--end-date',
            '2020-11-26']

        with patch.object(sys, 'argv', test_args):
            arg_dict = gd.cli("sim")
            [read_data, file_format, out_folder, no_raw, end_date,
             impute_dates, make_plot, moving_average, split_berlin,
             start_date] = [arg_dict["read_data"],
                            arg_dict["file_format"],
                            arg_dict["out_folder"],
                            arg_dict["no_raw"],
                            arg_dict["end_date"],
                            arg_dict["impute_dates"],
                            arg_dict["make_plot"],
                            arg_dict["moving_average"],
                            arg_dict["split_berlin"],
                            arg_dict["start_date"]]

            assert read_data == dd.defaultDict['read_data']
            assert file_format == 'json'
            assert out_folder == "some_folder"
            assert end_date == date(2020, 11, 26)
            assert start_date == date(2020, 11, 24)
            assert make_plot == True
            assert split_berlin == dd.defaultDict['split_berlin']
            assert moving_average == dd.defaultDict['moving_average']
            assert impute_dates == dd.defaultDict['impute_dates']
            assert no_raw == False

    def test_append_filename(self):
        test_moving_average = 2
        test_impute_dates = True
        testfilename1 = "FileName_ma2"
        testfilename2 = "FileName_all_dates"
        testfilename3 = "FileName_ma2"

        self.assertEqual(
            testfilename1, gd.append_filename(
                "FileName", dd.defaultDict['impute_dates'],
                test_moving_average))

        self.assertEqual(
            testfilename2, gd.append_filename(
                "FileName", test_impute_dates, dd.defaultDict
                ['moving_average']))

        self.assertEqual(
            testfilename3, gd.append_filename(
                "FileName", test_impute_dates, test_moving_average))

    def test_check_dir(self):

        self.assertFalse(os.path.exists(self.path))

        gd.check_dir(self.path)
        self.assertTrue(os.path.exists(self.path))

    test_string_json = (
        """[{"Date":1606348800000,"col2":"d1"},{"Date":1595721600000,"col2":"d2"}]""")
    test_string_json_timeasstring = (
        """[{"Date":"2020-11-26","col2":"d1"},{"Date":"2020-07-26","col2":"d2"}]""")

    def test_write_dataframe(self):

        gd.check_dir(self.path)

        d1 = date(2020, 11, 26)
        d1 = datetime.combine(d1, datetime.min.time())
        d2 = date(2020, 7, 26)
        d2 = datetime.combine(d2, datetime.min.time())

        d = {'Date': [d1, d2], 'col2': ["d1", "d2"]}
        df = pd.DataFrame(data=d)

        gd.write_dataframe(df, self.path, "test_json", 'json')

        file = "test_json.json"

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), [file])

        file_with_path = os.path.join(self.path, file)
        f = open(file_with_path, "r")
        fread = f.read()
        self.assertEqual(fread, self.test_string_json)

        d1_in_sec = (d1 - datetime(1970, 1, 1)).total_seconds()
        # date in milliseconds
        assert d1_in_sec*1000 == 1606348800000

        d2_in_sec = (d2 - datetime(1970, 1, 1)).total_seconds()
        assert d2_in_sec*1000 == 1595721600000

        gd.write_dataframe(
            df, self.path, "test_json_timeasstring", 'json_timeasstring')

        file2 = "test_json_timeasstring.json"

        self.assertEqual(len(os.listdir(self.path)), 2)
        self.assertEqual(os.listdir(self.path).sort(), [file, file2].sort())

        file_with_path = os.path.join(self.path, file2)
        f = open(file_with_path, "r")
        fread = f.read()
        self.assertEqual(fread, self.test_string_json_timeasstring)

        # TODO: Why does check of hdf5 not work?
        #gd.write_dataframe(df, self.path, "test_hdf", 'hdf5')
        # file = "test_hdf.h5"
        # self.assertEqual(len(os.listdir(self.path)), 3)
        # self.assertEqual(os.listdir(self.path), [file])

    def test_write_dataframe_error(self):

        gd.check_dir(self.path)

        d1 = date(2020, 11, 26)
        d1 = datetime.combine(d1, datetime.min.time())
        d2 = date(2020, 7, 26)
        d2 = datetime.combine(d2, datetime.min.time())

        d = {'Date': [d1, d2], 'col2': ["d1", "d2"]}
        df = pd.DataFrame(data=d)

        with self.assertRaises(SystemExit) as cm:
            gd.write_dataframe(df, self.path, "test_json", 'wrong')

        exit_string = "Error: The file format: " + 'wrong' + " does not exist. Use another one."
        self.assertEqual(cm.exception.code, exit_string)

    @patch('epidemiology.epidata.getDIVIData.get_divi_data')
    @patch('epidemiology.epidata.getRKIData.get_rki_data')
    @patch('epidemiology.epidata.getPopulationData.get_population_data')
    @patch('epidemiology.epidata.getPopulationData.get_age_population_data')
    @patch('epidemiology.epidata.getVaccinationData.get_vaccination_data')
    @patch('epidemiology.epidata.getRKIDatawithEstimations.get_rki_data_with_estimations')
    @patch('epidemiology.epidata.getJHData.get_jh_data')
    def test_call_functions(
            self, mock_jh, mock_rkiwe, mock_vaccination, mock_agep, mock_popul,
            mock_rki, mock_divi):

        arg_dict_all = {
            "read_data": dd.defaultDict['read_data'],
            "file_format": dd.defaultDict['file_format'],
            "out_folder": os.path.join(dd.defaultDict['out_folder']),
            'no_raw': dd.defaultDict["no_raw"]}

        arg_dict_vaccination = {
            **arg_dict_all,
            "end_date": dd.defaultDict['end_date'],
            "make_plot": dd.defaultDict['make_plot'],
            "moving_average": dd.defaultDict['moving_average'],
            "start_date": dd.defaultDict['start_date']}

        arg_dict_rki_est = {**arg_dict_all,
                            "make_plot": dd.defaultDict['make_plot']}

        arg_dict_rki = {
            **arg_dict_all,
            "impute_dates": dd.defaultDict['impute_dates'],
            "make_plot": dd.defaultDict['make_plot'],
            "moving_average": dd.defaultDict['moving_average'],
            "split_berlin": dd.defaultDict['split_berlin'],
            "rep_date": dd.defaultDict['rep_date']}

        arg_dict_divi = {
            **arg_dict_all, "end_date": dd.defaultDict['end_date'],
            "impute_dates": dd.defaultDict['impute_dates'],
            "moving_average": dd.defaultDict['moving_average'],
            "start_date": dd.defaultDict['start_date']}

        getVaccinationData.main()
        mock_vaccination.assert_called()
        mock_vaccination.assert_called_with(**arg_dict_vaccination)

        getPopulationData.main()
        mock_agep.assert_called()
        mock_agep.assert_called_with(**arg_dict_all)
        mock_popul.assert_called()
        mock_popul.assert_called_with(**arg_dict_all)

        getRKIData.main()
        mock_rki.assert_called()
        mock_rki.assert_called_with(**arg_dict_rki)

        getDIVIData.main()
        mock_divi.assert_called()
        mock_divi.assert_called_with(**arg_dict_divi)

        getRKIDatawithEstimations.main()
        mock_rkiwe.assert_called()
        mock_rkiwe.assert_called_with(**arg_dict_rki_est)

        getJHData.main()
        mock_jh.assert_called()
        mock_jh.assert_called_with(**arg_dict_all)


if __name__ == '__main__':
    unittest.main()
