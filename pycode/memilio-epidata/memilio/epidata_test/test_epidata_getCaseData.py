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
import json
import os
import unittest
from datetime import date, datetime
from unittest.mock import patch

import pandas as pd
from pyfakefs import fake_filesystem_unittest

from memilio.epidata import defaultDict as dd
from memilio.epidata import getCaseData as gcd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import progress_indicator


class TestGetCaseData(fake_filesystem_unittest.TestCase):
    path = '/home/Case_Data'

    # strings for read, download and update data
    # be careful: not completely realistic data
    # In the json file there is one data set for each county (ensuring that the
    # completeness check evaluates to true) that is combined with some addional
    # test data sets defined below

    # Get a file object with write permission.
    here = os.path.dirname(os.path.abspath(__file__))

    # load test data for read
    filename = os.path.join(
        here, 'test_data', 'TestSetCaseRead.json')
    # Load JSON file data to a python dict object.
    with open(filename) as file_object:
        dict_object = json.load(file_object)
    test_string_all_federal_states_and_counties_read = json.dumps(dict_object)[:-1] +\
        (""",{"Altersgruppe":"A60-A79","Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-08-11",\
    "IdLandkreis":1002,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-08-07","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":1,"IdBundesland":1},\
    {"Altersgruppe":"A60-A79","Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":1,"Meldedatum":"2020-03-24",\
    "IdLandkreis":1001,"NeuerFall":0,"NeuerTodesfall":0,"Refdatum":"2020-08-07","NeuGenesen":-9,"AnzahlGenesen":0,\
    "IstErkrankungsbeginn":1,"IdBundesland":1},\
    {"Altersgruppe":"A15-A34","Geschlecht":"W","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-08-13",\
    "IdLandkreis":2000,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-08-07","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":1,"IdBundesland":2},\
    {"Altersgruppe":"A60-A79","Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-08-07",\
    "IdLandkreis":3462,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2021-01-09","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":0,"IdBundesland":3},\
    {"Altersgruppe":"A15-A34","Geschlecht":"W","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2021-04-04",\
    "IdLandkreis":4011,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-08-07","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":1,"IdBundesland":4},\
    {"Altersgruppe":"A15-A34","Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-08-07",\
    "IdLandkreis":5111,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-06-22","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":0,"IdBundesland":5},\
    {"Altersgruppe":"A60-A79","Geschlecht":"W","AnzahlFall":1,"AnzahlTodesfall":1,"Meldedatum":"2020-04-21",\
    "IdLandkreis":6531,"NeuerFall":0,"NeuerTodesfall":0,"Refdatum":"2020-04-13","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":1,"IdBundesland":6},\
    {"Altersgruppe":"A35-A59","Geschlecht":"W","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-08-08",\
    "IdLandkreis":7231,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-08-08","NeuGenesen":-9,"AnzahlGenesen":0,\
    "IstErkrankungsbeginn":1,"IdBundesland":7},\
    {"Altersgruppe":"A15-A34","Geschlecht":"W","AnzahlFall":2,"AnzahlTodesfall":0,"Meldedatum":"2020-03-25",\
    "IdLandkreis":8118,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-03-25","NeuGenesen":0,"AnzahlGenesen":2,\
    "IstErkrankungsbeginn":0,"IdBundesland":8},\
    {"Altersgruppe":"A35-A59","Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-08-11",\
    "IdLandkreis":9186,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-08-11","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":1,"IdBundesland":9},\
    {"Altersgruppe":"A60-A79","Geschlecht":"W","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-06-10",\
    "IdLandkreis":10042,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-06-10","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":0,"IdBundesland":10},\
    {"Altersgruppe":"A60-A79","Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-06-04",\
    "IdLandkreis":11004,"NeuerFall":0,"NeuerTodesfall" :-9,"Refdatum":"2020-06-04","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":1,"IdBundesland":11},\
    {"Altersgruppe":"A60-A79","Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-06-04",\
    "IdLandkreis":11011,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-06-04","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":0,"IdBundesland":11},\
    {"Altersgruppe":"A35-A59","Geschlecht":"M","AnzahlFall":5,"AnzahlTodesfall":3,"Meldedatum":"2020-08-10",\
    "IdLandkreis":12066,"NeuerFall":0,"NeuerTodesfall":3,"Refdatum":"2020-08-10","NeuGenesen":0,"AnzahlGenesen":2,\
    "IstErkrankungsbeginn":1,"IdBundesland":12},\
    {"Altersgruppe": "A35-A59","Geschlecht":"M","AnzahlFall":2,"AnzahlTodesfall":0,"Meldedatum":"2021-01-12",\
    "IdLandkreis":13074,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-08-09","NeuGenesen":0,"AnzahlGenesen":2,\
    "IstErkrankungsbeginn":1,"IdBundesland":13},\
    {"Altersgruppe":"A05-A14","Geschlecht":"M","AnzahlFall":4,"AnzahlTodesfall":0,"Meldedatum":"2020-08-09",\
    "IdLandkreis":14511,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-12-10","NeuGenesen":0,"AnzahlGenesen":4,\
    "IstErkrankungsbeginn":0,"IdBundesland":14},\
    {"Altersgruppe":"A05-A14","Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-08-09",\
    "IdLandkreis":15003,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-08-09","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":1,"IdBundesland":15},\
    {"Altersgruppe":"A35-A59","Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-08-09",\
    "IdLandkreis":16061,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-08-09","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":1,"IdBundesland":16},\
    {"Altersgruppe":"A35-A59","Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":1,"Meldedatum":"2021-01-21",\
    "IdLandkreis":16061,"NeuerFall":0,"NeuerTodesfall":0,"Refdatum":"2021-01-21","NeuGenesen":0,"AnzahlGenesen":1,\
    "IstErkrankungsbeginn":1,"IdBundesland":16}]""")

    # load test data for download formatted as data from github
    # (https://github.com/robert-koch-institut/SARS-CoV-2-Infektionen_in_Deutschland)
    filename = os.path.join(
        here, 'test_data', 'TestSetCaseGithub.json')
    # Load JSON file data to a python dict object.
    with open(filename) as file_object:
        dict_object_github = json.load(file_object)

    test_string_all_federal_states_and_counties_github = json.dumps(
        dict_object_github)[:-1] + ("""]""")

    # load test data for download formatted as data from arcgis
    # (https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/66876b81065340a4a48710b062319336/about)
    filename = os.path.join(
        here, 'test_data', 'TestSetCaseArcgis.json')
    # Load JSON file data to a python dict object.
    with open(filename) as file_object:
        dict_object_arcgis = json.load(file_object)

    test_string_all_federal_states_and_counties_arcgis = json.dumps(
        dict_object_arcgis)[:-1] + ("""]""")

    string_not_all_states = (
        """[{"Altersgruppe":"A60-A79","Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"Meldedatum":"2020-08-11",\
        "IdLandkreis":1002,"NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020-08-07","NeuGenesen":0,"AnzahlGenesen":1,\
        "IstErkrankungsbeginn":1, "IdBundesland":1}]""")

    def setUp(self):
        self.setUpPyfakefs()
        progress_indicator.ProgressIndicator.disable_indicators(True)

    def write_case_data(self, out_folder):
        # write dataset for reading data
        case_file = "CaseDataFull.json"
        case_file_with_path = os.path.join(out_folder, case_file)
        with open(case_file_with_path, 'w') as f:
            f.write(self.test_string_all_federal_states_and_counties_read)

    def write_case_data_arcgis(self, out_folder):
        # write dataset from source for mocking download from arcgis
        case_file = "CaseDataArcgis.json"
        case_file_with_path = os.path.join(out_folder, case_file)
        with open(case_file_with_path, 'w') as f:
            f.write(self.test_string_all_federal_states_and_counties_arcgis)

    def write_case_data_not_all_states(self, out_folder):
        case_file = "CaseDataNotFull.json"
        case_file_with_path = os.path.join(out_folder, case_file)
        with open(case_file_with_path, 'w') as f:
            f.write(self.string_not_all_states)

    @patch('memilio.epidata.getDataIntoPandasDataFrame.get_file')
    def test_get_case_data(self, mock_file):
        # Test with downloading data
        read_data = True
        file_format = 'json_timeasstring'
        out_folder = self.path
        no_raw = False
        start_date = date(2020, 1, 1)
        end_date = date.today()
        impute_dates = False
        moving_average = 0
        make_plot = False
        split_berlin = False
        rep_date = False
        files = 'Plot'

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        self.write_case_data_arcgis(directory)
        self.write_case_data_not_all_states(directory)
        # check if expected files are written
        self.assertEqual(len(os.listdir(directory)), 2)

        # test case where all files are incomplete
        mock_file.return_value = pd.read_json(
            os.path.join(directory, "CaseDataNotFull.json"))
        with self.assertRaises(FileNotFoundError) as error:
            gcd.get_case_data(
                read_data=read_data, file_format=file_format,
                out_folder=out_folder, no_raw=no_raw, start_date=start_date,
                end_date=end_date, impute_dates=impute_dates,
                moving_average=moving_average, make_plot=make_plot,
                split_berlin=split_berlin, rep_date=rep_date, files=files)
        self.assertEqual(
            str(error.exception),
            "Something went wrong, dataframe is empty for csv and geojson!")

        self.assertEqual(mock_file.call_count, 3)

        # test case where csv files are incorrect
        mock_file.side_effect = [pd.DataFrame(),
                                 pd.DataFrame(),
                                 pd.read_json(
            os.path.join(
                directory,
                "CaseDataArcgis.json"))]

        gcd.get_case_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, start_date=start_date,
            end_date=end_date, impute_dates=impute_dates,
            moving_average=moving_average, make_plot=make_plot,
            split_berlin=split_berlin, rep_date=rep_date, files=files)

        # check if expected files are written
        # (5 written files + 2 source files: "CaseDataNotFull.json" and "CaseDataArcgis.json")
        self.assertEqual(len(os.listdir(directory)), 5+2)

    @patch('memilio.epidata.getDataIntoPandasDataFrame.get_file')
    def test_get_case_data_split_berlin(self, mock_file):
        # Test case with downloading data where first csv-source is incomplete and second one is used
        # and split_berlin = True
        read_data = False
        file_format = 'json_timeasstring'
        out_folder = self.path
        no_raw = False
        start_date = date(2020, 1, 1)
        end_date = date.today()
        impute_dates = False
        moving_average = 0
        make_plot = False
        split_berlin = True
        rep_date = False
        files = ['all_county', 'all_county_gender',
                 'all_state', 'all_germany', 'infected', 'deaths']

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # write file
        self.write_case_data_arcgis(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        # test case where first csv file is empty and second one is complete
        mock_file.side_effect = [pd.DataFrame(), pd.read_json(
            os.path.join(directory, "CaseDataArcgis.json"))]

        gcd.get_case_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, start_date=start_date,
            end_date=end_date, impute_dates=impute_dates,
            moving_average=moving_average, make_plot=make_plot,
            split_berlin=split_berlin, rep_date=rep_date, files=files)

        self.assertEqual(mock_file.call_count, 2)

        # check if expected files are written (7 files written + 1 source file "CaseDataArcgis.json")
        self.assertEqual(len(os.listdir(directory)), 7+1)

        # test output files (if all_germany is the same as without splitting Berlin)
        file = "cases_all_germany.json"
        f_read = os.path.join(directory, file)
        df = pd.read_json(f_read)

        file = 'cases_infected.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)

        file = 'cases_deaths.json'
        f_read = os.path.join(directory, file)
        df_deaths = pd.read_json(f_read)

        data_list = df.columns.values.tolist()
        self.assertEqual(
            data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        self.assertEqual(df[(df['Date'] == "2020-06-15")]['Confirmed'].item(),
                         df_infected[(df_infected['Date'] == "2020-06-15")]['Confirmed'].item())
        self.assertEqual(df[(df['Date'] == "2020-06-15")]
                         ['Confirmed'].item(), 3)
        self.assertEqual(
            df[(df['Date'] == "2021-01-06")]['Deaths'].item(),
            df_deaths[(df_deaths['Date'] == "2021-01-06")]['Deaths'].item())
        self.assertEqual(df[(df['Date'] == "2021-01-06")]['Deaths'].item(), 1)
        self.assertEqual(df[(df['Date'] == "2020-06-15")]
                         ["Recovered"].item(), 3)
        self.assertEqual(df[(df['Date'] == "2020-04-06")]
                         ['Confirmed'].item(), 2)
        self.assertEqual(df[(df['Date'] == "2021-03-31")]['Deaths'].item(), 2)
        self.assertEqual(df[(df['Date'] == "2020-04-06")]
                         ["Recovered"].item(), 2)

        # test files that should be different as in other cases
        file = 'cases_all_county_split_berlin.json'
        f_read = os.path.join(directory, file)
        df_county = pd.read_json(f_read)

        self.assertEqual(df_county.shape[0], 411)

        self.assertEqual(
            df_county[df_county['ID_County'] == 11004]
            ['Recovered'].shape[0],
            1)
        self.assertEqual(
            df_county[df_county['ID_County'] == 11011]
            ['Recovered'].shape[0],
            1)

        file = 'cases_all_county_gender_split_berlin.json'
        f_read = os.path.join(directory, file)
        df_gender = pd.read_json(f_read)

        self.assertEqual(
            df_gender[(df_gender['Date'] == "2020-12-28") &
                      (df_gender['Gender'] == "male")].shape[0], 1)

        # check if in state file the counties of Berlin are not splitted
        file = 'cases_all_state.json'
        f_read = os.path.join(directory, file)
        df_state = pd.read_json(f_read)
        self.assertEqual(df_state.shape[0], 286)

    @patch('memilio.epidata.getDataIntoPandasDataFrame.get_file',
           return_value=pd.read_json(
               test_string_all_federal_states_and_counties_read).copy())
    def test_get_case_data_moving_average(self, mock_file):

        read_data = True
        file_format = 'json_timeasstring'
        out_folder = self.path
        no_raw = False
        start_date = date(2020, 1, 1)
        end_date = date.today()
        impute_dates = False
        moving_average = 7
        make_plot = False
        split_berlin = False
        rep_date = False
        files = ['all_germany', 'infected', 'deaths', 'all_state']

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # write file
        self.write_case_data(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        gcd.get_case_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, start_date=start_date,
            end_date=end_date, impute_dates=impute_dates,
            moving_average=moving_average, make_plot=make_plot,
            split_berlin=split_berlin, rep_date=rep_date, files=files)

        # check if expected files are written
        self.assertEqual(len(os.listdir(directory)), 5)

        # test _ma files
        file = 'cases_all_germany_ma7.json'
        f_read = os.path.join(directory, file)
        df_ma = pd.read_json(f_read)

        data_list = df_ma.columns.values.tolist()
        self.assertEqual(
            data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        # test if 7 day average moving is calculated correctly
        self.assertAlmostEqual(
            df_ma[(df_ma['Date'] == "2020-08-07")]
            ['Confirmed'].item(),
            15 + 6 / 7)
        self.assertEqual(
            df_ma[(df_ma['Date'] == "2020-08-07")]['Deaths'].item(), 2)
        self.assertAlmostEqual(
            df_ma[(df_ma['Date'] == "2020-08-07")]
            ["Recovered"].item(),
            14 + 3 / 7)

        self.assertAlmostEqual(
            df_ma[(df_ma['Date'] == "2020-08-08")]
            ['Confirmed'].item(),
            18 + 6 / 7)
        self.assertAlmostEqual(
            df_ma[(df_ma['Date'] == "2020-08-08")]
            ['Deaths'].item(),
            2 + 4 / 7)
        self.assertAlmostEqual(
            df_ma[(df_ma['Date'] == "2020-08-08")]
            ["Recovered"].item(),
            16 + 5 / 7)

        self.assertEqual(
            df_ma[(df_ma['Date'] == "2020-08-11")]['Confirmed'].item(),
            27)
        self.assertAlmostEqual(
            df_ma[(df_ma['Date'] == "2020-08-11")]
            ['Deaths'].item(),
            4 + 1 / 7)
        self.assertAlmostEqual(
            df_ma[(df_ma['Date'] == "2020-08-11")]
            ["Recovered"].item(),
            22 + 6 / 7)

        self.assertEqual(
            df_ma[(df_ma['Date'] == "2020-08-20")]['Confirmed'].item(),
            30)
        self.assertEqual(
            df_ma[(df_ma['Date'] == "2020-08-20")]['Deaths'].item(), 5)
        self.assertEqual(
            df_ma[(df_ma['Date'] == "2020-08-20")]["Recovered"].item(),
            25)

        file = 'cases_infected_ma7.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)

        file = 'cases_deaths_ma7.json'
        f_read = os.path.join(directory, file)
        df_deaths = pd.read_json(f_read)

        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-11")]['Confirmed'].item(
        ), df_infected[(df_infected['Date'] == "2020-08-11")]['Confirmed'].item())
        self.assertEqual(
            df_ma[(df_ma['Date'] == "2020-08-11")]['Deaths'].item(),
            df_deaths[(df_deaths['Date'] == "2020-08-11")]['Deaths'].item())

        self.assertAlmostEqual(
            df_deaths[df_deaths['Date'] == "2020-04-13"]
            ['Deaths'].item(),
            4/7)
        self.assertAlmostEqual(
            df_deaths[df_deaths['Date'] == "2020-04-13"]['Deaths'].item(),
            df_ma[(df_ma['Date'] == "2020-04-13")]['Deaths'].item())
        self.assertAlmostEqual(
            df_deaths[df_deaths['Date'] == "2020-04-14"]['Deaths'].item(),
            df_ma[(df_ma['Date'] == "2020-04-14")]['Deaths'].item())
        self.assertAlmostEqual(
            df_ma[(df_ma['Date'] == "2020-04-14")]
            ['Deaths'].item(),
            5 / 7)

        file = 'cases_all_state_ma7.json'
        f_read = os.path.join(directory, file)
        df_state = pd.read_json(f_read)
        self.assertAlmostEqual(
            df_state
            [(df_state['Date'] == "2020-08-07") & (df_state['ID_State'] == 1)]
            ['Confirmed'].item(),
            1 + 1 / 7)
        self.assertAlmostEqual(
            df_state
            [(df_state['Date'] == "2020-08-08") & (df_state['ID_State'] == 1)]
            ['Confirmed'].item(),
            1 + 3 / 7)
        self.assertEqual(
            df_state
            [(df_state['Date'] == "2020-08-11") & (df_state['ID_State'] == 1)]
            ['Confirmed'].item(),
            2.0)
        self.assertEqual(
            df_state
            [(df_state['Date'] == "2020-08-20") & (df_state['ID_State'] == 1)]
            ['Confirmed'].item(),
            2.0)
        self.assertEqual(
            df_state
            [(df_state['Date'] == "2020-08-11") & (df_state['ID_State'] == 1)]
            ['Deaths'].item(),
            1.0)

    @patch('memilio.epidata.getDataIntoPandasDataFrame.get_file',
           return_value=pd.read_json(
               test_string_all_federal_states_and_counties_read).copy())
    def test_get_case_data_impute_dates(self, mock_file):
        read_data = True
        file_format = 'json_timeasstring'
        out_folder = self.path
        no_raw = False
        start_date = date(2020, 1, 1)
        end_date = date.today()
        impute_dates = True
        moving_average = 0
        make_plot = False
        split_berlin = False
        rep_date = False
        files = ['all_germany']

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        gcd.get_case_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, start_date=start_date,
            end_date=end_date, impute_dates=impute_dates,
            moving_average=moving_average, make_plot=make_plot,
            split_berlin=split_berlin, rep_date=rep_date, files=files)

        # check if expected files are written
        self.assertEqual(len(os.listdir(directory)), 2)

        self.assertTrue(
            'cases_all_germany_all_dates.json' in os.listdir(directory))

        # test _all_dates file
        file = 'cases_all_germany_all_dates.json'
        self.assertTrue(file in os.listdir(directory))
        f_read = os.path.join(directory, file)
        df_ad = pd.read_json(f_read)

        data_list = df_ad.columns.values.tolist()
        self.assertEqual(
            data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        # test if 7 day average moving is calculated correctly
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-06-10")]['Confirmed'].item(), 8)
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-06-10")]['Deaths'].item(), 1)
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-06-10")]["Recovered"].item(), 8)
        # Check an average date in between
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-07-08")]['Confirmed'].item(), 9)
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-07-20")]['Deaths'].item(), 1)
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-07-31")]["Recovered"].item(), 9)

        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-08-07")]['Confirmed'].item(),
            15)
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-08-07")]['Deaths'].item(), 2)
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-08-07")]["Recovered"].item(),
            14)

        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-08-08")]['Confirmed'].item(),
            16)
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-08-08")]['Deaths'].item(), 2)
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-08-08")]["Recovered"].item(),
            14)

        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-08-20")]['Confirmed'].item(),
            30)
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-08-20")]['Deaths'].item(), 5)
        self.assertEqual(
            df_ad[(df_ad['Date'] == "2020-08-20")]["Recovered"].item(),
            25)

    @patch('memilio.epidata.getDataIntoPandasDataFrame.get_file',
           return_value=pd.read_json(
               test_string_all_federal_states_and_counties_read).copy())
    def test_get_case_data_moving_average_and_split_berlin(self, mock_file):
        # test if split_berlin and moving_average = True are working together

        read_data = True
        file_format = 'json_timeasstring'
        out_folder = self.path
        no_raw = False
        start_date = date(2020, 1, 1)
        end_date = date.today()
        impute_dates = False
        moving_average = 7
        make_plot = False
        split_berlin = True
        rep_date = False
        files = ['all_county']

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # write file
        self.write_case_data(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        gcd.get_case_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, start_date=start_date,
            end_date=end_date, impute_dates=impute_dates,
            moving_average=moving_average, make_plot=make_plot,
            split_berlin=split_berlin, rep_date=rep_date, files=files)

        self.assertEqual(len(os.listdir(directory)), 2)

        file = 'cases_all_county_split_berlin_ma7.json'
        f_read = os.path.join(directory, file)
        df_county = pd.read_json(f_read)
        self.assertAlmostEqual(df_county[(df_county['ID_County'] == 11004) & (
            df_county['Date'] == '2020-06-04')]['Confirmed'].item(), 4/7)
        self.assertAlmostEqual(df_county[(df_county['ID_County'] == 11011) & (
            df_county['Date'] == '2020-06-04')]['Confirmed'].item(), 4/7)
        self.assertAlmostEqual(df_county[(df_county['ID_County'] == 11004) & (
            df_county['Date'] == '2020-06-09')]['Recovered'].item(), 1)
        self.assertEqual(df_county[(df_county['ID_County'] == 11011) & (
            df_county['Date'] == '2020-06-09')]['Recovered'].item(), 1)
        self.assertAlmostEqual(df_county[(df_county['ID_County'] == 11011) & (
            df_county['Date'] == '2020-06-09')]['Deaths'].item(), 0)

    def test_get_case_data_all_dates_and_split_berlin(self):
        # test if split_berlin and moving_average = True are working together
        read_data = True
        file_format = 'json_timeasstring'
        out_folder = self.path
        no_raw = False
        start_date = date(2020, 1, 1)
        end_date = date.today()
        impute_dates = True
        moving_average = 0
        make_plot = False
        split_berlin = True
        rep_date = False
        files = ['infected_county', 'all_county',
                 'all_county_age', 'all_county_gender']

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # write file
        self.write_case_data(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        gcd.get_case_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, start_date=start_date,
            end_date=end_date, impute_dates=impute_dates,
            moving_average=moving_average, make_plot=make_plot,
            split_berlin=split_berlin, rep_date=rep_date, files=files)

        self.assertEqual(len(os.listdir(directory)), 5)
        # many files are tested before, don't test them again
        files = [
            "cases_all_county_split_berlin_all_dates.json",
            "cases_infected_county_split_berlin_all_dates.json",
            "cases_all_county_gender_split_berlin_all_dates.json",
            "cases_all_county_age_split_berlin_all_dates.json"]
        for file in files:
            self.assertTrue(file in os.listdir(directory))

    @patch('memilio.epidata.getDataIntoPandasDataFrame.get_file')
    def test_no_raw(self, mock_file):
        read_data = False
        file_format = 'json_timeasstring'
        out_folder = self.path
        no_raw = True
        start_date = date(2020, 1, 1)
        end_date = date.today()
        impute_dates = False
        moving_average = 0
        make_plot = False
        split_berlin = False
        rep_date = False
        files = 'All'

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        mock_file.return_value = pd.read_json(
            self.test_string_all_federal_states_and_counties_github)

        gcd.get_case_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, start_date=start_date,
            end_date=end_date, impute_dates=impute_dates,
            moving_average=moving_average, make_plot=make_plot,
            split_berlin=split_berlin, rep_date=rep_date, files=files)

        mock_file.assert_called()

        # check if expected files are written
        # 13 is one less because CaseDataFull is not written
        self.assertEqual(len(os.listdir(directory)), 13)

        self.assertTrue("CaseDataFull.json" not in os.listdir(directory))

        # test output files
        file = "cases_all_germany.json"
        f_read = os.path.join(directory, file)
        df = pd.read_json(f_read)

        file = "cases_infected.json"
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)

        file = "cases_deaths.json"
        f_read = os.path.join(directory, file)
        df_deaths = pd.read_json(f_read)

        data_list = df.columns.values.tolist()
        self.assertEqual(
            data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        self.assertEqual(df[(df['Date'] == "2020-06-15")]['Confirmed'].item(),
                         df_infected[(df_infected['Date'] == "2020-06-15")]['Confirmed'].item())
        self.assertEqual(df[(df['Date'] == "2020-06-15")]
                         ['Confirmed'].item(), 3)
        self.assertEqual(
            df[(df['Date'] == "2021-01-06")]['Deaths'].item(),
            df_deaths[(df_deaths['Date'] == "2021-01-06")]['Deaths'].item())
        self.assertEqual(df[(df['Date'] == "2021-01-06")]['Deaths'].item(), 1)
        self.assertEqual(df[(df['Date'] == "2020-06-15")]
                         ["Recovered"].item(), 3)
        self.assertEqual(df[(df['Date'] == "2020-04-06")]
                         ['Confirmed'].item(), 2)
        self.assertEqual(df[(df['Date'] == "2021-03-31")]['Deaths'].item(), 2)
        self.assertEqual(df[(df['Date'] == "2020-04-06")]
                         ["Recovered"].item(), 2)

    def test_check_for_completeness(self):
        empty_df = pd.DataFrame()
        self.assertEqual(gcd.check_for_completeness(empty_df), False)

    @patch('memilio.epidata.getDataIntoPandasDataFrame.get_file')
    def test_rep_date(self, mock_file):

        mock_file.return_value = pd.read_json(
            self.test_string_all_federal_states_and_counties_github)

        read_data = False
        file_format = 'json_timeasstring'
        out_folder = self.path
        no_raw = False
        start_date = date(2020, 1, 1)
        end_date = date.today()
        impute_dates = False
        moving_average = 7
        make_plot = False
        split_berlin = False
        rep_date = True
        files = 'Plot'

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        gcd.get_case_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, start_date=start_date,
            end_date=end_date, impute_dates=impute_dates,
            moving_average=moving_average, make_plot=make_plot,
            split_berlin=split_berlin, rep_date=rep_date, files=files)

        self.assertEqual(len(os.listdir(directory)), 5)

    def test_get_case_data_timeframe(self):

        read_data = True
        file_format = 'json_timeasstring'
        out_folder = self.path
        no_raw = False
        start_date = date(2020, 12, 24)
        end_date = date(2021, 5, 17)
        impute_dates = False
        moving_average = 0
        make_plot = False
        split_berlin = False
        rep_date = False
        files = ['all_germany']

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # write file
        self.write_case_data(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        gcd.get_case_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, start_date=start_date,
            end_date=end_date, impute_dates=impute_dates,
            moving_average=moving_average, make_plot=make_plot,
            split_berlin=split_berlin, rep_date=rep_date, files=files)

        file = 'cases_all_germany.json'
        f_read = os.path.join(directory, file)
        df_germany_start_end_date = pd.read_json(f_read)

        # get case data without start date and end date
        gcd.get_case_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, impute_dates=impute_dates,
            moving_average=moving_average, make_plot=make_plot,
            split_berlin=split_berlin, rep_date=rep_date, files=files)

        f_read = os.path.join(directory, file)
        df_germany = pd.read_json(f_read)

        # extract dates which are between start and end date
        upperdate = datetime.strftime(end_date, '%Y-%m-%d')
        lowerdate = datetime.strftime(start_date, '%Y-%m-%d')
        df_germany = df_germany[df_germany[dd.EngEng['date']] <= upperdate]
        df_germany = df_germany[df_germany[dd.EngEng['date']] >= lowerdate]

        # dataframes should be equal
        self.assertEqual(len(df_germany_start_end_date), len(
            df_germany), "Dataframes don't have the same length.")
        self.assertEqual(
            list(df_germany_start_end_date['Confirmed']),
            list(df_germany['Confirmed']),
            "Dataframes don't have the same confirmed cases.")
        self.assertEqual(
            list(df_germany_start_end_date['Recovered']),
            list(df_germany['Recovered']),
            "Dataframes don't have the same recovered cases.")
        self.assertEqual(list(df_germany_start_end_date['Deaths']), list(
            df_germany['Deaths']), "Dataframes don't have the same death cases.")


if __name__ == '__main__':
    unittest.main()
