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
import json
import pandas as pd
import numpy as np

from epidemiology.epidata import getRKIData as grki
from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from unittest.mock import patch


class test_get_RKI_Data(fake_filesystem_unittest.TestCase):
    path = '/home/RKI_Data'

    # strings for read, download and update data
    # be careful: not completely realistic data
    # Get a file object with write permission.
    here = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(here, 'test_epidata_get_RKI_Data_data.json')
    file_object = open(filename, 'r')
    # Load JSON file data to a python dict object.
    dict_object = json.load(file_object)

    # Add aldo test data to new data with every county to reduce changes in tests to an minimum
    # With new data there are three additional cases with 1 case and 1 recovered for
    # 15.04.2020 (male, A0-A4), 15.03.2020 (male, A05-A14), 15.06.2020 (female, A15-A34)
    # All three cases are recovered
    test_string_all_federal_states_and_counties = json.dumps(dict_object)[:-1] +\
    (""",{"IdBundesland":1,"Bundesland":"Schleswig-Holstein","Landkreis":"SK Kiel",
    "Altersgruppe":"A60-A79" ,"Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"ObjectId":1,
    "Meldedatum":"2020\/08\/11 00:00:00+00", "IdLandkreis":1002,"Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0,
    "NeuerTodesfall":-9, "Refdatum":"2020\/08\/07 00:00:00+00","NeuGenesen":0,"AnzahlGenesen":1,
    "IstErkrankungsbeginn":1, "Altersgruppe2":"Nicht \\u00fcbermittelt"},
    {"IdBundesland":1,
    "Bundesland":"Schleswig-Holstein","Landkreis":"SK Flensburg","Altersgruppe":"A60-A79", "Geschlecht":"M",
    "AnzahlFall":1,"AnzahlTodesfall":1,"ObjectId":617,"Meldedatum":"2020\/03\/24 00:00:00+00", "IdLandkreis":1001,
    "Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0,"NeuerTodesfall":0, "Refdatum":"2020\/08\/07 00:00:00+00",
    "NeuGenesen":-9,"AnzahlGenesen":0,"IstErkrankungsbeginn":1, "Altersgruppe2":"Nicht \\u00fcbermittelt"},
    {"IdBundesland":2,"Bundesland":"Hamburg", "Landkreis":"SK Hamburg","Altersgruppe":"A15-A34","Geschlecht":"W",
    "AnzahlFall":1,"AnzahlTodesfall":0, "ObjectId":26489,"Meldedatum":"2020\/08\/13 00:00:00+00","IdLandkreis":2000,
    "Datenstand":"25.01.2021, 00:00 Uhr", "NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020\/08\/07 00:00:00+00",
    "NeuGenesen":0,"AnzahlGenesen":1, "IstErkrankungsbeginn":1,"Altersgruppe2":"Nicht \\u00fcbermittelt"},
    {"IdBundesland":3, "Bundesland":"Niedersachsen","Landkreis":"LK Wittmund","Altersgruppe":"A60-A79",
    "Geschlecht":"M","AnzahlFall":1, "AnzahlTodesfall":0,"ObjectId":121122,"Meldedatum":"2020\/08\/07 00:00:00+00",
    "IdLandkreis":3462, "Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0,"NeuerTodesfall":-9,
    "Refdatum":"2021\/01\/09 00:00:00+00", "NeuGenesen":0,"AnzahlGenesen":1,"IstErkrankungsbeginn":0,
    "Altersgruppe2":"Nicht \\u00fcbermittelt"}, {"IdBundesland":4,"Bundesland":"Bremen","Landkreis":"SK Bremen",
    "Altersgruppe":"A15-A34","Geschlecht":"W", "AnzahlFall":1,"AnzahlTodesfall":0,"ObjectId":123459,
    "Meldedatum":"2021\/04\/04 00:00:00+00","IdLandkreis":4011, "Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0,
    "NeuerTodesfall":-9,"Refdatum":"2020\/08\/07 00:00:00+00", "NeuGenesen":0,"AnzahlGenesen":1,
    "IstErkrankungsbeginn":1,"Altersgruppe2":"Nicht \\u00fcbermittelt"}, {"IdBundesland":5,
    "Bundesland":"Nordrhein-Westfalen","Landkreis":"SK D\\u00fcsseldorf","Altersgruppe":"A15-A34", "Geschlecht":"M",
    "AnzahlFall":1,"AnzahlTodesfall":0,"ObjectId":125996,"Meldedatum":"2020\/08\/07 00:00:00+00", "IdLandkreis":5111,
    "Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0,"NeuerTodesfall":-9, "Refdatum":"2020\/06\/22 00:00:00+00",
    "NeuGenesen":0,"AnzahlGenesen":1,"IstErkrankungsbeginn":0, "Altersgruppe2":"Nicht \\u00fcbermittelt"},
    {"IdBundesland":6,"Bundesland":"Hessen", "Landkreis":"LK Gie\\u00dfen", "Altersgruppe":"A60-A79",
    "Geschlecht":"W", "AnzahlFall":1,"AnzahlTodesfall":1, "ObjectId":408175, "Meldedatum":"2020\/04\/21 00:00:00+00",
    "IdLandkreis":6531, "Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0, "NeuerTodesfall":0,
    "Refdatum":"2020\/04\/13 00:00:00+00", "NeuGenesen":0, "AnzahlGenesen":1, "IstErkrankungsbeginn":1,
    "Altersgruppe2":"Nicht \\u00fcbermittelt"}, {"IdBundesland":7, "Bundesland":"Rheinland-Pfalz","Landkreis":\
    "LK Bernkastel-Wittlich","Altersgruppe":"A35-A59", "Geschlecht":"W", "AnzahlFall":1,"AnzahlTodesfall":0,
    "ObjectId":453169,"Meldedatum":"2020\/08\/08 00:00:00+00", "IdLandkreis":7231, "Datenstand":"25.01.2021,
    00:00 Uhr","NeuerFall":0,"NeuerTodesfall":-9,"Refdatum": "2020\/08\/08 00:00:00+00", "NeuGenesen":-9,
    "AnzahlGenesen":0,"IstErkrankungsbeginn":1, "Altersgruppe2":"Nicht \\u00fcbermittelt"}, {"IdBundesland":8,
    "Bundesland":"Baden-W\\u00fcrttemberg", "Landkreis":"LK Ludwigsburg","Altersgruppe":"A15-A34", "Geschlecht":"W",
    "AnzahlFall":2,"AnzahlTodesfall":0, "ObjectId":512433,"Meldedatum":"2020\/03\/25 00:00:00+00",
    "IdLandkreis":8118,"Datenstand":"25.01.2021, 00:00 Uhr", "NeuerFall":0,"NeuerTodesfall":-9,
    "Refdatum":"2020\/03\/25 00:00:00+00","NeuGenesen":0,"AnzahlGenesen":2, "IstErkrankungsbeginn":0,
    "Altersgruppe2":"Nicht \\u00fcbermittelt"}, {"IdBundesland":9, "Bundesland":"Bayern","Landkreis":\
    "LK Pfaffenhofen a.d.Ilm", "Altersgruppe":"A35-A59", "Geschlecht":"M", "AnzahlFall":1,"AnzahlTodesfall":0,
    "ObjectId":694772, "Meldedatum":"2020\/08\/11 00:00:00+00", "IdLandkreis":9186, "Datenstand":"25.01.2021, 00:00 Uhr"
    ,"NeuerFall":0, "NeuerTodesfall":-9,"Refdatum": "2020\/08\/11 00:00:00+00", "NeuGenesen":0,"AnzahlGenesen":1,
    "IstErkrankungsbeginn":1,"Altersgruppe2": "Nicht \\u00fcbermittelt"}, {"IdBundesland":10,"Bundesland":"Saarland",
    "Landkreis":"LK Merzig-Wadern","Altersgruppe":"A60-A79","Geschlecht":"W", "AnzahlFall":1,"AnzahlTodesfall":0,
    "ObjectId":853265,"Meldedatum":"2020\/06\/10 00:00:00+00","IdLandkreis":10042, "Datenstand":"25.01.2021,
    00:00 Uhr", "NeuerFall":0,"NeuerTodesfall":-9,"Refdatum":"2020\/06\/10 00:00:00+00", "NeuGenesen":0,
    "AnzahlGenesen":1,"IstErkrankungsbeginn":0,"Altersgruppe2":"Nicht \\u00fcbermittelt"}, {"IdBundesland":11,
    "Bundesland":"Berlin","Landkreis":"SK Berlin Charlottenburg-Wilmersdorf", "Altersgruppe": "A60-A79",
    "Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"ObjectId":879504,"Meldedatum": "2020\/06\/04 00:00:00+00",
    "IdLandkreis":11004,"Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0, "NeuerTodesfall" :-9,
    "Refdatum":"2020\/06\/04 00:00:00+00","NeuGenesen":0,"AnzahlGenesen":1, "IstErkrankungsbeginn":1,
    "Altersgruppe2": "Nicht \\u00fcbermittelt"}, {"IdBundesland":11, "Bundesland":"Berlin", "Landkreis":\
    "SK Berlin Lichtenberg","Altersgruppe":"A60-A79", "Geschlecht": "M","AnzahlFall":1, "AnzahlTodesfall":0,\
    "ObjectId":907757, "Meldedatum":"2020\/06\/04 00:00:00+00", "IdLandkreis":11011 , "Datenstand":\
    "25.01.2021, 00:00 Uhr","NeuerFall":0, "NeuerTodesfall":-9, "Refdatum":"2020\/06\/04 00:00:00+00", "NeuGenesen":0,\
    "AnzahlGenesen":1, "IstErkrankungsbeginn":0, "Altersgruppe2":"Nicht \\u00fcbermittelt"}, {"IdBundesland":12,
    "Bundesland":"Brandenburg","Landkreis":"LK Oberspreewald-Lausitz","Altersgruppe":"A35-A59", "Geschlecht":"M",
    "AnzahlFall":5,"AnzahlTodesfall":3,"ObjectId":941060,"Meldedatum":"2020\/08\/10 00:00:00+00",
    "IdLandkreis":12066, "Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0,"NeuerTodesfall":3,"Refdatum":
    "2020\/08\/10 00:00:00+00", "NeuGenesen":0,"AnzahlGenesen":2,"IstErkrankungsbeginn":1,"Altersgruppe2":
    "Nicht \\u00fcbermittelt"}, {"IdBundesland":13,"Bundesland":"Mecklenburg-Vorpommern","Landkreis":
    "LK Nordwestmecklenburg", "Altersgruppe": "A35-A59","Geschlecht":"M","AnzahlFall":2,"AnzahlTodesfall":0,
    "ObjectId":962502,"Meldedatum": "2021\/01\/12 00:00:00+00","IdLandkreis":13074,"Datenstand":"25.01.2021,
    00:00 Uhr","NeuerFall":0,"NeuerTodesfall" :-9, "Refdatum":"2020\/08\/09 00:00:00+00","NeuGenesen":0,
    "AnzahlGenesen":2,"IstErkrankungsbeginn":1, "Altersgruppe2": "Nicht \\u00fcbermittelt"}, {"IdBundesland":14,
    "Bundesland":"Sachsen","Landkreis":"SK Chemnitz", "Altersgruppe":"A05-A14","Geschlecht":"M", "AnzahlFall":4,
    "AnzahlTodesfall":0,"ObjectId":967744, "Meldedatum":"2020\/08\/09 00:00:00+00", "IdLandkreis":14511,
    "Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0, "NeuerTodesfall":-9,"Refdatum": "2020\/12\/10 00:00:00+00",
    "NeuGenesen":0,"AnzahlGenesen":4, "IstErkrankungsbeginn":0,"Altersgruppe2": "Nicht \\u00fcbermittelt"},
    {"IdBundesland":15, "Bundesland":"Sachsen-Anhalt","Landkreis":"SK Magdeburg","Altersgruppe":"A05-A14",
    "Geschlecht": "M", "AnzahlFall":1,"AnzahlTodesfall":0,"ObjectId":1032741,"Meldedatum":"2020\/08\/09 00:00:00+00",
    "IdLandkreis":15003,"Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0,"NeuerTodesfall":-9,
    "Refdatum":"2020\/08\/09 00:00:00+00","NeuGenesen":0,"AnzahlGenesen":1,"IstErkrankungsbeginn":1,
    "Altersgruppe2":"Nicht \\u00fcbermittelt"}, {"IdBundesland":16, "Bundesland":"Th\\u00fcringen","Landkreis":\
    "LK Eichsfeld","Altersgruppe":"A35-A59","Geschlecht" :"M", "AnzahlFall":1,"AnzahlTodesfall":0,"ObjectId":1066020,
    "Meldedatum":"2020\/08\/09 00:00:00+00", "IdLandkreis" :16061,"Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0,
    "NeuerTodesfall":-9, "Refdatum": "2020\/08\/09 00:00:00+00","NeuGenesen":0,"AnzahlGenesen":1,
    "IstErkrankungsbeginn":1, "Altersgruppe2": "Nicht \\u00fcbermittelt"}, {"IdBundesland":16,
    "Bundesland":"Th\\u00fcringen","Landkreis":"LK Eichsfeld","Altersgruppe":"A35-A59","Geschlecht" :"M",
    "AnzahlFall":1,"AnzahlTodesfall":1,"ObjectId":1066020, "Meldedatum":"2021\/01\/21 00:00:00+00", "IdLandkreis"
    :16061,"Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0, "NeuerTodesfall":0, "Refdatum": "2021\/01\/21
    00:00:00+00","NeuGenesen":0,"AnzahlGenesen":1, "IstErkrankungsbeginn":1, "Altersgruppe2":
    "Nicht \\u00fcbermittelt"}]""")

    string_not_all_states = ("""[{"IdBundesland":1,"Bundesland":"Schleswig-Holstein","Landkreis":"SK Kiel", 
    "Altersgruppe":"A60-A79" ,"Geschlecht":"M","AnzahlFall":1,"AnzahlTodesfall":0,"ObjectId":1, 
    "Meldedatum":"2020\/08\/11 00:00:00+00", "IdLandkreis":1002,"Datenstand":"25.01.2021, 00:00 Uhr","NeuerFall":0, 
    "NeuerTodesfall":-9, "Refdatum":"2020\/08\/07 00:00:00+00","NeuGenesen":0,"AnzahlGenesen":1, 
    "IstErkrankungsbeginn":1, "Altersgruppe2":"Nicht \\u00fcbermittelt"}]""")

    def setUp(self):
        self.setUpPyfakefs()
    
    def write_rki_data(self, out_folder):
        file_rki = "FullDataRKI.json"
        file_rki_with_path = os.path.join(out_folder, file_rki)
        with open(file_rki_with_path, 'w') as f:
            f.write(self.test_string_all_federal_states_and_counties)

    def write_rki_data_not_all_states(self, out_folder):
        file_rki = "notFullDataRKI.json"
        file_rki_with_path = os.path.join(out_folder, file_rki)
        with open(file_rki_with_path, 'w') as f:
            f.write(self.string_not_all_states)

    def test_get_rki_data_read(self):
        # Test without downloading data
        [read_data, file_format, out_folder, no_raw, impute_dates, make_plot,
         moving_average, split_berlin, rep_date] = [True, 'json_timeasstring', self.path,
                                          False, False, False, 0, False, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # Test case where file does not exist
        file = "FullDataRKI.json"
        file_with_path = os.path.join(directory, file)

        with self.assertRaises(SystemExit) as cm:
            grki.get_rki_data(read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average,
                              split_berlin,rep_date)

        self.assertEqual(cm.exception.code, "Error: The file: " + file_with_path +
                         " does not exist. Call program without -r flag to get it.")

        # Test case where file exists
        self.write_rki_data(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        grki.get_rki_data(read_data, file_format, out_folder, impute_dates, make_plot, moving_average, no_raw,
                          split_berlin, rep_date)

        # check if expected files are written
        self.assertEqual(len(os.listdir(directory)), 14)

        # test output files
        file = "all_germany_rki.json"
        f_read = os.path.join(directory, file)
        df = pd.read_json(f_read)

        file = 'infected_rki.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)

        file = 'deaths_rki.json'
        f_read = os.path.join(directory, file)
        df_deaths = pd.read_json(f_read)

        data_list = df.columns.values.tolist()
        self.assertEqual(data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(),
                         df_infected[(df_infected['Date'] == "2020-08-07")]['Confirmed'].item())
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(), 15)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(),
                         df_deaths[(df_deaths['Date'] == "2020-08-07")]['Deaths'].item())
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(), 2)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]["Recovered"].item(), 14)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Confirmed'].item(), 8)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Deaths'].item(), 1)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]["Recovered"].item(), 8)

        file = 'all_age_rki.json'
        f_read = os.path.join(directory, file)
        df_age = pd.read_json(f_read)
        self.assertEqual(
            df_age[(df_age['Date'] == "2020-08-07") & (df_age["Age_RKI"] == "A60-A79")]["Recovered"].item(), 6)
        self.assertEqual(
            df_age[(df_age['Date'] == "2020-08-07") & (df_age["Age_RKI"] == "A60-A79")]['Deaths'].item(), 2)
        self.assertEqual(
            df_age[(df_age['Date'] == "2020-08-07") & (df_age["Age_RKI"] == "A60-A79")]['Confirmed'].item(), 7)

        file = 'all_gender_rki.json'
        f_read = os.path.join(directory, file)
        df_gender = pd.read_json(f_read)
        self.assertEqual(
            df_gender[(df_gender['Date'] == "2020-08-07") & (df_gender["Gender"] == "male")]["Recovered"].item(), 7)
        self.assertEqual(
            df_gender[(df_gender['Date'] == "2020-08-07") & (df_gender["Gender"] == "female")]['Deaths'].item(), 1)
        self.assertEqual(
            df_gender[(df_gender['Date'] == "2020-08-07") & (df_gender["Gender"] == "male")]['Deaths'].item(), 1)
        self.assertEqual(
            df_gender[(df_gender['Date'] == "2020-08-07") & (df_gender["Gender"] == "female")]['Confirmed'].item(), 7)

        file = 'all_county_gender_rki.json'
        f_read = os.path.join(directory, file)
        df_gender = pd.read_json(f_read)
        self.assertEqual(df_gender.shape[0], 18+412)
        # checks if Berlins districts are concatenated

        #self.assertEqual(
        #    df_gender[(df_gender['County'] == "SK Berlin") & (df_gender['Gender'] == 'male')]['Confirmed'].shape[0],
        #    10)

        file = 'infected_county_rki.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)
        self.assertEqual(df_infected[(df_infected['County'] == "LK Ludwigsburg")]['Confirmed'].shape[0], 2)

        file = 'all_state_age_rki.json'
        f_read = os.path.join(directory, file)
        df_state = pd.read_json(f_read)
        # for every state one line + state 16 has two dates in string
        self.assertEqual(df_state.shape[0], 362)
        self.assertEqual(df_state[(df_state["ID_State"] == 1) &
                                  (df_state['Date'] == "2020-08-07")]['Confirmed'].item(), 2)

        file = 'infected_state_rki.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)
        self.assertEqual(df_state[(df_state["ID_State"] == 1) & (df_state['Date'] == "2020-08-07")]['Confirmed'].item(),
                         df_infected[(df_infected["ID_State"] == 1) &
                                     (df_infected['Date'] == "2020-08-07")]['Confirmed'].item())

    @patch('epidemiology.epidata.getRKIData.gd.loadGeojson')
    @patch('epidemiology.epidata.getRKIData.gd.loadCsv')
    def test_get_rki_data_dowload(self, mock_loadCsv, mock_loadGeojson):
        # Test with downloading data
        [read_data, file_format, out_folder, impute_dates, make_plot, moving_average, no_raw, split_berlin, rep_date] = \
            [False, 'json_timeasstring', self.path, False, False, False, 0, False, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        self.write_rki_data(directory)
        self.write_rki_data_not_all_states(directory)
        # check if expected files are written
        self.assertEqual(len(os.listdir(directory)), 2)

        # test case where all files are incomplete
        mock_loadCsv.return_value = pd.read_json(os.path.join(directory, "notFullDataRKI.json"))
        mock_loadGeojson.return_value = pd.read_json(os.path.join(directory, "notFullDataRKI.json"))
        with self.assertRaises(SystemExit) as cm:
            grki.get_rki_data(read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average,
                              split_berlin, rep_date)
        self.assertEqual(cm.exception.code, "Something went wrong, dataframe is empty for csv and geojson!")

        mock_loadGeojson.assert_called_once()
        mock_loadCsv.assert_called()

        # test case where csv files are incorrect
        mock_loadCsv.side_effect = [pd.DataFrame(), pd.read_json(os.path.join(directory, "notFullDataRKI.json"))]
        mock_loadGeojson.return_value = pd.read_json(os.path.join(directory, "FullDataRKI.json"))

        grki.get_rki_data(read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average,
                          split_berlin, rep_date)

        mock_loadGeojson.assert_called()
        mock_loadCsv.assert_called()

        # check if expected files are written
        # now 15 because file notFullDataRKI is written
        self.assertEqual(len(os.listdir(directory)), 15)

        # test output files
        file = "all_germany_rki.json"
        f_read = os.path.join(directory, file)
        df = pd.read_json(f_read)

        file = 'infected_rki.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)

        file = 'deaths_rki.json'
        f_read = os.path.join(directory, file)
        df_deaths = pd.read_json(f_read)

        data_list = df.columns.values.tolist()
        self.assertEqual(data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(),
                         df_infected[(df_infected['Date'] == "2020-08-07")]['Confirmed'].item())
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(), 15)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(),
                         df_deaths[(df_deaths['Date'] == "2020-08-07")]['Deaths'].item())
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(), 2)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]["Recovered"].item(), 14)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Confirmed'].item(), 8)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Deaths'].item(), 1)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]["Recovered"].item(), 8)

        # do not test all cases tested above because same input dataframe is used -> if test without downloading pass,
        # these files are correct too

    @patch('epidemiology.epidata.getRKIData.gd.loadGeojson')
    @patch('epidemiology.epidata.getRKIData.gd.loadCsv')
    def test_get_rki_data_dowload_split_berlin(self, mock_loadCsv, mock_loadGeojson):
        # Test case with downloading data where first csv-source is incomplete and second one is used
        # and split_berlin = True
        [read_data, file_format, out_folder, impute_dates, no_raw, make_plot, moving_average, split_berlin, rep_date] = \
            [False, 'json_timeasstring', self.path, False, False, False, 0, True, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # write file
        self.write_rki_data(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        # test case where first csv file is empty and second one is complete
        mock_loadCsv.side_effect = [pd.DataFrame(), pd.read_json(os.path.join(directory, "FullDataRKI.json"))]
        mock_loadGeojson.return_value = pd.DataFrame()

        grki.get_rki_data(read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average,
                          split_berlin, rep_date)

        mock_loadGeojson.assert_not_called()
        mock_loadCsv.assert_called()

        # check if expected files are written
        self.assertEqual(len(os.listdir(directory)), 14)

        # test output files (if all_germany is the same as without splitting Berlin)
        file = "all_germany_rki.json"
        f_read = os.path.join(directory, file)
        df = pd.read_json(f_read)

        file = 'infected_rki.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)

        file = 'deaths_rki.json'
        f_read = os.path.join(directory, file)
        df_deaths = pd.read_json(f_read)

        data_list = df.columns.values.tolist()
        self.assertEqual(data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(),
                         df_infected[(df_infected['Date'] == "2020-08-07")]['Confirmed'].item())
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(), 15)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(),
                         df_deaths[(df_deaths['Date'] == "2020-08-07")]['Deaths'].item())
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(), 2)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]["Recovered"].item(), 14)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Confirmed'].item(), 8)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Deaths'].item(), 1)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]["Recovered"].item(), 8)

        # test files that should be different as in other cases
        file = 'all_county_split_berlin_rki.json'
        f_read = os.path.join(directory, file)
        df_county = pd.read_json(f_read)

        file = 'infected_county_split_berlin_rki.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)

        self.assertEqual(df_county.shape[0], 19+412)

        self.assertEqual(
            df_county[df_county['County'] == "SK Berlin Charlottenburg-Wilmersdorf"]['Recovered'].shape[0], 2)
        self.assertEqual(df_county[df_county['County'] == "SK Berlin Lichtenberg"]['Recovered'].shape[0], 2)

        file = 'all_county_gender_split_berlin_rki.json'
        f_read = os.path.join(directory, file)
        df_gender = pd.read_json(f_read)

        self.assertEqual(df_gender[(df_gender['Date'] == "2020-06-04") & (df_gender['Gender'] == "male")].shape[0], 2)

        # check if in state file the counties of Berlin are not splitted
        file = 'all_state_rki.json'
        f_read = os.path.join(directory, file)
        df_state = pd.read_json(f_read)
        # last state has 2 different dates -> two rows
        self.assertEqual(df_state.shape[0], 304)

    def test_get_rki_data_read_moving_average(self):
        # Test without downloading data
        [read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average, split_berlin, rep_date] = \
            [True, 'json_timeasstring', self.path, False, False, False, 1, False, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # write file
        self.write_rki_data(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        grki.get_rki_data(read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average,
                          split_berlin, rep_date)

        # check if expected files are written
        self.assertEqual(len(os.listdir(directory)), 27)

        # test if normal file os the same
        file = "all_germany_rki.json"
        f_read = os.path.join(directory, file)
        df = pd.read_json(f_read)

        file = 'infected_rki.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)

        file = 'deaths_rki.json'
        f_read = os.path.join(directory, file)
        df_deaths = pd.read_json(f_read)

        data_list = df.columns.values.tolist()
        self.assertEqual(data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(),
                         df_infected[(df_infected['Date'] == "2020-08-07")]['Confirmed'].item())
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(), 15)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(),
                         df_deaths[(df_deaths['Date'] == "2020-08-07")]['Deaths'].item())
        # one deaths on 2020-04-13 + one on 2020-08-07
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(), 2)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]["Recovered"].item(), 14)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Confirmed'].item(), 8)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Deaths'].item(), 1)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]["Recovered"].item(), 8)

        self.assertEqual(df[(df['Date'] == "2020-08-10")]['Deaths'].item(), 5)

        # test _ma files
        file = 'all_germany_ma1_rki.json'
        f_read = os.path.join(directory, file)
        df_ma = pd.read_json(f_read)

        data_list = df_ma.columns.values.tolist()
        self.assertEqual(data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        # test if 7 day average moving is calculated correctly
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-07")]['Confirmed'].item(), 9 + 6)
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-07")]['Deaths'].item(), 1 + 1)
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-07")]["Recovered"].item(), 9 + 5)

        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-08")]['Confirmed'].item(), 9 +  6 + 1)
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-08")]['Deaths'].item(), 2)
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-08")]["Recovered"].item(), 9 + 5)

        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-11")]['Confirmed'].item(),
                               9 + 6   + 1   + 8   + 5   + 1  )
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-11")]['Deaths'].item(), 5)
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-11")]["Recovered"].item(),
                               9 + 5   + 8   + 2   + 1)

        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-20")]['Confirmed'].item(), 30)
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-20")]['Deaths'].item(), 5)
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-20")]["Recovered"].item(), 25)

        file = 'infected_ma1_rki.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)
        
        file = 'deaths_ma1_rki.json'
        f_read = os.path.join(directory, file)
        df_deaths = pd.read_json(f_read)

        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-11")]['Confirmed'].item(),
                               df_infected[(df_infected['Date'] == "2020-08-11")]['Confirmed'].item())
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-08-11")]['Deaths'].item(),
                               df_deaths[(df_deaths['Date'] == "2020-08-11")]['Deaths'].item())
        # Attention: deaths_rki_ma file and all_germany_rki_ma deaths-column are not identical in the first six days
        # after first death. This is the case because in all_germany file, zeros before the first death are included
        # in the calculation of the moving average and in deaths_rki-file first data are just cumulative deaths.
        self.assertEqual(df_deaths[df_deaths['Date'] == "2020-04-13"]['Deaths'].item(), 1.0)
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-04-13")]['Deaths'].item(), 1)
        self.assertNotEqual(df_deaths[df_deaths['Date'] == "2020-04-13"]['Deaths'].items(),
                            df_ma[(df_ma['Date'] == "2020-04-13")]['Deaths'].items())
        self.assertEqual(df_deaths[df_deaths['Date'] == "2020-04-14"]['Deaths'].item(), 1.0)
        self.assertEqual(df_ma[(df_ma['Date'] == "2020-04-14")]['Deaths'].item(), 1)

        file = 'all_state_ma1_rki.json'
        f_read = os.path.join(directory, file)
        df_state = pd.read_json(f_read)
        self.assertEqual(
            df_state[(df_state['Date'] == "2020-08-07") & (df_state['ID_State'] == 1)]['Confirmed'].item(),
            2  )
        self.assertEqual(
            df_state[(df_state['Date'] == "2020-08-08") & (df_state['ID_State'] == 1)]['Confirmed'].item(),
            2)
        self.assertEqual(
            df_state[(df_state['Date'] == "2020-08-11") & (df_state['ID_State'] == 1)]['Confirmed'].item(),
            2)
        self.assertEqual(
            df_state[(df_state['Date'] == "2020-08-20") & (df_state['ID_State'] == 1)]['Confirmed'].item(),
            2.0)
        self.assertEqual(
            df_state[(df_state['Date'] == "2020-08-11") & (df_state['ID_State'] == 1)]['Deaths'].item(),
            1)

    def test_get_rki_data_read_impute_dates(self):
        # Test without downloading data
        [read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average, split_berlin, rep_date] = \
            [True, 'json_timeasstring', self.path, False, True, False, 0,False, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # write file
        self.write_rki_data(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        grki.get_rki_data(read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average,
                          split_berlin, rep_date)

        # check if expected files are written
        self.assertEqual(len(os.listdir(directory)), 27)

        files = ['infected_all_dates_rki.json', 'deaths_all_dates_rki.json', 'all_state_all_dates_rki.json',
                 "infected_state_all_dates_rki.json", "all_state_all_dates_rki.json",
                 "infected_county_all_dates_rki.json", "all_county_all_dates_rki.json",
                 "all_gender_all_dates_rki.json", "all_state_gender_all_dates_rki.json",
                 "all_county_gender_all_dates_rki.json", "all_age_all_dates_rki.json",
                 "all_state_age_all_dates_rki.json", "all_county_age_all_dates_rki.json"]
        for file in files:
            self.assertTrue(file in os.listdir(directory))


        # test if normal file os the same
        file = "all_germany_rki.json"
        f_read = os.path.join(directory, file)
        df = pd.read_json(f_read)

        file = 'infected_rki.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)

        file = 'deaths_rki.json'
        f_read = os.path.join(directory, file)
        df_deaths = pd.read_json(f_read)

        data_list = df.columns.values.tolist()
        self.assertEqual(data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(),
                         df_infected[(df_infected['Date'] == "2020-08-07")]['Confirmed'].item())
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(), 15)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(),
                         df_deaths[(df_deaths['Date'] == "2020-08-07")]['Deaths'].item())
        # one deaths on 2020-04-13 + one on 2020-08-07
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(), 2)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]["Recovered"].item(), 14)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Confirmed'].item(), 8)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Deaths'].item(), 1)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]["Recovered"].item(), 8)

        self.assertEqual(df[(df['Date'] == "2020-08-10")]['Deaths'].item(), 5)

        # test _all_dates files
        file = 'all_germany_all_dates_rki.json'
        self.assertTrue(file in os.listdir(directory))
        f_read = os.path.join(directory, file)
        df_ad = pd.read_json(f_read)

        data_list = df_ad.columns.values.tolist()
        self.assertEqual(data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        # test if 7 day average moving is calculated correctly
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-06-10")]['Confirmed'].item(), 8)
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-06-10")]['Deaths'].item(), 1)
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-06-10")]["Recovered"].item(), 8)
        # Check an average date in between
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-07-08")]['Confirmed'].item(), 9)
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-07-20")]['Deaths'].item(), 1)
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-07-31")]["Recovered"].item(), 9)

        self.assertEqual(df_ad[(df_ad['Date'] == "2020-08-07")]['Confirmed'].item(), 15)
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-08-07")]['Deaths'].item(), 2)
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-08-07")]["Recovered"].item(), 14)

        self.assertEqual(df_ad[(df_ad['Date'] == "2020-08-08")]['Confirmed'].item(), 16)
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-08-08")]['Deaths'].item(), 2)
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-08-08")]["Recovered"].item(), 14)

        self.assertEqual(df_ad[(df_ad['Date'] == "2020-08-20")]['Confirmed'].item(), 30)
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-08-20")]['Deaths'].item(), 5)
        self.assertEqual(df_ad[(df_ad['Date'] == "2020-08-20")]["Recovered"].item(), 25)


    def test_get_rki_data_read_moving_average_and_split_berlin(self):
        # test if split_berlin and moving_average = True are working together
        [read_data, file_format, out_folder, no_raw, impute_dates, make_plot,
         moving_average, split_berlin, rep_date] = [True, 'json_timeasstring', self.path,
                                          False, False, False, 1, True, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # write file
        self.write_rki_data(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        grki.get_rki_data(read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average,
                          split_berlin, rep_date)

        # check if expected files are written (27  same number as with split_berlin=False)
        self.assertEqual(len(os.listdir(directory)), 27)
        # many files are tested before, don't test them again
        file = 'all_county_split_berlin_rki.json'
        f_read = os.path.join(directory, file)
        df_county = pd.read_json(f_read)
        self.assertEqual(df_county[(df_county['County'] == "SK Berlin Charlottenburg-Wilmersdorf") & (
                    df_county['Date'] == '2020-06-04')]['Confirmed'].item(),
                         1)
        self.assertEqual(
            df_county[(df_county['County'] == "SK Berlin Lichtenberg") & (df_county['Date'] == '2020-06-04')][
                'Confirmed'].item(), 1)

        file = 'all_county_split_berlin_ma1_rki.json'
        f_read = os.path.join(directory, file)
        df_county = pd.read_json(f_read)
        self.assertEqual(df_county[(df_county['County'] == "SK Berlin Charlottenburg-Wilmersdorf") & (
                    df_county['Date'] == '2020-06-04')]['Confirmed'].item(),
                         1  )
        self.assertEqual(
            df_county[(df_county['County'] == "SK Berlin Lichtenberg") & (df_county['Date'] == '2020-06-04')][
                'Confirmed'].item(), 1  )
        self.assertEqual(df_county[(df_county['County'] == "SK Berlin Charlottenburg-Wilmersdorf") & (
                    df_county['Date'] == '2020-06-09')]['Recovered'].item(),
                         1)
        self.assertEqual(
            df_county[(df_county['County'] == "SK Berlin Lichtenberg") & (df_county['Date'] == '2020-06-09')][
                'Recovered'].item(), 1)
        self.assertEqual(
            df_county[(df_county['County'] == "SK Berlin Lichtenberg") & (df_county['Date'] == '2020-06-09')][
                'Deaths'].item(), 0)

    def test_get_rki_data_read_all_dates_and_split_berlin(self):
        # test if split_berlin and moving_average = True are working together
        [read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average, split_berlin, rep_date] = \
            [True, 'json_timeasstring', self.path, False, True, False, 0, True, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # write file
        self.write_rki_data(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        grki.get_rki_data(read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average,
                          split_berlin, rep_date)

        # check if expected files are written (27  same number as with split_berlin=False)
        self.assertEqual(len(os.listdir(directory)), 27)
        # many files are tested before, don't test them again
        files = ['all_county_split_berlin_rki.json', 'all_county_split_berlin_all_dates_rki.json',
                 "infected_county_split_berlin_rki.json", "infected_county_split_berlin_all_dates_rki.json",
                 "all_county_gender_split_berlin_rki.json", "all_county_gender_split_berlin_all_dates_rki.json",
                 "all_county_age_split_berlin_rki.json", "all_county_age_split_berlin_all_dates_rki.json"]
        for file in files:
            self.assertTrue(file in os.listdir(directory))

    @patch('epidemiology.epidata.getRKIData.gd.loadCsv')
    def test_no_raw(self, mock_loadCsv):
        # Test with downloading data
        [read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average, split_berlin, rep_date] = \
            [False, 'json_timeasstring', self.path, True, False, False, 0, False, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # check if expected files are written
        mock_loadCsv.return_value = pd.read_json(self.test_string_all_federal_states_and_counties)

        grki.get_rki_data(read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average,
                          split_berlin, rep_date)

        mock_loadCsv.assert_called()

        # check if expected files are written
        # 13 is one less because FullDataRKI is not written
        self.assertEqual(len(os.listdir(directory)), 13)

        self.assertTrue("FullDataRKI.json" not in os.listdir(directory))

        # test output files
        file = "all_germany_rki.json"
        f_read = os.path.join(directory, file)
        df = pd.read_json(f_read)

        file = 'infected_rki.json'
        f_read = os.path.join(directory, file)
        df_infected = pd.read_json(f_read)

        file = 'deaths_rki.json'
        f_read = os.path.join(directory, file)
        df_deaths = pd.read_json(f_read)

        data_list = df.columns.values.tolist()
        self.assertEqual(data_list, ["Date", "Confirmed", "Deaths", "Recovered"])
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(),
                         df_infected[(df_infected['Date'] == "2020-08-07")]['Confirmed'].item())
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Confirmed'].item(), 15)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(),
                         df_deaths[(df_deaths['Date'] == "2020-08-07")]['Deaths'].item())
        self.assertEqual(df[(df['Date'] == "2020-08-07")]['Deaths'].item(), 2)
        self.assertEqual(df[(df['Date'] == "2020-08-07")]["Recovered"].item(), 14)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Confirmed'].item(), 8)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]['Deaths'].item(), 1)
        self.assertEqual(df[(df['Date'] == "2020-06-10")]["Recovered"].item(), 8)

    @patch('epidemiology.epidata.getRKIData.gd.cli')
    def test_main(self, mock_cli):


        mock_cli.return_value = {"read_data": True, "file_format": 'json_timeasstring', "out_folder": self.path,
                                 "impute_dates": False, "make_plot": False, "moving_average": 0,
                                 "split_berlin": False, "no_raw": False, "rep_date": False}

        out_folder = self.path
        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        # Test case where file does not exist
        file = "FullDataRKI.json"
        file_with_path = os.path.join(directory, file)
        # Test case where file exists
        self.write_rki_data(directory)
        # check if expected file is written
        self.assertEqual(len(os.listdir(directory)), 1)

        grki.main()
        # check if expected files are written
        self.assertEqual(len(os.listdir(directory)), 14)


    def test_check_for_completeness(self):
        empty_df = pd.DataFrame()
        self.assertEqual(grki.check_for_completeness(empty_df), False)
    
    @patch('epidemiology.epidata.getRKIData.gd.loadCsv')
    def test_rep_date(self, mocklcsv):

        mocklcsv.return_value = pd.read_json(self.test_string_all_federal_states_and_counties)

        [read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average, split_berlin, rep_date] = \
            [False, 'json_timeasstring', self.path, False, False, False, 7, False, True]
        
        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        grki.get_rki_data(read_data, file_format, out_folder, no_raw, impute_dates, make_plot, moving_average, split_berlin, rep_date)

        mocklcsv.assert_called()
        self.assertEqual(len(os.listdir(directory)), 27)

if __name__ == '__main__':
    unittest.main()
