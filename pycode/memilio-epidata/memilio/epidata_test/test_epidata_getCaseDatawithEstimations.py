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
import numpy as np
import io
from datetime import date, timedelta

import os
import pandas as pd

from memilio.epidata import getCaseDatawithEstimations as gcdwe
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from unittest.mock import patch, call


class TestGetCaseDatawithEstimations(fake_filesystem_unittest.TestCase):
    path = '/home/CaseEstimationData'

    # Notice data is not realistic
    str_whole_country_Germany_jh = \
        ("""[{"CountryRegion":"Germany","Date":"2020-01-22","Confirmed":0,"Recovered":0.0,"Deaths":0.0},\
{"CountryRegion":"Germany","Date":"2020-01-23","Confirmed":0,"Recovered":0.0,"Deaths":0.0},\
{"CountryRegion":"Germany","Date":"2020-01-24","Confirmed":0,"Recovered":0.0,"Deaths":0.0},\
{"CountryRegion":"Germany","Date":"2020-01-25","Confirmed":0,"Recovered":0.0,"Deaths":0.0},\
{"CountryRegion":"Germany","Date":"2020-01-26","Confirmed":0,"Recovered":1.0,"Deaths":1.0},\
{"CountryRegion":"Germany","Date":"2020-01-27","Confirmed":5,"Recovered":2.0,"Deaths":1.0},\
{"CountryRegion":"Germany","Date":"2020-01-28","Confirmed":6,"Recovered":3.0,"Deaths":1.0},\
{"CountryRegion":"Germany","Date":"2020-01-29","Confirmed":7,"Recovered":3.0,"Deaths":1.0},\
{"CountryRegion":"Germany","Date":"2020-01-30","Confirmed":8,"Recovered":3.0,"Deaths":1.0},\
{"CountryRegion":"Germany","Date":"2020-01-31","Confirmed":9,"Recovered":3.0,"Deaths":2.0}]""")

    str_cases_all_germany = (
        """[{"Date":1577836800000,"Confirmed":21,"Deaths":0,"Recovered":21},
{"Date":1577923200000,"Confirmed":23,"Deaths":0,"Recovered":23},
{"Date":1578009600000,"Confirmed":28,"Deaths":0,"Recovered":28},
{"Date":1578096000000,"Confirmed":32,"Deaths":0,"Recovered":32},
{"Date":1578268800000,"Confirmed":33,"Deaths":0,"Recovered":33},
{"Date":1578441600000,"Confirmed":34,"Deaths":0,"Recovered":34},
{"Date":1578528000000,"Confirmed":35,"Deaths":0,"Recovered":35},
{"Date":1578614400000,"Confirmed":37,"Deaths":0,"Recovered":37},
{"Date":1578700800000,"Confirmed":40,"Deaths":0,"Recovered":40},
{"Date":1578873600000,"Confirmed":41,"Deaths":0,"Recovered":41},
{"Date":1578960000000,"Confirmed":42,"Deaths":0,"Recovered":42},
{"Date":1579046400000,"Confirmed":43,"Deaths":0,"Recovered":43},
{"Date":1579132800000,"Confirmed":44,"Deaths":0,"Recovered":44},
{"Date":1579219200000,"Confirmed":47,"Deaths":1,"Recovered":46},
{"Date":1579305600000,"Confirmed":50,"Deaths":1,"Recovered":49},
{"Date":1579392000000,"Confirmed":51,"Deaths":1,"Recovered":50},
{"Date":1579478400000,"Confirmed":52,"Deaths":1,"Recovered":51},
{"Date":1579564800000,"Confirmed":53,"Deaths":1,"Recovered":52},
{"Date":1579651200000,"Confirmed":57,"Deaths":1,"Recovered":56},
{"Date":1579737600000,"Confirmed":59,"Deaths":1,"Recovered":58},
{"Date":1579824000000,"Confirmed":60,"Deaths":1,"Recovered":59},
{"Date":1579910400000,"Confirmed":62,"Deaths":1,"Recovered":61},
{"Date":1579996800000,"Confirmed":63,"Deaths":1,"Recovered":62},
{"Date":1580083200000,"Confirmed":68,"Deaths":1,"Recovered":67},
{"Date":1580169600000,"Confirmed":70,"Deaths":1,"Recovered":69},
{"Date":1580256000000,"Confirmed":73,"Deaths":1,"Recovered":72},
{"Date":1580342400000,"Confirmed":84,"Deaths":1,"Recovered":83},
{"Date":1580428800000,"Confirmed":96,"Deaths":1,"Recovered":95}]""")

    str_cases_all_gender = (
        """ [{"Gender": "female", "Date": 1577836800000, "Confirmed": 7, "Deaths": 0, "Recovered": 7},
{"Gender": "female", "Date": 1577923200000, "Confirmed": 9, "Deaths": 0, "Recovered": 9},
{"Gender": "female", "Date": 1578009600000, "Confirmed": 13, "Deaths": 0, "Recovered": 13},
{"Gender": "female", "Date": 1578096000000, "Confirmed": 16, "Deaths": 0, "Recovered": 16},
{"Gender": "female", "Date": 1578441600000, "Confirmed": 17, "Deaths": 0, "Recovered": 17},
{"Gender": "female", "Date": 1578614400000, "Confirmed": 18, "Deaths": 0, "Recovered": 18},
{"Gender": "female", "Date": 1578700800000, "Confirmed": 20, "Deaths": 0, "Recovered": 20},
{"Gender": "female", "Date": 1578873600000, "Confirmed": 21, "Deaths": 0, "Recovered": 21},
{"Gender": "female", "Date": 1579046400000, "Confirmed": 22, "Deaths": 0, "Recovered": 22},
{"Gender": "female", "Date": 1579132800000, "Confirmed": 23, "Deaths": 0, "Recovered": 23},
{"Gender": "female", "Date": 1579219200000, "Confirmed": 25, "Deaths": 1, "Recovered": 24},
{"Gender": "female", "Date": 1579392000000, "Confirmed": 26, "Deaths": 1, "Recovered": 25},
{"Gender": "female", "Date": 1579737600000, "Confirmed": 27, "Deaths": 1, "Recovered": 26},
{"Gender": "female", "Date": 1580083200000, "Confirmed": 30, "Deaths": 1, "Recovered": 29},
{"Gender": "female", "Date": 1580256000000, "Confirmed": 33, "Deaths": 1, "Recovered": 32},
{"Gender": "female", "Date": 1580342400000, "Confirmed": 38, "Deaths": 1, "Recovered": 37},
{"Gender": "female", "Date": 1580428800000, "Confirmed": 43, "Deaths": 1, "Recovered": 42},
{"Gender":"male","Date":1577836800000,"Confirmed":14,"Deaths":0,"Recovered":14},
{"Gender":"male","Date":1578009600000,"Confirmed":15,"Deaths":0,"Recovered":15},
{"Gender":"male","Date":1578096000000,"Confirmed":16,"Deaths":0,"Recovered":16},
{"Gender":"male","Date":1578268800000,"Confirmed":17,"Deaths":0,"Recovered":17},
{"Gender":"male","Date":1578528000000,"Confirmed":18,"Deaths":0,"Recovered":18},
{"Gender":"male","Date":1578614400000,"Confirmed":19,"Deaths":0,"Recovered":19},
{"Gender":"male","Date":1578700800000,"Confirmed":20,"Deaths":0,"Recovered":20},
{"Gender":"male","Date":1578960000000,"Confirmed":21,"Deaths":0,"Recovered":21},
{"Gender":"male","Date":1579219200000,"Confirmed":22,"Deaths":0,"Recovered":22},
{"Gender":"male","Date":1579305600000,"Confirmed":25,"Deaths":0,"Recovered":25},
{"Gender":"male","Date":1579478400000,"Confirmed":26,"Deaths":0,"Recovered":26},
{"Gender":"male","Date":1579564800000,"Confirmed":27,"Deaths":0,"Recovered":27},
{"Gender":"male","Date":1579651200000,"Confirmed":31,"Deaths":0,"Recovered":31},
{"Gender":"male","Date":1579737600000,"Confirmed":32,"Deaths":0,"Recovered":32},
{"Gender":"male","Date":1579824000000,"Confirmed":33,"Deaths":0,"Recovered":33},
{"Gender":"male","Date":1579910400000,"Confirmed":35,"Deaths":0,"Recovered":35},
{"Gender":"male","Date":1579996800000,"Confirmed":36,"Deaths":0,"Recovered":36},
{"Gender":"male","Date":1580083200000,"Confirmed":38,"Deaths":0,"Recovered":38},
{"Gender":"male","Date":1580169600000,"Confirmed":40,"Deaths":0,"Recovered":40},
{"Gender":"male","Date":1580342400000,"Confirmed":46,"Deaths":0,"Recovered":46},
{"Gender":"male","Date":1580428800000,"Confirmed":53,"Deaths":0,"Recovered":53}]""")

    str_age_germany = ("""[
{"Age_RKI":"A00-A04","Date":1579824000000,"Confirmed":1,"Deaths":0,"Recovered":1},
{"Age_RKI":"A00-A04","Date":1580256000000,"Confirmed":2,"Deaths":0,"Recovered":2},
{"Age_RKI":"A00-A04","Date":1580428800000,"Confirmed":3,"Deaths":0,"Recovered":3},
{"Age_RKI":"A05-A14","Date":1580256000000,"Confirmed":10,"Deaths":0,"Recovered":10},
{"Age_RKI":"A05-A14","Date":1580342400000,"Confirmed":11,"Deaths":0,"Recovered":11},
{"Age_RKI":"A05-A14","Date":1580428800000,"Confirmed":12,"Deaths":0,"Recovered":12},
{"Age_RKI":"A15-A34","Date":1579737600000,"Confirmed":42,"Deaths":0,"Recovered":42},
{"Age_RKI":"A15-A34","Date":1579824000000,"Confirmed":43,"Deaths":0,"Recovered":43},
{"Age_RKI":"A15-A34","Date":1579910400000,"Confirmed":44,"Deaths":0,"Recovered":44},
{"Age_RKI":"A15-A34","Date":1580083200000,"Confirmed":47,"Deaths":0,"Recovered":47},
{"Age_RKI":"A15-A34","Date":1580169600000,"Confirmed":48,"Deaths":0,"Recovered":48},
{"Age_RKI":"A15-A34","Date":1580256000000,"Confirmed":51,"Deaths":0,"Recovered":51},
{"Age_RKI":"A15-A34","Date":1580342400000,"Confirmed":56,"Deaths":0,"Recovered":56},
{"Age_RKI":"A15-A34","Date":1580428800000,"Confirmed":58,"Deaths":0,"Recovered":58},
{"Age_RKI":"A35-A59","Date":1579824000000,"Confirmed":69,"Deaths":0,"Recovered":69},
{"Age_RKI":"A35-A59","Date":1579910400000,"Confirmed":70,"Deaths":0,"Recovered":70},
{"Age_RKI":"A35-A59","Date":1579996800000,"Confirmed":74,"Deaths":0,"Recovered":74},
{"Age_RKI":"A35-A59","Date":1580083200000,"Confirmed":76,"Deaths":0,"Recovered":76},
{"Age_RKI":"A35-A59","Date":1580169600000,"Confirmed":78,"Deaths":0,"Recovered":78},
{"Age_RKI":"A35-A59","Date":1580256000000,"Confirmed":80,"Deaths":0,"Recovered":80},
{"Age_RKI":"A35-A59","Date":1580342400000,"Confirmed":81,"Deaths":0,"Recovered":81},
{"Age_RKI":"A35-A59","Date":1580428800000,"Confirmed":86,"Deaths":0,"Recovered":86},
{"Age_RKI":"A60-A79","Date":1579824000000,"Confirmed":31,"Deaths":3,"Recovered":28},
{"Age_RKI":"A60-A79","Date":1579910400000,"Confirmed":32,"Deaths":3,"Recovered":29},
{"Age_RKI":"A60-A79","Date":1580083200000,"Confirmed":33,"Deaths":3,"Recovered":30},
{"Age_RKI":"A60-A79","Date":1580169600000,"Confirmed":34,"Deaths":3,"Recovered":31},
{"Age_RKI":"A60-A79","Date":1580342400000,"Confirmed":35,"Deaths":3,"Recovered":32},
{"Age_RKI":"A80+","Date":1579737600000,"Confirmed":20,"Deaths":1,"Recovered":19},
{"Age_RKI":"A80+","Date":1579910400000,"Confirmed":21,"Deaths":1,"Recovered":20},
{"Age_RKI":"A80+","Date":1580083200000,"Confirmed":22,"Deaths":1,"Recovered":21},
{"Age_RKI":"A80+","Date":1580256000000,"Confirmed":23,"Deaths":1,"Recovered":22}]""")

    case_files_to_change = [
        "cases_all_germany", "cases_all_gender", "cases_all_age",
        "cases_all_state", "cases_all_state_gender", "cases_all_state_age",
        "cases_all_county", "cases_all_county_gender", "cases_all_county_age"]

    def setUp(self):
        self.setUpPyfakefs()

    def write_case_data(self, out_folder):

        for file_to_change in self.case_files_to_change:

            case_file = file_to_change + ".json"
            case_file_with_path = os.path.join(out_folder, case_file)

            if file_to_change == "cases_all_gender":
                with open(case_file_with_path, 'w') as f:
                    f.write(self.str_cases_all_gender)
            elif file_to_change == "cases_all_age":
                with open(case_file_with_path, 'w') as f:
                    f.write(self.str_age_germany)
            else:
                with open(case_file_with_path, 'w') as f:
                    f.write(self.str_cases_all_germany)

    def write_jh_data(self, out_folder):
        file_jh = "whole_country_Germany_jh.json"
        file_jh_with_path = os.path.join(out_folder, file_jh)

        with open(file_jh_with_path, 'w') as f:
            f.write(self.str_whole_country_Germany_jh)

    def write_weekly_deaths_xlsx_data(
            self, out_folder, file_name='Cases_deaths_weekly.xlsx'):
        sheet1 = pd.DataFrame(
            {'Sterbejahr': ['2020', '2020', '2020', '2020'],
             'Sterbewoche': ['1', '3', '10', '51'],
             'Anzahl verstorbene COVID-19 Fälle': ['0', '<4', '18', '3000']})
        sheet2 = pd.DataFrame(
            {'Sterbejahr': ['2020', '2020', '2020', '2020'],
             'Sterbewoche': ['1', '3', '10', '51'],
             'AG 0-9 Jahre': ['0', '<4', '30', '10'],
             'AG 10-19 Jahre': ['0', '<4', '30', '10'],
             'AG 20-29 Jahre': ['0', '<4', '30', '10'],
             'AG 30-39 Jahre': ['0', '<4', '30', '10'],
             'AG 40-49 Jahre': ['0', '<4', '30', '10'],
             'AG 50-59 Jahre': ['0', '<4', '30', '10'],
             'AG 60-69 Jahre': ['0', '<4', '30', '10'],
             'AG 70-79 Jahre': ['0', '<4', '30', '10'],
             'AG 80-89 Jahre': ['0', '<4', '30', '10'],
             'AG 90+ Jahre': ['0', '<4', '30', '10']})
        sheet3 = pd.DataFrame(
            {'Sterbejahr': ['2020', '2020', '2020', '2020'],
             'Sterbewoche': ['1', '3', '10', '51'],
             'Männer, AG 0-19 Jahre': ['0', '<4', '30', '10'],
             'Männer, AG 20-39 Jahre': ['0', '<4', '30', '10'],
             'Männer, AG 40-59 Jahre': ['0', '<4', '30', '10'],
             'Männer, AG 60-79 Jahre': ['0', '<4', '30', '10'],
             'Männer, AG 80+ Jahre': ['0', '<4', '30', '10'],
             'Frauen, AG 0-19 Jahre': ['0', '<4', '30', '10'],
             'Frauen, AG 20-39 Jahre': ['0', '<4', '30', '10'],
             'Frauen, AG 40-59 Jahre': ['0', '<4', '30', '10'],
             'Frauen, AG 60-79 Jahre': ['0', '<4', '30', '10'],
             'Frauen, AG 80+ Jahre': ['0', '<4', '30', '10']})

        income_sheets = {'COVID_Todesfälle': sheet1,
                         'COVID_Todesfälle_KW_AG10': sheet2,
                         'COVID_Todesfälle_KW_AG20_G': sheet3}
        path = os.path.join(out_folder, file_name)
        dummy = pd.ExcelWriter(path)

        for sheet_name in income_sheets.keys():
            income_sheets[sheet_name].to_excel(
                dummy, sheet_name=sheet_name, index=False)

        dummy.save()

    def test_get_case_data_with_estimations(self):

        [out_folder, read_data, make_plot, file_format, no_raw] \
            = [self.path, True, False, "json", False]

        # write files which should be read in by program

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        self.write_case_data(directory)
        self.write_jh_data(directory)
        self.write_weekly_deaths_xlsx_data(directory)

        # check if expected files are written
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(
            len(os.listdir(directory)),
            2 + len(self.case_files_to_change))

        gcdwe.get_case_data_with_estimations(
            out_folder, read_data, file_format, no_raw, make_plot)

        # check if expected files are written
        self.assertEqual(len(os.listdir(self.path)), 1)
        # 1 jh-file, 2*len(case file): original+estimated, 4 weekly deaths original+original&estimated+ageresolved+genderresolved
        self.assertEqual(
            len(os.listdir(directory)),
            1 + 2 * len(self.case_files_to_change) + 4)

        f_read = os.path.join(directory, "cases_all_germany_estimated.json")
        df = pd.read_json(f_read)

        confirmed = dd.EngEng['confirmed']
        recovered = dd.EngEng['recovered']
        deaths = dd.EngEng['deaths']
        date = dd.EngEng['date']
        recovered_estimated = recovered + "_estimated"
        deaths_estimated = deaths + "_estimated"

        data_list = df.columns.values.tolist()

        self.assertEqual(
            data_list,
            [date, confirmed, deaths, recovered, recovered_estimated,
             deaths_estimated])

        self.assertEqual(
            df[(df[date] == "2020-01-30")][recovered_estimated].item(),
            np.round(84 * 3. / 8.))
        self.assertEqual(
            df[(df[date] == "2020-01-31")][recovered_estimated].item(),
            np.round(96 * 3. / 9.))
        self.assertEqual(
            df[(df[date] == "2020-01-30")][deaths_estimated].item(),
            np.round(84 * 1. / 8.))
        self.assertEqual(
            df[(df[date] == "2020-01-31")][deaths_estimated].item(),
            np.round(96 * 2. / 9.))

        # gender specific data
        gender = dd.EngEng['gender']

        f_read = os.path.join(directory, "cases_all_gender_estimated.json")
        df = pd.read_json(f_read)

        data_list = df.columns.values.tolist()

        self.assertEqual(
            data_list,
            [gender, date, confirmed, deaths, recovered, recovered_estimated,
             deaths_estimated])

        self.assertEqual(
            df[(df[date] == "2020-01-28") & (df[gender] == "male")]
            [recovered_estimated].item(),
            np.round(40 * 3. / 6.))
        self.assertEqual(
            df[(df[date] == "2020-01-29") & (df[gender] == "female")]
            [recovered_estimated].item(),
            np.round(33 * 3. / 7.))
        self.assertEqual(
            df[(df[date] == "2020-01-31") & (df[gender] == "male")]
            [recovered_estimated].item(),
            np.round(53 * 3. / 9.))
        self.assertEqual(
            df[(df[date] == "2020-01-31") & (df[gender] == "female")]
            [recovered_estimated].item(),
            np.round(43 * 3. / 9.))
        self.assertEqual(
            df[(df[date] == "2020-01-28") & (df[gender] == "male")]
            [deaths_estimated].item(),
            np.round(40 * 1. / 6.))
        self.assertEqual(
            df[(df[date] == "2020-01-29") & (df[gender] == "female")]
            [deaths_estimated].item(),
            np.round(33 * 1. / 7.))
        self.assertEqual(
            df[(df[date] == "2020-01-31") & (df[gender] == "male")]
            [deaths_estimated].item(),
            np.round(53 * 2. / 9.))
        self.assertEqual(
            df[(df[date] == "2020-01-31") & (df[gender] == "female")]
            [deaths_estimated].item(),
            np.round(43 * 2. / 9.))

    def test_get_case_data_with_estimations_age_data(self):

        [out_folder, read_data, make_plot, file_format, no_raw] \
            = [self.path, True, False, "json", False]

        # write files which should be read in by program

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        self.write_case_data(directory)
        self.write_jh_data(directory)
        self.write_weekly_deaths_xlsx_data(directory)

        # check if expected files are written
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(
            len(os.listdir(directory)),
            2 + len(self.case_files_to_change))

        gcdwe.get_case_data_with_estimations(
            out_folder, read_data, file_format, no_raw, make_plot)

        # check if expected files are written
        self.assertEqual(len(os.listdir(self.path)), 1)
        # 1 jh-file, 2*len(case file): original+estimated, 4 weekly deaths original+original&estimated+ageresolved+genderresolved
        self.assertEqual(
            len(os.listdir(directory)),
            1 + 2 * len(self.case_files_to_change) + 4)

        f_read = os.path.join(directory, "cases_all_age_estimated.json")
        df = pd.read_json(f_read)

        confirmed = dd.EngEng['confirmed']
        recovered = dd.EngEng['recovered']
        deaths = dd.EngEng['deaths']
        date = dd.EngEng['date']
        recovered_estimated = recovered + "_estimated"
        deaths_estimated = deaths + "_estimated"

        ages = ["A0-A04", "A05-A14", "A15-A34", "A35-A59", "A60-A79", "A80+"]

        data_list = df.columns.values.tolist()

        self.assertEqual(
            data_list.sort(),
            [date, confirmed, deaths, recovered, recovered_estimated,
             deaths_estimated, "Age_RKI"].sort())

        age_values_rec = [3, 12, 58, 86]
        index = 0
        for value in ages:

            mask = (df['Age_RKI'] == value) & (df[date] == "2020-01-31")

            try:
                self.assertEqual(
                    df.loc[mask][recovered_estimated].item(),
                    np.round(age_values_rec[index] * 3. / 9.))
                self.assertEqual(
                    df.loc[mask][deaths_estimated].item(),
                    np.round(age_values_rec[index] * 2. / 9.))
            except ValueError:
                pass

            index = index + 1

    @patch('memilio.epidata.getCaseDatawithEstimations.gcd.get_case_data')
    @patch('memilio.epidata.getCaseDatawithEstimations.gjd.get_jh_data')
    @patch('memilio.epidata.getCaseDatawithEstimations.download_weekly_deaths_numbers')
    def test_get_case_data_with_estimations_download(
            self, mock_get_jh_data, mock_get_case_data,
            mock_download_weekly_deaths_numbers):

        [out_folder, read_data, make_plot, file_format, no_raw] \
            = [self.path, False, False, "json", False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        mock_get_case_data.side_effect = self.write_case_data(directory)
        mock_get_jh_data.side_effect = self.write_jh_data(directory)
        mock_download_weekly_deaths_numbers.side_effect = self.write_weekly_deaths_xlsx_data(
            directory)

        # write files which should be read in by program

        case_files_to_change = [
            "cases_all_germany", "cases_all_gender", "cases_all_age",
            "cases_all_state", "cases_all_state_gender", "cases_all_state_age",
            "cases_all_county", "cases_all_county_gender",
            "cases_all_county_age"]

        for file_to_change in case_files_to_change:

            case_file = file_to_change + ".json"
            case_file_with_path = os.path.join(directory, case_file)

            if file_to_change == "cases_all_gender":
                with open(case_file_with_path, 'w') as f:
                    f.write(self.str_cases_all_gender)
            else:
                with open(case_file_with_path, 'w') as f:
                    f.write(self.str_cases_all_germany)

        file_jh = "whole_country_Germany_jh.json"
        file_jh_with_path = os.path.join(directory, file_jh)

        with open(file_jh_with_path, 'w') as f:
            f.write(self.str_whole_country_Germany_jh)

        # check if expected files are written
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(
            len(os.listdir(directory)),
            2 + len(case_files_to_change))

        gcdwe.get_case_data_with_estimations(
            out_folder, read_data, file_format, no_raw, make_plot)

        # check if expected files are written
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(
            len(os.listdir(directory)),
            1 + 2 * len(case_files_to_change) + 4)

        confirmed = dd.EngEng['confirmed']
        recovered = dd.EngEng['recovered']
        deaths = dd.EngEng['deaths']
        date = dd.EngEng['date']
        recovered_estimated = recovered + "_estimated"
        deaths_estimated = deaths + "_estimated"

        f_read = os.path.join(directory, "cases_all_germany_estimated.json")
        df = pd.read_json(f_read)

        data_list = df.columns.values.tolist()

        self.assertEqual(
            data_list,
            [date, confirmed, deaths, recovered, recovered_estimated,
             deaths_estimated])

        self.assertEqual(
            df[(df[date] == "2020-01-30")][recovered_estimated].item(),
            np.round(84 * 3. / 8.))
        self.assertEqual(
            df[(df[date] == "2020-01-31")][recovered_estimated].item(),
            np.round(96 * 3. / 9.))
        self.assertEqual(
            df[(df[date] == "2020-01-30")][deaths_estimated].item(),
            np.round(84 * 1. / 8.))
        self.assertEqual(
            df[(df[date] == "2020-01-31")][deaths_estimated].item(),
            np.round(96 * 2. / 9.))

        # gender specific data
        gender = dd.EngEng['gender']

        f_read = os.path.join(directory, "cases_all_gender_estimated.json")
        df = pd.read_json(f_read)

        data_list = df.columns.values.tolist()

        self.assertEqual(
            data_list,
            [gender, date, confirmed, deaths, recovered, recovered_estimated,
             deaths_estimated])

        self.assertEqual(
            df[(df[date] == "2020-01-28") & (df[gender] == "male")]
            [recovered_estimated].item(),
            np.round(40 * 3. / 6.))
        self.assertEqual(
            df[(df[date] == "2020-01-29") & (df[gender] == "female")]
            [recovered_estimated].item(),
            np.round(33 * 3. / 7.))
        self.assertEqual(
            df[(df[date] == "2020-01-31") & (df[gender] == "male")]
            [recovered_estimated].item(),
            np.round(53 * 3. / 9.))
        self.assertEqual(
            df[(df[date] == "2020-01-31") & (df[gender] == "female")]
            [recovered_estimated].item(),
            np.round(43 * 3. / 9.))
        self.assertEqual(
            df[(df[date] == "2020-01-28") & (df[gender] == "male")]
            [deaths_estimated].item(),
            np.round(40 * 1. / 6.))
        self.assertEqual(
            df[(df[date] == "2020-01-29") & (df[gender] == "female")]
            [deaths_estimated].item(),
            np.round(33 * 1. / 7.))
        self.assertEqual(
            df[(df[date] == "2020-01-31") & (df[gender] == "male")]
            [deaths_estimated].item(),
            np.round(53 * 2. / 9.))
        self.assertEqual(
            df[(df[date] == "2020-01-31") & (df[gender] == "female")]
            [deaths_estimated].item(),
            np.round(43 * 2. / 9.))

    def test_download_weekly(self):
        directory = os.path.join(self.path, 'Germany/')
        gd.check_dir(directory)

        self.write_weekly_deaths_xlsx_data(
            directory, file_name='Cases_deaths_weekly_fake.xlsx')
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(len(os.listdir(directory)), 1)

        with patch('requests.get') as mock_request:
            df = gd.loadExcel(
                'Cases_deaths_weekly_fake', apiUrl=directory,
                extension='.xlsx',
                param_dict={"sheet_name": 'COVID_Todesfälle', "header": 0,
                            "engine": 'openpyxl'})
            towrite = io.BytesIO()
            df.to_excel(towrite, index=False)
            towrite.seek(0)
            mock_request.return_value.content = towrite.read()
            gcdwe.download_weekly_deaths_numbers(directory)
        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(len(os.listdir(directory)), 2)

        df_real_deaths_per_week = gd.loadExcel(
            'Cases_deaths_weekly', apiUrl=directory, extension='.xlsx',
            param_dict={"sheet_name": 0, "header": 0, "engine": 'openpyxl'})
        self.assertEqual(df_real_deaths_per_week.shape, (4, 3))
        self.assertEqual(pd.to_numeric(
            df_real_deaths_per_week['Sterbejahr'])[0], 2020)

    @patch('builtins.print')
    def test_except_non_existing_file(self, mock_print):

        [read_data, make_plot, file_format, out_folder, no_raw] \
            = [True, False, "json", self.path, False]

        # write files which should be read in by program

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)
        self.write_jh_data(directory)

        gcdwe.get_case_data_with_estimations(
            out_folder, read_data, file_format, no_raw, make_plot)

        # print is called 9 times, because no file exists
        self.assertEqual(len(mock_print.mock_calls), 9)


if __name__ == '__main__':
    unittest.main()
