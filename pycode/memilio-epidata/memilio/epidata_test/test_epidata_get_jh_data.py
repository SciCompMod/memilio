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
import os
import unittest
from datetime import date, datetime
from unittest.mock import patch

import pandas as pd
from pyfakefs import fake_filesystem_unittest

from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gD
from memilio.epidata import getJHData as gJHD


class TestGetJHData(fake_filesystem_unittest.TestCase):
    path = '/home/JHData'

    str_FullData_JohnHopkins = \
        (r"""[{"Date":"2020-01-22","Country\/Region":"China","Province\/State":"Anhui",
        "Confirmed":1,"Recovered":0.0,"Deaths":0},
        {"Date":"2020-02-05","Country\/Region":"Germany","Province\/State":null,
        "Confirmed":12,"Recovered":0.0,"Deaths":0},
        {"Date":"2020-02-07","Country\/Region":"Germany","Province\/State":null,
        "Confirmed":13,"Recovered":0.0,"Deaths":0},
        {"Date":"2021-01-04","Country\/Region":"Germany","Province\/State":null,
        "Confirmed":1796216,"Recovered":1445442.0,"Deaths":35748},
        {"Date":"2020-02-05","Country\/Region":"France","Province\/State":"Guadeloupe",
        "Confirmed":20,"Recovered":30,"Deaths":10},
        {"Date":"2020-09-26","Country\/Region":"France","Province\/State":"Guadeloupe",
        "Confirmed":4487,"Recovered":2199.0,"Deaths":42},
        {"Date":"2020-09-26","Country\/Region":"France","Province\/State":"French Polynesia",
        "Confirmed":1579,"Recovered":1335.0,"Deaths":6},
        {"Date":"2020-09-26","Country\/Region":"France","Province\/State":"French Guiana",
        "Confirmed":9863,"Recovered":9500.0,"Deaths":65},
        {"Date":"2020-09-26","Country\/Region":"France","Province\/State":"Martinique",
        "Confirmed":1290,"Recovered":98.0,"Deaths":20},
        {"Date":"2020-09-26","Country\/Region":"France","Province\/State":"Mayotte",
        "Confirmed":3541,"Recovered":2964.0,"Deaths":40},
        {"Date":"2020-02-26","Country\/Region":"Spain","Province\/State":null,
        "Confirmed":13,"Recovered":2.0,"Deaths":0},
        {"Date":"2020-03-02","Country\/Region":"Spain","Province\/State":null,
        "Confirmed":120,"Recovered":2.0,"Deaths":0},
        {"Date":"2020-02-02","Country\/Region":"Italy","Province\/State":null,
        "Confirmed":2,"Recovered":0.0,"Deaths":0},
        {"Date":"2020-02-03","Country\/Region":"Italy","Province\/State":null,
        "Confirmed":2,"Recovered":0.0,"Deaths":0},
        {"Date":"2020-03-31","Country\/Region":"US","Province\/State":null,
        "Confirmed":192301,"Recovered":7024.0,"Deaths":5278},
        {"Date":"2020-04-01","Country\/Region":"US","Province\/State":null,
        "Confirmed":224544,"Recovered":8474.0,"Deaths":6539},
        {"Date":"2020-02-20","Country\/Region":"Korea, South","Province\/State":null,
        "Confirmed":104,"Recovered":16.0,"Deaths":1},
        {"Date":"2020-02-21","Country\/Region":"Korea, South","Province\/State":null,
        "Confirmed":204,"Recovered":16.0,"Deaths":2},
        {"Date":"2020-01-22","Country\/Region":"Afghanistan","Province\/State":null,
        "Confirmed":0,"Recovered":0.0,"Deaths":0},
        {"Date":"2020-01-23","Country\/Region":"Afghanistan","Province\/State":null,
        "Confirmed":0,"Recovered":0.0,"Deaths":0}]""")

    def setUp(self):
        self.setUpPyfakefs()

    def write_jh_data(self, out_folder):
        file_jh = "FullData_JohnHopkins.json"
        file_jh_with_path = os.path.join(out_folder, file_jh)

        with open(file_jh_with_path, 'w') as f:
            f.write(self.str_FullData_JohnHopkins)

    @patch('builtins.input', return_value='n')
    def test_get_JH_Data(self, mockin):
        # Test without downloading data
        [read_data, file_format, out_folder, no_raw] \
            = [True, "json", self.path, False]

        gD.check_dir(out_folder)

        # Test case where file does not exist
        file = "FullData_JohnHopkins.json"
        file_with_path = os.path.join(out_folder, file)

        with self.assertRaises(FileNotFoundError) as error:
            gJHD.get_jh_data(read_data=read_data, file_format=file_format,
                             out_folder=out_folder, no_raw=no_raw)
        self.assertEqual(str(error.exception),
                         "Error: The file from " + file_with_path +
                         " does not exist. Call program without -r "
                         "flag to get it.")

        # Test case where file exists

        # write files which should be read in by program
        self.write_jh_data(out_folder)

        # check if expected file is written
        self.assertEqual(len(os.listdir(self.path)), 1)

        gJHD.get_jh_data(read_data=read_data, file_format=file_format,
                         out_folder=out_folder, no_raw=no_raw)

        # check if expected files are written
        # 7 country-folders+3 all countries-files
        self.assertEqual(len(os.listdir(self.path)), 3 + 7)
        # check if files are written in folders
        directory_ger = os.path.join(out_folder, 'Germany/')
        directory_es = os.path.join(out_folder, 'Spain/')
        directory_fr = os.path.join(out_folder, 'France/')
        directory_it = os.path.join(out_folder, 'Italy/')
        directory_us = os.path.join(out_folder, 'US/')
        directory_rok = os.path.join(out_folder, 'SouthKorea/')
        directory_prc = os.path.join(out_folder, 'China/')
        self.assertEqual(len(os.listdir(directory_ger)), 1)
        self.assertEqual(len(os.listdir(directory_es)), 1)
        self.assertEqual(len(os.listdir(directory_fr)), 1)
        self.assertEqual(len(os.listdir(directory_it)), 1)
        self.assertEqual(len(os.listdir(directory_us)), 1)
        self.assertEqual(len(os.listdir(directory_rok)), 1)
        self.assertEqual(len(os.listdir(directory_prc)), 1)

        # test whole Germany file
        f_read = os.path.join(directory_ger, "whole_country_Germany_jh.json")
        df = pd.read_json(f_read)

        data_list = df.columns.values.tolist()

        self.assertEqual(
            data_list,
            ["CountryRegion", "Date", "Confirmed", "Recovered", "Deaths"])
        self.assertEqual(df[(df["Date"] == "2021-01-04")]
                         ["Recovered"].item(), 1445442.0)

        # test whole France file (country with Provinces)
        f_read = os.path.join(directory_fr, "whole_country_France_jh.json")
        df = pd.read_json(f_read)

        data_list = df.columns.values.tolist()
        self.assertEqual(
            data_list,
            ["CountryRegion", "Date", "Confirmed", "Recovered", "Deaths"])
        # check if data is added up correctly
        self.assertEqual(
            df[(df["Date"] == "2020-09-26")]["Recovered"].item(),
            2199.0 + 1335.0 + 9500.0 + 98.0 + 2964.0)
        self.assertEqual(df[(df["Date"] == "2020-02-05")]
                         ["Recovered"].item(), 30)

        # test all_countries_jh file
        f_read = os.path.join(self.path, "all_countries_jh.json")
        df = pd.read_json(f_read)

        data_list = df.columns.values.tolist()

        self.assertEqual(
            data_list,
            ["CountryRegion", "Date", "Confirmed", "Recovered", "Deaths"])
        self.assertEqual('Germany' in df["CountryRegion"].values, True)
        self.assertEqual(df[(df["CountryRegion"] == 'Germany') & (
            df["Date"] == "2021-01-04")]["Recovered"].item(), 1445442.0)
        # check if Korea, South is renamed
        self.assertEqual('SouthKorea' in df["CountryRegion"].values, True)
        self.assertEqual('Afghanistan' in df["CountryRegion"].values, True)
        self.assertEqual('France' in df["CountryRegion"].values, True)
        self.assertEqual(
            df[(df["CountryRegion"] == 'France') &
               (df["Date"] == "2020-09-26")]["Recovered"].item(),
            2199.0 + 1335.0 + 9500.0 + 98.0 + 2964.0)
        self.assertEqual(
            df[(df["Date"] == "2020-02-05") &
               (df["CountryRegion"] == 'France')]["Recovered"].item(),
            30)
        self.assertEqual(
            df[(df["CountryRegion"] == 'France') &
               (df["Date"] == "2020-09-26")]["Recovered"].item(),
            2199.0 + 1335.0 + 9500.0 + 98.0 + 2964.0)
        self.assertEqual(
            df[(df["CountryRegion"] == 'France') &
               (df["Date"] == "2020-09-26")]["Confirmed"].item(),
            4487 + 1579 + 9863 + 1290 + 3541)
        self.assertEqual(
            df[(df["CountryRegion"] == 'France') &
               (df["Date"] == "2020-09-26")]["Deaths"].item(),
            42 + 6 + 65 + 20 + 40)

        # test all_provincestate_jh file
        f_read = os.path.join(self.path, "all_provincestate_jh.json")
        df = pd.read_json(f_read)

        data_list = df.columns.values.tolist()

        self.assertEqual(
            data_list,
            ["CountryRegion", "ProvinceState", "Date", "Confirmed",
             "Recovered", "Deaths"])
        self.assertEqual('Germany' in df["CountryRegion"].values, False)
        self.assertEqual('China' in df["CountryRegion"].values, True)
        self.assertEqual(df[(df["CountryRegion"] == 'China') & (
            df["Date"] == "2020-01-22")].shape[0], 1)
        self.assertEqual('France' in df["CountryRegion"].values, True)
        self.assertEqual(df[(df["CountryRegion"] == 'France') & (
            df["Date"] == "2020-09-26")].shape[0], 5)
        self.assertEqual(
            df
            [(df["CountryRegion"] == 'France') & (df["Date"] == "2020-09-26") &
             (df['ProvinceState'] == 'Martinique')]["Deaths"].item(),
            20)

    def test_get_JH_Data_extract_subframe_between_timeframe(self):

        [read_data, file_format, out_folder, no_raw] \
            = [True, "json", self.path, False]

        gD.check_dir(out_folder)

        # write files which should be read in by program
        self.write_jh_data(out_folder)

        # define start and end date
        test_start_date = date(2020, 12, 22)
        test_end_date = date(2021, 4, 22)

        # read in JH data with start and end date
        gJHD.get_jh_data(
            read_data=read_data, file_format=file_format,
            out_folder=out_folder, no_raw=no_raw, start_date=test_start_date,
            end_date=test_end_date)

        # read in german data
        directory_ger = os.path.join(out_folder, 'Germany/')
        f_read = os.path.join(directory_ger, "whole_country_Germany_jh.json")
        df_test_start_end_date = pd.read_json(f_read)

        # do the same without the date constraints
        gJHD.get_jh_data(read_data=read_data, file_format=file_format,
                         out_folder=out_folder, no_raw=no_raw)
        directory_ger = os.path.join(out_folder, 'Germany/')
        f_read = os.path.join(directory_ger, "whole_country_Germany_jh.json")
        df_test = pd.read_json(f_read)

        # extract dates which are between start and end date
        upperdate = datetime.strftime(test_end_date, '%Y-%m-%d')
        lowerdate = datetime.strftime(test_start_date, '%Y-%m-%d')
        df_test = df_test[df_test[dd.EngEng['date']] <= upperdate]
        df_test = df_test[df_test[dd.EngEng['date']] >= lowerdate]

        # compare if they are the same
        self.assertEqual(len(df_test), len(df_test_start_end_date),
                         "Dataframes don't have the same length.")
        self.assertEqual(
            list(df_test['Confirmed']),
            list(df_test_start_end_date['Confirmed']),
            "Dataframes don't have the same confirmed cases.")
        self.assertEqual(
            list(df_test['Recovered']),
            list(df_test_start_end_date['Recovered']),
            "Dataframes don't have the same recovered cases.")
        self.assertEqual(
            list(df_test['Deaths']),
            list(df_test_start_end_date['Deaths']),
            "Dataframes don't have the same death cases.")


if __name__ == '__main__':
    unittest.main()
