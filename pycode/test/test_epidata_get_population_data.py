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
from logging import exception
from sys import exec_prefix
import unittest
from pyfakefs import fake_filesystem_unittest
from collections import namedtuple

import os
import pandas as pd
import numpy as np

from epidemiology.epidata import getPopulationData as gpd
from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from unittest.mock import patch


class Test_getPopulationData(fake_filesystem_unittest.TestCase):

    path = '/home/x'

    Data = namedtuple("Data", "filename item columns_wanted filename_out")

    d1 = Data("FullDataB", '5dc2fc92850241c3be3d704aa0945d9c_2', ["LAN_ew_RS", 'LAN_ew_GEN', 'LAN_ew_EWZ'],
                  "PopulStates")
    d2 = Data("FullDataL", 'b2e6d8854d9744ca88144d30bef06a76_1', ['RS', 'GEN', 'EWZ'], "PopulCounties")

    test_string1 = ("""[\
{"FID":1,"LAN_ew_RS":1,"LAN_ew_AGS":1,"LAN_ew_SDV_RS":10020000000,"LAN_ew_GEN":"Schleswig-Holstein",\
"LAN_ew_BEZ":"Land","LAN_ew_IBZ":20,"LAN_ew_BEM":"--","LAN_ew_SN_L":1,"LAN_ew_SN_R":0,"LAN_ew_SN_K":0,"LAN_ew_SN_V1":0,\
"LAN_ew_SN_V2":0,"LAN_ew_SN_G":0,"LAN_ew_FK_S3":0,"LAN_ew_NUTS":"DEF","LAN_ew_WSK":"2012\/02\/01 00:00:00",\
"LAN_ew_EWZ":2889821,"LAN_ew_KFL":15804.35,"SHAPE_Length":20.9191264621,"SHAPE_Area":2.1595456768},\
{"FID":2,"LAN_ew_RS":2,"LAN_ew_AGS":2,"LAN_ew_SDV_RS":20000000000,"LAN_ew_GEN":"Hamburg",\
"LAN_ew_BEZ":"Freie und Hansestadt","LAN_ew_IBZ":22,"LAN_ew_BEM":"--","LAN_ew_SN_L":2,"LAN_ew_SN_R":0,"LAN_ew_SN_K":0,\
"LAN_ew_SN_V1":0,"LAN_ew_SN_V2":0,"LAN_ew_SN_G":0,"LAN_ew_FK_S3":0,"LAN_ew_NUTS":"DE6",\
"LAN_ew_WSK":"1974\/01\/01 00:00:00","LAN_ew_EWZ":1830584,"LAN_ew_KFL":755.09,"SHAPE_Length":3.1095198283,\
"SHAPE_Area":0.1001785991},\
{"FID":3,"LAN_ew_RS":3,"LAN_ew_AGS":3,"LAN_ew_SDV_RS":32410001001,"LAN_ew_GEN":"Niedersachsen","LAN_ew_BEZ":"Land",\
"LAN_ew_IBZ":20,"LAN_ew_BEM":"--","LAN_ew_SN_L":3,"LAN_ew_SN_R":0,"LAN_ew_SN_K":0,"LAN_ew_SN_V1":0,"LAN_ew_SN_V2":0,\
"LAN_ew_SN_G":0,"LAN_ew_FK_S3":0,"LAN_ew_NUTS":"DE9","LAN_ew_WSK":"2015\/01\/01 00:00:00","LAN_ew_EWZ":7962775,\
"LAN_ew_KFL":47709.8,"SHAPE_Length":29.8156067698,"SHAPE_Area":6.3454724588}]""")

    test_string1r = ("""[\
{"ID_State":1,"State":"Schleswig-Holstein","Population":2889821},\
{"ID_State":2,"State":"Hamburg","Population":1830584},\
{"ID_State":3,"State":"Niedersachsen","Population":7962775}]""")

    test_string2 = """[\
{"FID":1,"RS":1001,"AGS":1001,"SDV_RS":10010000000,"GEN":"Flensburg","BEZ":"Kreisfreie Stadt","IBZ":40,"BEM":"--",\
"SN_L":1,"SN_R":0,"SN_K":1,"SN_V1":0,"SN_V2":0,"SN_G":0,"FK_S3":"R","NUTS":"DEF01","WSK":"2008\/01\/01 00:00:00",\
"EWZ":88519,"KFL":56.73,"Kennziffer":1001,"EWZ_18":null,"SHAPE_Length":0.5247234366,"SHAPE_Area":0.0068727541},\
{"FID":2,"RS":1002,"AGS":1002,"SDV_RS":10020000000,"GEN":"Kiel","BEZ":"Kreisfreie Stadt","IBZ":40,"BEM":"--",\
"SN_L":1,"SN_R":0,"SN_K":2,"SN_V1":0,"SN_V2":0,"SN_G":0,"FK_S3":"R","NUTS":"DEF02","WSK":"2006\/01\/01 00:00:00",\
"EWZ":247943,"KFL":118.65,"Kennziffer":1002,"EWZ_18":null,"SHAPE_Length":1.2755450552,"SHAPE_Area":0.0155057123}]"""
#{"FID":3,"RS":1003,"AGS":1003,"SDV_RS":10030000000,"GEN":"L\u00fcbeck","BEZ":"Kreisfreie Stadt","IBZ":40,"BEM":"--",
#"SN_L":1,"SN_R":0,"SN_K":3,"SN_V1":0,"SN_V2":0,"SN_G":0,"FK_S3":"R","NUTS":"DEF03","WSK":"2006\/02\/01 00:00:00",
#"EWZ":216318,"KFL":214.19,"Kennziffer":1003,"EWZ_18":null,"SHAPE_Length":1.8350372077,"SHAPE_Area":0.0289309207}

    test_string2r = """[{"ID_County":1001,"County":"Flensburg","Population":88519},\
{"ID_County":1002,"County":"Kiel","Population":247943}]"""
#{"ID_County":1003,"County":"L\u00fcbeck","Population":216318}

    test_old_counties = np.zeros((18,2))
    test_old_counties[:,0] = [3152, 3156, 13056, 13002, 13055, 13052, 13051, 13053, 13061, 13005, 13057, 13006, 13058,
                              13059, 13062, 13001, 13054, 13060]
    test_old_counties[:,1] = np.arange(len(test_old_counties))

    test_new_counties = np.zeros((7,2))
    test_new_counties[:,0] = [3159, 13071, 13072, 13073, 13074, 13075, 13076]
    test_new_counties[:,1] = [1,14,13,27,23, 42, 33]

    data = np.zeros((5, 30))
    data[:,0] = np.arange(1,6)
    for i in range(len(data)):
        data[i, 3] = 22*(i+1)
        data[i, 4] = 11*(i+1)
        data[i,5:-2] = 1*(i+1)
        data[i, 16] = 11*(i+1)

    test_zensus = pd.DataFrame(data, columns=["FID", "DES", "Name", "EWZ", "Gesamt_Maennlich", 'M_Unter_3', 'M_3_bis_5',
                                              'M_6_bis_14', 'M_15_bis_17', 'M_18_bis_24','M_25_bis_29', 'M_30_bis_39',
                                              'M_40_bis_49', 'M_50_bis_64','M_65_bis_74', 'M_75_und_aelter', "Gesamt_Weiblich",
                                              'W_Unter_3', 'W_3_bis_5', 'W_6_bis_14', 'W_15_bis_17', 'W_18_bis_24',
                                              'W_25_bis_29', 'W_30_bis_39', 'W_40_bis_49', 'W_50_bis_64',
                                              'W_65_bis_74', 'W_75_und_aelter', 'SHAPE_Length', 'SHAPE_Area'])
    test_zensus["DES"] = "Kreis"
    test_zensus["Name"] = ["Hogwarts", "Narnia", "MittelErde", "Westeros", "Wakanda"]

    data = np.zeros((5, 3))
    data[:, 0] = [1001, 1002, 1003, 1004, 1005]
    data[:, 2] = [(x+1)*22/1000 for x in range(len(data))]

    test_reg_key = pd.DataFrame(data, columns=['AGS', 'NAME', 'Zensus_EWZ'])
    test_reg_key['NAME'] = ["Hogwarts", "Narnia", "MittelErde", "Westeros", "Wakanda"]

    data = np.zeros((5,2))
    data[:, 0] = [1001, 1002, 1003, 1004, 1005]
    data[:, 1] = [(x+1)*44 for x in range(len(data))]

    test_counties = pd.DataFrame(data, columns=['Schlüssel-nummer', 'Bevölkerung2)'])

    columns = ['ID_County', 'Total', '<3 years', '3-5 years', '6-14 years', '15-17 years', '18-24 years',
               '25-29 years', '30-39 years', '40-49 years', '50-64 years',
               '65-74 years', '>74 years']

    data = np.zeros((5, len(columns)))
    for i in range(len(data)):
        data[i, 0] = 1001 + i
        data[i, 1] = 2*22*(i+1)
        data[i, 2:] = 2*2*(i+1)
    test_population_result = pd.DataFrame(data, columns=columns)
    test_population_result = test_population_result.astype('int64')

    data = np.zeros((5, len(columns)))
    for i in range(len(data)):
        data[i, 0] = 1001 + i
        data[i, 1] = 44 * (i+1)
        data[i, 2:] = 4 * (i+1)

    test_current_population_result = pd.DataFrame(data, columns=columns)
    test_current_population_result = test_current_population_result.astype('int64')

    def setUp(self):
            self.setUpPyfakefs()

    @patch('epidemiology.epidata.getPopulationData.gd.loadCsv')
    def test_gpd_download_data1(self, mock_loadCSV):

        mock_loadCSV.return_value = pd.read_json(self.test_string1)

        [read_data, file_format, out_folder, no_raw] = [False, "json", self.path, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        gpd.get_one_data_set(read_data, file_format, no_raw, directory, self.d1)

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ["Germany"])

        self.assertEqual(len(os.listdir(directory)), 2)
        self.assertEqual(os.listdir(directory), ["FullDataB.json", "PopulStates.json"])

        f_path = os.path.join(directory, "FullDataB.json")
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_string1)

        f_path = os.path.join(directory, "PopulStates.json")
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_string1r)

    @patch('epidemiology.epidata.getPopulationData.gd.loadCsv')
    def test_gpd_download_data1_no_raw(self, mock_loadCSV):

        mock_loadCSV.return_value = pd.read_json(self.test_string1)

        [read_data, file_format, out_folder, no_raw] = [False, "json", self.path, True]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        gpd.get_one_data_set(read_data, file_format, no_raw, directory, self.d1)

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ["Germany"])

        self.assertEqual(len(os.listdir(directory)), 1)
        self.assertEqual(os.listdir(directory), ["PopulStates.json"])

    def test_gpd_read_data1(self):

        [read_data, file_format, out_folder, no_raw] = [True,  "json", self.path, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        with open(os.path.join(directory, "FullDataB.json"), 'w') as f:
            f.write(self.test_string1)

        gpd.get_one_data_set(read_data, file_format, no_raw, directory, self.d1)

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ["Germany"])

        self.assertEqual(len(os.listdir(directory)), 2)
        self.assertEqual(os.listdir(directory), ["FullDataB.json", "PopulStates.json"])

        f_path = os.path.join(directory, "PopulStates.json")
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_string1r)

    @patch('epidemiology.epidata.getPopulationData.gd.loadCsv')
    def test_gpd_download_data2(self, mock_loadCSV):

        mock_loadCSV.return_value = pd.read_json(self.test_string2)

        [read_data, file_format, out_folder, no_raw] = [False, "json", self.path, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        d = self.d2

        gpd.get_one_data_set(read_data, file_format, no_raw, directory, d)

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ["Germany"])

        self.assertEqual(len(os.listdir(directory)), 2)
        self.assertEqual(os.listdir(directory), [d.filename + ".json", d.filename_out + ".json"])

        f_path = os.path.join(directory, d.filename + ".json")
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_string2)

        f_path = os.path.join(directory, d.filename_out + ".json")
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_string2r)

    def test_gpd_read_data2(self):

        [read_data, file_format, out_folder, no_raw] = [True, "json", self.path, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        d = self.d2

        with open(os.path.join(directory, d.filename + ".json" ), 'w') as f:
            f.write(self.test_string2)

        gpd.get_one_data_set(read_data, file_format, no_raw, directory, d)

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ["Germany"])

        self.assertEqual(len(os.listdir(directory)), 2)
        self.assertEqual(os.listdir(directory), [d.filename + ".json", d.filename_out + ".json"])

        f_path = os.path.join(directory, d.filename_out + ".json")
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_string2r)

    def test_gpd_read_data_full(self):

        [read_data, file_format, out_folder, no_raw] = [True, "json", self.path, False]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        file1 = "FullDataB.json"
        file2 = "FullDataL.json"

        file_out1 = "PopulStates.json"
        file_out2 = "PopulCounties.json"

        with open(os.path.join(directory, file1), 'w') as f:
            f.write(self.test_string1)
        with open(os.path.join(directory, file2), 'w') as f:
            f.write(self.test_string2)

        gpd.get_population_data(read_data, file_format, out_folder)

        self.assertEqual(len(os.listdir(directory)), 3)
        self.assertEqual(os.listdir(directory), [file1, file2, file_out1])

        f_path = os.path.join(directory, file_out1)
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_string1r)

    def test_get_new_counties(self):
        test = gpd.get_new_counties(self.test_old_counties)
        self.assertTrue(np.array_equal(test, self.test_new_counties))

    @patch('epidemiology.epidata.getPopulationData.load_age_population_data', return_value=(test_counties, test_reg_key, test_zensus))
    def test_get_age_population(self, mock_data):

        gpd.get_age_population_data(False, 'json', self.path)

        test_df = pd.read_json(os.path.join(self.path, 'Germany/', 'county_population.json'))
        pd.testing.assert_frame_equal(test_df, self.test_population_result)

        test_df = pd.read_json(os.path.join(self.path, 'Germany/', 'county_current_population.json'))
        pd.testing.assert_frame_equal(test_df, self.test_current_population_result)

    @patch('pandas.read_excel', return_value=test_counties)
    @patch('pandas.read_excel', return_value=test_reg_key)
    @patch('epidemiology.epidata.getDataIntoPandasDataFrame.loadCsv', return_value=test_zensus)
    def test_load_age_population_data(self, mock_read_excel1, mock_read_excel2, mock_read_csv):

        directory = os.path.join(self.path, 'Germany/')

        counties_write, reg_key_write, zensus_write = gpd.load_age_population_data(self.path)
        self.assertEqual(len(os.listdir(directory)), 3)

        counties_read, reg_key_read, zensus_read = gpd.load_age_population_data(self.path)

        pd.testing.assert_frame_equal(counties_read, counties_write, check_dtype=False)
        pd.testing.assert_frame_equal(reg_key_read, reg_key_write, check_dtype=False)
        pd.testing.assert_frame_equal(zensus_read, zensus_write, check_dtype=False)

    @patch('epidemiology.epidata.getPopulationData.gd.loadCsv')
    @patch('epidemiology.epidata.getPopulationData.pandas.read_json')
    @patch('epidemiology.epidata.getPopulationData.gd.loadExcel')
    def test_errors(self, mocklexcel, mockrjson, mocklcsv):
        mockrjson.side_effect = ValueError
        [read_data, file_format, out_folder, no_raw] = [True, "json", self.path, True]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)
        with self.assertRaises(SystemExit) as cm:
            gpd.get_one_data_set(read_data, file_format, no_raw, directory, self.d1)
        file_in = os.path.join(self.path, "Germany/FullDataB.json")
        exit_string = "Error: The file: " + file_in + " does not exist. "\
            "Call program without -r flag to get it."
        self.assertEqual(cm.exception.code, exit_string)

        mocklexcel.side_effect = ValueError
        with self.assertRaises(SystemExit) as cm:
            gpd.load_age_population_data(self.path)
        exit_string = "Error: The counties file does not exist."
        self.assertEqual(cm.exception.code, exit_string)

        mocklexcel.side_effect = None
        mocklexcel.return_value = self.test_zensus.copy()
        mocklcsv.side_effect = ValueError
        with self.assertRaises(SystemExit) as cm:
            gpd.load_age_population_data(self.path)
        exit_string = "Error: The zensus file does not exist."
        self.assertEqual(cm.exception.code, exit_string)


if __name__ == '__main__':
    unittest.main()