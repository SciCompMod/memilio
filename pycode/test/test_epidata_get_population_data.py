import unittest
from pyfakefs import fake_filesystem_unittest
from collections import namedtuple

import os
import pandas as pd

from epidemiology.epidata import getPopulationData as gpd
from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from unittest.mock import patch, mock_open


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

    def setUp(self):
        self.setUpPyfakefs()

    @patch('epidemiology.epidata.getPopulationData.gd.loadCsv')
    def test_gpd_download_data1(self, mock_loadCSV):

        mock_loadCSV.return_value = pd.read_json(self.test_string1)

        [read_data, make_plot, out_form, out_folder] = [False, False, "json", self.path]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        gpd.get_one_data_set(read_data, out_form, directory, self.d1)

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

    def test_gpd_read_data1(self):

        [read_data, make_plot, out_form, out_folder] = [True, False, "json", self.path]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        with open(os.path.join(directory, "FullDataB.json"), 'w') as f:
            f.write(self.test_string1)

        gpd.get_one_data_set(read_data, out_form, directory, self.d1)

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

        [read_data, make_plot, out_form, out_folder] = [False, False, "json", self.path]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        d = self.d2

        gpd.get_one_data_set(read_data, out_form, directory, d)

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

        [read_data, make_plot, out_form, out_folder] = [True, False, "json", self.path]

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        d = self.d2

        with open(os.path.join(directory, d.filename + ".json" ), 'w') as f:
            f.write(self.test_string2)

        gpd.get_one_data_set(read_data, out_form, directory, d)

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ["Germany"])

        self.assertEqual(len(os.listdir(directory)), 2)
        self.assertEqual(os.listdir(directory), [d.filename + ".json", d.filename_out + ".json"])

        f_path = os.path.join(directory, d.filename_out + ".json")
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_string2r)

    def test_gpd_read_data_full(self):

        [read_data, update_data, make_plot, out_form, out_folder] = [True, False, False, "json", self.path]

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


        gpd.get_population_data(read_data, update_data, make_plot, out_form, out_folder)

        #self.assertEqual(len(os.listdir(directory)), 4)
        self.assertEqual(len(os.listdir(directory)), 3)
        #self.assertEqual(os.listdir(directory), [file1, file2, file_out1, file_out2])
        self.assertEqual(os.listdir(directory), [file1, file2, file_out1])

        f_path = os.path.join(directory, file_out1)
        f = open(f_path, "r")
        self.assertEqual(f.read(), self.test_string1r)

        #f_path = os.path.join(directory, file_out2)
        #f = open(f_path, "r")
        #self.assertEqual(f.read(), self.test_string2r)


# TODO: How to test hdf5 export?

if __name__ == '__main__':
    unittest.main()