import unittest
from pyfakefs import fake_filesystem_unittest

import os
import pandas as pd
from io import StringIO
import numpy as np

from epidemiology.epidata import getSpainData as gsp
from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from unittest.mock import patch


class test_get_spain_data(fake_filesystem_unittest.TestCase):
    path = '/home/SpainData'

    string_age = ("""fecha,rango_edad,sexo,casos_confirmados,hospitalizados,ingresos_uci,fallecidos
2020-03-23,0-9,ambos,129,34,1,0
2020-03-23,10-19,ambos,221,15,0,1
2020-03-23,20-29,ambos,1285,183,8,4
2020-03-23,30-39,ambos,2208,365,15,3""")

    string_state = ("""Fecha,cod_ine,CCAA,Casos,PCR+,TestAc+,Hospitalizados,UCI,Fallecidos
2020-05-20,10,C. Valenciana,,10987,3810,5747,730,1383
2020-05-20,11,Extremadura,,3042,1001,1780,110,505
2020-05-20,12,Galicia,,9077,2026,2943,334,608
2020-05-20,13,Madrid,,67049,5293,42497,3617,8931
2020-05-20,19,Melilla,,121,13,44,3,2
2020-05-20,14,Murcia,,1570,1039,680,112,149
2020-05-20,15,Navarra,,5195,3157,2048,136,506
2020-05-20,16,País Vasco,,13421,5358,7032,578,1483
2020-05-20,17,La Rioja,,4033,1395,1504,91,354""")

    def setUp(self):
        self.setUpPyfakefs()

    @patch('epidemiology.epidata.getSpainData.gd.loadCsv')
    def test_download_data_state(self, mock_csv):
        [read_data, out_form, out_folder, no_raw] = \
            [False, 'json', self.path, False]

        directory = os.path.join(out_folder, 'Spain/')
        gd.check_dir(directory)

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ['Spain'])

        mock_csv.return_value = pd.read_csv(StringIO(self.string_state))

        gsp. _get_spain_all_state(read_data, out_form, no_raw, directory)

        self.assertEqual(len(os.listdir(directory)), 2)
        self.assertEqual(os.listdir(directory).sort(), ["spain_all_state.json", "raw_spain_all_state.json"].sort())

    @patch('epidemiology.epidata.getSpainData.gd.loadCsv')
    def test_download_data_state_no_raw(self, mock_csv):
        [read_data, out_form, out_folder, no_raw] = \
            [False, 'json', self.path, True]

        directory = os.path.join(out_folder, 'Spain/')
        gd.check_dir(directory)

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ['Spain'])

        mock_csv.return_value = pd.read_csv(StringIO(self.string_state))

        gsp._get_spain_all_state(read_data, out_form, no_raw, directory)

        self.assertEqual(len(os.listdir(directory)), 1)
        self.assertEqual(os.listdir(directory).sort(), ["spain_all_state.json"].sort())

    @patch('epidemiology.epidata.getSpainData.gd.loadCsv')
    def test_download_data_age(self, mock_csv):
        [read_data, out_form, out_folder, no_raw] = \
            [False, 'json', self.path, False]

        directory = os.path.join(out_folder, 'Spain/')
        gd.check_dir(directory)

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ['Spain'])

        mock_csv.return_value = pd.read_csv(StringIO(self.string_age))

        gsp. _get_spain_all_age(read_data, out_form, no_raw, directory)

        self.assertEqual(len(os.listdir(directory)), 3)
        self.assertEqual(os.listdir(directory).sort(), ["spain.json", "spain_all_age.json", "raw_spain_all_age.json"].sort())

    @patch('epidemiology.epidata.getSpainData.gd.loadCsv')
    def test_download_data_age_no_raw(self, mock_csv):
        [read_data, out_form, out_folder, no_raw] = \
            [False, 'json', self.path, True]

        directory = os.path.join(out_folder, 'Spain/')
        gd.check_dir(directory)

        self.assertEqual(len(os.listdir(self.path)), 1)
        self.assertEqual(os.listdir(self.path), ['Spain'])

        mock_csv.return_value = pd.read_csv(StringIO(self.string_age))

        gsp._get_spain_all_age(read_data, out_form, no_raw, directory)

        self.assertEqual(len(os.listdir(directory)), 2)
        self.assertEqual(os.listdir(directory).sort(), ["spain.json", "spain_all_age.json"].sort())


if __name__ == '__main__':
    unittest.main()


