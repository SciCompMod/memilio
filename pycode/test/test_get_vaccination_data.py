######################################################################
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
######################################################################
import sys

import unittest
from unittest.mock import patch
from pyfakefs import fake_filesystem_unittest

import pandas as pd
import numpy as np

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import getCommuterMobility as gcm
from epidemiology.epidata import getVaccinationData as gvd
from epidemiology.epidata import defaultDict as dd


class TestGetVaccinationData(fake_filesystem_unittest.TestCase):
    maxDiff = None

    path = '/home/VaccinationData'

    col_names_vacc_data = [
        'Impfdatum', 'LandkreisId_Impfort', 'Altersgruppe', 'Impfschutz',
        'Anzahl']
    df_vacc_data = pd.DataFrame(columns=col_names_vacc_data)
    CountyMerging = {
        # Different districts to Berlin; reporting differs according to source
        11000: [11001, 11002, 11003, 11004, 11005, 11006, 11007, 11008, 11009,
                11010, 11011, 11012],
        # Wartburgkreis and Eisenach to Wartburgkreis (decision from July 1, 2021)
        16063: [16063, 16056]
    }
    counties = sorted(set(dd.County.keys()))
    for i in CountyMerging[11000]:
        counties.remove(i)

    for county in counties:
        vacc_data = [
            ('2020-12-27', str(county), '12-17', 1, 10),
            ('2020-12-27', str(county), '12-17', 2, 15),
            ('2020-12-27', str(county), '12-17', 3, 72),
            ('2020-12-27', str(county), '18-59', 1, 2),
            ('2020-12-27', str(county), '18-59', 2, 3),
            ('2020-12-27', str(county), '18-59', 3, 222),
            ('2020-12-27', str(county), '60+', 1, 22),
            ('2020-12-27', str(county), '60+', 2, 332),
            ('2020-12-27', str(county), '60+', 3, 76)
        ]
        df_to_append = pd.DataFrame(
            vacc_data, columns=col_names_vacc_data)
        df_vacc_data = df_vacc_data.append(
            df_to_append, ignore_index=True)

    df_vacc_data.astype({'LandkreisId_Impfort': 'string'}).dtypes
    df_vacc_data.astype({'Altersgruppe': 'string'}).dtypes
    df_vacc_data.astype({'Impfschutz': int}).dtypes
    df_vacc_data.astype({'Anzahl': int}).dtypes

    def setUp(self):
        self.setUpPyfakefs()

    @patch('epidemiology.epidata.getVaccinationData.download_vaccination_data',
           return_value=df_vacc_data)
    def test_get_standart_vaccination_data_no_errors_with_plots(
            self, mockv):
        gvd.get_vaccination_data(out_folder=self.path)

    def test_download_vaccination_data(self):
        df = gvd.download_vaccination_data()
        # Normally this should'nt be empty.
        self.assertFalse(
            df.empty,
            "Vaccination Data is empty. Should'nt be.")

    @patch('epidemiology.epidata.getVaccinationData.pd.read_csv',
           side_effect=ImportError())
    def test_download_not_working(self, mock_download):
        df = gvd.download_vaccination_data()
        self.assertTrue(df.empty, "Vaccination Data is empty.")

    def test_split_column_based_on_values(self):
        col_names_vacc_data = [
            'Impfdatum', 'LandkreisId_Impfort', 'Altersgruppe', 'Impfschutz',
            'Anzahl']

        vacc_data = [
            ('2020-12-27', str(10001), '12-17', 1, 10),
            ('2020-12-27', str(10001), '12-17', 2, 15),
            ('2020-12-27', str(10001), '12-17', 3, 72),
            ('2020-12-27', str(10001), '18-59', 1, 2),
            ('2020-12-27', str(10001), '18-59', 2, 3),
            ('2020-12-27', str(10001), '18-59', 3, 222),
            ('2020-12-27', str(10001), '60+', 1, 22),
            ('2020-12-27', str(10001), '60+', 2, 332),
            ('2020-12-27', str(10001), '60+', 3, 76)
        ]
        df_to_split = pd.DataFrame(
            vacc_data, columns=col_names_vacc_data)

        column_ident = 'Impfschutz'
        column_vals_name = 'Anzahl'
        new_col_labels = ['Vacc_partially', 'Vacc_completed', 'Vacc_refreshed']
        df_split = gvd.split_column_based_on_values(
            df_to_split, column_ident, column_vals_name, new_col_labels)

        new_column_names1 = [
            'Impfdatum',
            'LandkreisId_Impfort',
            'Altersgruppe',
            'Vacc_partially']
        new_column_names2 = [
            'Impfdatum',
            'LandkreisId_Impfort',
            'Altersgruppe',
            'Vacc_completed']
        new_column_names3 = [
            'Impfdatum',
            'LandkreisId_Impfort',
            'Altersgruppe',
            'Vacc_refreshed']
        new_column_names = []
        new_column_names.append(new_column_names1)
        new_column_names.append(new_column_names2)
        new_column_names.append(new_column_names3)

        indent_unique = df_to_split[column_ident].unique()
        number_unique_idents_df = []
        for i in range(0, len(indent_unique)):
            number_unique_idents_df.append(
                len(df_to_split[df_to_split[column_ident] == indent_unique[i]]))

        self.assertEqual(
            len(df_split),
            3, 'Did not divide dataframe in three sperate ones.')
        for i in range(0, len(new_col_labels)):
            self.assertEqual(
                number_unique_idents_df[i],
                len(df_split[i]),
                'Length of Dataframe is wrong.')
            column_names_i = df_split[i].columns
            for j in range(0, len(column_names_i)):
                self.assertEqual(
                    new_column_names[i][j],
                    column_names_i[j],
                    'Column names do not match.')

    def test_split_column_based_on_values_ident_length_wrong(self):
        col_names_vacc_data = [
            'Impfdatum', 'LandkreisId_Impfort', 'Altersgruppe', 'Impfschutz',
            'Anzahl']

        vacc_data = [
            ('2020-12-27', str(10001), '12-17', 1, 10),
            ('2020-12-27', str(10001), '12-17', 2, 15),
            ('2020-12-27', str(10001), '12-17', 3, 72),
            ('2020-12-27', str(10001), '18-59', 1, 2),
            ('2020-12-27', str(10001), '18-59', 2, 3),
            ('2020-12-27', str(10001), '18-59', 3, 222),
            ('2020-12-27', str(10001), '60+', 1, 22),
            ('2020-12-27', str(10001), '60+', 2, 332),
            ('2020-12-27', str(10001), '60+', 3, 76)
        ]
        df_to_split = pd.DataFrame(
            vacc_data, columns=col_names_vacc_data)

        column_ident = 'Impfschutz'
        column_vals_name = 'Anzahl'
        # wrong amount for data
        new_col_labels = ['Vacc_partially', 'Vacc_completed']
        df_split = gvd.split_column_based_on_values(
            df_to_split, column_ident, column_vals_name, new_col_labels)

        new_column_names1 = [
            'Impfdatum',
            'LandkreisId_Impfort',
            'Altersgruppe',
            'new_column_1']
        new_column_names2 = [
            'Impfdatum',
            'LandkreisId_Impfort',
            'Altersgruppe',
            'new_column_2']
        new_column_names3 = [
            'Impfdatum',
            'LandkreisId_Impfort',
            'Altersgruppe',
            'new_column_3']
        new_column_names = []
        new_column_names.append(new_column_names1)
        new_column_names.append(new_column_names2)
        new_column_names.append(new_column_names3)

        indent_unique = df_to_split[column_ident].unique()
        number_unique_idents_df = []
        for i in range(0, len(indent_unique)):
            number_unique_idents_df.append(
                len(df_to_split[df_to_split[column_ident] == indent_unique[i]]))

        self.assertEqual(
            len(df_split),
            3, 'Did not divide dataframe in three sperate ones.')
        for i in range(0, len(new_col_labels)):
            self.assertEqual(
                number_unique_idents_df[i],
                len(df_split[i]),
                'Length of Dataframe is wrong.')
            column_names_i = df_split[i].columns
            for j in range(0, len(column_names_i)):
                self.assertEqual(
                    new_column_names[i][j],
                    column_names_i[j],
                    'Column names do not match.')

    def test_intervall_mapping(self):
        lower_bounds1 = np.array(
            [0, 3, 6, 15, 18, 25, 30, 40, 50, 65, 74, 100])
        upper_bounds1 = np.array(
            [0, 3, 5, 6, 12, 15, 18, 25, 30, 35, 40, 50, 60, 65, 74, 80, 100])
        lower_bounds2 = np.array([0, 10, 20, 70, 100])
        upper_bounds2 = np.array([0, 5, 15, 20, 50, 60, 70, 85, 90, 100])

        test_map1 = [
            [[1, 0]],
            [[2 / 3, 1], [1 / 3, 2]],
            [[2 / 3, 3], [1 / 3, 4]],
            [[1, 5]],
            [[1, 6]],
            [[1, 7]],
            [[0.5, 8], [0.5, 9]],
            [[1, 10]],
            [[2/3, 11], [1/3, 12]],
            [[1, 13]],
            [[6/26, 14], [20/26, 15]]]

        test_map2 = [
            [[1 / 2, 0], [1 / 2, 1]],
            [[1.0, 1], [0.0, 2]],
            [[0.6, 3], [0.2, 4], [0.2, 5]],
            [[1/2, 6], [1/6, 7], [1/3, 8]]]

        map_bounds1 = gvd.create_intervals_mapping(
            lower_bounds1, upper_bounds1)
        map_bounds2 = gvd.create_intervals_mapping(
            lower_bounds2, upper_bounds2)

        for i in range(1, len(test_map1)):
            for j in range(1, len(test_map1[i])):
                self.assertTrue(
                    np.allclose(
                        np.array(test_map1[i][j]),
                        np.array(map_bounds1[i][j]),
                        rtol=1e-05),
                    "Not the same Arrays")

        for i in range(1, len(test_map2)):
            for j in range(1, len(test_map2[i])):
                self.assertTrue(
                    np.allclose(
                        np.array(test_map2[i][j]),
                        np.array(map_bounds2[i][j]),
                        rtol=1e-05),
                    "Not the same Arrays")


if __name__ == '__main__':
    unittest.main()
