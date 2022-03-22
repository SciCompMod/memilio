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
import unittest
from unittest.mock import patch, call
from pyfakefs import fake_filesystem_unittest

import pandas as pd
import numpy as np

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import getCommuterMobility as gcm
from memilio.epidata import getVaccinationData as gvd
from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger


class TestGetVaccinationData(fake_filesystem_unittest.TestCase):
    maxDiff = None

    path = '/home/VaccinationData'

    col_names_vacc_data = [
        'Impfdatum', 'LandkreisId_Impfort', 'Altersgruppe', 'Impfschutz',
        'Anzahl']
    df_vacc_data = pd.DataFrame(columns=col_names_vacc_data)

    counties = geoger.get_county_ids(merge_eisenach=False)

    for county in counties:
        vacc_data = [
            ('2020-12-27', str(county), '05-11', 1, 3),
            ('2020-12-27', str(county), '05-11', 2, 2),
            ('2020-12-27', str(county), '05-11', 3, 1),
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

    df_vacc_data = df_vacc_data.astype(
        {'LandkreisId_Impfort': 'string', 'Altersgruppe': "string",
         'Impfschutz': int, 'Anzahl': int})

    df_vacc_data_altern = pd.DataFrame(columns=col_names_vacc_data)
    for i in range(len(counties)):
        vacc_data_altern = [
            ('2020-12-27', str(counties[i]), '02-03', 2, i),
            ('2020-12-27', str(counties[i]), '04-10', 2, i),
            ('2020-12-27', str(counties[i]), '11-17', 2, i),
            ('2020-12-27', str(counties[i]), '18-55', 2, i),
            ('2020-12-27', str(counties[i]), '56+', 1, i),
            ('2020-12-27', str(counties[i]), '56+', 2, i),
            ('2020-12-27', str(counties[i]), '56+', 3, i)
        ]
        df_to_append = pd.DataFrame(
            vacc_data_altern, columns=col_names_vacc_data)
        df_vacc_data_altern = df_vacc_data_altern.append(
            df_to_append, ignore_index=True)

    df_vacc_data_altern = df_vacc_data_altern.astype(
        {'LandkreisId_Impfort': 'string', 'Altersgruppe': "string",
         'Impfschutz': int, 'Anzahl': int})

    def setUp(self):
        self.setUpPyfakefs()

    @patch('memilio.epidata.getVaccinationData.download_vaccination_data',
           return_value=df_vacc_data_altern)
    def test_get_vaccination_data_alternative_ages_no_errors_with_plots(
            self, mockv):
        gvd.get_vaccination_data(out_folder=self.path)

    @patch('memilio.epidata.getVaccinationData.download_vaccination_data',
           return_value=df_vacc_data)
    def test_get_standard_vaccination_data_no_errors_with_plots(
            self, mockv):
        gvd.get_vaccination_data(out_folder=self.path)

    @patch('memilio.epidata.getVaccinationData.download_vaccination_data',
           return_value=df_vacc_data)
    def test_get_standard_vaccination_data_no_errors_with_plots(
            self, mockv):
        gvd.get_vaccination_data(out_folder=self.path, sanitize_data=1)

    @patch('memilio.epidata.getVaccinationData.download_vaccination_data',
           return_value=df_vacc_data)
    def test_get_standard_vaccination_data_no_errors_with_plots(
            self, mockv):
        gvd.get_vaccination_data(out_folder=self.path, sanitize_data=2)

    @patch('memilio.epidata.getVaccinationData.download_vaccination_data',
           return_value=df_vacc_data)
    def test_get_standard_vaccination_data_no_errors_with_plots(
            self, mockv):
        gvd.get_vaccination_data(out_folder=self.path, sanitize_data=3)

    @patch('memilio.epidata.getVaccinationData.download_vaccination_data',
           return_value=df_vacc_data)
    def test_get_standard_vaccination_data_no_errors_with_plots(
            self, mockv):
        gvd.get_vaccination_data(out_folder=self.path, sanitize_data=0)

    def test_download_vaccination_data(self):
        df = gvd.download_vaccination_data()
        self.assertFalse(
            df.empty,
            "Vaccination Data is empty. Shouldn't be.")

        # test agegroups in raw dataframe
        agegr_arr = df['Altersgruppe'].unique()
        self.assertIn('05-11', agegr_arr)
        self.assertIn('12-17', agegr_arr)
        self.assertIn('18-59', agegr_arr)
        self.assertIn('60+', agegr_arr)

        if 'u' in agegr_arr:
            self.assertEqual(len(agegr_arr), 5)
        if not 'u' in agegr_arr:
            self.assertEqual(len(agegr_arr), 4)

    @patch('builtins.print')
    @patch('memilio.epidata.getVaccinationData.pd.read_csv',
           side_effect=ImportError())
    def test_download_not_working(self, mock_download, mock_print):
        with self.assertRaises(ImportError):
            df = gvd.download_vaccination_data()
        expected_call = call(
            "Error in reading csv while downloading vaccination data.")
        mock_print.assert_has_calls([expected_call])

    def test_split_column_based_on_values(self):
        col_names_vacc_data = [
            'Impfdatum', 'LandkreisId_Impfort', 'Altersgruppe', 'Impfschutz',
            'Anzahl']

        vacc_data = [
            ('2020-12-27', '1001', '05-11', 1, 3),
            ('2020-12-27', '1001', '05-11', 2, 2),
            ('2020-12-27', '1001', '05-11', 3, 1),
            ('2020-12-27', '1001', '12-17', 1, 10),
            ('2020-12-27', '1001', '12-17', 2, 15),
            ('2020-12-27', '1001', '12-17', 3, 72),
            ('2020-12-27', '1001', '18-59', 1, 2),
            ('2020-12-27', '1001', '18-59', 2, 3),
            ('2020-12-27', '1001', '18-59', 3, 222),
            ('2020-12-27', '1001', '60+', 1, 22),
            ('2020-12-27', '1001', '60+', 2, 332),
            ('2020-12-27', '1001', '60+', 3, 76)
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
            3, 'Did not divide dataframe in three separate ones.')
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
            ('2020-12-27', '1001', '05-11', 1, 3),
            ('2020-12-27', '1001', '05-11', 2, 2),
            ('2020-12-27', '1001', '05-11', 3, 1),
            ('2020-12-27', '1001', '12-17', 1, 10),
            ('2020-12-27', '1001', '12-17', 2, 15),
            ('2020-12-27', '1001', '12-17', 3, 72),
            ('2020-12-27', '1001', '18-59', 1, 2),
            ('2020-12-27', '1001', '18-59', 2, 3),
            ('2020-12-27', '1001', '18-59', 3, 222),
            ('2020-12-27', '1001', '60+', 1, 22),
            ('2020-12-27', '1001', '60+', 2, 332),
            ('2020-12-27', '1001', '60+', 3, 76)
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
            3, 'Did not divide dataframe in three separate ones.')
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

    def test_interval_mapping(self):
        # testing refinement
        from_lower_bounds1 = [0, 3, 6, 15, 18, 25, 30, 40, 50, 65, 74, 100]
        to_lower_bounds1 = [0, 3, 5, 6, 12, 15, 18, 25, 30, 35, 40, 50, 60, 65,
                            74, 80, 100]
        # testing if the lists of bounds are not contained in one another
        from_lower_bounds2 = [0, 10, 20, 70, 100]
        to_lower_bounds2 = [0, 5, 15, 20, 50, 60, 70, 85, 90, 100]
        # testing with different boundaries
        from_lower_bounds3 = [5, 20, 30, 80, 85, 90]
        to_lower_bounds3 = [0, 15, 20, 60, 100]
        # testing error handling for invalid combination of boundaries
        from_lower_bounds4 = [0, 10, 100]
        to_lower_bounds4 = [10, 20, 90]

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
            [[0.5, 1], [0.5, 2]],
            [[0.6, 3], [0.2, 4], [0.2, 5]],
            [[1/2, 6], [1/6, 7], [1/3, 8]]]

        test_map3 = [
            [[2/3, 0], [1/3, 1]],
            [[1, 2]],
            [[3/5, 2], [2/5, 3]],
            [[1, 3]],
            [[1, 3]]]

        map_bounds1 = gvd.create_intervals_mapping(
            from_lower_bounds1, to_lower_bounds1)
        map_bounds2 = gvd.create_intervals_mapping(
            from_lower_bounds2, to_lower_bounds2)
        map_bounds3 = gvd.create_intervals_mapping(
            from_lower_bounds3, to_lower_bounds3)
        with self.assertRaises(ValueError):
            gvd.create_intervals_mapping(from_lower_bounds4, to_lower_bounds4)

        for test_map, calculated_map in zip(test_map1, map_bounds1):
            for test_val, calculated_val in zip(test_map, calculated_map):
                self.assertEqual(test_val[1], calculated_val[1])
                self.assertAlmostEqual(test_val[0], calculated_val[0])

        for test_map, calculated_map in zip(test_map2, map_bounds2):
            for test_val, calculated_val in zip(test_map, calculated_map):
                self.assertEqual(test_val[1], calculated_val[1])
                self.assertAlmostEqual(test_val[0], calculated_val[0])

        for test_map, calculated_map in zip(test_map3, map_bounds3):
            for test_val, calculated_val in zip(test_map, calculated_map):
                self.assertEqual(test_val[1], calculated_val[1])
                self.assertAlmostEqual(test_val[0], calculated_val[0])

if __name__ == '__main__':
    unittest.main()
