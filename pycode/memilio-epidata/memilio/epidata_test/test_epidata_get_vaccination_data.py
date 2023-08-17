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
from unittest.mock import patch

import os
import json
import pandas as pd
from pyfakefs import fake_filesystem_unittest

from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import getVaccinationData as gvd
from memilio.epidata import modifyDataframeSeries as mdfs
from memilio.epidata import progress_indicator


class TestGetVaccinationData(fake_filesystem_unittest.TestCase):
    maxDiff = None

    path = '/home/VaccinationData'

    here = os.path.dirname(os.path.abspath(__file__))

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
            ('2020-12-27', str(county), '60+', 3, 76),
            ('2020-12-27', str(county), '60+', 4, 1)
        ]
        df_to_concatenate = pd.DataFrame(
            vacc_data, columns=col_names_vacc_data)
        df_vacc_data = pd.concat(
            [df_vacc_data, df_to_concatenate], ignore_index=True)

    df_vacc_data = df_vacc_data.astype(
        {'LandkreisId_Impfort': 'string', 'Altersgruppe': "string",
         'Impfschutz': int, 'Anzahl': int})

    df_vacc_data_altern = pd.DataFrame(columns=col_names_vacc_data)
    for i in range(len(counties)):
        vacc_data_altern = [
            ('2020-12-27', str(counties[i]), '02-03', 1, 3),
            ('2020-12-27', str(counties[i]), '04-10', 1, 10),
            ('2020-12-27', str(counties[i]), '11-17', 1, 2),
            ('2020-12-27', str(counties[i]), '18-55', 1, 5),
            ('2020-12-27', str(counties[i]), '56+', 1, 7),
            ('2020-12-27', str(counties[i]), '02-03', 2, 17),
            ('2020-12-27', str(counties[i]), '04-10', 2, 36),
            ('2020-12-27', str(counties[i]), '11-17', 2, 40),
            ('2020-12-27', str(counties[i]), '18-55', 2, 6),
            ('2020-12-27', str(counties[i]), '56+', 2, 11),
            ('2020-12-27', str(counties[i]), '02-03', 3, 20),
            ('2020-12-27', str(counties[i]), '04-10', 3, 40),
            ('2020-12-27', str(counties[i]), '11-17', 3, 40),
            ('2020-12-27', str(counties[i]), '18-55', 3, 6),
            ('2020-12-27', str(counties[i]), '56+', 3, 57),
            ('2020-12-27', str(counties[i]), '56+', 4, 1)
        ]
        df_to_concatenate = pd.DataFrame(
            vacc_data_altern, columns=col_names_vacc_data)
        df_vacc_data_altern = pd.concat(
            [df_vacc_data_altern, df_to_concatenate], ignore_index=True)

    df_vacc_data_altern = df_vacc_data_altern.astype(
        {'LandkreisId_Impfort': 'string', 'Altersgruppe': "string",
         'Impfschutz': int, 'Anzahl': int})

    filename = os.path.join(
        here, 'test_data', 'TestSetPopulationFinal.json')
    with open(filename) as file_object:
        df_pop = pd.DataFrame(json.load(file_object))

    def setUp(self):
        self.setUpPyfakefs()
        progress_indicator.ProgressIndicator.disable_indicators(True)

    @patch('memilio.epidata.getVaccinationData.download_vaccination_data',
           return_value=df_vacc_data_altern)
    @patch('memilio.epidata.getPopulationData.get_population_data', return_value=df_pop)
    @patch('builtins.input', return_value='y')
    def test_get_vaccination_data_alternative_ages(self, mockin, mockp, mockv):
        gvd.get_vaccination_data(out_folder=self.path)

    @patch('memilio.epidata.getVaccinationData.download_vaccination_data',
           return_value=df_vacc_data)
    @patch('memilio.epidata.getPopulationData.get_population_data', return_value=df_pop)
    @patch('builtins.input', return_value='y')
    def test_get_standard_vaccination_sanitize_3(self, mockin, mockp, mockv):
        gvd.get_vaccination_data(out_folder=self.path, sanitize_data=3)

    @patch('memilio.epidata.getVaccinationData.pd.read_csv',
           return_value=df_vacc_data_altern)
    def test_sanity_checks(self, mockv):
        # test empty dataframe
        df_empty = pd.DataFrame()
        with self.assertRaises(gd.DataError) as error:
            gvd.sanity_checks(df_empty)
        error_message = "Download of Vaccination Data failed. File is empty."
        self.assertEqual(str(error.exception), error_message)

        # test wrong number of data categories
        df_wrong_categories = pd.DataFrame(
            {"Impfdatum": ["2022-01-12"],
             "LandkreisId_Impfort": ["05754"],
             "Altersgruppe": ["01-59"],
             "Impfschutz": [1],
             "Anzahl": [10000],
             "xyz": ['Error']})
        with self.assertRaises(gd.DataError) as error:
            gvd.sanity_checks(df_wrong_categories)
        error_message = "Error: Number of data categories changed."
        self.assertEqual(str(error.exception), error_message)

        # test wrong columns names
        df_wrong_column_names = pd.DataFrame(
            {"Impfdatum": ["2022-01-12"],
             "LandkreisId_Impfort": ["05754"],
             "Altersgruppe": ["01-59"],
             "Impfschutz": [1],
             "xyz": ['Error']})
        with self.assertRaises(gd.DataError) as error:
            gvd.sanity_checks(df_wrong_column_names)
        error_message = "Error: Data categories have changed."
        self.assertEqual(str(error.exception), error_message)

        # test no errors
        # this test should not raise any errors
        df_no_errors = pd.DataFrame(
            {"Impfdatum":
                ["2022-01-12", "2022-01-12", "2022-01-12", "2022-01-12",
                 "2022-01-12"],
                "LandkreisId_Impfort": ['05754', '1', '2', '3', '4'],
                "Altersgruppe": ["01-59", "01-59", "01-59", "01-59", "01-59"],
                "Impfschutz": [1, 1, 2, 3, 1],
                "Anzahl": [10000, 1, 2, 3, 4]})
        gvd.sanity_checks(df_no_errors)

    def test_sanitizing_based_on_regions(self):
        to_county_map = {0: [1001, 1002, 2001], 1: [6000], 2: [6005, 6006]}
        age_groups = ['0-1', '2-3', '4-10', '11+']
        data = pd.DataFrame({
            'Date': ['2022-05-11'for i in range(0, 24)],
            'ID_State': sorted(4*[1, 1, 2, 6, 6, 6]),
            'ID_County': sorted(4*[1001, 1002, 2001, 6000, 6005, 6006]),
            'Age_RKI': 6*age_groups,
            'vacc_1': [0, 0, 0, 0, 2, 5, 7, 9, 2, 4, 6, 8, 4, 4, 4, 4, 1, 6, 1, 2, 0, 0, 5, 17],
            'vacc_2': [0, 1, 0, 2, 1, 4, 3, 2, 1, 1, 6, 4, 4, 4, 4, 1, 2, 1, 2, 0, 0, 5, 4, 0]
        })
        population = pd.DataFrame({
            'ID_County': [1001, 1002, 2001, 6000, 6005, 6006],
            '0-1': [100, 200, 150, 200, 100, 100],
            '2-3': [400, 200, 300, 150, 150, 400],
            '4-10': [2000, 3000, 2000, 5000, 1000, 2800],
            '11+': [30000, 20000, 35000, 100000, 0, 20000]
        })
        test_1 = data[data['Age_RKI'] ==
                      '0-1'].drop(['vacc_1', 'vacc_2'], axis=1)
        test_1['vacc_1'] = [4/4.5, 8/4.5, 6/4.5, 4.0, 0.5, 0.5]
        test_1['vacc_2'] = [2/4.5, 4/4.5, 3/4.5, 4.0, 1.0, 1.0]

        test_2 = data[data['Age_RKI'] ==
                      '2-3'].drop(['vacc_1', 'vacc_2'], axis=1)
        test_2['vacc_1'] = [4.0, 2.0, 3.0, 4.0, 9/5.5, 6*4/5.5]
        test_2['vacc_2'] = [8/3, 8/6, 2.0, 4.0, 9/5.5, 6*4/5.5]

        test_3 = data[data['Age_RKI'] ==
                      '4-10'].drop(['vacc_1', 'vacc_2'], axis=1)
        test_3['vacc_1'] = [26/7, 39/7, 26/7, 4.0, 6/3.8, 6*2.8/3.8]
        test_3['vacc_2'] = [18/7, 27/7, 18/7, 4.0, 6/3.8, 6*2.8/3.8]

        test_4 = data[data['Age_RKI'] ==
                      '11+'].drop(['vacc_1', 'vacc_2'], axis=1)
        test_4['vacc_1'] = [17*3/8.5, 17*2/8.5, 17*3.5/8.5, 4.0, 0.0, 19.0]
        test_4['vacc_2'] = [8*3/8.5, 8*2/8.5, 8*3.5/8.5, 1.0, 0.0, 0.0]
        returndata = gvd.sanitizing_average_regions(data, to_county_map, age_groups, [
            'vacc_1', 'vacc_2'], population)
        pd.testing.assert_frame_equal(
            returndata[returndata['Age_RKI'] == '0-1'].reset_index(drop=True),
            test_1.reset_index(drop=True), check_dtype=False)
        pd.testing.assert_frame_equal(
            returndata[returndata['Age_RKI'] == '2-3'].reset_index(drop=True),
            test_2.reset_index(drop=True), check_dtype=False)
        pd.testing.assert_frame_equal(
            returndata[returndata['Age_RKI'] == '4-10'].reset_index(drop=True),
            test_3.reset_index(drop=True), check_dtype=False)
        pd.testing.assert_frame_equal(
            returndata[returndata['Age_RKI'] == '11+'].reset_index(drop=True),
            test_4.reset_index(drop=True), check_dtype=False)

    def test_extrapolate_age_groups_vaccinations(self):
        unique_age_groups_old = ['00-04', '05-11', '12-17', '18+']
        unique_age_groups_new = ['00-05', '06-16', '17-19', '20+']
        column_names = ['vacc_1', 'vacc_2']
        age_old_to_all_ages_indices = [[0], [1, 2], [3, 4], [5, 6]]
        min_all_ages = ['0', '5', '6', '12', '17', '18', '20']
        all_ages_to_age_new_share = [
            [[1, 0]],
            [[1, 0]],
            [[1, 1]],
            [[1, 1]],
            [[1, 2]],
            [[1, 2]],
            [[1, 3]]]

        df_data = pd.DataFrame({'Date': ['2022-05-11' for i in range(
            0, 24)] +
            ['2022-05-12' for i in range(0, 24)] +
            ['2022-05-13' for i in range(0, 24)],
            'ID_State': 3 * sorted(
            4 * [1, 1, 2, 6, 6, 6]),
            'ID_County': 3 *
            sorted(
            4 *
            [1001, 1002, 2001, 6000, 6005, 6006]),
            'County': 3 *
            sorted(
            4 *
            ['test_county_1', 'test_county_2',
             'test_county_3', 'test_county_4',
             'test_county_5', 'test_county_6']),
            'Age_RKI': 18 *
            ['00-04', '05-11', '12-17', '18+'],
            'vacc_1':
            [0, 0, 0, 0, 2, 5, 7, 9, 2, 4, 6, 8, 4, 4, 4,
             4, 1, 6, 1, 2, 0, 0, 5, 17, 8, 3, 7, 2, 12,
             5, 7, 9, 2, 4, 4, 8, 1, 5, 1, 2, 1, 6, 3, 2,
             10, 2, 5, 13, 10, 10, 12, 8, 12, 9, 9, 20, 12,
             18, 0, 14, 12, 11, 9, 13, 12, 12, 1, 2, 10, 8,
             5, 7],
            'vacc_2':
            [0, 1, 0, 2, 1, 4, 3, 2, 1, 1, 6, 4, 4, 4, 4,
             1, 2, 1, 2, 0, 0, 5, 4, 0, 5, 2, 3, 1, 0, 9,
             7, 8, 4, 3, 6, 4, 4, 1, 8, 3, 8, 1, 2, 7, 6,
             7, 3, 1, 7, 1, 0, 2, 1, 4, 3, 15, 1, 1, 16,
             4, 4, 17, 4, 11, 5, 1, 2, 8, 6, 5, 4, 3]})

        population_all_ages = pd.DataFrame({
            'ID_County': [1001, 1002, 2001, 6000, 6005, 6006],
            '0': [100, 200, 150, 200, 100, 100],
            '5': [400, 200, 300, 150, 150, 400],
            '6': [2000, 3000, 2000, 5000, 1000, 2800],
            '12': [200, 40, 500, 700, 100, 600],
            '17': [1000, 600, 350, 750, 400, 2200],
            '18': [1600, 400, 800, 900, 1100, 2000],
            '20': [30000, 20000, 35000, 100000, 0, 20000]})

        # return should be population ratio in new age group * vaccination
        # for all old age groups
        test_1 = pd.DataFrame({'vacc_1': [0, 2+200/3200*5, 2+300/2300*4, 4+150/5150*4, 1+150/1150*6, 0], 'vacc_2': [
                              0+400/2400*1, 1+200/3200*4, 1+300/2300*1, 4+150/5150*4, 2+150/1150*1, 400/3200*5]})

        test_2 = pd.DataFrame(
            {
                'vacc_1':
                [0, 3000 / 3200 * 5 + 40 / 640 * 7, 2000 / 2300 * 4 + 500 / 850 * 6,
                 5000 / 5150 * 4 + 700 / 1450 * 4, 1000 / 1150 * 6 + 100 / 500 * 1,
                 600 / 2800 * 5],
                'vacc_2':
                [2000 / 2400 * 1, 3000 / 3200 * 4 + 40 / 640 * 3, 2000 / 2300 * 1 +
                    500 / 850 * 6, 5000 / 5150 * 4 + 700 / 1450 * 4, 1000 / 1150 * 1 +
                    100 / 500 * 2, 2800 / 3200 * 5 + 600 / 2800 * 4]})

        test_3 = pd.DataFrame({
            'vacc_1':
            [0, 7 * 600 / 640 + 9 * 400 / 20400, 6 * 350 /
             850 + 8 * 800 / 35800, 4 * 750 / 1450 + 4 * 900 /
             100900, 400 / 500 + 2, 5 * 2200 / 2800 + 17 *
             2000 / 22000],
            'vacc_2':
            [2 * 1600 / 31600, 3 * 600 / 640 + 2 * 400 /
             20400, 6 * 350 / 850 + 4 * 800 / 35800, 4 * 750
             / 1450 + 1 * 900 / 100900, 2 * 400 / 500, 4 *
             2200 / 2800]})

        test_4 = pd.DataFrame({'vacc_1': [0, 9*20000/20400, 8*35000/35800, 4*100000/100900, 0, 17*20000/22000],
                              'vacc_2': [2*30000/31600, 2*20000/20400, 4*35000/35800, 1*100000/100900, 0, 0]})

        returndata = gvd.extrapolate_age_groups_vaccinations(
            df_data.copy(), population_all_ages, unique_age_groups_old,
            unique_age_groups_new, column_names, age_old_to_all_ages_indices,
            min_all_ages, all_ages_to_age_new_share)

        # test if number of all vaccinations remains the same
        self.assertEqual(returndata['vacc_1'].sum(), df_data['vacc_1'].sum())
        self.assertEqual(returndata['vacc_2'].sum(), df_data['vacc_2'].sum())

        # test if one day is calculated correctly
        returndata_one_day = returndata[returndata['Date'] == '2022-05-11']
        # test each age group seperately
        pd.testing.assert_frame_equal(
            returndata_one_day
            [returndata_one_day['Age_RKI'] == '00-05']
            [column_names].reset_index(drop=True),
            test_1)
        pd.testing.assert_frame_equal(
            returndata_one_day
            [returndata_one_day['Age_RKI'] == '06-16']
            [column_names].reset_index(drop=True),
            test_2)
        pd.testing.assert_frame_equal(
            returndata_one_day
            [returndata_one_day['Age_RKI'] == '17-19']
            [column_names].reset_index(drop=True),
            test_3)
        pd.testing.assert_frame_equal(
            returndata_one_day
            [returndata_one_day['Age_RKI'] == '20+']
            [column_names].reset_index(drop=True),
            test_4)


if __name__ == '__main__':
    unittest.main()
