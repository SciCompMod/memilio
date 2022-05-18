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

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import getVaccinationData as gvd
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
            ('2020-12-27', str(county), '60+', 3, 76),
            ('2020-12-27', str(county), '18-59', 4, 1)
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
            ('2020-12-27', str(counties[i]), '18-59', 4, 1)
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
    def test_get_standard_vaccination_sanitize_3(
            self, mockv):
        gvd.get_vaccination_data(out_folder=self.path, sanitize_data=3)

    @patch('memilio.epidata.getVaccinationData.download_vaccination_data',
           return_value=df_vacc_data)
    def test_get_standard_vaccination_sanitize_0(
            self, mockv):
        gvd.get_vaccination_data(out_folder=self.path, sanitize_data=0)

    @patch('memilio.epidata.getVaccinationData.pd.read_csv',
           return_value=df_vacc_data_altern)
    def test_sanity_checks(self, mockv):

        # test empty dataframe
        df = pd.DataFrame()
        with self.assertRaises(gd.DataError) as error:
            gvd.sanity_checks(df)
        error_message = "Download of Vaccination Data failed. File is empty."
        self.assertEqual(str(error.exception), error_message)

        # test more agegroups than expected
        df = pd.DataFrame(
            {'Altersgruppe': ['02-04', '05-11', '12-17', '18-59', '60+'],
             'b': [4, 5, 6, 7, 8]})
        with self.assertRaises(gd.DataError) as error:
            gvd.sanity_checks(df)
        error_message = "Number of agegroups has changed. Please report this as an issue."
        self.assertEqual(str(error.exception), error_message)

        # test different agegroups 
        df = pd.DataFrame(
            {'Altersgruppe': ['02-04', '05-12', '13-17', '18+'],
             'b': [4, 5, 6, 7]})
        with self.assertRaises(gd.DataError) as error:
            gvd.sanity_checks(df)
        error_message = "Agegroups have changed. Please report this as an issue."
        self.assertEqual(str(error.exception), error_message)

        # test if sanity checks are executed
        with self.assertRaises(gd.DataError) as error:
            gvd.download_vaccination_data()
        error_message = "Number of agegroups has changed. Please report this as an issue."
        self.assertEqual(str(error.exception), error_message)


        # Following test should not raise any errors
        
        # test actual agegroups
        df = pd.DataFrame(
            {'Altersgruppe': ['05-11', '12-17', '18-59', '60+'],
             'b': [4, 5, 6, 7]})
        gvd.sanity_checks(df)

        # test actual agegroups with unknown
        df = pd.DataFrame(
            {'Altersgruppe': ['05-11', '12-17', '18-59', '60+', 'u'],
             'b': [4, 5, 6, 7, 8]})
        gvd.sanity_checks(df)

        # test uniqueness
        df = pd.DataFrame({
            'Altersgruppe':
            ['05-11', '60+', '12-17', '05-11', 'u', '12-17',
             '18-59', '18-59', '60+', 'u', 'u', 'u', 'u'],
            'b': [1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3]})
        gvd.sanity_checks(df)


    @patch('builtins.print')
    @patch('memilio.epidata.getVaccinationData.pd.read_csv',
           side_effect=ImportError())
    def test_download_not_working(self, mock_download, mock_print):
        with self.assertRaises(ImportError):
            df = gvd.download_vaccination_data()
        expected_call = call(
            "Error in reading csv while downloading vaccination data.")
        mock_print.assert_has_calls([expected_call])


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

    def test_sanitizing_based_on_regions(self):
        to_county_map={0:[1001,1002,2001],1:[6000],2:[6005,6006]}
        age_groups = ['0-1','2-3','4-10','11+']
        data = pd.DataFrame({
            'Date':['2022-05-11'for i in range(0,24)],
             'ID_State':sorted(4*[1,1,2,6,6,6]),
             'ID_County':sorted(4*[1001,1002,2001,6000,6005,6006]),
             'Age_RKI':6*age_groups,
             'vacc_1':[0,0,0,0, 2,5,7,9, 2,4,6,8, 4,4,4,4, 1,6,1,2, 0,0,5,17],
             'vacc_2':[0,1,0,2, 1,4,3,2, 1,1,6,4, 4,4,4,1, 2,1,2,0, 0,5,4,0]
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
        test_3['vacc_2'] = [18/7, 27/7, 18/7 ,4.0, 6/3.8, 6*2.8/3.8]

        test_4 = data[data['Age_RKI'] ==
                      '11+'].drop(['vacc_1', 'vacc_2'], axis=1)
        test_4['vacc_1'] = [17*3/8.5, 17*2/8.5, 17*3.5/8.5, 4.0, 0.0, 19.0]
        test_4['vacc_2'] = [8*3/8.5, 8*2/8.5, 8*3.5/8.5, 1.0, 0.0, 0.0]
        returndata = gvd.sanitizing_based_on_regions(data, to_county_map, age_groups, [
                                                     'vacc_1', 'vacc_2'], population)
        pd.testing.assert_frame_equal(
            returndata[returndata['Age_RKI'] == '0-1'],
            test_1, check_dtype=False)
        pd.testing.assert_frame_equal(
            returndata[returndata['Age_RKI'] == '2-3'],
            test_2, check_dtype=False)
        pd.testing.assert_frame_equal(
            returndata[returndata['Age_RKI'] == '4-10'],
            test_3, check_dtype=False)
        pd.testing.assert_frame_equal(
            returndata[returndata['Age_RKI'] == '11+'],
            test_4, check_dtype=False)
if __name__ == '__main__':
    unittest.main()
