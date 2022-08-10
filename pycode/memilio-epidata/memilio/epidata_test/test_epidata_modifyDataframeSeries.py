#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Sascha Korf, Patrick Lenz
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
import pandas as pd
import numpy as np
from pyfakefs import fake_filesystem_unittest
from datetime import date
from memilio.epidata import modifyDataframeSeries as mdfs

class Test_modifyDataframeSeries(fake_filesystem_unittest.TestCase):
    test_df1 = pd.DataFrame(
        {
         'Date':
         ['2021-01-06', '2021-01-06', '2021-01-06', '2021-01-07', '2021-01-07',
          '2021-01-07', '2021-01-08', '2021-01-08', '2021-01-08', '2021-01-09',
          '2021-01-09', '2021-01-09', '2021-01-10', '2021-01-10',
          '2021-01-10'],
         'test_col1': [12, 3, 6, 0, 3, 1, 4, 7, 11, 15, 19, 19, 27, 13, 5],
         'test_col2': ['a', 'x', 't', 'a', 'b', 'a', 'x', 't', 'a', 'b', 'a', 'x', 't', 'a', 'b'],
         'test_col3': [1, 0, 1, 9, 4, 3, 2, 1, 1, 1, 0, 6, 5, 3, 1],
         'ID': [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]})
    test_df2 = pd.DataFrame(
        {
         'Date':
         ['2021-01-06', '2021-01-06', '2021-01-06', '2021-01-07', '2021-01-07',
          '2021-01-07', '2021-01-08', '2021-01-08', '2021-01-08', '2021-01-09',
          '2021-01-09', '2021-01-09', '2021-01-13', '2021-01-13',
          '2021-01-13'],
         'test_col1': [12, 3, 6, 0, 3, 1, 4, 7, 11, 15, 19, 19, 27, 13, 5],
         'test_col2': ['a', 'x', 't', 'a', 'b', 'a', 'x', 't', 'a', 'b', 'a', 'x', 't', 'a', 'b'],
         'test_col3': [1, 0, 1, 9, 4, 3, 2, 1, 1, 1, 0, 6, 5, 3, 1],
         'ID': [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]})
    df_dates = pd.DataFrame({
        'Date':[
            '2022-01-01', '2022-01-01', '2022-01-03', '2022-01-03',
            '2022-01-04', '2022-01-06'],
        'col_int': [1, 10, 3, 3, 4, 6],
        'col_str': ['a', 'a', 'b', 'c', 'd', 'e']})
    df_dates_unsorted = pd.DataFrame({
        'Date':[
            '2022-01-06', '2022-01-03', '2022-01-01', '2022-01-06',
            '2022-01-03', '2022-01-04'],
        'col_int': [60, 3, 1, 6, 3, 4],
        'col_str': ['a', 'b', 'a', 'e', 'c', 'd']})
    df_dates_result = pd.DataFrame({
        'Date': ['2022-01-03', '2022-01-03', '2022-01-04'],
        'col_int': [3, 3, 4], 'col_str': ['b', 'c', 'd']})
    int_map = [(10*i, i) for i in range(10)]
    str_map = [('A','a'), ('X','x'), ('T','t'), ('B','b')]
    df_str_map_col = pd.Series(
        data = ['A', 'X', 'T', 'A', 'B', 'A', 'X', 'T', 'A', 'B', 'A', 'X', 'T',
                'A', 'B'],
        name = 'inserted_col')

    date_df = pd.DataFrame(
        {
        'Date':
        ['2021-09-08', '2021-09-08', '2021-09-09', '2021-09-09', '2021-09-10',
         '2021-09-10', '2021-09-11', '2021-09-11', '2021-09-12', '2021-09-12',
         '2021-09-13', '2021-09-13', '2021-09-14', '2021-09-14', '2021-09-14',
         '2021-09-15'],
        'test_col1':
        [16, 52, 111, 7, 432, 126, 74, 175, 208, 33, 79, 16, 11, 27, 5, 15],
        'test_col2':
        [13, 34, 63, 5, 220, 53, 38, 79, 111, 15, 53, 8, 7, 13, 2, 9],
        'test_col_str':
        ['SH', 'H', 'N', 'B', 'NW', 'H', 'RP', 'BW', 'B', 'S', 'B', 'B',
         'MV', 'S', 'S', 'T']})
    test_date_df1 = pd.DataFrame({
        'Date': ['2021-09-12', '2021-09-12', '2021-09-13', '2021-09-13', '2021-09-14', '2021-09-14', '2021-09-14'],
        'test_col1': [208, 33, 79, 16, 11, 27, 5], 
        'test_col2': [111, 15, 53, 8, 7, 13, 2],
        'test_col_str': ['B', 'S', 'B', 'B', 'MV', 'S', 'S']})
    test_date_df2 = pd.DataFrame({
        'Date': ['2021-09-09', '2021-09-09'],
        'test_col1': [111, 7], 
        'test_col2': [63, 5],
        'test_col_str': ['N', 'B']})

    def setUp(self):
        self.setUpPyfakefs()

    def test_impute_and_reduce_df(self):

        group_by_cols = {'ID': sorted(set(self.test_df1['ID'].unique()))}
        mod_cols = ['test_col1', 'test_col3']

        # test impute forward and fill dates with moving average = 3
        df = mdfs.impute_and_reduce_df(
            self.test_df1, group_by_cols, mod_cols, impute='forward',
            moving_average=3, min_date='2021-01-05', max_date='2021-01-11',
            start_w_firstval=False)
        # test that values at new start days are zero since start_w_firstval=False
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-05") & (df['ID'] == 2.0)]['test_col1'].item(), (0+0+3)/3)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-05") & (df['ID'] == 2.0)]['test_col3'].item(), (0+0+0)/3)
        # test that the values at first original date are obtained by the value itself plus the value right of it divided by 3
        # (6 + 1) / 3 = 2 + 2 / 3
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 3.0)]['test_col1'].item(), 2 + 1 / 3)
        # (3 + 1) / 3 = 1 + 1 / 3
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 3.0)]['test_col3'].item(), 1 + 1 / 3)
        
        # check that last entry of original frame is computed by the value left of it and twice its value since impute='forward'
        # (15 + 27 + 27) / 3 = 23
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 1.0)]['test_col1'].item(), 23)
        # (1 + 5 + 5) / 3 = 3 + 2 / 3
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 1.0)]['test_col3'].item(), 3 + 2/ 3)
        
        # check that new date gets imputed the value the column had the day before because impute='forward'
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-11") & (df['ID'] == 3.0)]['test_col3'].item(), 1)

        # test impute zeros with moving average = 3
        df = mdfs.impute_and_reduce_df(
            self.test_df1, group_by_cols, mod_cols, impute='zeros',
            moving_average=3, min_date='2021-01-05', max_date='2021-01-11',
            start_w_firstval=False)
        # test that values at new start days are zero since start_w_firstval=False
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-05") & (df['ID'] == 2.0)]['test_col1'].item(), (0+0+3)/3)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-05") & (df['ID'] == 2.0)]['test_col3'].item(), (0+0+0)/3)

        # test that the values at first original date are obtained by the value itself plus the value right of it divided by 3
        # (6 + 1) / 3 = 2 + 2 / 3
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 3.0)]['test_col1'].item(), 2 + 1 / 3)
        # (3 + 1) / 3 = 1 + 1 / 3
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 3.0)]['test_col3'].item(), 1 + 1 / 3)
        
        # check that last entry of original frame is computed by the value left of it and the value itself because impute = "zeros"
        # (15 + 27) / 3 = 14
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 1.0)]['test_col1'].item(), 14)
        # (1 + 5) / 3 = 2
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 1.0)]['test_col3'].item(), 2)

        # check that new date gets imputed 0 because impute = "zeros"
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-11") & (df['ID'] == 3.0)]['test_col3'].item(), (1+0+0)/3)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-11") & (df['ID'] == 3.0)]['test_col1'].item(), (5+0+0)/3)

        # test fill missing dates moving average = 4
        # if moving average is an even number it always should calculate with one more earlier date
        df = mdfs.impute_and_reduce_df(
            self.test_df2, group_by_cols, mod_cols, impute='forward',
            moving_average=4, min_date='2021-01-06', max_date='2021-01-13',
            start_w_firstval=False)
        # test that the values at first original date arent changed since there is no value left of it
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 2.0)]['test_col3'].item(), (0+4+0+0)/4)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 2.0)]['test_col1'].item(), (3+3+0+0)/4)
        # test that a value is computed by two values left of it, the value itself and the value right of it divided by 4
        # (6 + 1 + 11 + 19) / 4 = 9 + 1 / 4
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-08") & (df['ID'] == 3.0)]['test_col1'].item(), 9 + 1 / 4)
        # (1 + 3 + 1 + 6) = 2 + 3 / 4
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-08") & (df['ID'] == 3.0)]['test_col3'].item(), 2 + 3 / 4)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-13") & (df['ID'] == 3.0)]['test_col3'].item(), (2*1+2*6)/4)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-13") & (df['ID'] == 3.0)]['test_col1'].item(), (2*5+2*19)/4)
        # test that the first of three missing dates is computed by three times the value left of it plus the penultimate value devided by 4
        # (19 + 19 + 19 + 11) / 4 = 17
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 3.0)]['test_col1'].item(), 17)
        # (6 + 6 + 6 + 1) / 4 = 4 + 3 / 4
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 3.0)]['test_col3'].item(), 4 + 3/4)

        # test mod_cols = ['test_col1']. test_col3 should not be modified
        mod_cols = ['test_col1']
        df = mdfs.impute_and_reduce_df(
            self.test_df2, group_by_cols, mod_cols, impute='forward',
            moving_average=4, min_date='2021-01-06', max_date='2021-01-13',
            start_w_firstval=False)
        # test same tests as in the previous test with moving average = 4
        # 'test_col1' should be same same as in the previous test
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 2.0)]['test_col1'].item(), (3+3+0+0)/4)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-08") & (df['ID'] == 3.0)]['test_col1'].item(), 9 + 1 / 4)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-13") & (df['ID'] == 3.0)]['test_col1'].item(), (2*5+2*19)/4)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 3.0)]['test_col1'].item(), 17)
        # 'test_col3' should not be changed
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 2.0)]['test_col3'].item(), 0)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-08") & (df['ID'] == 3.0)]['test_col3'].item(), 1)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-13") & (df['ID'] == 3.0)]['test_col3'].item(), 1)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 3.0)]['test_col3'].item(), 6)

        mod_cols = ['test_col1', 'test_col3']
        # test start date higher than end date
        # empty dataframe should be returned
        df = mdfs.impute_and_reduce_df(
            self.test_df1, group_by_cols, mod_cols, impute='forward',
            moving_average=4, min_date='2021-01-13', max_date='2021-01-06',
            start_w_firstval=False)
        edf = pd.DataFrame()
        self.assertEqual(len(edf),len(df))

        # test start_w_firstval = True
        df = mdfs.impute_and_reduce_df(
            self.test_df1, group_by_cols, mod_cols, impute='forward',
            moving_average=3, min_date='2021-01-05', max_date='2021-01-11',
            start_w_firstval=True)
        # test that values at new start days are the same as in the first original date since start_w_firstval=True
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-05") & (df['ID'] == 2.0)]['test_col1'].item(), 3)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-05") & (df['ID'] == 2.0)]['test_col3'].item(), 0)
        # test that the values at first original date are twice the value itself plus the value right of it divided by 3 since start_w_firstval=True
        # (6 + 6 + 1) / 3 = 4 + 2 / 3
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 3.0)]['test_col1'].item(), 4 + 1 / 3)
        # (1+ 1 + 3) / 3 = 1 + 2 / 3
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 3.0)]['test_col3'].item(), 1 + 2 / 3)

    def test_split_column_based_on_values(self):
        col_names_vacc_data = [
            'Impfdatum', 'LandkreisId_Impfort', 'Altersgruppe', 'Impfschutz',
            'Anzahl']

        vacc_data = [
            ('2020-12-27', '1001', '05-11', 1, 3),
            ('2020-12-27', '1001', '05-11', 2, 2),
            ('2020-12-27', '1001', '05-11', 3, 1),
            ('2020-12-27', '1001', '05-11', 4, 1),
            ('2020-12-27', '1001', '12-17', 1, 10),
            ('2020-12-27', '1001', '12-17', 2, 15),
            ('2020-12-27', '1001', '12-17', 3, 72),
            ('2020-12-27', '1001', '12-17', 4, 20),
            ('2020-12-27', '1001', '18-59', 1, 2),
            ('2020-12-27', '1001', '18-59', 2, 3),
            ('2020-12-27', '1001', '18-59', 3, 222),
            ('2020-12-27', '1001', '18-59', 4, 12),
            ('2020-12-27', '1001', '60+', 1, 22),
            ('2020-12-27', '1001', '60+', 2, 332),
            ('2020-12-27', '1001', '60+', 3, 76),
            ('2020-12-27', '1001', '60+', 4, 2)
        ]
        df_to_split = pd.DataFrame(
            vacc_data, columns=col_names_vacc_data)

        test_df = pd.DataFrame(
            {'Vacc_partially': [3, 10, 2, 22],
             'Vacc_completed': [2, 15, 3, 332],
             'Vacc_refreshed': [1, 72, 222, 76],
             'Vacc_refreshed_2': [1, 20, 12, 2]})

        groupby_list = ['Impfdatum', 'LandkreisId_Impfort', 'Altersgruppe']
        column_ident = 'Impfschutz'
        column_vals_name = 'Anzahl'
        new_col_labels = ['Vacc_partially', 'Vacc_completed', 'Vacc_refreshed', 'Vacc_refreshed_2']
        test_labels = new_col_labels.copy()
        returned_column_labels, df_split = mdfs.split_column_based_on_values(
            df_to_split, column_ident, column_vals_name, groupby_list, new_col_labels, compute_cumsum=False)

        self.assertEqual(test_labels, returned_column_labels)
        pd.testing.assert_frame_equal(df_split[test_labels],test_df)
        new_column_names = [
            'Impfdatum', 'LandkreisId_Impfort', 'Altersgruppe',
            'Vacc_partially', 'Vacc_completed', 'Vacc_refreshed',
            'Vacc_refreshed_2']

        for i in range(0, len(new_column_names)):
            self.assertIn(df_split.columns[i], new_column_names)
            self.assertIn(new_column_names[i], df_split.columns)

        # test wrong length of new columns

        vacc_data_length_wrong = [
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
        df_to_split_length_wrong = pd.DataFrame(
            vacc_data_length_wrong, columns=col_names_vacc_data)

        test_df_length_wrong = pd.DataFrame(
            {'Vacc_partially': [3, 10, 2, 22],
             'Vacc_completed': [2, 15, 3, 332],
             'Vacc_completed_2': [1, 72, 222, 76]})

        # wrong amount for data
        new_col_labels_length_wrong = ['Vacc_partially', 'Vacc_completed']
        test_labels_length_wrong = new_col_labels_length_wrong.copy()
        returned_column_labels_length_wrong, df_split_length_wrong = mdfs.split_column_based_on_values(
            df_to_split_length_wrong, column_ident, column_vals_name, groupby_list, new_col_labels_length_wrong, compute_cumsum=False)

        new_column_names_length_wrong = [
            'Impfdatum', 'LandkreisId_Impfort', 'Altersgruppe',
            'Vacc_partially', 'Vacc_completed', 'Vacc_completed_2']

        col_labels = ['Vacc_partially', 'Vacc_completed', 'Vacc_completed_2']

        self.assertNotEqual(test_labels_length_wrong, returned_column_labels_length_wrong)
        pd.testing.assert_frame_equal(df_split_length_wrong[col_labels],test_df_length_wrong)

        for i in range(0, len(new_column_names_length_wrong)):
            self.assertIn(df_split_length_wrong.columns[i], new_column_names_length_wrong)
            self.assertIn(new_column_names_length_wrong[i], df_split_length_wrong.columns)

    def test_split_column_based_on_values_compute_cumsum(self):
        df_to_split = pd.DataFrame({'Date': 
            3*(['2022-05-11' for i in range(0, 16)] +
            ['2022-05-12' for i in range(0, 16)] +
            ['2022-05-13' for i in range(0, 16)]),
            'ID_County': 9*sorted(8*[1000+i for i in range (0,2)]),
            'Age_RKI': 48*['18+'] + 48*['12-17']+48*['0-11'],
            'Impfschutz':36*[1,2,3,4],
            'Anzahl':[i for i in range(0,144)]
            })
        groupby_list = ['Date', 'ID_County', 'Age_RKI']
        column_ident = 'Impfschutz'
        column_vals_name = 'Anzahl'
        new_col_labels = ['Vacc_partially', 'Vacc_completed', 'Vacc_refreshed', 'Vacc_refreshed_2']
        test_labels = new_col_labels.copy()
        returned_column_labels, df_split = mdfs.split_column_based_on_values(
            df_to_split, column_ident, column_vals_name, groupby_list, new_col_labels, compute_cumsum=True)
        # test returned labels
        self.assertEqual(test_labels, returned_column_labels)
        # test if cumulative sum is correct
        df_county1 = df_split[df_split.ID_County == 1000]
        df_county2 = df_split[df_split.ID_County == 1001]
        # second dataframe should have + 16 at first date, then cumulativ sum
        for column in returned_column_labels:
            self.assertEqual(list(df_county1[column].values + np.array(3*[16]+ 3*[32]+ 3*[48])), list(df_county2[column].values))
        # test return for county 1000 in age group 18+
        return_1 = pd.DataFrame({'Vacc_partially':[4,40,108],'Vacc_completed':[6,44,114], 'Vacc_refreshed':[8,48,120], 'Vacc_refreshed_2':[10,52,126]})
        pd.testing.assert_frame_equal(df_county1[df_county1.Age_RKI=='18+'][returned_column_labels].reset_index(drop=True), return_1)
        # test return for county 1001 in age group 0-11
        return_2 = pd.DataFrame({'Vacc_partially':[212,456,732],'Vacc_completed':[214,460,738], 'Vacc_refreshed':[216,464,744], 'Vacc_refreshed_2':[218,468,750]})
        pd.testing.assert_frame_equal(df_county2[df_county2.Age_RKI=='0-11'][returned_column_labels].reset_index(drop=True), return_2)

    def test_extract_subframe_based_on_dates(self):
        test_df = self.df_dates.copy()
        # test with start and end value in dataframe
        extracted_df = mdfs.extract_subframe_based_on_dates(
            self.df_dates, date(2022, 1, 3), date(2022, 1, 5))
        pd.testing.assert_frame_equal(extracted_df, self.df_dates_result)
        # check that input frame is not manipulated in the function
        pd.testing.assert_frame_equal(self.df_dates, test_df)
        # test with start and end value not in dataframe
        extracted_df = mdfs.extract_subframe_based_on_dates(
            self.df_dates, date(2021, 1, 2), date(2023, 1, 5))
        pd.testing.assert_frame_equal(extracted_df, self.df_dates)
        pd.testing.assert_frame_equal(self.df_dates, test_df)
        # test with unsorted dataframe
        extracted_df = mdfs.extract_subframe_based_on_dates(
            self.df_dates_unsorted, date(2022, 1, 3), date(2022, 1, 5))
        pd.testing.assert_frame_equal(extracted_df, self.df_dates_result)


    def test_insert_column_by_map(self):
        old_cols = self.test_df1.columns.to_list()

        #test with integer mapping
        df = mdfs.insert_column_by_map(
            self.test_df1, 'test_col3', 'inserted_col', self.int_map)
        new_cols = df.columns.to_list()
        exp_cols = ['Date', 'test_col1', 'test_col2', 'test_col3',
            'inserted_col', 'ID']
        self.assertEqual(new_cols, exp_cols)
        pd.testing.assert_frame_equal(df[old_cols], self.test_df1)
        exp_new_col = (10*self.test_df1['test_col3']).rename('inserted_col')
        pd.testing.assert_series_equal(df['inserted_col'], exp_new_col)

        # test with string mapping
        df = mdfs.insert_column_by_map(
            self.test_df1, 'test_col2', 'inserted_col', self.str_map)
        new_cols = df.columns.to_list()
        exp_cols = ['Date', 'test_col1', 'test_col2', 'inserted_col',
            'test_col3', 'ID']
        self.assertEqual(new_cols, exp_cols)
        pd.testing.assert_frame_equal(df[old_cols], self.test_df1)
        pd.testing.assert_series_equal(df['inserted_col'], self.df_str_map_col)

    def test_extract_subframe_based_on_dates_single_date(self):
        # test if only dates from 2021-09-09 are returned
        extracted_df = mdfs.extract_subframe_based_on_dates(
            self.date_df, date(2021, 9, 9),
            date(2021, 9, 9))
        pd.testing.assert_frame_equal(self.test_date_df2, extracted_df)
    
    def test_extract_subframe_based_on_dates_multiple_dates(self):
        # test if only dates from 2021-09-12 to 2021-09-14 are returned
        extracted_df = mdfs.extract_subframe_based_on_dates(
            self.date_df, date(2021, 9, 12),
            date(2021, 9, 14))
        pd.testing.assert_frame_equal(self.test_date_df1, extracted_df)

if __name__ == '__main__':
    unittest.main()
