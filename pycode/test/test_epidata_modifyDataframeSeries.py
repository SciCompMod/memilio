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
from pyfakefs import fake_filesystem_unittest
from epidemiology.epidata import modifyDataframeSeries as mDfS

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


    def setUp(self):
        self.setUpPyfakefs()

    def test_impute_and_reduce_df(self):

        group_by_cols = {'ID': sorted(set(self.test_df1['ID'].unique()))}
        mod_cols = ['test_col1', 'test_col3']

        # test impute forward and fill dates with moving average = 3
        df = mDfS.impute_and_reduce_df(
            self.test_df1, group_by_cols, mod_cols, impute='forward',
            moving_average=3, min_date='2021-01-05', max_date='2021-01-11',
            start_w_firstval=False)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 3.0)]['test_col3'].item(), 1 + 1 / 3)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 3.0)]['test_col1'].item(), 2 + 1 / 3)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-05") & (df['ID'] == 2.0)]['test_col3'].item(), 0)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-05") & (df['ID'] == 2.0)]['test_col1'].item(), 0)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-11") & (df['ID'] == 3.0)]['test_col3'].item(), 1)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 1.0)]['test_col1'].item(), 23)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 1.0)]['test_col3'].item(), 3+2/3)

        # test impute zeros with moving average = 3
        df = mDfS.impute_and_reduce_df(
            self.test_df1, group_by_cols, mod_cols, impute='zeros',
            moving_average=3, min_date='2021-01-05', max_date='2021-01-11',
            start_w_firstval=False)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 3.0)]['test_col3'].item(), 1 + 1 / 3)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 3.0)]['test_col1'].item(), 2 + 1 / 3)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-05") & (df['ID'] == 2.0)]['test_col3'].item(), 0)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-05") & (df['ID'] == 2.0)]['test_col1'].item(), 0)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-11") & (df['ID'] == 3.0)]['test_col3'].item(), 0)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-11") & (df['ID'] == 3.0)]['test_col1'].item(), 0)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 1.0)]['test_col1'].item(), 14)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 1.0)]['test_col3'].item(), 2)

        # test fill missing dates moving average = 2
        # if moving average is an even number it always should calculate with one more earlier date
        df = mDfS.impute_and_reduce_df(
            self.test_df2, group_by_cols, mod_cols, impute='forward',
            moving_average=4, min_date='2021-01-06', max_date='2021-01-13',
            start_w_firstval=False)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-08") & (df['ID'] == 3.0)]['test_col3'].item(), 2 + 3 / 4)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-08") & (df['ID'] == 3.0)]['test_col1'].item(), 9 + 1 / 4)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 2.0)]['test_col3'].item(), 0)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-06") & (df['ID'] == 2.0)]['test_col1'].item(), 3)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-13") & (df['ID'] == 3.0)]['test_col3'].item(), 1)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-13") & (df['ID'] == 3.0)]['test_col1'].item(), 5)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 1.0)]['test_col1'].item(), 12 + 1/4)
        self.assertAlmostEqual(df[(df['Date'] == "2021-01-10") & (df['ID'] == 1.0)]['test_col3'].item(), 1 + 1/4)

        # test start date higher than end date
        # empty dataframe should be returned
        df = mDfS.impute_and_reduce_df(
            self.test_df1, group_by_cols, mod_cols, impute='forward',
            moving_average=4, min_date='2021-01-13', max_date='2021-01-06',
            start_w_firstval=False)
        edf = pd.DataFrame()
        self.assertEqual(len(edf),len(df))



if __name__ == '__main__':
    unittest.main()
