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

from datetime import datetime
import unittest
import numpy as np
import pandas as pd
import pandas.testing as pd_testing
from unittest.mock import patch
from pyfakefs import fake_filesystem_unittest
from epidemiology.epidata import modifyDataframeSeries as mDfS


class Test_modifyDataframeSeries(fake_filesystem_unittest.TestCase):
    filldates_df = pd.DataFrame(
        {'Date': ['2021-01-05', '2021-01-06', '2021-01-08'],
         'test_col1': [11, 3, 6],
         'test_col2': ['a', 'x', 't']})
    test_df1 = pd.DataFrame(
        {'Date': 
        [datetime(2021, 1, 5),
         datetime(2021, 1, 6),
         datetime(2021, 1, 7),
         datetime(2021, 1, 8)],
         'test_col1': [11.0, 3.0, 3.0, 6.0],
         'test_col2': ['a', 'x', 'x', 't']})
    test_df2 = pd.DataFrame({
        'Date':
        [datetime(2021, 1, 4),
         datetime(2021, 1, 5),
         datetime(2021, 1, 6),
         datetime(2021, 1, 7),
         datetime(2021, 1, 8)],
        'test_col1': [11.0, 11.0, 3.0, 3.0, 6.0],
        'test_col2': ['a', 'a', 'x', 'x', 't']})
    test_df3 = pd.DataFrame({
        'Date':
        [datetime(2021, 1, 4),
         datetime(2021, 1, 5),
         datetime(2021, 1, 6),
         datetime(2021, 1, 7),
         datetime(2021, 1, 8)],
        'test_col1': [11.0, 11.0, 3.0, 3.0, 6.0],
        'test_col2': ['a', 'a', 'x', 'x', 't']})

    def setUp(self):
        self.setUpPyfakefs()

    def test_impute_and_reduce_df(self):
        group_by_cols = {}
        mod_cols = ''
        # test fill date
        test_df_new = mDfS.impute_and_reduce_df(
            self.filldates_df, group_by_cols, mod_cols, impute='forward',
            moving_average=0, min_date='2021-01-05', max_date='2021-01-08',
            start_w_firstval=False)
        pd_testing.assert_frame_equal(test_df_new, self.test_df1)

        # test start_w_firstval
        test_df_new = mDfS.impute_and_reduce_df(
            self.filldates_df, group_by_cols, mod_cols, impute='forward',
            moving_average=0, min_date='2021-01-04', max_date='2021-01-08',
            start_w_firstval=True)
        pd_testing.assert_frame_equal(test_df_new, self.test_df2)
        test_df_new = mDfS.impute_and_reduce_df(
            self.filldates_df, group_by_cols, mod_cols, impute='forward',
            moving_average=0, min_date='2021-01-04', max_date='2021-01-08',
            start_w_firstval=False)
        pd_testing.assert_frame_equal(test_df_new, self.test_df3)


if __name__ == '__main__':
    unittest.main()
