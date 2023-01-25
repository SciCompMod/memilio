#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Patrick Lenz
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
import os
import unittest
from unittest.mock import MagicMock, patch

from pyfakefs import fake_filesystem_unittest

from memilio.epidata import customPlot


class Test_customPlot(fake_filesystem_unittest.TestCase):

    path = '/home/figures/'

    def setUp(self):
        self.setUpPyfakefs()

    @patch('memilio.epidata.customPlot.plt')
    def test_plot_list(self, mock_plt):

        mock_plt.subplots.return_value = (MagicMock(), MagicMock())

        xvals = [i for i in range(100)]
        yvals = [xvals for i in range(3)]

        customPlot.plot_multiple_series(xvals, yvals, ['yvals' + str(i)
                                                       for i in range(len(yvals))], title='Test', xlabel='Date',
                                        ylabel='Number of Test', xticks_idx=[0, 17, 47, 66, 99], linewidth=2,
                                        loc_legend='upper right', fig_size=(9, 6), fig_name='Test', path_rel=self.path, dpi=50,
                                        outercolor=[0.3, 0.5, 0.3], innercolor=[1, 1, 1])

        mock_plt.subplots.assert_called_once_with(
            figsize=(9, 6), facecolor=[0.3, 0.5, 0.3])
        mock_plt.savefig.assert_called_once_with(
            self.path + 'Test' + '.png', bbox_inches='tight', dpi=50)


if __name__ == '__main__':
    unittest.main()
