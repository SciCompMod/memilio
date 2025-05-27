#############################################################################
# Copyright (C) 2020-2024 MEmilio
#
# Authors: Sascha Korf
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
from unittest.mock import patch, MagicMock
import numpy as np
import pandas as pd

import memilio.plot.plotAbmInfectionStates as abm

class TestPlotAbmInfectionStates(unittest.TestCase):

    @patch('memilio.plot.plotAbmInfectionStates.h5py.File')
    def test_load_h5_results(self, mock_h5file):
        mock_group = {'Time': np.arange(10), 'Total': np.ones((10, 8))}
        mock_h5file().__enter__().get.return_value = {'0': mock_group}
        mock_h5file().__enter__().items.return_value = [('0', mock_group)]
        mock_h5file().__enter__().__getitem__.return_value = mock_group
        mock_h5file().__enter__().items.return_value = [('Time', np.arange(10)), ('Total', np.ones((10, 8)))]
        with patch('memilio.plot.plotAbmInfectionStates.h5py.File', mock_h5file):
            result = abm.load_h5_results('dummy_path', 'p50')
            assert 'Time' in result
            assert 'Total' in result
            np.testing.assert_array_equal(result['Time'], np.arange(10))
            np.testing.assert_array_equal(result['Total'], np.ones((10, 8)))

    @patch('memilio.plot.plotAbmInfectionStates.load_h5_results')
    @patch('memilio.plot.plotAbmInfectionStates.matplotlib')
    @patch('memilio.plot.plotAbmInfectionStates.gaussian_filter1d', side_effect=lambda x, sigma, mode: x)
    @patch('memilio.plot.plotAbmInfectionStates.pd.DataFrame')
    def test_plot_infections_loc_types_average(self, mock_df, mock_gauss, mock_matplotlib, mock_load):
        mock_load.return_value = {'Time': np.arange(48), 'Total': np.ones((48, 7))}
        mock_df.return_value.rolling.return_value.sum.return_value.to_numpy.return_value = np.ones((48, 1))
        mock_matplotlib.colormaps.get_cmap.return_value.colors = [(1,0,0)]*7

        # Patch plt.gca().plot to a MagicMock
        with patch.object(abm.plt, 'gca') as mock_gca:
            mock_ax = MagicMock()
            mock_gca.return_value = mock_ax
            abm.plot_infections_loc_types_average('dummy_path')
            assert mock_ax.plot.called
            assert mock_ax.set_xticks.called
            assert mock_ax.set_xticklabels.called

    @patch('memilio.plot.plotAbmInfectionStates.load_h5_results')
    @patch('memilio.plot.plotAbmInfectionStates.plot_infection_states')
    @patch('memilio.plot.plotAbmInfectionStates.plot_infection_states_individual')
    def test_plot_infection_states_results(self, mock_indiv, mock_states, mock_load):
        mock_load.side_effect = [
            {'Time': np.arange(10), 'Total': np.ones((10, 8)), 'Group1': np.ones((10, 8)), 'Group2': np.ones((10, 8)), 'Group3': np.ones((10, 8)), 'Group4': np.ones((10, 8)), 'Group5': np.ones((10, 8)), 'Group6': np.ones((10, 8)), 'Total': np.ones((10, 8))},
            {'Time': np.arange(10), 'Total': np.ones((10, 8)), 'Group1': np.ones((10, 8)), 'Group2': np.ones((10, 8)), 'Group3': np.ones((10, 8)), 'Group4': np.ones((10, 8)), 'Group5': np.ones((10, 8)), 'Group6': np.ones((10, 8)), 'Total': np.ones((10, 8))},
            {'Time': np.arange(10), 'Total': np.ones((10, 8)), 'Group1': np.ones((10, 8)), 'Group2': np.ones((10, 8)), 'Group3': np.ones((10, 8)), 'Group4': np.ones((10, 8)), 'Group5': np.ones((10, 8)), 'Group6': np.ones((10, 8)), 'Total': np.ones((10, 8))}
        ]
        abm.plot_infection_states_results('dummy_path')
        assert mock_indiv.called
        assert mock_states.called

    @patch('memilio.plot.plotAbmInfectionStates.matplotlib')
    def test_plot_infection_states(self, mock_matplotlib):
        x = np.arange(10)
        y50 = np.ones((10, 8))
        y25 = np.zeros((10, 8))
        y75 = np.ones((10, 8))*2
        mock_matplotlib.colormaps.get_cmap.return_value.colors = [(1,0,0)]*8

        # Patch plt.gca().plot and fill_between
        with patch.object(abm.plt, 'gca') as mock_gca:
            mock_ax = MagicMock()
            mock_gca.return_value = mock_ax
            abm.plot_infection_states(x, y50, y25, y75)
            assert mock_ax.plot.called
            assert mock_ax.fill_between.called
            assert mock_ax.set_xticks.called
            assert mock_ax.set_xticklabels.called

    @patch('memilio.plot.plotAbmInfectionStates.matplotlib')
    def test_plot_infection_states_individual(self, mock_matplotlib):
        x = np.arange(10)
        group_data = np.ones((10, 8))
        p50_bs = {g: group_data for g in ['Group1', 'Group2', 'Group3', 'Group4', 'Group5', 'Group6', 'Total']}
        p25_bs = {g: group_data for g in ['Group1', 'Group2', 'Group3', 'Group4', 'Group5', 'Group6', 'Total']}
        p75_bs = {g: group_data for g in ['Group1', 'Group2', 'Group3', 'Group4', 'Group5', 'Group6', 'Total']}
        mock_matplotlib.colormaps.get_cmap.return_value.colors = [(1,0,0)]*8

        # Patch plt.subplots to return a grid of MagicMock axes (as np.array with dtype=object)
        with patch.object(abm.plt, 'subplots') as mock_subplots:
            fig_mock = MagicMock()
            ax_mock = np.empty((6, 7), dtype=object)
            for i in range(6):
                for j in range(7):
                    ax_mock[i, j] = MagicMock()
            mock_subplots.return_value = (fig_mock, ax_mock)
            abm.plot_infection_states_individual(x, p50_bs, p25_bs, p75_bs)
            # Check that at least one ax's plot was called
            assert any(ax_mock[i, j].plot.called for i in range(6) for j in range(7))
            assert fig_mock.suptitle.called

    def test__format_x_axis(self):
        with patch('memilio.plot.plotAbmInfectionStates.plt') as mock_plt:
            abm._format_x_axis(np.arange(10), '2021-03-01', 2)
            assert mock_plt.gca.called
            assert mock_plt.gcf.called

if __name__ == '__main__':
    unittest.main()