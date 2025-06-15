#############################################################################
# Copyright (C) 2020-2025 MEmilio
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
        mock_load.return_value = {
            'Time': np.arange(48), 'Total': np.ones((48, 7))}
        mock_df.return_value.rolling.return_value.sum.return_value.to_numpy.return_value = np.ones(
            (48, 1))
        mock_matplotlib.colormaps.get_cmap.return_value.colors = [(1, 0, 0)]*7

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
        test_data = {
            'Time': np.arange(10),
            'Total': np.ones((10, 8)),
            'Group1': np.ones((10, 8)),
            'Group2': np.ones((10, 8)),
            'Group3': np.ones((10, 8)),
            'Group4': np.ones((10, 8)),
            'Group5': np.ones((10, 8)),
            'Group6': np.ones((10, 8))
        }
        mock_load.side_effect = [test_data, test_data, test_data]

        abm.plot_infection_states_results('dummy_path')

        # Verify functions are called with correct arguments
        self.assertEqual(mock_load.call_count, 3,
                         "load_h5_results should be called 3 times (p25, p50, p75)")

        # Check that load_h5_results was called with correct percentiles
        expected_calls = [
            unittest.mock.call('dummy_path', 'p25'),
            unittest.mock.call('dummy_path', 'p50'),
            unittest.mock.call('dummy_path', 'p75')
        ]
        mock_load.assert_has_calls(expected_calls, any_order=True)

        # Verify plotting functions were called with the loaded data
        mock_states.assert_called_once()
        mock_indiv.assert_called_once()

        # Check that plot_infection_states was called with correct data structure
        states_call_args = mock_states.call_args
        self.assertIsNotNone(states_call_args)
        self.assertEqual(len(states_call_args[0]), 4)  # x, y50, y25, y75

        # Check that plot_infection_states_individual was called with correct data structure
        indiv_call_args = mock_indiv.call_args
        self.assertIsNotNone(indiv_call_args)
        # x, p50_bs, p25_bs, p75_bs
        self.assertEqual(len(indiv_call_args[0]), 4)

    @patch('memilio.plot.plotAbmInfectionStates.matplotlib')
    def test_plot_infection_states(self, mock_matplotlib):
        x = np.arange(10)
        y50 = np.ones((10, 8))
        y25 = np.zeros((10, 8))
        y75 = np.ones((10, 8))*2
        y05 = np.ones((10, 8))*-1
        y95 = np.ones((10, 8))*3
        mock_matplotlib.colormaps.get_cmap.return_value.colors = [(1, 0, 0)]*8

        # Patch plt.gca().plot and fill_between
        with patch.object(abm.plt, 'gca') as mock_gca:
            mock_ax = MagicMock()
            mock_gca.return_value = mock_ax

            abm.plot_infection_states(
                x, y50, y25, y75,
                start_date='2021-03-01',
                colormap='Set1',
                xtick_step=2,
                y05=y05,
                y95=y95,
                show_90=True
            )

            # Verify plot was called with correct data
            self.assertTrue(mock_ax.plot.called, "Plot should be called")
            plot_calls = mock_ax.plot.call_args_list
            self.assertEqual(len(plot_calls), 8,
                             "Should plot 8 infection states")

            # Verify fill_between was called for confidence intervals
            self.assertTrue(mock_ax.fill_between.called,
                            "fill_between should be called for confidence intervals")
            fill_calls = mock_ax.fill_between.call_args_list
            # Should have calls for both 50% and 90% confidence intervals if show_90=True
            self.assertGreater(len(fill_calls), 0)

            # Verify axis formatting with correct parameters
            mock_ax.set_xticks.assert_called_once()
            xticks_call = mock_ax.set_xticks.call_args[0][0]
            expected_ticks = np.arange(0, len(x), 2)  # xtick_step=2
            np.testing.assert_array_equal(xticks_call, expected_ticks)

            # Verify xticklabels are set correctly
            mock_ax.set_xticklabels.assert_called_once()
            xticklabels_call = mock_ax.set_xticklabels.call_args[0][0]
            self.assertEqual(len(xticklabels_call), len(expected_ticks))

            # Verify that start_date is used in label formatting
            if len(xticklabels_call) > 0:
                # Labels should contain date strings when start_date is provided
                self.assertTrue(all('2021' in str(label)
                                for label in xticklabels_call))

    @patch('memilio.plot.plotAbmInfectionStates.matplotlib')
    def test_plot_infection_states_individual(self, mock_matplotlib):
        x = np.arange(10)
        group_data = np.ones((10, 8))
        groups = ['Group1', 'Group2', 'Group3',
                  'Group4', 'Group5', 'Group6', 'Total']
        p50_bs = {g: group_data for g in groups}
        p25_bs = {g: group_data for g in groups}
        p75_bs = {g: group_data for g in groups}
        p05_bs = {g: group_data*-1 for g in groups}
        p95_bs = {g: group_data*3 for g in groups}
        mock_matplotlib.colormaps.get_cmap.return_value.colors = [(1, 0, 0)]*8

        # Patch plt.subplots to return a grid of MagicMock axes
        with patch.object(abm.plt, 'subplots') as mock_subplots:
            fig_mock = MagicMock()
            ax_mock = MagicMock()  # Let the function determine the actual structure
            mock_subplots.return_value = (fig_mock, ax_mock)

            abm.plot_infection_states_by_age_group(
                x, p50_bs, p25_bs, p75_bs,
                colormap='Set1',
                p05_bs=p05_bs,
                p95_bs=p95_bs,
                show90=True
            )

            # Verify subplot was called (without assuming specific dimensions)
            mock_subplots.assert_called_once()
            subplot_call = mock_subplots.call_args

            # Verify that subplots was called with reasonable dimensions
            if subplot_call and len(subplot_call[0]) >= 2:
                rows, cols = subplot_call[0][:2]
                self.assertGreater(rows, 0, "Should have at least one row")
                self.assertGreater(cols, 0, "Should have at least one column")

                # The dimensions should be related to our data structure
                total_subplots = rows * cols
                expected_min_subplots = len(groups) * group_data.shape[1]
                self.assertGreaterEqual(total_subplots, expected_min_subplots,
                                        f"Should have at least {expected_min_subplots} subplots for {len(groups)} groups and {group_data.shape[1]} infection states")

        # Verify figure title is set
        fig_mock.suptitle.assert_called_once()

    def test__format_x_axis(self):
        test_x = np.arange(10)
        test_start_date = '2021-03-01'
        test_xtick_step = 2

        with patch('memilio.plot.plotAbmInfectionStates.plt') as mock_plt:
            mock_ax = MagicMock()
            mock_plt.gca.return_value = mock_ax

            abm._format_x_axis(test_x, test_start_date, test_xtick_step)

            # Verify that gca was called to get current axis
            mock_plt.gca.assert_called_once()

            # Verify that gcf was called to get current figure
            mock_plt.gcf.assert_called_once()

            # Verify axis formatting methods were called
            self.assertTrue(mock_ax.set_xticks.called,
                            "set_xticks should be called")
            self.assertTrue(mock_ax.set_xticklabels.called,
                            "set_xticklabels should be called")

            # Verify correct tick positions
            xticks_call = mock_ax.set_xticks.call_args
            if xticks_call and xticks_call[0]:
                tick_positions = xticks_call[0][0]
                expected_positions = np.arange(0, len(test_x), test_xtick_step)
                np.testing.assert_array_equal(
                    tick_positions, expected_positions)

            # Verify tick labels are date strings
            xticklabels_call = mock_ax.set_xticklabels.call_args
            if xticklabels_call and xticklabels_call[0]:
                tick_labels = xticklabels_call[0][0]
                self.assertIsInstance(tick_labels, (list, np.ndarray))
                if len(tick_labels) > 0:
                    # Should contain date information when start_date is provided
                    self.assertTrue(any('2021' in str(label)
                                    for label in tick_labels))


if __name__ == '__main__':
    unittest.main()
