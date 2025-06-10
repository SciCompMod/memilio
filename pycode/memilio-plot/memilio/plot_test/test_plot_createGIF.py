#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Henrik Zunker
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
import pandas as pd
import numpy as np
from memilio.plot import createGIF


class TestCreateGif(unittest.TestCase):
    """ """

    @patch('memilio.plot.createGIF.pm.extract_data')
    @patch('pandas.read_json')
    @patch('memilio.plot.createGIF.pm.scale_dataframe_relative')
    @patch('memilio.plot.createGIF.pm.plot_map')
    def test_create_plot_map(self, mock_plot_map, mock_scale_dataframe_relative, mock_read_json, mock_extract_data):
        """

        :param mock_plot_map: 
        :param mock_scale_dataframe_relative: 
        :param mock_read_json: 
        :param mock_extract_data: 

        """
        # Mock the return values of the functions
        mock_extract_data.return_value = pd.DataFrame(
            {'Region': [0, 1], 'Value': [1, 2]})
        mock_read_json.return_value = pd.DataFrame(
            {'Region': [0, 1], 'Population': [100, 200]})
        mock_scale_dataframe_relative.return_value = pd.DataFrame(
            {'Region': [0, 1], 'Value': [0.01, 0.01]})

        # Call the function with test parameters
        createGIF.create_plot_map(day=1, filename='test', files_input={'file1': 'path1'}, output_path='output', compartments=[
                                  'compartment1'], file_format='h5', relative=True, age_groups={0: '0-4', 1: '5-14', 2: '15-34', 3: '35-59', 4: '60-79', 5: '80+'})

        # Assert that the mocked functions were called with the correct parameters
        mock_extract_data.assert_called()
        mock_read_json.assert_called()
        mock_scale_dataframe_relative.assert_called()

        assert mock_plot_map.called

        # Get the arguments passed to the mock and assert that they are correct
        args, kwargs = mock_plot_map.call_args

        pd.testing.assert_frame_equal(args[0], pd.DataFrame(
            {'Region': [0, 1], 'Value 0': [0.01, 0.01]}))
        np.testing.assert_array_equal(
            kwargs['scale_colors'], np.array([0.01, 0.01]))
        self.assertEqual(kwargs['legend'], ['', ''])
        self.assertEqual(kwargs['title'], 'Synthetic data (relative) day  1')
        self.assertEqual(kwargs['plot_colorbar'], True)
        self.assertEqual(kwargs['output_path'], 'output')
        self.assertEqual(kwargs['fig_name'], 'test')
        self.assertEqual(kwargs['dpi'], 300)
        self.assertEqual(kwargs['outercolor'], [
                         205 / 255, 238 / 255, 251 / 255])

    @patch('memilio.plot.createGIF.pm.extract_time_steps')
    @patch('memilio.plot.createGIF.create_plot_map')
    @patch('tempfile.TemporaryDirectory')
    @patch('imageio.v2.imread')
    @patch('imageio.mimsave')
    def test_create_gif_map_plot(self, mock_mimsave, mock_imread, mock_tempdir, mock_create_plot_map, mock_extract_time_steps):
        """

        :param mock_mimsave: 
        :param mock_imread: 
        :param mock_tempdir: 
        :param mock_create_plot_map: 
        :param mock_extract_time_steps: 

        """
        # Mock the return values of the functions
        mock_extract_time_steps.return_value = 10
        mock_tempdir.return_value.__enter__.return_value = 'tempdir'
        mock_imread.return_value = 'image'

        # Call the function with test parameters
        createGIF.create_gif_map_plot(input_data='input', output_dir='output', compartments=[
                                      'compartment1'], filename='test')

        # Assert that the mocked functions were called with the correct parameters
        mock_extract_time_steps.assert_called_with(
            'input/Results', file_format='h5')
        mock_create_plot_map.assert_called()
        mock_tempdir.assert_called()
        mock_imread.assert_called_with('tempdir/test.png')
        mock_mimsave.assert_called_with(
            'output/test.gif', ['image']*10, duration=0.2, loop=0)

        # Get the arguments passed to the mock and assert that they are correct
        args, kwargs = mock_create_plot_map.call_args
        self.assertEqual(args[0], 9)
        self.assertEqual(args[1], 'test')
        self.assertEqual(args[2], {'Data set': 'input/Results'})
        self.assertEqual(args[3], 'tempdir')
        self.assertEqual(args[4], ['compartment1'])
        self.assertEqual(args[5], 'h5')
        self.assertEqual(args[6], True)
        self.assertEqual(args[7], {
                         0: '0-4', 1: '5-14', 2: '15-34', 3: '35-59', 4: '60-79', 5: '80+'})


if __name__ == '__main__':
    unittest.main()
