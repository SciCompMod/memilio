import unittest
from unittest.mock import patch, MagicMock
import pandas as pd
import numpy as np
from memilio.plot import createGIF


class TestCreateGif(unittest.TestCase):

    @patch('memilio.plot.createGIF.pm.extract_data')
    @patch('memilio.plot.createGIF.gpd.get_population_data')
    @patch('memilio.plot.createGIF.pm.scale_dataframe_relative')
    @patch('memilio.plot.createGIF.pm.plot_map')
    def test_create_plot_map(self, mock_plot_map, mock_scale_dataframe_relative, mock_get_population_data, mock_extract_data):
        # Mock the return values of the functions
        mock_extract_data.return_value = pd.DataFrame(
            {'Region': ['A', 'B'], 'Value': [1, 2]})
        mock_get_population_data.return_value = pd.DataFrame(
            {'Region': ['A', 'B'], 'Population': [100, 200]})
        mock_scale_dataframe_relative.return_value = pd.DataFrame(
            {'Region': ['A', 'B'], 'Value': [0.01, 0.01]})

        # Call the function with test parameters
        createGIF.create_plot_map(day=1, filename='test', files_input={'file1': 'path1'}, output_path='output', compartments=[
                                  'compartment1'], file_format='h5', relative=True, age_groups={0: '0-4', 1: '5-14', 2: '15-34', 3: '35-59', 4: '60-79', 5: '80+'})

        # Assert that the mocked functions were called with the correct parameters
        mock_extract_data.assert_called()
        mock_get_population_data.assert_called()
        mock_scale_dataframe_relative.assert_called()
        mock_plot_map.assert_called_with(pd.DataFrame({'Region': ['A', 'B'], 'Value': [0.01, 0.01]}), scale_colors=np.array([0.01, 0.01]), legend=[
                                         '', ''], title='Synthetic data (relative) day  1', plot_colorbar=True, output_path='output', fig_name='test', dpi=300, outercolor=[205 / 255, 238 / 255, 251 / 255])

    @patch('memilio.plot.createGIF.pm.extract_time_steps')
    @patch('memilio.plot.createGIF.create_plot_map')
    @patch('memilio.plot.createGIF.imageio.v2.imread')
    @patch('memilio.plot.createGIF.imageio.mimsave')
    @patch('memilio.plot.createGIF.tempfile.TemporaryDirectory')
    def test_createGIF(self, mock_temp_dir, mock_mimsave, mock_imread, mock_create_plot_map, mock_extract_time_steps):
        # Mock the return values of the functions
        mock_extract_time_steps.return_value = 5
        mock_temp_dir.return_value.__enter__.return_value = 'temp_dir'
        mock_imread.return_value = 'image'

        # Call the function with test parameters
        createGIF.createGIF(files_input={'file1': 'path1'}, output_dir='output', filename='test', compartments=[
                            'compartment1'], file_format='h5', relative=True, age_groups={0: '0-4', 1: '5-14', 2: '15-34', 3: '35-59', 4: '60-79', 5: '80+'})

        # Assert that the mocked functions were called with the correct parameters
        mock_extract_time_steps.assert_called_with('path1', file_format='h5')
        mock_create_plot_map.assert_called()
        mock_imread.assert_called_with('temp_dir/filename.png')
        mock_mimsave.assert_called_with(
            'output/test.gif', ['image', 'image', 'image', 'image', 'image'], duration=10, loop=0)


if __name__ == '__main__':
    unittest.main()
