#############################################################################
# Copyright (C) 2020-2025 MEmilio
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
from pyfakefs import fake_filesystem_unittest

from memilio.surrogatemodel.ode_secir_groups import (data_generation, model,
                                                     network_architectures)
from unittest.mock import patch
import os
import unittest

import numpy as np
import logging

# suppress all autograph warnings from tensorflow

logging.getLogger("tensorflow").setLevel(logging.ERROR)


class TestSurrogatemodelOdeSecirGroups(fake_filesystem_unittest.TestCase):
    """ """

    path = '/home/'

    def setUp(self):
        """ """
        self.setUpPyfakefs()

    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getMinimumMatrix',
           return_value=0.1 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    def test_simulation_run(self, mock_baseline, mock_minimum):
        """

        :param mock_baseline: 
        :param mock_minimum: 

        """
        days_1 = 10
        days_2 = 30
        days_3 = 50

        damping_date = 5
        population = [5256.0, 10551, 32368.5,
                      43637.833333333336, 22874.066666666666, 8473.6]

        simulation_1 = data_generation.run_secir_groups_simulation(
            days_1, damping_date, population)
        simulation_2 = data_generation.run_secir_groups_simulation(
            days_2, damping_date, population)
        simulation_3 = data_generation.run_secir_groups_simulation(
            days_3, damping_date, population)

        # result length
        self.assertEqual(len(simulation_1[0]), days_1+1)
        self.assertEqual(len(simulation_2[0]), days_2+1)
        self.assertEqual(len(simulation_3[0]), days_3+1)

        # damping
        self.assertEqual(
            simulation_1[1].size,
            len(population) * len(population))
        self.assertEqual(
            simulation_2[1].size,
            len(population) * len(population))
        self.assertEqual(
            simulation_3[1].size,
            len(population) * len(population))

        self.assertLessEqual(simulation_1[1].max(), 1)
        self.assertGreaterEqual(simulation_1[1].max(), 0)

        self.assertLessEqual(simulation_2[1].max(), 1)
        self.assertGreaterEqual(simulation_2[1].max(), 0)

        self.assertLessEqual(simulation_3[1].max(), 1)
        self.assertGreaterEqual(simulation_3[1].max(), 0)

    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getMinimumMatrix',
           return_value=0.1 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.get_population',
           return_value=[[5256.0, 10551, 32368.5,
                         43637.833333333336, 22874.066666666666, 8473.6]])
    def test_data_generation_runs(
            self, mock_population, mock_baseline, mock_minimum):
        """

        :param mock_population: 
        :param mock_baseline: 
        :param mock_minimum: 

        """

        input_width_1 = 1
        input_width_2 = 5

        label_width_1 = 1
        label_width_2 = 10

        num_runs_1 = 1
        num_runs_2 = 5

        data_1 = data_generation.generate_data(
            num_runs_1, self.path, "", input_width_1, label_width_1,
            save_data=False)
        self.assertEqual(len(data_1['inputs']), num_runs_1)
        self.assertEqual(len(data_1['inputs'][0]), input_width_1)
        self.assertEqual(len(data_1['inputs'][0][0]), 48)
        self.assertEqual(len(data_1['labels']), num_runs_1)
        self.assertEqual(len(data_1['labels'][0]), label_width_1)
        self.assertEqual(len(data_1['labels'][0][0]), 48)

        data_2 = data_generation.generate_data(
            num_runs_2, self.path, "", input_width_2, label_width_2,
            save_data=False)
        self.assertEqual(len(data_2['inputs']), num_runs_2)
        self.assertEqual(len(data_2['inputs'][0]), input_width_2)
        self.assertEqual(len(data_2['inputs'][0][0]), 48)
        self.assertEqual(len(data_2['labels']), num_runs_2)
        self.assertEqual(len(data_2['labels'][0]), label_width_2)
        self.assertEqual(len(data_2['labels'][0][0]), 48)

    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getMinimumMatrix',
           return_value=0.1 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.get_population',
           return_value=[[5256.0, 10551, 32368.5,
                         43637.833333333336, 22874.066666666666, 8473.6]])
    def test_data_generation_save(
            self, mock_population, mock_baseline, mock_minimum):
        """

        :param mock_population: 
        :param mock_baseline: 
        :param mock_minimum: 

        """

        input_width = 2
        label_width = 3
        num_runs = 1

        data_generation.generate_data(num_runs, self.path, "", input_width,
                                      label_width)
        self.assertEqual(len(os.listdir(self.path)), 1)

        self.assertEqual(os.listdir(self.path),
                         ['data_secir_groups.pickle'])

    def test_structures_networks(self):
        """ """

        model_mlp_multi_input_single_output = network_architectures.mlp_multi_input_single_output()
        self.assertEqual(len(model_mlp_multi_input_single_output.layers), 5)
        input_zero = np.zeros((1, 5, 8))
        output_zeros = model_mlp_multi_input_single_output(input_zero)
        self.assertEqual(output_zeros.shape[0], 1)
        self.assertEqual(output_zeros.shape[1], 1)
        self.assertEqual(output_zeros.shape[2], 48)

        label_width = 30
        model_mlp_multi_input_multi_output = network_architectures.mlp_multi_input_multi_output(
            label_width)
        self.assertEqual(len(model_mlp_multi_input_multi_output.layers), 5)
        input_zero = np.zeros((1, 1, 8))
        output_zeros = model_mlp_multi_input_multi_output(input_zero)
        self.assertEqual(output_zeros.shape[0], 1)
        self.assertEqual(output_zeros.shape[1], 30)
        self.assertEqual(output_zeros.shape[2], 48)

        label_width = 5
        model_cnn = network_architectures.cnn_multi_input_multi_output(
            label_width, num_age_groups=1)
        self.assertEqual(len(model_cnn.layers), 7)
        input_zero = np.zeros((1, label_width, 41))
        output_zeros = model_cnn(input_zero)
        self.assertEqual(output_zeros.shape[0], 1)
        self.assertEqual(output_zeros.shape[1], 5)
        self.assertEqual(output_zeros.shape[2], 8)

        model_lstm_multi = network_architectures.lstm_multi_input_multi_output(
            label_width)
        self.assertEqual(len(model_lstm_multi.layers), 3)
        input_zero = np.zeros((1, 5, 8))
        output_zeros = model_lstm_multi(input_zero)
        self.assertEqual(output_zeros.shape[0], 1)
        self.assertEqual(output_zeros.shape[1], 5)
        self.assertEqual(output_zeros.shape[2], 48)

    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getMinimumMatrix',
           return_value=0.1 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.get_population',
           return_value=[[5256.0, 10551, 32368.5,
                         43637.833333333336, 22874.066666666666, 8473.6]])
    def test_network_fit(self, mock_population, mock_baseline, mock_minimum):
        """

        :param mock_population: 
        :param mock_baseline: 
        :param mock_minimum: 

        """
        max_epochs = 5
        # models with single output
        model_mlp_multi_input_single_output = network_architectures.mlp_multi_input_single_output()

        # no existing dataset
        with self.assertRaises(FileNotFoundError) as error:
            model.network_fit(
                self.path, model=model_mlp_multi_input_single_output,
                modeltype='classic', max_epochs=max_epochs)
        error_message = "[Errno 2] No such file or directory in the fake filesystem: '" + \
            self.path + "data_secir_groups.pickle'"
        self.assertEqual(str(error.exception), error_message)

        # generate dataset with single output
        input_width = 5
        label_width = 1
        num_runs = 5
        data_generation.generate_data(num_runs, self.path, "", input_width,
                                      label_width)

        # test different network_architectures
        mlp_output = model.network_fit(
            self.path, model=model_mlp_multi_input_single_output,
            modeltype='classic', max_epochs=max_epochs, plot=False)

        self.assertEqual(
            len(mlp_output.history['val_loss']), max_epochs)

        # generate dataset with multiple output
        input_width = 5
        label_width = 10
        num_runs = 5
        data_generation.generate_data(num_runs, self.path, "", input_width,
                                      label_width)

        # models with multiple outputs
        model_mlp_multi_input_multi_output = model.network_fit(
            self.path,
            model=network_architectures.mlp_multi_input_multi_output(
                label_width),
            modeltype='classic', max_epochs=max_epochs, plot=False)
        self.assertEqual(
            model_mlp_multi_input_multi_output.model.output_shape[1],
            label_width)
        self.assertEqual(
            len(model_mlp_multi_input_multi_output.history['val_loss']),
            max_epochs)

        model_lstm_multi_output = model.network_fit(
            self.path,
            model=network_architectures.lstm_multi_input_multi_output(
                label_width),
            modeltype='timeseries', max_epochs=max_epochs, plot=False)
        self.assertEqual(
            model_lstm_multi_output.model.output_shape[1], label_width)
        self.assertEqual(
            len(model_lstm_multi_output.history['val_loss']), max_epochs)

        cnn_output = model.network_fit(
            self.path,
            model=network_architectures.cnn_multi_input_multi_output(
                label_width),
            modeltype='timeseries', max_epochs=max_epochs, plot=False)
        self.assertEqual(
            cnn_output.model.output_shape[1], label_width)
        self.assertEqual(
            len(cnn_output.history['val_loss']), max_epochs)


if __name__ == '__main__':
    unittest.main()
