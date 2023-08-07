#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

from memilio.surrogatemodel.ode_secir_simple import (data_generation, model,
                                                     network_architectures)
import os
import unittest

import numpy as np
import logging

# suppress all autograph warnings from tensorflow

logging.getLogger("tensorflow").setLevel(logging.ERROR)


class TestSurrogatemodelOdeSecirSimple(fake_filesystem_unittest.TestCase):

    path = '/home/'

    def setUp(self):
        self.setUpPyfakefs()

    def test_simulation_run(self):

        days_1 = 10
        days_2 = 30
        days_3 = 50

        simulation_1 = data_generation.run_secir_simulation(days_1)
        simulation_2 = data_generation.run_secir_simulation(days_2)
        simulation_3 = data_generation.run_secir_simulation(days_3)

        self.assertEqual(len(simulation_1), days_1+1)
        self.assertEqual(len(simulation_2), days_2+1)
        self.assertEqual(len(simulation_3), days_3+1)

    def test_data_generation_runs(self):

        input_width_1 = 1
        input_width_2 = 3

        label_width_1 = 1
        label_width_2 = 10

        num_runs_1 = 1
        num_runs_2 = 5

        data_1 = data_generation.generate_data(
            num_runs_1, self.path, input_width_1, label_width_1,
            save_data=False)
        self.assertEqual(len(data_1['inputs']), num_runs_1)
        self.assertEqual(len(data_1['inputs'][0]), input_width_1)
        self.assertEqual(len(data_1['inputs'][0][0]), 8)
        self.assertEqual(len(data_1['labels']), num_runs_1)
        self.assertEqual(len(data_1['labels'][0]), label_width_1)
        self.assertEqual(len(data_1['labels'][0][0]), 8)

        data_2 = data_generation.generate_data(
            num_runs_2, self.path, input_width_2, label_width_2,
            save_data=False)
        self.assertEqual(len(data_2['inputs']), num_runs_2)
        self.assertEqual(len(data_2['inputs'][0]), input_width_2)
        self.assertEqual(len(data_2['inputs'][0][0]), 8)
        self.assertEqual(len(data_2['labels']), num_runs_2)
        self.assertEqual(len(data_2['labels'][0]), label_width_2)
        self.assertEqual(len(data_2['labels'][0][0]), 8)

    def test_data_generation_save(self):

        input_width = 2
        label_width = 3
        num_runs = 1

        data_generation.generate_data(num_runs, self.path, input_width,
                                      label_width)
        self.assertEqual(len(os.listdir(self.path)), 1)

        self.assertEqual(os.listdir(self.path),
                         ['data_secir_simple.pickle'])

    def test_structures_networks(self):

        model_mlp_multi_input_single_output = network_architectures.mlp_multi_input_single_output()
        self.assertEqual(len(model_mlp_multi_input_single_output.layers), 5)
        input_zero = np.zeros((1, 5, 8))
        output_zeros = model_mlp_multi_input_single_output(input_zero)
        self.assertEqual(output_zeros.shape[0], 1)
        self.assertEqual(output_zeros.shape[1], 1)
        self.assertEqual(output_zeros.shape[2], 8)

        model_lstm_single = network_architectures.lstm_network_multi_input_single_output()
        self.assertEqual(len(model_lstm_single.layers), 2)
        input_zero = np.zeros((1, 1, 8))
        output_zeros = model_lstm_single(input_zero)
        self.assertEqual(output_zeros.shape[0], 1)
        self.assertEqual(output_zeros.shape[1], 1)
        self.assertEqual(output_zeros.shape[2], 8)

        label_width = 5
        model_cnn = network_architectures.cnn_multi_input_multi_output(
            label_width)
        self.assertEqual(len(model_cnn.layers), 4)
        input_zero = np.zeros((1, 5, 8))
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
        self.assertEqual(output_zeros.shape[2], 8)

    def test_network_fit(self):

        max_epochs = 5

        # models with single output
        model_mlp_multi_input_single_output = network_architectures.mlp_multi_input_single_output()
        model_lstm_multi_input_single_output = network_architectures.lstm_network_multi_input_single_output()

        # no existing dataset
        with self.assertRaises(FileNotFoundError) as error:
            model.network_fit(
                self.path, model=model_mlp_multi_input_single_output,
                max_epochs=max_epochs)
        error_message = "[Errno 2] No such file or directory in the fake filesystem: '" + \
            self.path + "data_secir_simple.pickle'"
        self.assertEqual(str(error.exception), error_message)

        # generate dataset with single output
        input_width = 5
        label_width = 1
        num_runs = 15
        data_generation.generate_data(num_runs, self.path, input_width,
                                      label_width)

        # test different network_architectures
        mlp_output = model.network_fit(
            self.path, model=model_mlp_multi_input_single_output,
            max_epochs=max_epochs, plot=False)
        self.assertEqual(
            len(mlp_output.history['val_loss']), max_epochs)
        lstm_single_output = model.network_fit(
            self.path, model=model_lstm_multi_input_single_output,
            max_epochs=max_epochs, plot=False)
        self.assertEqual(
            len(lstm_single_output.history['val_loss']), max_epochs)

        # generate dataset with multiple output
        input_width = 5
        label_width = 10
        num_runs = 25
        data_generation.generate_data(num_runs, self.path, input_width,
                                      label_width)

        # models with multiple outputs
        model_cnn = network_architectures.cnn_multi_input_multi_output(
            label_width)
        model_lstm = network_architectures.lstm_multi_input_multi_output(
            label_width)

        cnn_output = model.network_fit(
            self.path, model=model_cnn, max_epochs=max_epochs, plot=False)
        self.assertEqual(
            cnn_output.model.output_shape[1], label_width)
        self.assertEqual(
            len(cnn_output.history['val_loss']), max_epochs)

        lstm_output = model.network_fit(
            self.path, model=model_lstm, max_epochs=max_epochs, plot=False)
        self.assertEqual(
            lstm_output.model.output_shape[1], label_width)
        self.assertEqual(
            len(lstm_output.history['val_loss']), max_epochs)


if __name__ == '__main__':
    unittest.main()
