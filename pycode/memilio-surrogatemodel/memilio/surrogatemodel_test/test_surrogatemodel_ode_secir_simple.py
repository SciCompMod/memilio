#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Manuel Heger
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
from memilio.surrogatemodel.ode_secir_simple import (data_generation,
                                                     network_architectures)
from memilio.surrogatemodel.ode_secir_simple import model as md
import os
import unittest
from unittest.mock import patch

import numpy as np
import tensorflow as tf
import logging

# suppress all autograph warnings from tensorflow

logging.getLogger("tensorflow").setLevel(logging.ERROR)


class TestSurrogatemodelOdeSecirSimple(fake_filesystem_unittest.TestCase):
    """ """

    path = '/home/'

    def setUp(self):
        """ """
        self.setUpPyfakefs()

    def test_simulation_run(self):
        """ Checking if simulated data for one sample using run_secir_simple_simulation has the desired length"""

        days_1 = 10
        days_2 = 30
        days_3 = 50

        simulation_1 = data_generation.run_secir_simple_simulation(days_1)
        simulation_2 = data_generation.run_secir_simple_simulation(days_2)
        simulation_3 = data_generation.run_secir_simple_simulation(days_3)

        self.assertEqual(len(simulation_1), days_1+1)
        self.assertEqual(len(simulation_2), days_2+1)
        self.assertEqual(len(simulation_3), days_3+1)

    def test_data_generation_runs(self):
        """ Testing if data_generation produces data with desired dimensions """

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
        """ Checking if saving results as pickle-file is possible """

        input_width = 2
        label_width = 3
        num_runs = 1

        data_generation.generate_data(num_runs, self.path, input_width,
                                      label_width)
        self.assertEqual(len(os.listdir(self.path)), 1)

        self.assertEqual(os.listdir(self.path),
                         ['data_secir_simple_3days_1.pickle'])

    def test_mlp_multi_single(self):
        """" testing initialization of MLP architectures works correctly """
        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_single_output(
                0, 1, 1
            )
        error_message = "Output dimension must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_single_output(
                1, -1, 1
            )
        error_message = "Number of layers must be at least 0, here -1"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_single_output(
                1, 1, 0
            )
        error_message = "Number of neurons per layer must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        model = network_architectures.mlp_multi_input_single_output(3, 4, 84)
        # Layers: Flatten layer + hidden layers + flatten layer + output layer
        self.assertEqual(len(model.layers), 7)
        input_zero = np.zeros((1, 5, 8))
        output_zeros = model(input_zero)
        # Output has shape (1, num_outputs)
        self.assertEqual(output_zeros.shape[0], 1)
        self.assertEqual(output_zeros.shape[1], 3)

    def test_mlp_multi_multi(self):
        """" testing initialization of MLP architectures works correctly """
        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_multi_output(
                0, 12, 12, 12
            )
        error_message = "label width has to be a positive integer"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_multi_output(
                12, 0, 12, 12
            )
        error_message = "Output dimension must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_multi_output(
                12, 12, -1, 12
            )
        error_message = "Number of layers must be at least 0, here -1"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_multi_output(
                12, 12, 12, 0
            )
        error_message = "Number of neurons per layer must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        model = network_architectures.mlp_multi_input_multi_output(
            12, 3, 4, 84)
        # Layers: flatten layers + hidden layers + output layer + reshape layer
        self.assertEqual(len(model.layers), 7)
        input_zero = np.zeros((2, 5, 8))
        output_zeros = model(input_zero)
        # Output has shape (Number of samples, label width, number of outputs)
        self.assertEqual(output_zeros.shape[0], 2)
        self.assertEqual(output_zeros.shape[1], 12)
        self.assertEqual(output_zeros.shape[2], 3)

    def test_lstm_multi_single(self):
        """" testing initialization of LSTM architectures works correctly """
        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_network_multi_input_single_output(
                0, 12, 12, 12
            )
        error_message = "Output dimension must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_network_multi_input_single_output(
                12, 0, 12, 12
            )
        error_message = "Internal dimension must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_network_multi_input_single_output(
                12, 12, -1, 12
            )
        error_message = "Number of layers must be at least 0, here -1"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_network_multi_input_single_output(
                12, 12, 12, 0
            )
        error_message = "Number of neurons per layer must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)
        model = network_architectures.lstm_network_multi_input_single_output(
            num_outputs=5, num_hidden_layers=4)
        # model: LSTM-layer + hidden layers + output layer + reshape layer, here 1+4+1+1
        self.assertEqual(len(model.layers), 7)
        input_zero = np.zeros((1, 5, 3))
        output_zeros = model(input_zero)
        # Dimensionality output (1, Number of outputs), here (1,5)
        self.assertEqual(output_zeros.shape[0], 1)
        self.assertEqual(output_zeros.shape[1], 5)

    def test_lstm_multi_multi(self):
        """" testing initialization of LSTM architectures works correctly """
        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_multi_input_multi_output(
                0, 12, 12, 12, 12
            )
        error_message = "label width has to be a positive integer"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_multi_input_multi_output(
                12, 0, 12, 12, 12
            )
        error_message = "Output dimension must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_multi_input_multi_output(
                12, 12, 0, 12, 12
            )
        error_message = "Internal dimension must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_multi_input_multi_output(
                12, 12, 12, -1, 12
            )
        error_message = "Number of hidden layers must be at least 0, here -1"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_multi_input_multi_output(
                12, 12, 12, 12, 0
            )
        error_message = "Number of neurons per layer must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        model = network_architectures.lstm_multi_input_multi_output(
            3, 84, 4, 12)
        # model: LSTM-layer + hidden layers + output layer + reshape layer, here 1+12+1+1
        self.assertEqual(len(model.layers), 15)
        input_zero = np.zeros((1, 5, 8))
        output_zeros = model(input_zero)
        # Dimensionality output (Number of samples, label width, Number of outputs per day), here (1,3,84)
        self.assertEqual(output_zeros.shape[0], 1)
        self.assertEqual(output_zeros.shape[1], 3)
        self.assertEqual(output_zeros.shape[2], 84)

    def test_cnn_multi_multi(self):
        """" testing initialization of CNN architectures works correctly """
        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                label_width=0, conv_size=3, num_outputs=8, num_filters=12,
                num_hidden_layers=12, num_neurons_per_layer=12
            )
        error_message = "label width has to be a positive integer"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                label_width=2, conv_size=1, num_outputs=8, num_filters=12,
                num_hidden_layers=12, num_neurons_per_layer=12
            )
        error_message = "Size of the convolution kernel has to be larger than 1, here 1"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                label_width=3, conv_size=3, num_outputs=-23, num_filters=12,
                num_hidden_layers=12, num_neurons_per_layer=12
            )
        error_message = "Output dimension must be at least 1, here -23"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                label_width=3, conv_size=3, num_outputs=8, num_filters=-2,
                num_hidden_layers=12, num_neurons_per_layer=12
            )
        error_message = "Number of filters must be at least 1, here -2"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                label_width=3, conv_size=3, num_outputs=8, num_filters=2,
                num_hidden_layers=-23, num_neurons_per_layer=12
            )
        error_message = "Number of hidden layers must be at least 0, here -23"
        self.assertEqual(str(error.exception), error_message)

        model = network_architectures.cnn_multi_input_multi_output(
            label_width=3, conv_size=3, num_outputs=42, num_filters=2,
            num_hidden_layers=23, num_neurons_per_layer=12
        )
        # lambda layer + conv layer + hidden layers + output layer + reshape layer -> 1+1+23+1+1 = 27
        self.assertEqual(len(model.layers), 27)
        input_zero = np.zeros((1, 5, 8))
        output_zeros = model(input_zero)
        # Dimensionality output (Number of samples, label width, Number of outputs per day), here (1,3,42)
        self.assertEqual(output_zeros.shape[0], 1)
        self.assertEqual(output_zeros.shape[1], 3)
        self.assertEqual(output_zeros.shape[2], 42)

    def test_initialize_model(self):
        """" Testing wheter initialization via initialize_model works correctly """
        # Deleting any automatically added numbering like ("Dense_01" -> "Dense")
        def normalize_config(config):
            config.pop('name', None)
            for layer in config["layers"]:
                layer["config"].pop("name", None)
            return config

        label_width = 8
        num_outputs = 32
        hidden_layers = [2]
        neurons_in_hidden_layer = [32]
        activation_function = ['relu']
        models = ["LSTM", "Dense", "CNN", "MLP"]

        model_lstm = network_architectures.lstm_multi_input_multi_output(
            label_width=label_width, num_outputs=num_outputs, num_hidden_layers=2,
            num_neurons_per_layer=32, activation='relu'
        )
        model_cnn = network_architectures.cnn_multi_input_multi_output(
            label_width=label_width, num_outputs=num_outputs, num_hidden_layers=2,
            num_neurons_per_layer=32, activation='relu'
        )
        model_mlp = network_architectures.mlp_multi_input_multi_output(
            label_width=label_width, num_outputs=num_outputs, num_hidden_layers=2,
            num_neurons_per_layer=32, activation='relu'
        )

        model_parameters = [(label_width, num_outputs, layer, neuron_number, activation, modelname)
                            for layer in hidden_layers for neuron_number in neurons_in_hidden_layer
                            for activation in activation_function for modelname in models]
        for parameter in model_parameters:
            _, _, _, _, _, modelname = parameter
            if modelname == "MLP":
                with self.assertRaises(ValueError) as error:
                    model = md.initialize_model(parameter)
                error_message = "name_architecture must be one of 'Dense', 'LSTM' or 'CNN'"
                self.assertEqual(str(error.exception), error_message)
            else:
                model = md.initialize_model(parameter)
                if modelname == "Dense":
                    self.assertDictEqual(
                        normalize_config(model.get_config()),
                        normalize_config(model_mlp.get_config()))
                if modelname == "LSTM":
                    self.assertDictEqual(
                        normalize_config(model.get_config()),
                        normalize_config(model_lstm.get_config()))
                if modelname == "CNN":
                    self.assertDictEqual(
                        normalize_config(model.get_config()),
                        normalize_config(model_cnn.get_config()))

    def test_network_fit(self):
        """ Test for training of network """
        max_epochs = 1

        early_stop = 100
        loss = tf.keras.losses.MeanAbsolutePercentageError()
        optimizer = 'Adam'
        metric = [tf.keras.metrics.MeanAbsoluteError(),
                  tf.keras.metrics.MeanAbsolutePercentageError()]

        training_parameter = (
            early_stop, max_epochs, loss, optimizer, metric)

        # generate dataset with multiple outputs
        label_width = 10

        inputs = tf.ones([5, 5, 8])
        labels = tf.ones([5, 10, 8])

        model_cnn = network_architectures.cnn_multi_input_multi_output(
            label_width=label_width, num_filters=1, num_neurons_per_layer=1,
            num_hidden_layers=0)
        model_lstm = network_architectures.lstm_multi_input_multi_output(
            label_width=label_width, internal_dimension=1, num_neurons_per_layer=1,
            num_hidden_layers=0)

        cnn_output = md.network_fit(
            model=model_cnn, inputs=inputs, labels=labels,
            training_parameter=training_parameter, plot=False)
        self.assertEqual(
            cnn_output.model.output_shape[1], label_width)
        self.assertEqual(
            len(cnn_output.history['val_loss']), max_epochs)

        lstm_output = md.network_fit(
            model=model_lstm, inputs=inputs, labels=labels,
            training_parameter=training_parameter, plot=False)
        self.assertEqual(
            lstm_output.model.output_shape[1], label_width)
        self.assertEqual(
            len(lstm_output.history['val_loss']), max_epochs)


if __name__ == '__main__':
    unittest.main()
