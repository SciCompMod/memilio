#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Manuel Heger, Henrik Zunker
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
                                                     network_architectures, dampings)
from unittest.mock import patch
import os
import unittest

import numpy as np
import tensorflow as tf
import logging

# suppress all autograph warnings from tensorflow

logging.getLogger("tensorflow").setLevel(logging.ERROR)


class TestSurrogatemodelOdeSecirGroups(fake_filesystem_unittest.TestCase):
    """ """

    path = '/home/'

    def setUp(self):
        """ """
        self.setUpPyfakefs()

    def test_dampings_active(self):
        # length of generated sequence
        factors1 = dampings.calc_factors_active(10)
        factors2 = dampings.calc_factors_active(18)

        self.assertEqual(len(factors1), 10)
        self.assertEqual(len(factors2), 18)

        with self.assertRaises(ValueError) as error:
            dampings.dampings_active(10, 200)
        error_message = "Number of dampings must be smaller than total number of days!"
        self.assertEqual(str(error.exception), error_message)

        # Combine days and factors
        days1, factors1 = dampings.dampings_active(30, 5)
        days2, factors2 = dampings.dampings_active(40, 10)

        self.assertEqual(len(days1), 5)
        self.assertEqual(len(factors1), len(days1))
        self.assertEqual(len(days2), 10)
        self.assertEqual(len(days2), len(factors2))

    def test_dampings_classic(self):
        days1 = dampings.generate_dampings_withshadowdamp(4, 30, 3, 4, 1)
        days2 = dampings.generate_dampings_withshadowdamp(10, 90, 2, 10, 1)
        self.assertEqual(len(days1[0]), 4)
        self.assertEqual(len(days2[0]), 10)

        with self.assertRaises(ValueError) as error:
            dampings.dampings_classic(10, 200)
        error_message = "Invalid input: It's not possible to generate this number of damping" \
            "in the desired time interval."
        self.assertEqual(str(error.exception), error_message)

        # Combine days and factors
        days1, factors1 = dampings.dampings_classic(30, 5)
        days2, factors2 = dampings.dampings_classic(40, 5)

        self.assertEqual(len(days1), 5)
        self.assertEqual(len(factors1), len(days1))
        self.assertEqual(len(days2), 5)
        self.assertEqual(len(days2), len(factors2))

    def test_dampings_random(self):
        with self.assertRaises(ValueError) as error:
            dampings.dampings_random(10, 200)
        error_message = "Invalid input: It's not possible to generate this number of damping" \
            "in the desired time interval."
        self.assertEqual(str(error.exception), error_message)

        days1, factors1 = dampings.dampings_random(30, 5)
        days2, factors2 = dampings.dampings_random(40, 5)

        self.assertEqual(len(days1), 5)
        self.assertEqual(len(factors1), len(days1))
        self.assertEqual(len(days2), 5)
        self.assertEqual(len(days2), len(factors2))

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

        damping_days = [3, 5]
        damping_factors = [0.3, 0.4]
        population = [5256.0, 10551, 32368.5,
                      43637.833333333336, 22874.066666666666, 8473.6]

        simulation_1 = data_generation.run_secir_groups_simulation(
            days_1, damping_days, damping_factors,  population)
        simulation_2 = data_generation.run_secir_groups_simulation(
            days_2, damping_days, damping_factors, population)
        simulation_3 = data_generation.run_secir_groups_simulation(
            days_3, damping_days, damping_factors, population)

        # result length
        self.assertEqual(len(simulation_1[0]), days_1+1)
        self.assertEqual(len(simulation_2[0]), days_2+1)
        self.assertEqual(len(simulation_3[0]), days_3+1)

        # damping
        self.assertEqual(
            simulation_1[1][0].size,
            len(population) * len(population))
        self.assertEqual(
            simulation_2[1][0].size,
            len(population) * len(population))
        self.assertEqual(
            simulation_3[1][0].size,
            len(population) * len(population))

        self.assertLessEqual(simulation_1[1][0].max(), 1)
        self.assertGreaterEqual(simulation_1[1][0].max(), 0)

        self.assertLessEqual(simulation_2[1][0].max(), 1)
        self.assertGreaterEqual(simulation_2[1][0].max(), 0)

        self.assertLessEqual(simulation_3[1][0].max(), 1)
        self.assertGreaterEqual(simulation_3[1][0].max(), 0)

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
        input_width_1 = 10
        input_width_2 = 5
        label_width_1 = 20
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
        : param mock_population:
        : param mock_baseline:
        : param mock_minimum:
        """
        input_width = 5
        label_width = 30
        num_runs = 1

        data_generation.generate_data(num_runs, self.path, "", input_width,
                                      label_width, damping_method="random")
        self.assertEqual(len(os.listdir(self.path)), 1)

        self.assertEqual(os.listdir(self.path),
                         ['data_secir_groups_30days_1_random.pickle'])


# # Testing network_architectures.py

    def test_mlp_multi_single(self):
        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_single_output(
                0, 12, 12, 12
            )
        error_message = "Number of age groups has to be positive, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_single_output(
                12, 0, 12, 12
            )
        error_message = "Output dimension must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_single_output(
                12, 12, -1, 12
            )
        error_message = "Number of layers must be at least 0, here -1"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_single_output(
                12, 12, 12, 0
            )
        error_message = "Number of neurons per layer must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        model = network_architectures.mlp_multi_input_single_output(
            7, 3, 2, 12
        )
        self.assertEqual(len(model.layers), 5)
        input_zero = np.zeros((1, 5, 8))
        output_zeros = model(input_zero)
        self.assertEqual(output_zeros.shape[0], 1)
        self.assertEqual(output_zeros.shape[1], 21)

    def test_mlp_multi_multi(self):
        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_multi_output(
                0, 12, 12, 12, 12
            )
        error_message = "label width has to be a positive integer"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_multi_output(
                12, 0, 12, 12, 12
            )
        error_message = "Number of age groups must be positive, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_multi_output(
                12, 12, 0, 12, 12
            )
        error_message = "Output dimension must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_multi_output(
                12, 12, 12, -1, 12
            )
        error_message = "Number of layers must be at least 0, here -1"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.mlp_multi_input_multi_output(
                12, 12, 12, 12, 0
            )
        error_message = "Number of neurons per layer must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        model = network_architectures.mlp_multi_input_multi_output(
            12, 3, 7, 2, 16)
        self.assertEqual(len(model.layers), 5)
        input_zero = np.zeros((5, 3, 7))
        output_zeros = model(input_zero)
        # Number of time series
        self.assertEqual(output_zeros.shape[0], 5)
        # Number of days per predicted time series
        self.assertEqual(output_zeros.shape[1], 12)
        # Size of feature space per day
        self.assertEqual(output_zeros.shape[2], 21)

    def test_cnn_multi_multi(self):
        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                0, 1, 1, 2, 1, 1, 1
            )
        error_message = "label width has to be a positive integer"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                1, 0, 1, 2, 1, 1, 1
            )
        error_message = "Number of age groups has to be positive, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                1, 1, 0, 2, 1, 1, 1
            )
        error_message = "Output dimension must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                1, 1, 1, 0, 1, 1, 1
            )
        error_message = "Size of the convolution kernel has to be larger than 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                1, 1, 1, 2, 0, 1, 1
            )
        error_message = "Number of filters must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                1, 1, 1, 2, 1, -1, 1
            )
        error_message = "Number of hidden layers must be at least 0, here -1"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.cnn_multi_input_multi_output(
                1, 1, 1, 2, 1, 1, 0
            )
        error_message = "Number of neurons per layer must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        model = network_architectures.cnn_multi_input_multi_output(
            21, 4, 4, 3, 256, 2)
        self.assertEqual(len(model.layers), 6)
        input_zero = np.zeros((12, 5, 7))
        output_zeros = model(input_zero)
        # Number of time series
        self.assertEqual(output_zeros.shape[0], 12)
        # length of one time series
        self.assertEqual(output_zeros.shape[1], 21)
        # Dimension of one time step (16 = 4*4)
        self.assertEqual(output_zeros.shape[2], 16)

    def test_lstm_multi_multi(self):
        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_multi_input_multi_output(
                0, 1, 1, 1, 1, 1
            )
        error_message = "label width has to be a positive integer"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_multi_input_multi_output(
                1, 0, 1, 1, 1, 1
            )
        error_message = "Number of age groups has to be positive, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_multi_input_multi_output(
                1, 1, 0, 1, 1, 1
            )
        error_message = "Output dimension must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_multi_input_multi_output(
                1, 1, 1, 0, 1, 1
            )
        error_message = "Internal dimension must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_multi_input_multi_output(
                1, 1, 1, 1, -1, 1
            )
        error_message = "Number of hidden layers must be at least 0, here -1"
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            network_architectures.lstm_multi_input_multi_output(
                1, 1, 1, 1, 1, 0
            )
        error_message = "Number of neurons per layer must be at least 1, here 0"
        self.assertEqual(str(error.exception), error_message)

        model = network_architectures.lstm_multi_input_multi_output(
            21, 4, 4, 12, 3, 12)
        self.assertEqual(len(model.layers), 6)
        input_zero = np.zeros((12, 5, 7))
        output_zeros = model(input_zero)
        # Number of time series
        self.assertEqual(output_zeros.shape[0], 12)
        # length of one time series
        self.assertEqual(output_zeros.shape[1], 21)
        # Dimension of one time step (16 = 4*4)
        self.assertEqual(output_zeros.shape[2], 16)

# Testing model.py
    def test_calc_split_index(self):
        with self.assertRaises(ValueError) as error:
            model.calc_split_index(
                10, 0.9, 0.1, 0.1
            )
        error_message = "Summed data set shares are greater than 1. Please adjust the values."
        self.assertEqual(str(error.exception), error_message)
        split_index = model.calc_split_index(10, 0.7, 0.1, 0.2)
        self.assertEqual(split_index, [7, 1, 2])

    def test_flat_input(self):
        A = np.zeros((12, 12))
        a1 = [A for _ in np.arange(5)]
        a2 = np.zeros((5, 144))
        a1_flatten = model.flat_input(a1)
        b = a2 == a1_flatten
        self.assertTrue(np.asarray(b).all())

    def test_prepare_data_classic(self):
        data = {
            "inputs": np.zeros((10, 5, 48)),
            "labels": np.zeros((10, 25, 48)),
            "contact_matrices": [[np.zeros((6, 6)) for _ in np.arange(2)] for _ in np.arange(10)],
            "damping_days": [[1, 2] for _ in np.arange(10)]
        }

        data_new = model.prepare_data_classic(data)
        a = np.zeros(5*48 + 2*36)
        b = np.array([1, 2])
        res = tf.cast(np.hstack((a, b)).reshape(
            (1, 5*48 + 2*36+2)), tf.float32)
        self.assertTrue(np.asarray(res == data_new["train_inputs"][0]).all())
        self.assertEqual(data_new["train_inputs"].shape, (7, 5*48 + 2*36 + 2))
        self.assertEqual(data_new["valid_labels"].shape, (2, 25, 48))

    def test_prod_time_series(self):
        A = np.ones((2, 3))
        obj = [i*A for i in np.arange(3)]
        ts = model.prod_time_series(obj, 3, 2)
        res = [[0*A, 0*A], [A, A], [2*A, 2*A]]
        bl = (ts == tf.reshape(
            tf.stack(res), [3, 2, -1]))
        self.assertTrue(tf.math.reduce_all(bl))

    def test_prepare_data_timeseries(self):
        data = {
            "inputs": np.zeros((10, 5, 48)),
            "labels": np.zeros((10, 25, 48)),
            "contact_matrices": [[np.zeros((6, 6)) for _ in np.arange(2)] for _ in np.arange(10)],
            "damping_days": [[1, 2] for _ in np.arange(10)]
        }
        data_new = model.prepare_data_timeseries(data)
        a = np.zeros(48 + 2*36)
        b = np.array([1, 2])
        res = tf.cast(np.hstack((a, b)).reshape(
            (1, 48 + 2*36+2)), tf.float32)
        res = tf.cast([res for _ in np.arange(5)], tf.float32)
        self.assertTrue(tf.math.reduce_all(
            data_new["train_inputs"].shape == (7, 5, 48 + 2*36 + 2)))
        self.assertTrue(tf.math.reduce_all(res == data_new["train_inputs"][0]))

    def test_initialize_model(self):
        # Helper function to normalize the .getconfig() output
        def normalize_config(config):
            config.pop('name', None)
            for layer in config["layers"]:
                layer["config"].pop("name", None)
            return config

        label_width = 8
        number_age_groups = 6
        num_outputs = 30
        hidden_layers = 2
        neurons_in_hidden_layer = 32
        activation_function = "relu"
        models = ["LSTM", "Dense", "CNN", "Heinz"]

        # Generating valid models
        model_lstm = network_architectures.lstm_multi_input_multi_output(
            label_width, number_age_groups, num_outputs, 32, hidden_layers,
            neurons_in_hidden_layer, activation_function
        )
        model_cnn = network_architectures.cnn_multi_input_multi_output(
            label_width=label_width, num_age_groups=number_age_groups, num_outputs=num_outputs,
            num_hidden_layers=hidden_layers, num_neurons_per_layer=neurons_in_hidden_layer, activation=activation_function
        )
        model_mlp = network_architectures.mlp_multi_input_multi_output(
            label_width=label_width, num_age_groups=number_age_groups, num_outputs=num_outputs,
            num_hidden_layers=hidden_layers, num_neurons_per_layer=neurons_in_hidden_layer, activation=activation_function
        )

        model_parameters = [(label_width, number_age_groups, num_outputs, hidden_layers,
                             neurons_in_hidden_layer, activation_function, modelname) for modelname in models]

        for parameter in model_parameters:
            _, _, _, _, _, _, modelname = parameter
            if modelname == "Heinz":
                with self.assertRaises(ValueError) as error:
                    md = model.initialize_model(parameter)
                error_message = "name_architecture must be one of 'Dense', 'LSTM' or 'CNN'"
                self.assertEqual(str(error.exception), error_message)
            else:
                md = model.initialize_model(parameter)
                if modelname == "Dense":
                    self.assertDictEqual(
                        normalize_config(md.get_config()),
                        normalize_config(model_mlp.get_config())
                    )
                if modelname == "LSTM":
                    self.assertDictEqual(
                        normalize_config(md.get_config()),
                        normalize_config(model_lstm.get_config())
                    )
                if modelname == "CNN":
                    self.assertDictEqual(
                        normalize_config(md.get_config()),
                        normalize_config(model_cnn.get_config())
                    )

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
        max_epochs = 1
        early_stop = 100
        loss = tf.keras.losses.MeanAbsolutePercentageError()
        optimizer = "Adam"
        metric = [tf.keras.metrics.MeanAbsoluteError(),
                  tf.keras.metrics.MeanAbsolutePercentageError()]

        training_parameter = (
            early_stop, max_epochs, loss, optimizer, metric)

        md = network_architectures.mlp_multi_input_multi_output(
            label_width=10
        )
        # no existing dataset
        with self.assertRaises(FileNotFoundError) as error:
            model.network_fit(
                path=self.path, model=md,
                modeltype='classic', training_parameter=training_parameter, filename="data_secir_groups.pickle")
        error_message = "[Errno 2] No such file or directory in the fake filesystem: '" + \
            self.path + "data_secir_groups.pickle'"
        self.assertEqual(str(error.exception), error_message)

        # generate dataset with multiple output
        input_width = 5
        label_width = 10
        num_runs = 5
        data_generation.generate_data(num_runs, self.path, "", input_width,
                                      label_width, max_number_damping=1)

        # models with multiple outputs
        model_mlp_multi_input_multi_output = model.network_fit(
            model=network_architectures.mlp_multi_input_multi_output(
                label_width),
            training_parameter=training_parameter,
            path=self.path,
            filename="data_secir_groups_10days_5_active.pickle",
            modeltype='classic',  plot_stats=False)
        self.assertEqual(
            model_mlp_multi_input_multi_output.model.output_shape[1],
            label_width)
        self.assertEqual(
            len(model_mlp_multi_input_multi_output.history['val_loss']),
            max_epochs)

        model_lstm_multi_output = model.network_fit(
            model=network_architectures.lstm_multi_input_multi_output(
                label_width),
            modeltype='timeseries',
            training_parameter=training_parameter,
            path=self.path,
            filename="data_secir_groups_10days_5_active.pickle",
            plot_stats=False)
        self.assertEqual(
            model_lstm_multi_output.model.output_shape[1], label_width)
        self.assertEqual(
            len(model_lstm_multi_output.history['val_loss']), max_epochs)

        cnn_output = model.network_fit(
            model=network_architectures.cnn_multi_input_multi_output(
                label_width),
            modeltype='timeseries',
            training_parameter=training_parameter,
            path=self.path,
            filename="data_secir_groups_10days_5_active.pickle",
            plot_stats=False)
        self.assertEqual(
            cnn_output.model.output_shape[1], label_width)
        self.assertEqual(
            len(cnn_output.history['val_loss']), max_epochs)

    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getMinimumMatrix',
           return_value=0.1 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.get_population',
           return_value=[[5256.0, 10551, 32368.5,
                         43637.833333333336, 22874.066666666666, 8473.6]])
    def test_save_model(
            self, mock_population, mock_baseline, mock_minimum):
        """
        : param mock_population:
        : param mock_baseline:
        : param mock_minimum:
        """
        input_width = 5
        label_width = 10
        num_runs = 2

        data_generation.generate_data(num_runs, self.path, "", input_width,
                                      label_width)
        max_epochs = 1
        early_stop = 100
        loss = tf.keras.losses.MeanAbsolutePercentageError()
        optimizer = "Adam"
        metric = [tf.keras.metrics.MeanAbsoluteError(),
                  tf.keras.metrics.MeanAbsolutePercentageError()]

        training_parameter = (
            early_stop, max_epochs, loss, optimizer, metric)

        model_mlp_multi_input_multi_output = model.network_fit(
            model=network_architectures.mlp_multi_input_multi_output(
                label_width),
            training_parameter=training_parameter,
            path=self.path,
            filename="data_secir_groups_10days_2_active.pickle",
            modeltype='classic',  plot_stats=False)

        model.save_model(model_mlp_multi_input_multi_output.model,
                         self.path, "mlp_multi_multi")

        self.assertEqual(len(os.listdir(self.path)), 2)

        self.assertEqual(os.listdir(self.path),
                         ['data_secir_groups_10days_2_active.pickle', 'mlp_multi_multi.keras'])

    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getMinimumMatrix',
           return_value=0.1 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.getBaselineMatrix',
           return_value=0.6 * np.ones((6, 6)))
    @patch('memilio.surrogatemodel.ode_secir_groups.data_generation.get_population',
           return_value=[[5256.0, 10551, 32368.5,
                         43637.833333333336, 22874.066666666666, 8473.6]])
    def test_load_model(
            self, mock_population, mock_baseline, mock_minimum):
        """
        : param mock_population:
        : param mock_baseline:
        : param mock_minimum:
        """
        input_width = 5
        label_width = 10
        num_runs = 2

        data_generation.generate_data(num_runs, self.path, "", input_width,
                                      label_width)
        max_epochs = 1
        early_stop = 100
        loss = tf.keras.losses.MeanAbsolutePercentageError()
        optimizer = "Adam"
        metric = [tf.keras.metrics.MeanAbsoluteError(),
                  tf.keras.metrics.MeanAbsolutePercentageError()]

        training_parameter = (
            early_stop, max_epochs, loss, optimizer, metric)

        mlp1 = model.network_fit(
            model=network_architectures.mlp_multi_input_multi_output(
                label_width),
            training_parameter=training_parameter,
            path=self.path,
            filename="data_secir_groups_10days_2_active.pickle",
            modeltype='classic',  plot_stats=False)

        model.save_model(mlp1.model,
                         self.path, "mlp_multi_multi")

        path = os.path.join(self.path, "mlp_multi_multi.keras")
        mlp2 = model.load_model(path)

        weights1 = mlp1.model.get_weights()
        weights2 = mlp2.get_weights()
        for w1, w2 in zip(weights1, weights2):
            np.testing.assert_allclose(w1, w2)


if __name__ == '__main__':
    unittest.main()
