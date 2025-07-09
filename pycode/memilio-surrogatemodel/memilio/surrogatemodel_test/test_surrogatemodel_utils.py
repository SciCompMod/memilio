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
import memilio.epidata.modifyDataframeSeries
from memilio.surrogatemodel.ode_secir_groups import (data_generation, model,
                                                     network_architectures)
import memilio.surrogatemodel.utils.helper_functions as utils
from memilio.surrogatemodel.utils import dampings, grid_search
from sklearn.preprocessing import FunctionTransformer

from unittest.mock import patch
import pandas as pd
import os
import unittest
import pickle

import numpy as np
import tensorflow as tf
import logging

# suppress all autograph warnings from tensorflow

logging.getLogger("tensorflow").setLevel(logging.ERROR)


class TestSurrogatemodelUtils(fake_filesystem_unittest.TestCase):
    """ """

    path = '/home/'

    def setUp(self):
        """ """
        self.setUpPyfakefs()

    def test_interpolate_age_groups(self):
        entry = {'ID_County': 1000, 'Population': 11000, '<3 years': 1000, '3-5 years': 1000, '6-14 years': 1000, '15-17 years': 1000,
                 '18-24 years': 1000, '25-29 years': 1000, '30-39 years': 1000, '40-49 years': 1000, '50-64 years': 1000, '65-74 years': 1000, '>74 years': 1000}
        res = [1666.6666666666665, 1333.3333333333335, 3500.0,
               2166.6666666666665, 1533.3333333333335, 800.0]
        interpol_entry = utils.interpolate_age_groups(
            entry, ['0-4', '5-14', '15-34', '35-59', '60-79', '>79'])
        self.assertEqual(res, interpol_entry)

    def test_calc_dist_days(self):
        with self.assertRaises(ValueError) as error:
            dampings.calc_dist_days(10, 2, 10, 2)
        error_message = "It's not possible to arrange this number of dampings in the desired interval with the given minimal distance."
        self.assertEqual(str(error.exception), error_message)

        with self.assertRaises(ValueError) as error:
            dampings.calc_dist_days(102, 2, 50, 3)
        error_message = "It's not possible to arrange this number of dampings in the desired interval with the given minimal distance."
        self.assertEqual(str(error.exception), error_message)

        self.assertEqual(dampings.calc_dist_days(102, 2, 50, 1), 2)
        self.assertEqual(dampings.calc_dist_days(30, 2, 10, 1), 2)

    def test_generate_dampings(self):
        with self.assertRaises(ValueError) as error:
            dampings.generate_dampings(100, 2, "gelb")
        error_message = "The method argument has to be one of the following: 'classic', 'active' or 'random'."
        self.assertEqual(str(error.exception), error_message)

    def test_dampings_active(self):
        # length of generated sequence
        factors1 = dampings.calc_factors_active(10)
        factors2 = dampings.calc_factors_active(18)
        # testing if damping factors are in the desired interval
        self.assertTrue(np.all(np.log(1-np.array(factors1)) <= 2))
        self.assertTrue(np.all(np.log(1-np.array(factors2)) <= 2))
        self.assertEqual(len(factors1), 10)
        self.assertEqual(len(factors2), 18)

        with self.assertRaises(ValueError) as error:
            dampings.dampings_active(10, 200, 3)
        error_message = "It's not possible to arrange this number of dampings in the desired interval with the given minimal distance."
        self.assertEqual(str(error.exception), error_message)

        # Combine days and factors
        days1, factors1 = dampings.dampings_active(30, 5, 4)
        days2, factors2 = dampings.dampings_active(40, 10, 4)

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
        error_message = "It's not possible to arrange this number of dampings in the desired interval with the given minimal distance."
        self.assertEqual(str(error.exception), error_message)

        # Combine days and factors
        days1, factors1 = dampings.dampings_classic(30, 5)
        days2, factors2 = dampings.dampings_classic(40, 5)

        self.assertEqual(len(days1), 5)
        self.assertEqual(len(factors1), len(days1))
        self.assertEqual(len(days2), 5)
        self.assertEqual(len(days2), len(factors2))
        # testing if damping factors are in the desired interval
        self.assertTrue(np.all(np.array(factors1) <= 2))
        self.assertTrue(np.all(np.array(factors2) <= 2))

    def test_dampings_random(self):
        with self.assertRaises(ValueError) as error:
            dampings.dampings_random(10, 200)
        error_message = "It's not possible to arrange this number of dampings in the desired interval with the given minimal distance."
        self.assertEqual(str(error.exception), error_message)

        days1, factors1 = dampings.dampings_random(30, 5)
        days2, factors2 = dampings.dampings_random(40, 5)

        self.assertEqual(len(days1), 5)
        self.assertEqual(len(factors1), len(days1))
        self.assertEqual(len(days2), 5)
        self.assertEqual(len(days2), len(factors2))
        # testing if damping factors are in the desired interval
        self.assertTrue(np.all(np.array(factors1) <= 2))
        self.assertTrue(np.all(np.array(factors2) <= 2))

    # Testing helper_functions.py
    def test_calc_split_index(self):
        with self.assertRaises(ValueError) as error:
            utils.calc_split_index(
                10, 0.9, 0.1, 0.1
            )
        error_message = "Summed data set shares are greater than 1. Please adjust the values."
        self.assertEqual(str(error.exception), error_message)
        split_index = utils.calc_split_index(10, 0.7, 0.1, 0.2)
        self.assertEqual(split_index, [7, 1, 2])

    def test_remove_confirmed_compartments(self):
        A = np.asarray([[1, 2, 3, 4, 5]])

        A_new = utils.remove_confirmed_compartments(A, [1, 2])
        self.assertTrue(np.all(np.asarray([[1, 4, 5]]) == A_new))

    def test_normalize_simulation_data(self):
        transformer = FunctionTransformer(np.log1p, validate=True)
        A = [[[1, 2], [1, 2]],
             [[3, 4], [3, 4]]
             ]
        res = [[[0.69314718, 1.09861229], [0.69314718, 1.09861229],
                [1.38629436, 1.60943791], [1.38629436, 1.60943791]]]
        output = utils.normalize_simulation_data(A, transformer, 1, 1, 2)
        for w1, w2 in zip(res, output):
            np.testing.assert_allclose(w1, w2)

    def test_flat_input(self):
        A = np.zeros((12, 12))
        a1 = [A for _ in np.arange(5)]
        a2 = np.zeros((5, 144))
        a1_flatten = utils.flat_input(a1)
        b = a2 == a1_flatten
        self.assertTrue(np.asarray(b).all())

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
                                      label_width, number_dampings=2)
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
            filename="data_secir_groups_10days_2_random.pickle",
            modeltype='classic',  plot_stats=False)

        utils.save_model(model_mlp_multi_input_multi_output.model,
                         self.path, "mlp_multi_multi")

        self.assertEqual(len(os.listdir(self.path)), 2)

        self.assertEqual(os.listdir(self.path),
                         ['data_secir_groups_10days_2_random.pickle', 'mlp_multi_multi.keras'])

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
                                      label_width, number_dampings=2)
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
            filename="data_secir_groups_10days_2_random.pickle",
            modeltype='classic',  plot_stats=False)

        utils.save_model(mlp1.model, self.path, "mlp_multi_multi")
        path_file = os.path.join(self.path, "mlp_multi_multi.keras")
        if os.path.isfile(path_file):
            print("File !")

        mlp2 = utils.load_model(path_file)

        weights1 = mlp1.model.get_weights()
        weights2 = mlp2.get_weights()
        for w1, w2 in zip(weights1, weights2):
            np.testing.assert_allclose(w1, w2)

    def test_train_and_evaluate_model(self):
        """ Testing evaluation procedure """
        inputs = tf.ones([5, 5, 8])
        labels = tf.ones([5, 10, 8])
        max_epochs = 1
        early_stop = 100
        loss = tf.keras.losses.MeanAbsolutePercentageError()
        optimizer = 'Adam'
        metric = [tf.keras.metrics.MeanAbsoluteError(),
                  tf.keras.metrics.MeanAbsolutePercentageError()]
        model_parameter = (
            10, 8, 1, 5, "relu", "Dense"
        )

        training_parameter = (
            early_stop, max_epochs, loss, optimizer, metric)

        res = grid_search.train_and_evaluate_model(
            model_parameter, inputs, labels, training_parameter, False)
        self.assertEqual(len(res["val_losses"][0]), 5)
        self.assertEqual(len(res["val_losses"][0][0]), max_epochs)
        self.assertEqual(res["optimizer"], "Adam")

    @patch('memilio.surrogatemodel.utils.grid_search.train_and_evaluate_model')
    def test_perform_grid_search_simple(self, mock_train_and_evaluate_model):
        """ Testing grid search for simple model along a sequence of possible parameter values"""
        mock_train_and_evaluate_model.return_value = {
            "model": "._.",
            "activation": "relu",
            "optimizer": "Eva",
            "mean_train_loss_kfold": 42,
            "mean_val_loss_kfold": 12,
            "training_time": 1,
            "train_losses": [3],
            "val_losses": [6]
        }
        filename_df = "dataframe_optimizer"
        os.mkdir(self.path)

        # General grid search parameters for the training process:
        early_stop = [100]
        max_epochs = [1]
        losses = [tf.keras.losses.MeanAbsolutePercentageError()]
        optimizers = ['Adam', 'Adam']
        metrics = [[tf.keras.metrics.MeanAbsoluteError(),
                    tf.keras.metrics.MeanAbsolutePercentageError()]]

        # Define grid search parameters for the architecture
        label_width = 3
        num_outputs = 2
        hidden_layers = [1, 2]
        neurons_in_hidden_layer = [32]
        activation_function = ['relu']
        models = ["Dense"]

        # Collecting parameters
        training_parameters = [(early, epochs, loss, optimizer, metric)
                               for early in early_stop for epochs in max_epochs for loss in losses
                               for optimizer in optimizers for metric in metrics]

        model_parameters = [(label_width, num_outputs, layer, neuron_number, activation, modelname)
                            for layer in hidden_layers for neuron_number in neurons_in_hidden_layer
                            for activation in activation_function for modelname in models]

        # generate dataset with multiple output
        inputs = tf.ones([5, 1, 2])
        labels = tf.ones([5, 3, 2])

        grid_search.perform_grid_search(
            model_parameters, inputs, labels, training_parameters,
            filename_df, self.path
        )

        # testing saving the results
        self.assertEqual(len(os.listdir(self.path)), 1)

        self.assertEqual(os.listdir(os.path.join(self.path,  'secir_simple_grid_search')),
                         [filename_df+".pickle"])

        # testing size of the output
        path_name = os.path.join(self.path,  'secir_simple_grid_search')
        with open(os.path.join(path_name, filename_df + ".pickle"), "rb") as f:
            d = pickle.load(f)
        df = pd.DataFrame(d)
        self.assertEqual(len(df['model']), 4)

    @patch('memilio.surrogatemodel.utils.grid_search.train_and_evaluate_model')
    def test_perform_grid_search_groups(self, mock_train_and_evaluate_model):
        """ Testing grid search for groups model along a sequence of possible parameter values"""
        mock_train_and_evaluate_model.return_value = {
            "model": "._.",
            "activation": "relu",
            "optimizer": "Eva",
            "mean_train_loss_kfold": 42,
            "mean_val_loss_kfold": 12,
            "training_time": 1,
            "train_losses": [3],
            "val_losses": [6]
        }
        filename_df = "dataframe_optimizer"
        os.mkdir(self.path)

        # General grid search parameters for the training process:
        early_stop = [100]
        max_epochs = [1]
        losses = [tf.keras.losses.MeanAbsolutePercentageError()]
        optimizers = ['Adam', 'Adam']
        metrics = [[tf.keras.metrics.MeanAbsoluteError(),
                    tf.keras.metrics.MeanAbsolutePercentageError()]]

        # Define grid search parameters for the architecture
        label_width = 3
        num_outputs = 2
        num_age_groups = 1
        hidden_layers = [1, 2]
        neurons_in_hidden_layer = [32]
        activation_function = ['relu']
        models = ["Dense"]

        # Collecting parameters
        training_parameters = [(early, epochs, loss, optimizer, metric)
                               for early in early_stop for epochs in max_epochs for loss in losses
                               for optimizer in optimizers for metric in metrics]

        model_parameters = [(label_width, num_age_groups, num_outputs, layer, neuron_number, activation, modelname)
                            for layer in hidden_layers for neuron_number in neurons_in_hidden_layer
                            for activation in activation_function for modelname in models]

        # generate dataset with multiple output
        inputs = tf.ones([5, 1, 2])
        labels = tf.ones([5, 3, 2])

        grid_search.perform_grid_search(
            model_parameters, inputs, labels, training_parameters,
            filename_df, self.path, "groups"
        )

        # testing saving the results
        self.assertEqual(len(os.listdir(self.path)), 1)

        self.assertEqual(os.listdir(os.path.join(self.path,  'secir_groups_grid_search')),
                         [filename_df+".pickle"])

        # testing size of the output
        path_name = os.path.join(self.path,  'secir_groups_grid_search')
        with open(os.path.join(path_name, filename_df + ".pickle"), "rb") as f:
            d = pickle.load(f)
        df = pd.DataFrame(d)
        self.assertEqual(len(df['model']), 4)


if __name__ == '__main__':
    unittest.main()
