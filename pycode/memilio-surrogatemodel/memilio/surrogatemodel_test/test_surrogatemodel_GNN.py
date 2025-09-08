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
import memilio.surrogatemodel.GNN.network_architectures as gnn_arch

from unittest.mock import patch
import os
import unittest
import pickle
import numpy as np
import pandas as pd
import tensorflow as tf
import logging

# suppress all autograph warnings from tensorflow

logging.getLogger("tensorflow").setLevel(logging.ERROR)


class TestSurrogatemodelGNN(fake_filesystem_unittest.TestCase):
    path = "/home/"

    def setUp(self):
        self.setUpPyfakefs()

    def test_generate_model_class(self):
        from memilio.surrogatemodel.GNN.network_architectures import generate_model_class

        # Test parameters
        layer_types = [
            lambda: tf.keras.layers.Dense(10, activation="relu"),
        ]
        num_layers = [3]
        num_output = 2

        # Generate the model class
        ModelClass = generate_model_class(
            "TestModel", layer_types, num_layers, num_output)
        # Check if the generated class is a subclass of tf.keras.Model
        self.assertTrue(issubclass(ModelClass, tf.keras.Model))

        # Instantiate the model
        model = ModelClass()

        # Check if the model has the expected number of layers
        expected_num_layers = num_layers[0] + 1  # +1 for the output layer
        self.assertEqual(len(model.layers), expected_num_layers)
        self.assertIsInstance(model.layers[-1], tf.keras.layers.Dense)
        self.assertEqual(model.layers[-1].units, num_output)
        self.assertEqual(
            model.layers[-1].activation.__name__, "relu")

    def test_get_model(self):
        from memilio.surrogatemodel.GNN.network_architectures import get_model

        # Test parameters
        layer_type = "GCNConv"
        num_layers = 2
        num_channels = 16
        activation = "relu"
        num_output = 3

        # Generate the model
        model = get_model(layer_type, num_layers,
                          num_channels, activation, num_output)

        # Check if the model is an instance of tf.keras.Model
        self.assertIsInstance(model, tf.keras.Model)

        # Check if the model has the expected number of layers
        expected_num_layers = num_layers + 1  # +1 for the output layer
        self.assertEqual(len(model.layers), expected_num_layers)

        layer_type = "MonvConv"
        with self.assertRaises(ValueError) as error:
            model = get_model(layer_type, num_layers,
                              num_channels, activation, num_output)
        self.assertEqual(str(
            error.exception), "Unsupported layer_type: MonvConv. "
            "Supported types are 'ARMAConv', 'GCSConv', 'GATConv', 'GCNConv', 'APPNPConv'.")
        layer_type = "GATConv"
        num_layers = 0
        with self.assertRaises(ValueError) as error:
            model = get_model(layer_type, num_layers,
                              num_channels, activation, num_output)
        self.assertEqual(str(
            error.exception), "num_layers must be at least 1.")

        num_layers = 2
        num_output = 0
        with self.assertRaises(ValueError) as error:
            model = get_model(layer_type, num_layers,
                              num_channels, activation, num_output)
        self.assertEqual(str(
            error.exception), "num_output must be at least 1.")

        num_output = 2
        num_channels = 0
        with self.assertRaises(ValueError) as error:
            model = get_model(layer_type, num_layers,
                              num_channels, activation, num_output)
        self.assertEqual(str(
            error.exception), "num_channels must be at least 1.")

        num_channels = 16
        activation = 5
        with self.assertRaises(ValueError) as error:
            model = get_model(layer_type, num_layers,
                              num_channels, activation, num_output)
        self.assertEqual(str(
            error.exception), "activation must be a string representing the activation function.")

    def test_create_dataset(self):
        from memilio.surrogatemodel.GNN.evaluate_and_train import create_dataset

        # Create dummy data in the fake filesystem for testing
        num_samples = 10
        num_nodes = 5
        num_node_features = 3
        output_dim = 4
        X = np.random.rand(num_samples, 1,
                           num_node_features, num_nodes).astype(np.float32)
        A = np.random.randint(0, 2, (num_nodes, num_nodes)).astype(np.float32)
        y = np.random.rand(num_samples, 1,
                           output_dim,  num_nodes).astype(np.float32)
        data = {"inputs": X, "labels": y}
        path_cases_dir = os.path.join(self.path, "cases")
        self.fs.create_dir(path_cases_dir)
        path_cases = os.path.join(path_cases_dir, "cases.pickle")
        with open(path_cases, 'wb') as f:
            pickle.dump(data, f)
        path_mobility = os.path.join(self.path, "mobility")
        mobility_file = os.path.join(
            path_mobility, "commuter_mobility_2022.txt")
        self.fs.create_file(mobility_file)
        with open(mobility_file, 'w') as f:
            np.savetxt(f, A, delimiter=" ")
        # Create dataset
        dataset = create_dataset(
            path_cases, path_mobility, number_of_nodes=num_nodes)
        self.assertEqual(len(dataset), num_samples)
        for graph in dataset:
            self.assertEqual(graph.x.shape, (num_nodes, num_node_features))
            self.assertEqual(dataset.a.shape, (num_nodes, num_nodes))
            self.assertEqual(graph.y.shape, (num_nodes, output_dim))
        # Clean up
        self.fs.remove_object(path_cases)
        self.fs.remove_object(os.path.join(
            path_mobility, "commuter_mobility_2022.txt"))
        self.fs.remove_object(path_mobility)

    def test_train_step(self):
        from memilio.surrogatemodel.GNN.evaluate_and_train import train_step
        from tensorflow.keras.losses import MeanAbsolutePercentageError

        # Create a simple model for testing
        model = tf.keras.Sequential([
            tf.keras.layers.Dense(10, activation='relu'),
            tf.keras.layers.Dense(2)
        ])
        optimizer = tf.keras.optimizers.Adam()
        loss_fn = MeanAbsolutePercentageError()

        # Create dummy data
        batch_size = 4
        num_nodes = 5
        num_node_features = 3
        output_dim = 2
        inputs = np.random.rand(batch_size, num_nodes,
                                num_node_features).astype(np.float32)
        y = np.random.rand(batch_size, num_nodes,
                           output_dim).astype(np.float32)

        # Perform a training step
        loss, acc = train_step(inputs, y, loss_fn, model, optimizer)

        # Check if the loss is a scalar tensor
        self.assertIsInstance(loss, tf.Tensor)
        self.assertEqual(loss.shape, ())
        self.assertIsInstance(acc, tf.Tensor)
        self.assertEqual(acc.shape, ())
        self.assertGreaterEqual(acc.numpy(), 0)
        self.assertGreaterEqual(loss.numpy(), 0)

    def test_evaluate(self):
        from memilio.surrogatemodel.GNN.evaluate_and_train import (
            evaluate, create_dataset, MixedLoader)
        from memilio.surrogatemodel.GNN.network_architectures import get_model
        from tensorflow.keras.losses import MeanAbsolutePercentageError

        # Create a simple model for testing
        model = get_model(
            layer_type="GCNConv", num_layers=2, num_channels=16, activation="relu", num_output=4)
        loss_fn = MeanAbsolutePercentageError()

        # Create dummy data in the fake filesystem for testing
        num_samples = 10
        num_nodes = 5
        num_node_features = 3
        output_dim = 4
        X = np.random.rand(num_samples, 1,
                           num_node_features, num_nodes).astype(np.float32)
        A = np.random.randint(0, 2, (num_nodes, num_nodes)).astype(np.float32)
        y = np.random.rand(num_samples, 1,
                           output_dim,  num_nodes).astype(np.float32)
        data = {"inputs": X, "labels": y}
        path_cases_dir = os.path.join(self.path, "cases")
        self.fs.create_dir(path_cases_dir)
        path_cases = os.path.join(path_cases_dir, "cases.pickle")
        with open(path_cases, 'wb') as f:
            pickle.dump(data, f)
        path_mobility = os.path.join(self.path, "mobility")
        mobility_file = os.path.join(
            path_mobility, "commuter_mobility_2022.txt")
        self.fs.create_file(mobility_file)
        with open(mobility_file, 'w') as f:
            np.savetxt(f, A, delimiter=" ")
        # Create dataset
        dataset = create_dataset(
            path_cases, path_mobility, number_of_nodes=num_nodes)
        loader = MixedLoader(dataset, batch_size=2, epochs=1)
        res = evaluate(loader, model, loss_fn)

        self.assertEqual(len(res), 2)
        self.assertGreaterEqual(res[0], 0)
        self.assertGreaterEqual(res[1], 0)

        # Test with retransformation
        loader = MixedLoader(dataset, batch_size=2, epochs=1)
        res = evaluate(loader, model, loss_fn, True)

        self.assertEqual(len(res), 2)
        self.assertGreaterEqual(res[0], 0)
        self.assertGreaterEqual(res[1], 0)
        # Clean up
        self.fs.remove_object(path_cases)
        self.fs.remove_object(os.path.join(
            path_mobility, "commuter_mobility_2022.txt"))
        self.fs.remove_object(path_mobility)

    def test_train_and_evaluate(self):
        from memilio.surrogatemodel.GNN.evaluate_and_train import (
            train_and_evaluate, create_dataset, MixedLoader)
        from memilio.surrogatemodel.GNN.network_architectures import get_model
        from tensorflow.keras.losses import MeanAbsolutePercentageError

        number_of_epochs = 2
        # Create a simple model for testing
        model = get_model(
            layer_type="GCNConv", num_layers=2, num_channels=16, activation="relu", num_output=4)

        # Create dummy data in the fake filesystem for testing
        num_samples = 20
        num_nodes = 5
        num_node_features = 3
        output_dim = 4
        X = np.random.rand(num_samples, 1,
                           num_node_features, num_nodes).astype(np.float32)
        A = np.random.randint(0, 2, (num_nodes, num_nodes)).astype(np.float32)
        y = np.random.rand(num_samples, 1,
                           output_dim,  num_nodes).astype(np.float32)
        data = {"inputs": X, "labels": y}
        path_cases_dir = os.path.join(self.path, "cases")
        self.fs.create_dir(path_cases_dir)
        path_cases = os.path.join(path_cases_dir, "cases.pickle")
        with open(path_cases, 'wb') as f:
            pickle.dump(data, f)
        path_mobility = os.path.join(self.path, "mobility")
        mobility_file = os.path.join(
            path_mobility, "commuter_mobility_2022.txt")
        self.fs.create_file(mobility_file)
        with open(mobility_file, 'w') as f:
            np.savetxt(f, A, delimiter=" ")
        # Create dataset
        dataset = create_dataset(
            path_cases, path_mobility, number_of_nodes=num_nodes)

        res = train_and_evaluate(
            dataset,
            batch_size=2,
            epochs=number_of_epochs,
            model=model,
            loss_fn=MeanAbsolutePercentageError(),
            optimizer=tf.keras.optimizers.Adam(),
            es_patience=100)

        self.assertEqual(len(res["train_losses"][0][0]), number_of_epochs)

    @patch("os.path.realpath", return_value="/home/")
    def test_perform_grid_search(self, mock_realpath):
        from memilio.surrogatemodel.GNN.evaluate_and_train import (
            create_dataset)
        from memilio.surrogatemodel.GNN.grid_search import perform_grid_search
        from memilio.surrogatemodel.GNN.network_architectures import get_model
        from tensorflow.keras.losses import MeanAbsolutePercentageError

        # Create dummy data in the fake filesystem for testing
        num_samples = 20
        num_nodes = 5
        num_node_features = 3
        output_dim = 4
        X = np.random.rand(num_samples, 1,
                           num_node_features, num_nodes).astype(np.float32)
        A = np.random.randint(0, 2, (num_nodes, num_nodes)).astype(np.float32)
        y = np.random.rand(num_samples, 1,
                           output_dim,  num_nodes).astype(np.float32)
        data = {"inputs": X, "labels": y}
        path_cases_dir = os.path.join(self.path, "cases")
        self.fs.create_dir(path_cases_dir)
        path_cases = os.path.join(path_cases_dir, "cases.pickle")
        with open(path_cases, 'wb') as f:
            pickle.dump(data, f)
        path_mobility = os.path.join(self.path, "mobility")
        mobility_file = os.path.join(
            path_mobility, "commuter_mobility_2022.txt")
        self.fs.create_file(mobility_file)
        with open(mobility_file, 'w') as f:
            np.savetxt(f, A, delimiter=" ")
        # Create dataset
        dataset = create_dataset(
            path_cases, path_mobility, number_of_nodes=num_nodes)
        layers = ["GCNConv", "GATConv"]
        num_layers = [1, 2]
        activation = "relu"
        num_channels = 2

        model_parameters = [
            (layer, n_layer, num_channels, activation) for layer in layers for n_layer in num_layers
        ]
        batch_size = 2
        loss_function = MeanAbsolutePercentageError()
        optimizer = tf.keras.optimizers.Adam()
        es_patience = 5
        max_epochs = 2
        training_parameters = [batch_size, loss_function,
                               optimizer, es_patience, max_epochs]
        # Perform grid search
        perform_grid_search(model_parameters, training_parameters,
                            dataset)

        # Check if the results file is created
        results_file = os.path.join(
            self.path, "saves", "grid_search_results.csv")
        self.assertTrue(os.path.exists(results_file))

        # Check if the results file has the expected number of rows
        df_results = pd.read_csv(results_file)
        print(df_results.columns)
        self.assertEqual(len(df_results), len(model_parameters))
        self.assertEqual(len(df_results.columns), 10)
        # Clean up
        self.fs.remove_object(path_cases)
        self.fs.remove_object(os.path.join(
            path_mobility, "commuter_mobility_2022.txt"))
        self.fs.remove_object(path_mobility)
        self.fs.remove_object(os.path.join(self.path, "saves"))


if __name__ == '__main__':
    unittest.main()
