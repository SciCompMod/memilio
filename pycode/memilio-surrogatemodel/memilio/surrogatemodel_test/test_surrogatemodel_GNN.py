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
from unittest.mock import patch
import os
import unittest
import pickle
import numpy as np
import pandas as pd
import tensorflow as tf
import logging
import spektral

import memilio.surrogatemodel.GNN.network_architectures as gnn_arch
import memilio.surrogatemodel.GNN.GNN_utils as utils

from memilio.surrogatemodel.GNN.evaluate_and_train import (
    create_dataset, train_and_evaluate, evaluate, train_step, MixedLoader)
from memilio.surrogatemodel.GNN.grid_search import perform_grid_search
from memilio.surrogatemodel.GNN.network_architectures import generate_model_class, get_model
from tensorflow.keras.losses import MeanAbsolutePercentageError
# suppress all autograph warnings from tensorflow

logging.getLogger("tensorflow").setLevel(logging.ERROR)


class TestSurrogatemodelGNN(fake_filesystem_unittest.TestCase):
    path = "/home/"

    def create_dummy_data(
            self, num_samples, num_nodes, num_node_features, output_dim):
        """
        Create dummy data for testing.

        :param num_samples: Number of samples in the dataset.
        :param num_nodes: Number of nodes in each graph.
        :param num_node_features: Number of features per node.
        :param output_dim: Number of output dimensions per node.
        :return: A dictionary containing inputs, adjacency matrix, and labels.
        """
        # Shape should be (num_samples, input_width, num_nodes, num_features)
        X = np.random.rand(num_samples, 1, num_nodes,
                           num_node_features).astype(np.float32)
        A = np.random.randint(0, 2, (num_nodes, num_nodes)).astype(np.float32)
        # Shape should be (num_samples, label_width, num_nodes, label_features)
        y = np.random.rand(num_samples, 1, num_nodes,
                           output_dim).astype(np.float32)
        return {"inputs": X, "adjacency": A, "labels": y}

    def setup_fake_filesystem(self, fs, path, data):
        """
        Save dummy data to the fake file system.

        :param fs: The fake file system object.
        :param path: The base path in the fake file system.
        :param data: The dummy data dictionary containing inputs, adjacency, and labels.
        :return: Paths to the cases and mobility files.
        """
        path_cases_dir = os.path.join(path, "cases")
        fs.create_dir(path_cases_dir)
        path_cases = os.path.join(path_cases_dir, "cases.pickle")
        with open(path_cases, 'wb') as f:
            pickle.dump({"inputs": data["inputs"],
                        "labels": data["labels"]}, f)

        path_mobility = os.path.join(path, "mobility")
        mobility_file = os.path.join(
            path_mobility, "commuter_mobility_2022.txt")
        fs.create_dir(path_mobility)
        fs.create_file(mobility_file)
        with open(mobility_file, 'w') as f:
            np.savetxt(f, data["adjacency"], delimiter=" ")

        return path_cases, path_mobility

    def setUp(self):
        self.setUpPyfakefs()

    def test_generate_model_class(self):

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

        # Test with invalid parameters
        layer_types = [
            lambda: tf.keras.layers.Dense(10, activation="relu"),
        ]
        num_layers = [0]
        num_output = 2
        with self.assertRaises(ValueError) as error:
            ModelClass = generate_model_class(
                "TestModel", layer_types, num_layers, num_output)
            model = ModelClass()
        self.assertEqual(str(
            error.exception), "All values in num_repeat must be at least 1.")

        num_layers = [3, 2]
        with self.assertRaises(ValueError) as error:
            ModelClass = generate_model_class(
                "TestModel", layer_types, num_layers, num_output)
            model = ModelClass()
        self.assertEqual(
            str(error.exception),
            "layer_types and num_repeat must have the same length. "
            "Got 1 and 2.")

        # Test with multiple layer types
        layer_types = [
            lambda: tf.keras.layers.Dense(10, activation="relu"),
            lambda: tf.keras.layers.Dense(20, activation="relu"),
        ]
        num_repeat = [2, 3]
        num_output = 4

        ModelClass = generate_model_class(
            "TestModel", layer_types, num_repeat, num_output)
        model = ModelClass()

        # Check the number of layers
        self.assertEqual(len(model.layer_seq), sum(num_repeat))
        self.assertEqual(model.output_layer.units, num_output)

    def test_get_model(self):

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

        # Check handling of invalid parameters
        # Test with invalid layer type
        layer_type = "MonvConv"
        with self.assertRaises(ValueError) as error:
            model = get_model(layer_type, num_layers,
                              num_channels, activation, num_output)
        self.assertEqual(
            str(error.exception),
            "Unsupported layer_type: 'MonvConv'. "
            "Supported types are: ARMAConv, GCSConv, GATConv, GCNConv, APPNPConv")
        # Test with invalud num_layers
        layer_type = "GATConv"
        num_layers = 0
        with self.assertRaises(ValueError) as error:
            model = get_model(layer_type, num_layers,
                              num_channels, activation, num_output)
        self.assertEqual(str(
            error.exception), "num_layers must be at least 1, got 0.")
        # Test with invalid num_output
        num_layers = 2
        num_output = 0
        with self.assertRaises(ValueError) as error:
            model = get_model(layer_type, num_layers,
                              num_channels, activation, num_output)
        self.assertEqual(str(
            error.exception), "num_output must be at least 1, got 0.")
        # Test with invalid num_channels
        num_output = 2
        num_channels = 0
        with self.assertRaises(ValueError) as error:
            model = get_model(layer_type, num_layers,
                              num_channels, activation, num_output)
        self.assertEqual(str(
            error.exception), "num_channels must be at least 1, got 0.")
        # Test with invalid activation
        num_channels = 16
        activation = 5
        with self.assertRaises(ValueError) as error:
            model = get_model(layer_type, num_layers,
                              num_channels, activation, num_output)
        self.assertEqual(
            str(error.exception),
            "activation must be a string, got int.")

    def test_create_dataset(self):

        # Create dummy data in the fake filesystem for testing
        num_samples = 10
        num_nodes = 5
        num_node_features = 3
        output_dim = 4
        data = self.create_dummy_data(
            num_samples, num_nodes, num_node_features, output_dim)
        # Save dummy data to the fake file system
        path_cases, path_mobility = self.setup_fake_filesystem(
            self.fs, self.path, data)

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

        # Create a simple model for testing
        model = get_model(
            layer_type="GCNConv", num_layers=2, num_channels=16,
            activation="relu", num_output=2)
        optimizer = tf.keras.optimizers.Adam()
        loss_fn = MeanAbsolutePercentageError()

        # Create dummy data
        num_samples = 10
        num_nodes = 5
        num_node_features = 3
        output_dim = 2
        data = self.create_dummy_data(
            num_samples, num_nodes, num_node_features, output_dim)
        # Save dummy data to the fake file system
        path_cases, path_mobility = self.setup_fake_filesystem(
            self.fs, self.path, data)
        # Create dataset
        dataset = create_dataset(
            path_cases, path_mobility, number_of_nodes=num_nodes)
        # Build the model by calling it on a batch of data
        loader = MixedLoader(dataset, batch_size=4, epochs=1)
        inputs, y = next(loader)
        model(inputs)

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

        # Create a simple model for testing
        model = get_model(
            layer_type="GCNConv", num_layers=2, num_channels=16,
            activation="relu", num_output=4)
        loss_fn = MeanAbsolutePercentageError()

        # Create dummy data in the fake filesystem for testing
        num_samples = 10
        num_nodes = 5
        num_node_features = 3
        output_dim = 4
        data = self.create_dummy_data(
            num_samples, num_nodes, num_node_features, output_dim)
        # Save dummy data to the fake file system
        path_cases, path_mobility = self.setup_fake_filesystem(
            self.fs, self.path, data)
        # Create dataset
        dataset = create_dataset(
            path_cases, path_mobility, number_of_nodes=num_nodes)
        # Build the model by calling it on a batch of data
        loader = MixedLoader(dataset, batch_size=2, epochs=1)
        inputs, _ = loader.__next__()
        model(inputs)
        # Redefine the loader
        loader = MixedLoader(dataset, batch_size=2, epochs=1)
        res = evaluate(loader, model, loss_fn)
        # Check if the result is a tuple of (loss, accuracy)
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

    @patch("os.path.realpath", return_value="/home/")
    def test_train_and_evaluate(self, mock_realpath):

        number_of_epochs = 2
        # Create a simple model for testing
        model = get_model(
            layer_type="GCNConv", num_layers=2, num_channels=16,
            activation="relu", num_output=4)

        # Create dummy data in the fake filesystem for testing
        num_samples = 20
        num_nodes = 5
        num_node_features = 3
        output_dim = 4
        data = self.create_dummy_data(
            num_samples, num_nodes, num_node_features, output_dim)
        # Save dummy data to the fake file system
        path_cases, path_mobility = self.setup_fake_filesystem(
            self.fs, self.path, data)
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

        self.assertEqual(len(res["train_losses"][0]), number_of_epochs)
        self.assertEqual(len(res["val_losses"][0]), number_of_epochs)
        self.assertGreater(res["mean_test_loss"], 0)

        # Testing with saving the results
        res = train_and_evaluate(
            dataset,
            batch_size=2,
            epochs=number_of_epochs,
            model=model,
            loss_fn=MeanAbsolutePercentageError(),
            optimizer=tf.keras.optimizers.Adam(),
            es_patience=100,
            save_dir=self.path)
        save_results_path = os.path.join(self.path, "model_evaluations_paper")
        save_model_path = os.path.join(self.path, "saved_weights")
        self.assertTrue(os.path.exists(save_results_path))
        self.assertTrue(os.path.exists(save_model_path))

        file_path_df = save_results_path+"/model.csv"
        df = pd.read_csv(file_path_df)
        self.assertEqual(len(df), 1)
        for item in [
            "train_loss", "val_loss", "test_loss",
            "test_loss_orig", "training_time",
                "loss_history", "val_loss_history"]:
            self.assertIn(item, df.columns)

        file_path_model = save_model_path+"/model.pickle"
        with open(file_path_model, 'rb') as f:
            weights_loaded = pickle.load(f)
        weights = model.get_weights()
        for w1, w2 in zip(weights_loaded, weights):
            np.testing.assert_array_equal(w1, w2)
        # Clean up
        self.fs.remove_object(path_cases)
        self.fs.remove_object(os.path.join(
            path_mobility, "commuter_mobility_2022.txt"))
        self.fs.remove_object(path_mobility)
        self.fs.remove_object(save_results_path)
        self.fs.remove_object(save_model_path)

    def test_perform_grid_search(self):

        # Create dummy data in the fake filesystem for testing
        num_samples = 10
        num_nodes = 4
        num_node_features = 3
        output_dim = 4
        data = self.create_dummy_data(
            num_samples, num_nodes, num_node_features, output_dim)
        # Save dummy data to the fake file system
        path_cases, path_mobility = self.setup_fake_filesystem(
            self.fs, self.path, data)
        # Create dataset
        dataset = create_dataset(
            path_cases, path_mobility, number_of_nodes=num_nodes)

        # Define model parameters for grid search
        layers = ["GCNConv"]
        num_layers = [1]
        num_channels = [8]
        activations = ["relu"]
        parameter_grid = [(layer, n_layer, channel, activation)
                          for layer in layers for n_layer in num_layers
                          for channel in num_channels
                          for activation in activations]
        batch_size = 2
        es_patience = 5
        max_epochs = 2
        # Perform grid search with explicit save_dir to avoid os.path.realpath issues
        perform_grid_search(dataset, parameter_grid, self.path,
                            batch_size=batch_size, max_epochs=max_epochs,
                            es_patience=es_patience)

        # Check if the results file is created
        results_file = os.path.join(
            self.path, "saves", "grid_search_results.csv")
        self.assertTrue(os.path.exists(results_file))

        # Check if the results file has the expected number of rows
        df_results = pd.read_csv(results_file)
        self.assertEqual(len(df_results), len(parameter_grid))
        self.assertEqual(len(df_results.columns), 12)

    def test_scale_data_valid_data(self):
        """Test utils.scale_data with valid input and label data."""
        data = {
            # 10 samples, 1 day, 5 nodes, 8 groups
            "inputs": np.random.rand(10, 1, 8, 5),
            "labels": np.random.rand(10, 1, 8, 5)
        }

        scaled_inputs, scaled_labels = utils.scale_data(data, True)

        # Check that the scaled data is not equal to the original data
        assert not np.allclose(
            data["inputs"].transpose(0, 3, 1, 2), scaled_inputs)
        assert not np.allclose(
            data["labels"].transpose(0, 3, 1, 2), scaled_labels)

        # Check that the scaled data is log-transformed
        assert np.allclose(scaled_inputs, np.log1p(
            data["inputs"]).transpose(0, 3, 1, 2))
        assert np.allclose(scaled_labels, np.log1p(
            data["labels"]).transpose(0, 3, 1, 2))

    def test_scale_data_invalid_data(self):
        """Test utils.scale_data with invalid (non-numeric) data."""
        data = {
            "inputs": np.array([["a", "b"], ["c", "d"]]),  # Non-numeric data
            "labels": np.array([["e", "f"], ["g", "h"]])
        }

        with self.assertRaises(ValueError):
            utils.scale_data(data)


if __name__ == '__main__':
    unittest.main()
