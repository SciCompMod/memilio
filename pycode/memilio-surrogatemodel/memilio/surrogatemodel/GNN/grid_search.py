#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Agatha Schmidt, Henrik Zunker, Manuel Heger
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

import os
from typing import Optional

import pandas as pd
import tensorflow as tf
from tensorflow.keras.losses import MeanAbsolutePercentageError as MeanAbsolutePercentageError
from tensorflow.keras.metrics import MeanAbsolutePercentageError as mean_absolute_percentage_error

from tensorflow.keras.optimizers import Adam
from spektral.data import MixedLoader
from memilio.surrogatemodel.GNN.network_architectures import get_model
from memilio.surrogatemodel.GNN.evaluate_and_train import (
    train_and_evaluate, create_dataset)

# Define the parameter grid for grid search
layer_types = [
    "ARMAConv",
    "GCSConv",
    "GATConv",
    "GCNConv",
    "APPNPConv"
]

numbers_layers = [2, 3, 4, 5, 6, 7]
numbers_channels = [2, 3, 4, 5, 6, 7]
activation_functions = ["elu", "relu", "tanh", "sigmoid"]
model_parameters = [
    (layer, num_layers, num_channels, activation)
    for layer in layer_types
    for num_layers in numbers_layers
    for num_channels in numbers_channels for activation in
    activation_functions]

# Fix training parameters
batch_size = 32
loss_function = MeanAbsolutePercentageError()
optimizer = Adam()
es_patience = 30
max_epochs = 100

training_parameters = [batch_size, loss_function,
                       optimizer, es_patience, max_epochs]


def perform_grid_search(
        model_parameters, training_parameters, data,
        save_dir: Optional[str] = None):
    """Perform grid search over the specified model parameters. 
    Trains and evaluates models with different configurations and stores the results in a CSV file.
    Results are saved in 'grid_search_results.csv' in the 'saves' directory.
    :param model_parameters: List of tuples containing model parameters (layer_type, num_layers, num_channels, activation).
    :param training_parameters: List containing training parameters (batch_size, loss_function, optimizer, es_patience, max_epochs).
    :param data: The dataset to be used for training and evaluation.
    :return: None
    """
    # Setting hyprparameter
    output_dim = data[0].y.shape[-1]
    batch_size, loss_function, optimizer, es_patience, max_epochs = training_parameters
    # Create a DataFrame to store the results
    df_results = pd.DataFrame(
        columns=['model', 'optimizer',
                 'number_of_hidden_layers',
                 'number_of_channels', 'activation',
                                       'mean_train_loss',
                                       'mean_validation_loss', 'training_time',
                                       'train_losses', 'val_losses'])

    for param in model_parameters:
        layer_type, num_layers, num_channels, activation = param
        print(
            f"Training model with {layer_type}, {num_layers} layers, {num_channels} channels, activation: {activation}")

        # Create a model instance
        model = get_model(
            layer_type=layer_type,
            num_layers=num_layers,
            num_channels=num_channels,
            activation=activation,
            num_output=output_dim
        )
        optimizer = Adam()

        loader = MixedLoader(data, batch_size=batch_size, epochs=1)
        inputs, _ = loader.__next__()
        model(inputs)  # Build the model by calling it on a batch of data

        # Initialize optimizer variables
        optimizer.build(model.trainable_variables)
        model.compile(
            optimizer=optimizer,
            loss=MeanAbsolutePercentageError(),
            metrics=[mean_absolute_percentage_error()]
        )

        results = train_and_evaluate(
            data, batch_size, epochs=max_epochs, model=model,
            loss_fn=loss_function, optimizer=optimizer,
            es_patience=es_patience, save_name="", save_dir=save_dir)

        df_results.loc[len(df_results.index)] = [
            layer_type, optimizer.__class__.__name__, num_layers, num_channels, activation,
            results["mean_train_loss"],
            results["mean_val_loss"],
            results["training_time"],
            results["train_losses"],
            results["val_losses"]
        ]
        # Save intermediate results to avoid data loss
        # If save_dir is provided, use it. Otherwise store under a local 'saves' folder next to this module
        base_dir = save_dir if (
            save_dir and len(save_dir) > 0) else os.path.realpath(
            os.path.dirname(__file__))
        saves_dir = os.path.join(base_dir, 'saves')
        os.makedirs(saves_dir, exist_ok=True)
        df_results.to_csv(os.path.join(
            saves_dir, 'grid_search_results.csv'), index=False)
        print(f"Saved intermediate results to {saves_dir}")
        # Clear session to free memory after each iteration
        tf.keras.backend.clear_session()


if __name__ == "__main__":
    # Generate the Dataset
    path_cases = "/localdata1/hege_mn/memilio/saves/GNN_data_30days_3dampings_classic5.pickle"
    path_mobility = '/localdata1/hege_mn/memilio/data/Germany/mobility'
    data = create_dataset(path_cases, path_mobility)
    perform_grid_search(model_parameters, data)
