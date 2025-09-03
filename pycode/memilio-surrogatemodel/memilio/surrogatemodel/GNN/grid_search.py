import os
import pickle

import pandas as pd
import numpy as np
from sklearn.preprocessing import FunctionTransformer
import tensorflow as tf
from tensorflow.keras.layers import Dense
from tensorflow.keras.losses import MeanAbsolutePercentageError
from tensorflow.keras.metrics import mean_absolute_percentage_error
from tensorflow.keras.models import Model
import tensorflow.keras.initializers as initializers

from tensorflow.keras.optimizers import Adam, Nadam, RMSprop, SGD, Adagrad
import spektral.layers as spektral_layers
import spektral.utils.convolution as spektral_convolution
from spektral.data import MixedLoader

from sklearn.model_selection import KFold

from memilio.surrogatemodel.GNN.network_architectures import get_model
from memilio.surrogatemodel.GNN.evaluate_and_train import (
    train_and_evaluate, create_dataset)

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

output_dim = 1440

batch_size = 32
loss_function = MeanAbsolutePercentageError()
optimizer = Adam()
es_patience = 30
max_epochs = 100

model_parameters = [(layer, num_layers, num_channels, activation) for layer in layer_types
                    for num_layers in numbers_layers
                    for num_channels in numbers_channels
                    for activation in activation_functions]


def perform_grid_search(model_parameters, data):
    """Perform grid search over the specified model parameters. 
    Trains and evaluates models with different configurations and stores the results in a CSV file.
    """
    # Create a DataFrame to store the results
    df_results = pd.DataFrame(columns=['model', 'optimizer', 'number_of_hidden_layers', 'number_of_channels', 'activation',
                                       'mean_train_loss',
                                       'mean_validation_loss', 'training_time', 'train_losses', 'val_losses'])

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
            metrics=[mean_absolute_percentage_error]
        )

        results = train_and_evaluate(
            data, batch_size, epochs=max_epochs, model=model,
            loss_fn=loss_function, optimizer=optimizer,
            es_patience=es_patience, save_name="")

        df_results.loc[len(df_results.index)] = [
            layer_type, optimizer.__class__.__name__, num_layers, num_channels, activation,
            results["mean_train_loss"],
            results["mean_val_loss"],
            results["training_time"],
            results["train_losses"],
            results["val_losses"]
        ]

        path = os.path.dirname(os.path.realpath(__file__))

        file_path_df = os.path.join(os.path.join(
            os.path.dirname(
                os.path.realpath(os.path.dirname(os.path.realpath(path)))), 'saves'))
        if not os.path.isdir(file_path_df):
            os.mkdir(file_path_df)

        df_results.to_csv(os.path.join(
            file_path_df, 'grid_search_results.csv'))
        print(f"Saved intermediate results to {file_path_df}")
        tf.keras.backend.clear_session()


if __name__ == "__main__":
    # Generate the Dataset
    path_cases = "/localdata1/hege_mn/memilio/saves/GNN_data_30days_3dampings_classic5.pickle"
    path_mobility = '/localdata1/hege_mn/memilio/data/Germany/mobility'
    data = create_dataset(path_cases, path_mobility)
    perform_grid_search(model_parameters, data)
