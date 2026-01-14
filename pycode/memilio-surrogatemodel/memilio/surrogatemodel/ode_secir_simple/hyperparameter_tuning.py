#############################################################################
# Copyright (C) 2020-2026 MEmilio
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

import os
import tensorflow as tf
import pickle
import pandas as pd
import time
from sklearn.model_selection import KFold
import numpy as np
import memilio.surrogatemodel.utils.grid_search as grid_search
import memilio.surrogatemodel.ode_secir_simple.model as md

# Setting random seed
np.random.seed(42)

label_width = 30
num_outputs = 8

# Name of the file containing training data
filename = "data_secir_simple_30days_10k.pickle"

# Saving grid search results
filename_df = "dataframe_optimizer_paper.csv"

# General grid search parameters for the training process:
early_stop = [100]
max_epochs = [200]
losses = [tf.keras.losses.MeanAbsolutePercentageError()]
optimizers = ['Adam', 'AdamW', 'Nadam', 'SGD', 'Adagrad', 'RMSProp']
metrics = [[tf.keras.metrics.MeanAbsoluteError(
), tf.keras.metrics.MeanAbsolutePercentageError()]]

# Define grid search parameters for the architecture
hidden_layers = [1, 2]
neurons_in_hidden_layer = [32]
activation_function = ['relu', 'elu', 'softmax',  'sigmoid', 'linear', 'tanh']
models = ["LSTM"]

# Collecting parameters
training_parameters = [(early, epochs, loss, optimizer, metric)
                       for early in early_stop for epochs in max_epochs for loss in losses
                       for optimizer in optimizers for metric in metrics]

model_parameters = [(label_width, num_outputs, layer, neuron_number, activation, modelname)
                    for layer in hidden_layers for neuron_number in neurons_in_hidden_layer
                    for activation in activation_function for modelname in models]

# Loading data
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'data')

if not os.path.isfile(os.path.join(path_data, filename)):
    raise ValueError(f"No dataset found in path: {path_data}")

with open(os.path.join(path_data, filename), 'rb') as file:
    data = pickle.load(file)

# Split the data: 80% will be used for grid search with cross-validation, and the remaining 20% is withheld for testing
data_splitted = md.split_data(data['inputs'], data['labels'], 0.8, 0, 0.2)
inputs_grid_search = data_splitted['train_inputs']
labels_grid_search = data_splitted['train_labels']
inputs_withhold = data_splitted['test_inputs']
labels_withhold = data_splitted['test_labels']


start_hyper = time.perf_counter()

# Performing grid search
grid_search.perform_grid_search(
    model_parameters, inputs_grid_search, labels_grid_search, training_parameters, filename_df)

elapsed_hyper = time.perf_counter() - start_hyper

print(
    "Time for hyperparameter testing: {:.4f} minutes".format(
        elapsed_hyper / 60))
print(
    "Time for hyperparameter testing: {:.4f} hours".format(
        elapsed_hyper / 60 / 60))
