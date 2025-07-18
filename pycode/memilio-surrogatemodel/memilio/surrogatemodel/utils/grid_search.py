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

import os
import tensorflow as tf
import pickle

import pandas as pd
import time
from sklearn.model_selection import KFold
import numpy as np
import memilio.surrogatemodel.ode_secir_simple.network_architectures as architectures_simple
import memilio.surrogatemodel.ode_secir_groups.network_architectures as architectures_groups
import memilio.surrogatemodel.ode_secir_simple.model as md_simple
import memilio.surrogatemodel.ode_secir_groups.model as md_groups

# Function to train and evaluate the model using cross-validation


def train_and_evaluate_model(param, inputs, labels, training_parameter, show_results=False, modeltype="simple", n_kfold=5):
    """ Training and evaluating a model with given architecture using 5-fold cross validation, returning a dictionary with the main training statistics. 

    :param param: tuple of parameters describing the model architecture, it should be of the form 
            (num_days_per_output, num_outputs, num_hidden_layers, neurons_per_layer, name_activation, name_architecture)
    :param inputs: training inputs 
    :param labels: training output labels 
    :param training_parameter: tuple of parameters used for the training process, it should be of the form
        (early_stop, max_epochs, loss, optimizer, metrics), where loss is a loss-function implemented in keras, optimizer is the name of the used optimizer, 
        metrics is a list of used training metrics, e.g. [tf.keras.metrics.MeanAbsoluteError(), tf.keras.metrics.MeanAbsolutePercentageError()]
    :param show_results:  Boolean, whether or not the evaluation results are printed. 
    :param modeltype: String, specifying, which type of model is going to be tested. The possible values are "simple" - refering to the surrogate for the SECIR-model 
        without age_resolution, "groups" - for the age resolved SECIR model. 
    :param n_kfold: number of partizions used to cross-validate
    :returns: a dictionary of training statistics of the form 
        {"model", "activation","optimizer","mean_train_loss_kfold","mean_val_loss_kfold","training_time", "train_losses", "val_losses"}

    """
    # Unpacking parameters
    activation, modelname = param[-2:]
    early_stop, max_epochs, loss, optimizer, metrics = training_parameter
    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')

    # Preparing K-Fold Cross-Validation
    kf = KFold(n_splits=n_kfold)
    train_losses = []
    val_losses = []

    losses_history_all = []
    val_losses_history_all = []

    start = time.perf_counter()
    for train_idx, val_idx in kf.split(inputs):
        # Clearing any information about the previous model
        tf.keras.backend.clear_session()

        # Gather training and validation data based on the fold
        train_inputs = tf.gather(inputs, indices=train_idx)
        train_labels = tf.gather(labels, indices=train_idx)
        valid_inputs = tf.gather(inputs, indices=val_idx)
        valid_labels = tf.gather(labels, indices=val_idx)

        # Initializing model
        if modeltype == "simple":
            model = md_simple.initialize_model(param)
        elif modeltype == "groups":
            model = md_groups.initialize_model(param)
        else:
            raise ValueError(modeltype+" is not known.")
        # Compile the model
        model.compile(loss=loss,
                      optimizer=optimizer,
                      metrics=metrics)

        # Train the model
        history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                            validation_data=(valid_inputs, valid_labels),
                            callbacks=[early_stopping])

        train_losses.append(np.min(history.history['loss']))
        val_losses.append(np.min(history.history['val_loss']))
        losses_history_all.append(history.history['loss'])
        val_losses_history_all.append(history.history['val_loss'])

    elapsed = time.perf_counter() - start

    # Print out the results
    if show_results:
        print(f"Best train losses: {train_losses}")
        print(f"Best validation losses: {val_losses}")
        print("--------------------------------------------")
        print(f"K-Fold Train Score: {np.mean(train_losses)}")
        print(f"K-Fold Validation Score: {np.mean(val_losses)}")
        print(f"Time for training: {elapsed:.4f} seconds")
        print(f"Time for training: {elapsed / 60:.4f} minutes")

    # After cross-validation, we can test on the withhold dataset (outside of the loop) ?
    return {
        "model": modelname,
        "activation": activation,
        "optimizer": optimizer,
        "mean_train_loss_kfold": np.mean(train_losses),
        "mean_val_loss_kfold": np.mean(val_losses),
        "training_time": elapsed/60,
        "train_losses": [losses_history_all],
        "val_losses": [val_losses_history_all]
    }


def perform_grid_search(model_parameters, inputs, labels, training_parameters, filename_df, path=None,  modeltype="simple"):
    """ Performing grid search for a given set of model parameters

    The results are stored in directory 'secir_simple_grid_search', each row has the form 
    ['model', 'optimizer', 'number_of_hidden_layers', 'number_of_neurons', 'activation',
                                    'mean_test_MAPE', 'kfold_train', 'kfold_val',
                                    'kfold_test', 'training_time', 'train_losses', 'val_losses']

    :param model_parameters: List of tuples of model parameters, each entry should be a tuple of the form 
        (num_days_per_output, num_outputs, num_hidden_layers, neurons_per_layer, name_activation, name_architecture)
    :param inputs: training input data 
    :param labels: training label data 
    :param training_parameters: List of tuples of parameters used for the training process, each should be of the form
        (early_stop, max_epochs, loss, optimizer, metrics), where loss is a loss-function implemented in keras, optimizer is the name of the used optimizer, 
        metrics is a list of used training metrics, e.g. [tf.keras.metrics.MeanAbsoluteError(), tf.keras.metrics.MeanAbsolutePercentageError()]
    :param filename_df: String, giving name of the file, where the data is stored, actual filename is given by filename_df + ".pickle"
    :param path: String representing the path, where dataframe should be stored
    :param modeltype: String, specifying, which type of model is going to be tested. The possible values are "simple" - refering to the surrogate for the SECIR-model 
        without age_resolution, "groups" - for the age resolved SECIR model. 
    """
    # Create a DataFrame to store the results
    df_results = pd.DataFrame(columns=['model', 'optimizer', 'number_of_hidden_layers', 'number_of_neurons', 'activation',
                                       'mean_test_MAPE', 'kfold_train', 'kfold_val',
                                       'kfold_test', 'training_time', 'train_losses', 'val_losses'])

    # Iterate the different model architectures and save the training results
    for param in model_parameters:
        for training_parameter in training_parameters:
            layer, neuron_number, activation, modelname = param[-4:]
            results = train_and_evaluate_model(
                param, inputs, labels, training_parameter)
            df_results.loc[len(df_results.index)] = [
                # Placeholder for test score
                modelname, results["optimizer"], layer, neuron_number, activation, np.nan,
                results["mean_train_loss_kfold"],
                results["mean_val_loss_kfold"],
                np.nan,  # Placeholder for final test score
                results["training_time"],
                results["train_losses"],
                results["val_losses"]
            ]

    # Save the results in file
    folder_name = 'secir_' + modeltype + '_grid_search'

    if path is None:
        path = os.path.dirname(os.path.realpath(__file__))
        file_path = os.path.join(os.path.dirname(os.path.realpath(path)),
                                 folder_name)
    else:
        file_path = os.path.join(path, folder_name)

    if not os.path.isdir(file_path):
        os.mkdir(file_path)

    file_path = os.path.join(file_path, filename_df + ".pickle")

    with open(file_path, "wb") as f:
        pickle.dump(df_results, f)
