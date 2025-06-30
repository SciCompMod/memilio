#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Agatha Schmidt, Henrik Zunker
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
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf

from memilio.simulation.osecir import InfectionState
import memilio.surrogatemodel.ode_secir_simple.network_architectures as architectures
from memilio.surrogatemodel.surrogate_utils import (calc_split_index)


def initialize_model(parameters):
    """ Initialize model from given list of parameters

    :param parameters: tuple of parameters describing the model architecture, it should be of the form 
            (num_days_per_output, num_outputs, num_hidden_layers, neurons_per_layer, name_activation, name_architecture)
    :returns: tensor flow keras model with the given architecture 
     """
    label_width, num_outputs, layer, neuron_number, activation, modelname = parameters

    if modelname == "Dense":
        return architectures.mlp_multi_input_multi_output(
            label_width=label_width,
            num_outputs=num_outputs,
            num_hidden_layers=layer,
            num_neurons_per_layer=neuron_number,
            activation=activation
        )
    elif modelname == "LSTM":
        return architectures.lstm_multi_input_multi_output(
            label_width=label_width,
            num_outputs=num_outputs,
            num_hidden_layers=layer,
            num_neurons_per_layer=neuron_number,
            activation=activation
        )
    elif modelname == "CNN":
        return architectures.cnn_multi_input_multi_output(
            label_width=label_width,
            num_outputs=num_outputs,
            num_hidden_layers=layer,
            num_neurons_per_layer=neuron_number,
            activation=activation
        )
    else:
        raise ValueError(
            "name_architecture must be one of 'Dense', 'LSTM' or 'CNN'")


def network_fit(model, inputs, labels, training_parameter, plot=True):
    """ Training and evaluation of a given model under given training parameters.

    :param model: Keras sequential model.
    :param inputs: Input data used for training  
    :param labels: Output data (labels) used for training
    :param training_parameter: tuple of parameters used for the training process, it should be of the form
        (early_stop, max_epochs, loss, optimizer, metrics), where loss is a loss-function implemented in keras, optimizer is the name of the used optimizer, 
        metrics is a list of used training metrics, e.g. [tf.keras.metrics.MeanAbsoluteError(), tf.keras.metrics.MeanAbsolutePercentageError()]
    :param plot: Default: True] 
    :returns: trained model with training history (output of model.fit())
    """
    early_stop, max_epochs, loss, optimizer, metrics = training_parameter
    data_splitted = split_data(inputs, labels)

    train_inputs = data_splitted['train_inputs']
    train_labels = data_splitted['train_labels']
    valid_inputs = data_splitted['valid_inputs']
    valid_labels = data_splitted['valid_labels']
    test_inputs = data_splitted['test_inputs']
    test_labels = data_splitted['test_labels']

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')

    model.compile(
        loss=loss,
        optimizer=optimizer,
        metrics=metrics)

    history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                        validation_data=(valid_inputs, valid_labels),
                        callbacks=[early_stopping])

    if (plot):
        plot_losses(history)
        plot_compartment_prediction_model(
            test_inputs, test_labels, model=model,
            plot_compartment='InfectedSymptoms', max_subplots=3)
        df = get_test_statistic(test_inputs, test_labels, model)
        print(df)

    return history


def split_data(inputs, labels, split_train=0.7,
               split_valid=0.2, split_test=0.1):
    """ Split data set in training, validation and testing data sets.

    :param inputs: input dataset
    :param labels: label dataset
    :param split_train: Share of training data sets. (Default value = 0.7)
    :param split_valid: Share of validation data sets. (Default value = 0.2)
    :param split_test: Share of testing data sets. (Default value = 0.1)
    :returns: dictionary of the form {'train_inputs', 'train_labels', 'valid_inputs',
        'valid_labels', 'test_inputs', 'test_labels'}
    """

    if inputs.shape[0] != labels.shape[0] or inputs.shape[2] != labels.shape[2]:
        raise ValueError(
            "Number of batches or features different for input and labels")

    n = inputs.shape[0]
    splitting = calc_split_index(n, split_train, split_valid, split_test)

    inputs_train, inputs_valid, inputs_test = tf.split(
        inputs, splitting, 0)
    labels_train, labels_valid, labels_test = tf.split(
        labels, splitting, 0)

    data = {
        'train_inputs': inputs_train,
        'train_labels': labels_train,
        'valid_inputs': inputs_valid,
        'valid_labels': labels_valid,
        'test_inputs': inputs_test,
        'test_labels': labels_test
    }

    return data


def plot_compartment_prediction_model(
        inputs, labels, model=None, plot_compartment='InfectedSymptoms',
        max_subplots=8):
    """ Plot prediction of the model and label for one compartment.

    If model is none, we just plot the inputs and labels for the selected compartment without any predictions.
    Plots are saved in directory ./plots/evaluation_secir_simple_{compartment}.png

    :param inputs: test inputs for model prediction.
    :param labels: test labels.
    :param model: Default: None] trained model. 
    :param plot_col: string name of compartment to be plotted, 
    :param max_subplots: Default: 8] Number of the simulation runs to be plotted and compared against. 
    :param plot_compartment: possible compartment names "Susceptible", "Exposed",
        "InfectedNoSymptoms", "InfectedSevere", "InfectedCritical", "Recovered", "Dead" 
    """

    input_width = inputs.shape[1]
    label_width = labels.shape[1]

    plt.figure(figsize=(12, 8))
    plot_compartment_index = 0
    for compartment in InfectionState.values():
        if compartment.name == plot_compartment:
            break
        plot_compartment_index += 1
    if plot_compartment_index == len(InfectionState.values()):
        raise ValueError('Compartment name given could not be found.')
    max_n = min(max_subplots, inputs.shape[0])

    for n in range(max_n):
        plt.subplot(max_n, 1, n+1)
        plt.ylabel(plot_compartment)

        input_array = inputs[n].numpy()
        label_array = labels[n].numpy()
        plt.plot(
            np.arange(0, input_width),
            input_array[:, plot_compartment_index],
            label='Inputs', marker='.', zorder=-10)
        plt.scatter(
            np.arange(input_width, input_width + label_width),
            label_array[:, plot_compartment_index],
            edgecolors='k', label='Labels', c='#2ca02c', s=64)

        if model is not None:
            input_series = tf.expand_dims(inputs[n], axis=0)
            pred = model(input_series)
            pred = pred.numpy()
            plt.scatter(np.arange(input_width, input_width+pred.shape[-2]),
                        pred[0, :, plot_compartment_index],
                        marker='X', edgecolors='k', label='Predictions',
                        c='#ff7f0e', s=64)

    plt.xlabel('days')
    if os.path.isdir("plots") == False:
        os.mkdir("plots")
    plt.savefig('plots/evaluation_secir_simple_' + plot_compartment + '.png')


def plot_losses(history):
    """ Plotting the losses of the model training.

    :param history: model training history.

    """
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    if os.path.isdir("plots") == False:
        os.mkdir("plots")
    plt.savefig('plots/losses_plot.png')
    plt.show()


def get_test_statistic(test_inputs, test_labels, model):
    """ Calculates the mean absolute percentage error based on the test dataset.

    :param test_inputs: inputs from test data.
    :param test_labels: labels (output) from test data.
    :param model: trained model.
    :returns: Dataframe with the corresponding MAPEs of the different compartments 
    """

    pred = model(test_inputs)
    pred = pred.numpy()
    test_labels = np.array(test_labels)

    diff = pred - test_labels
    relative_err = (abs(diff))/abs(test_labels)
    # reshape [batch, time, features] -> [features, time * batch]
    relative_err_transformed = relative_err.transpose(2, 0, 1).reshape(8, -1)
    relative_err_means_percentage = relative_err_transformed.mean(axis=1) * 100
    compartments = [str(compartment).split('.')[1]
                    for compartment in InfectionState.values()]

    # Excluding combined compartments
    compartments.pop(5)
    compartments.pop(3)
    mean_percentage = pd.DataFrame(
        data=relative_err_means_percentage,
        index=[compartments],
        columns=['Percentage Error'])

    return mean_percentage


if __name__ == "__main__":
    label_width = 30
    num_outputs = 8
    # General parameters for the training process:
    early_stop = [100]
    max_epochs = [200]
    losses = [tf.keras.losses.MeanAbsolutePercentageError()]
    optimizers = ['AdamW']
    metrics = [[tf.keras.metrics.MeanAbsoluteError(),
                tf.keras.metrics.MeanAbsolutePercentageError()]]

    # Define parameters for the architecture
    hidden_layers = [3]
    neurons_in_hidden_layer = [32]
    activation_function = ['relu']
    models = ["LSTM"]

    # Collecting parameters
    training_parameters = [(early, epochs, loss, optimizer, metric)
                           for early in early_stop for epochs in max_epochs for loss in losses
                           for optimizer in optimizers for metric in metrics]

    model_parameters = [(label_width, num_outputs, layer, neuron_number, activation, modelname)
                        for layer in hidden_layers for neuron_number in neurons_in_hidden_layer
                        for activation in activation_function for modelname in models]

    # getting training data
    filename = "data_secir_simple_30days_10k.pickle"
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')

    if not os.path.isfile(os.path.join(path_data, filename)):
        raise ValueError(f"No dataset found in path: {path_data}")

    with open(os.path.join(path_data, filename), 'rb') as file:
        data = pickle.load(file)

    input_data = data['inputs']
    label_data = data['labels']
    model = initialize_model(model_parameters[0])
    model_output = network_fit(
        model, input_data, label_data, training_parameters[0])
    plot_compartment_prediction_model(input_data, label_data, model)
