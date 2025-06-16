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
from memilio.surrogatemodel.ode_secir_groups import network_architectures
from memilio.simulation.osecir import InfectionState
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf


def plot_compartment_prediction_model(
        inputs, labels, modeltype,  model=None,
        plot_compartment='InfectedSymptoms', max_subplots=8):
    """ Plot prediction of the model and label for one compartment. The average of all age groups is plotted.

    If model is none, we just plot the inputs and labels for the selected compartment without any predictions.

    :param inputs: test inputs for model prediction.
    :param labels: test labels.
    :param modeltype: type of model. Can be 'classic' or 'timeseries'
    :param model: trained model. (Default value = None) 
    :param plot_compartment:  (Default value = 'InfectedSymptoms')
    :param max_subplots: Number of the simulation runs to be plotted and compared against. (Default value = 8)
    :returns: No return 
    """
    num_groups = 6
    num_compartments = 8
    if modeltype == 'classic':
        # to get the input_width, we first subtract the damping date and the contact matrix entries.
        # Next, we divide by the number of age groups * the number of compartments
        input_width = int(
            (inputs.shape[1] - (1 + num_groups * num_groups)) / (num_groups * num_compartments))
    elif modeltype == 'timeseries':
        input_width = int(inputs.shape[1])
    else:
        ValueError("Modeltype "+modeltype + " not known.")

    label_width = int(labels.shape[1])

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

        if modeltype == 'classic':
            input_plot = input_array[:(
                input_width*num_groups*num_compartments)]
            input_plot = input_plot.reshape(
                input_width, num_groups*num_compartments)

            mean_per_day_input = []
            for i in input_plot:
                x = i[plot_compartment_index::8]
                mean_per_day_input.append(x.mean())

            plt.plot(
                np.arange(0, input_width),
                mean_per_day_input,
                label='Inputs', marker='.', zorder=-10)

        elif modeltype == 'timeseries':
            mean_per_day_input = []
            for i in input_array:
                # Inputs has the dimensions [num_runs, input_width, features].
                # The features consist of the compartment data, contact matrices and the damping day.
                # Here, we want to get the mean of the plot_compartment over all age groups. Therefore,
                # we subtract the damping day and the contact matrix entries.
                x = i[plot_compartment_index: inputs.shape[2] -
                      (1 + num_groups * num_groups):8]
                mean_per_day_input.append(x.mean())

            plt.plot(
                np.arange(0, input_width),
                mean_per_day_input,
                label='Inputs', marker='.', zorder=-10)

        mean_per_day = []
        for i in label_array:
            x = i[plot_compartment_index::8]
            mean_per_day.append(x.mean())
        plt.scatter(
            np.arange(input_width, input_width + label_width),
            mean_per_day,
            edgecolors='k', label='Labels', c='#2ca02c', s=64)

        if model is not None:
            input_series = tf.expand_dims(inputs[n], axis=0)
            pred = model(input_series)
            pred = pred.numpy()
            pred = pred.reshape((30, 48))

            mean_per_day_pred = []
            for i in pred:
                x = i[plot_compartment_index::8]
                mean_per_day_pred.append(x.mean())

            plt.scatter(np.arange(input_width, input_width+pred.shape[-2]),
                        # pred[0, :, plot_compartment_index],
                        mean_per_day_pred,
                        marker='X', edgecolors='k', label='Predictions',
                        c='#ff7f0e', s=64)

    plt.xlabel('days')
    plt.legend()
    if os.path.isdir("plots") == False:
        os.mkdir("plots")
    plt.savefig('plots/evaluation_secir_groups_' + plot_compartment + '.png')


####################
# Helper functions #
####################
def calc_split_index(n, split_train=0.7,
                     split_valid=0.2, split_test=0.1):
    """
    Calculating the indixes for a split_train:split_valid:split_test decomposition of a set with size n 

    It must hold split_train + split_valid + split_test = 1

    :param n: integer value 
    :param split_train: value between 0 and 1
    :param split_valid: value between 0 and 1
    :param split_test: value between 0 and 1
    :returns: a list of the form [i_train, i_valid, i_test]
    """
    if split_train + split_valid + split_test > 1 + 1e-10:
        raise ValueError(
            "Summed data set shares are greater than 1. Please adjust the values.")
    n_train = int(n * split_train)
    n_valid = int(n * split_valid)
    n_test = n - n_train - n_valid

    return [n_train, n_valid, n_test]


def flat_input(input):
    """ Flatten input dimension

    :param input: input array of size (n,k,l)
    :returns: reshaped array of size (n, k*l)

    """
    dim = tf.reduce_prod(tf.shape(input)[1:])
    return tf.reshape(input, [-1, dim])


def prepare_data_classic(data):
    """
    Transforming data to be processable by "classic" network, simply flattening and concatenating for each data instance.

    :param data: dictionary produced by data_generation
    :returns: dictionary with entries {
        "train_inputs", "train_labels", "valid_inputs",
        "valid_labels", "test_inputs", "test_labels"}
    """
    # Getting number of samples
    n = data["inputs"].shape[0]

    # Calculate split inidces
    split_indices = calc_split_index(n)

    # Flattening all inputs
    inputs = flat_input(data["inputs"])
    labels = data["labels"]
    cmatrices = flat_input(data["contact_matrices"])
    dampdays = flat_input(data["damping_days"])

    aggregated_inputs = tf.concat(
        [tf.cast(inputs, tf.float32),
         tf.cast(cmatrices, tf.float32),
         tf.cast(dampdays, tf.float32)],
        axis=1, name='concat')

    # Splitting the data
    labels_train, labels_valid, labels_test = tf.split(
        labels, split_indices, 0)
    inputs_train, inputs_valid, inputs_test = tf.split(
        aggregated_inputs, split_indices, 0)

    return {
        "train_inputs": inputs_train,
        "train_labels": labels_train,
        "valid_inputs": inputs_valid,
        "valid_labels": labels_valid,
        "test_inputs":  inputs_test,
        "test_labels":  labels_test
    }


def prod_time_series(obj, n, length_input):
    """
    Repeating static informations to fit into a time series framework 

    :param obj: an array of objects of shape (n, shape_rest), which should be repeated
    :param n: total number of samples 
    :param length_input: number of days observed per input 
    :returns: a tensor of shape [n, length_input, -1], where for each sample the static object is repeated length_input times
    """
    new_obj = []
    for i in obj:
        new_obj.extend([i for _ in range(length_input)])

    new_obj = tf.reshape(
        tf.stack(new_obj),
        [n, length_input, -1])
    return new_obj


def prepare_data_timeseries(data):
    """
    Transforming data to be processable by "time_series" network, simply repeating static values, flattening and concatenating for each data instance.

    :param data: dictionary produces by data_generation
    :returns: dictionary with entries {
        "train_inputs", "train_labels", "valid_inputs",
        "valid_labels", "test_inputs", "test_labels"
    }
    """
    # Getting the number of samples
    n = data["inputs"].shape[0]

    # number of days per input sample
    input_width = data["inputs"][0].shape[0]

    # Reshaping the matrix input
    cmatrices = flat_input(tf.stack(data["contact_matrices"]))
    dampdays = flat_input(tf.stack(data["damping_days"]))

    # Repeat data (contact matrix and damping day) to produce time series
    cmatrices_repeated = prod_time_series(cmatrices, n, input_width)
    dampdays_repeated = prod_time_series(dampdays, n, input_width)

    # Calculate split indices
    split_indices = calc_split_index(n)

    # Splitting the data
    compinputs_train, compinputs_valid, compinputs_test = tf.split(
        data["inputs"], split_indices, 0)
    labels_train, labels_valid, labels_test = tf.split(
        data["labels"], split_indices, 0)
    cmatrices_train, cmatrices_valid, cmatrices_test = tf.split(
        cmatrices_repeated, split_indices, 0)
    dampdays_train, dampdays_valid, dampdays_test = tf.split(
        dampdays_repeated, split_indices, 0)

    # Combining the ingredients to one input object
    inputs_train = tf.concat(
        [tf.cast(compinputs_train, tf.float32),
         tf.cast(cmatrices_train, tf.float32),
         tf.cast(dampdays_train, tf.float32)],
        axis=2, name='concat')
    inputs_valid = tf.concat(
        [tf.cast(compinputs_valid, tf.float32),
         tf.cast(cmatrices_valid, tf.float32),
         tf.cast(dampdays_valid, tf.float32)],
        axis=2, name='concat')
    inputs_test = tf.concat(
        [tf.cast(compinputs_test, tf.float32),
         tf.cast(cmatrices_test, tf.float32),
         tf.cast(dampdays_test, tf.float32)],
        axis=2, name='concat')

    return {
        "train_inputs": inputs_train,
        "train_labels": labels_train,
        "valid_inputs": inputs_valid,
        "valid_labels": labels_valid,
        "test_inputs":  inputs_test,
        "test_labels":  labels_test
    }


#########################################
# Initialization and Training of Models #
#########################################

def initialize_model(parameters):
    """ Initialize model from given list of parameters 

    :param parameters: tuple of parameters describing the model architecture, it should be of the form 
                (number_of_output_days, number_age_groups, number_compartments, 
                        hidden_layers, neurons_in_hidden_layer, activation_function, modelname)
    :returns: tensor flow keras model with the given architecture, see network_architectures.py 
    """
    label_width, number_age_groups, number_compartments, hidden_layers, neurons_in_hidden_layer, activation_function, modelname = parameters

    if modelname == "Dense":
        return network_architectures.mlp_multi_input_multi_output(
            label_width=label_width,
            num_age_groups=number_age_groups,
            num_outputs=number_compartments,
            num_hidden_layers=hidden_layers,
            num_neurons_per_layer=neurons_in_hidden_layer,
            activation=activation_function
        )
    elif modelname == "LSTM":
        return network_architectures.lstm_multi_input_multi_output(
            label_width=label_width,
            num_age_groups=number_age_groups,
            num_outputs=number_compartments,
            num_hidden_layers=hidden_layers,
            num_neurons_per_layer=neurons_in_hidden_layer,
            activation=activation_function
        )
    elif modelname == "CNN":
        return network_architectures.cnn_multi_input_multi_output(
            label_width=label_width,
            num_age_groups=number_age_groups,
            num_outputs=number_compartments,
            num_hidden_layers=hidden_layers,
            num_neurons_per_layer=neurons_in_hidden_layer,
            activation=activation_function
        )
    else:
        raise ValueError(
            "name_architecture must be one of 'Dense', 'LSTM' or 'CNN'"
        )


def network_fit(
        model, modeltype, training_parameter, path, filename='data_secir_groups_30days_10k_active.pickle', plot_stats=True):
    """ Training and evaluation of a given model with mean squared error loss and Adam optimizer using the mean absolute error as a metric.

    :param model: Keras sequential model.
    :param modeltype: type of model. Can be 'classic' or 'timeseries'. Data preparation is made based on the modeltype.
    :param training_parameter: tuple of parameters used for the training process, it should be of the form
        (early_stop, max_epochs, loss, optimizer, metrics), where loss is a loss-function implemented in keras, optimizer is the name of the used optimizer, 
        metrics is a list of used training metrics, e.g. [tf.keras.metrics.MeanAbsoluteError(), tf.keras.metrics.MeanAbsolutePercentageError()]
    :param path: path of the dataset.
    :param filename: name of the file containing the data 
    :param plot_stats:  (Default value = True)
    :returns: training history as returned by the keras fit() method. 
    """
    # Unpacking training parameters
    early_stop, max_epochs, loss, optimizer, metrics = training_parameter

    # Getting data and loading it
    if not os.path.isfile(os.path.join(path, filename)):
        ValueError("no dataset found in path: " + path)
    file = open(os.path.join(
        path, filename), 'rb')

    data = pickle.load(file)

    # preprocessing the data
    if modeltype == 'classic':
        data_prep = prepare_data_classic(data)

    elif modeltype == 'timeseries':
        data_prep = prepare_data_timeseries(data)

    else:
        raise ValueError("modeltype must be either classic or timeseries!")

    # Setting up the training parameters
    batch_size = 32
    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')
    model.compile(
        loss=loss,
        optimizer=optimizer,
        metrics=metrics)

    history = model.fit(data_prep["train_inputs"],
                        data_prep["train_labels"],
                        epochs=max_epochs,
                        validation_data=(
                            data_prep["valid_inputs"],
                            data_prep["valid_labels"]),
                        batch_size=batch_size,
                        callbacks=[early_stopping])

    if (plot_stats):
        plot_losses(history)
        plot_compartment_prediction_model(
            data_prep["test_inputs"], data_prep["test_labels"], modeltype, model=model,
            plot_compartment='InfectedSymptoms', max_subplots=3)
        df = get_test_statistic(
            data_prep["test_inputs"], data_prep["test_labels"], model)
        print(df)
    return history


def save_model(model, path, modelname):
    """
    Saving a trained model. 

    :param model: trained tensorflow keras model 
    :param path: path where the model should be stored 
    :param modelname: the name of the model 
    """
    if not os.path.isdir(path):
        os.mkdir(path)
    path_to_file = os.path.join(path, modelname + ".keras")
    model.save(path_to_file)
    print("Model successfully saved")


def load_model(path):
    """
    Loading a trained model. 

    :param path: path to the .keras file containing the desired model
    :returns: trained tf.keras model 
    """
    if not os.path.isfile(path):
        raise FileExistsError(
            "There is no .keras model stored at the given directory.")
    return tf.keras.models.load_model(path)


#####################
# Plots etc.
#####################

def plot_losses(history):
    """ Plots the losses of the model training.

    :param history: model training history.
    :returns: No return 
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
    :returns: dataframe containing the MAPE for the different compartments  
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
    compartments = [x for x in compartments if x !=
                    'InfectedNoSymptomsConfirmed' and x != 'InfectedSymptomsConfirmed']
    mean_percentage = pd.DataFrame(
        data=relative_err_means_percentage,
        index=compartments,
        columns=['Percentage Error'])

    return mean_percentage


def get_input_dim_lstm(path_to_file):
    """ Extract the dimension of the input data

    :param path_to_file: path to the data

    """
    file = open(path_to_file, 'rb')

    data = pickle.load(file)
    input_dim = data['inputs'].shape[2] + np.asarray(
        data['contact_matrix']).shape[1] * np.asarray(data['contact_matrix']).shape[2]+1

    return input_dim


if __name__ == "__main__":
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')

    # Defining model parameters
    label_width = 30
    number_age_groups = 6
    number_compartments = 8
    hidden_layers = 3
    neurons_in_hidden_layer = 512
    activation_function = 'relu'
    modelname = "Dense"
    modeltype = "timeseries"  # or "classic"

    model_parameters = (label_width, number_age_groups, number_compartments,
                        hidden_layers, neurons_in_hidden_layer, activation_function, modelname)

    # Defining training parameters
    early_stop = 100
    max_epochs = 200
    loss = tf.keras.losses.MeanAbsolutePercentageError()
    optimizer = "AdamW"
    metrics = [tf.keras.metrics.MeanAbsoluteError()]
    training_parameters = (early_stop, max_epochs, loss, optimizer, metrics)

    # input_dim = get_input_dim_lstm(path_data) -> Warum?
    model = initialize_model(model_parameters)

    model_output = network_fit(
        model=model, modeltype=modeltype,
        training_parameter=training_parameters, path=path_data)
