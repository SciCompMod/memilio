#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

# from memilio.surrogatemodel.ode_secir_groups import network_architectures
# from memilio.simulation.osecir import InfectionState
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf


def network_fit(
        path, model, modeltype, filename, modelname,  max_epochs=30, early_stop=100, plot=True):
    """! Training and evaluation of a given model with mean squared error loss and Adam optimizer using the mean absolute error as a metric.

    @param path path of the dataset. 
    @param model Keras sequential model.
    @param modeltype type of model. Can be 'classic' or 'timeseries'. Data preparation is made based on the modeltype.
    @param max_epochs int maximum number of epochs in training. 
    @param early_stop Integer that forces an early stop of training if the given number of epochs does not give a significant reduction of validation loss. 

    """

    file = open(path, 'rb')
    data = pickle.load(file)
    data_splitted = split_data(data['inputs'], data['labels'])

    if modeltype == 'classic':

        train_inputs_compartments = flat_input(data_splitted["train_inputs"])
        train_labels = (data_splitted["train_labels"])
        valid_inputs_compartments = flat_input(data_splitted["valid_inputs"])
        valid_labels = (data_splitted["valid_labels"])
        test_inputs_compartments = flat_input(data_splitted["test_inputs"])
        test_labels = (data_splitted["test_labels"])

        contact_matrices = split_contact_matrices(
            tf.stack(data["contact_matrix"]))
        contact_matrices_train = flat_input(contact_matrices['train'])
        contact_matrices_valid = flat_input(contact_matrices['valid'])
        contact_matrices_test = flat_input(contact_matrices['test'])

        damping_days = data['damping_day']
        damping_days_splitted = split_damping_days(damping_days)
        damping_days_train = damping_days_splitted['train']
        damping_days_valid = damping_days_splitted['valid']
        damping_days_test = damping_days_splitted['test']

        train_inputs = tf.concat(
            [tf.cast(train_inputs_compartments, tf.float32),
             tf.cast(contact_matrices_train, tf.float32),
             tf.cast(damping_days_train, tf.float32)],
            axis=1, name='concat')
        valid_inputs = tf.concat(
            [tf.cast(valid_inputs_compartments, tf.float32),
             tf.cast(contact_matrices_valid, tf.float32),
             tf.cast(damping_days_valid, tf.float32)],
            axis=1, name='concat')
        test_inputs = tf.concat(
            [tf.cast(test_inputs_compartments, tf.float32),
             tf.cast(contact_matrices_test, tf.float32),
             tf.cast(damping_days_test, tf.float32)],
            axis=1, name='concat')

    elif modeltype == 'timeseries':

        train_inputs_compartments = (data_splitted["train_inputs"])
        train_labels = (data_splitted["train_labels"])
        valid_inputs_compartments = (data_splitted["valid_inputs"])
        valid_labels = (data_splitted["valid_labels"])
        test_inputs_compartments = (data_splitted["test_inputs"])
        test_labels = (data_splitted["test_labels"])

        contact_matrices = split_contact_matrices(
            tf.stack(data["contact_matrix"]))
        contact_matrices_train = flat_input(contact_matrices['train'])
        contact_matrices_valid = flat_input(contact_matrices['valid'])
        contact_matrices_test = flat_input(contact_matrices['test'])

        # damping days
        n = np.array(data['damping_day']).shape[0]
        train_days = data['damping_day'][:int(n*0.64)]
        valid_days = data['damping_day'][int(n*0.64):int(n*0.8)]
        test_days = data['damping_day'][int(n*0.8):]

        # damping coefficient

        n = np.array(data['damping_coeff']).shape[0]
        train_coeff = np.asarray(data['damping_coeff']).transpose()[
            0][0][:int(n*0.64)]
        valid_coeff = np.asarray(data['damping_coeff']).transpose()[
            0][0][int(n*0.64):int(n*0.8)]
        test_coeff = np.asarray(data['damping_coeff']).transpose()[
            0][0][int(n*0.8):]

        # concatenate the compartment data with contact matrices and damping days
        # to receive complete input data
        new_contact_train = []
        for i in contact_matrices_train:
            new_contact_train.extend([i for j in range(5)])

        new_contact_train = tf.reshape(
            tf.stack(new_contact_train),
            [train_inputs_compartments.shape[0],
             5, np.asarray(new_contact_train).shape[1]])

        new_damping_days_train = []
        for i in train_days:
            new_damping_days_train.extend([i for j in range(5)])
        new_damping_days_train = tf.reshape(
            tf.stack(new_damping_days_train),
            [train_inputs_compartments.shape[0],
             5, 1])

        new_damping_coeffs_train = []
        for i in train_coeff:
            new_damping_coeffs_train.extend([i for j in range(5)])
        new_damping_coeffs_train = tf.reshape(
            tf.stack(new_damping_coeffs_train),
            [train_inputs_compartments.shape[0],
             5, 1])

        train_inputs = tf.concat(
            (tf.cast(train_inputs_compartments, tf.float16),
             tf.cast(new_contact_train, tf.float16),
             tf.cast(new_damping_days_train, tf.float16),
             tf.cast(new_damping_coeffs_train, tf.float16)),
            axis=2)

        new_contact_test = []
        for i in contact_matrices_test:
            new_contact_test.extend([i for j in range(5)])

        new_contact_test = tf.reshape(tf.stack(new_contact_test), [
            contact_matrices_test.shape[0], 5, contact_matrices_test.shape[1]])

        new_damping_days_test = []
        for i in test_days:
            new_damping_days_test.extend([i for j in range(5)])
        new_damping_days_test = tf.reshape(
            tf.stack(new_damping_days_test),
            [test_inputs_compartments.shape[0],
             5, 1])

        new_damping_coeff_test = []
        for i in test_coeff:
            new_damping_coeff_test.extend([i for j in range(5)])
        new_damping_coeff_test = tf.reshape(
            tf.stack(new_damping_coeff_test),
            [test_inputs_compartments.shape[0],
             5, 1])

        test_inputs = tf.concat(
            (tf.cast(test_inputs_compartments, tf.float16),
             tf.cast(new_contact_test, tf.float16),
             tf.cast(new_damping_days_test, tf.float16),
             tf.cast(new_damping_coeff_test, tf.float16)),
            axis=2)

        new_contact_val = []
        for i in contact_matrices_valid:
            new_contact_val.extend([i for j in range(5)])

        new_contact_val = tf.reshape(
            tf.stack(new_contact_val),
            [contact_matrices_valid.shape[0],
             5, contact_matrices_valid.shape[1]])

        new_damping_days_valid = []
        for i in valid_days:
            new_damping_days_valid.extend([i for j in range(5)])
        new_damping_days_valid = tf.reshape(
            tf.stack(new_damping_days_valid),
            [valid_inputs_compartments.shape[0],
             5, 1])

        new_damping_coeff_valid = []
        for i in valid_coeff:
            new_damping_coeff_valid.extend([i for j in range(5)])
        new_damping_coeff_valid = tf.reshape(
            tf.stack(new_damping_coeff_valid),
            [valid_inputs_compartments.shape[0],
             5, 1])

        valid_inputs = tf.concat(
            (tf.cast(valid_inputs_compartments, tf.float16),
             tf.cast(new_contact_val, tf.float16),
             tf.cast(new_damping_days_valid, tf.float16),
             tf.cast(new_damping_coeff_valid, tf.float16)),
            axis=2)

    batch_size = 32

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')

    model.compile(
        loss=tf.keras.losses.MeanAbsolutePercentageError(),
        optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
        metrics=[tf.keras.metrics.MeanSquaredError(), tf.keras.metrics.MeanAbsolutePercentageError()])

    history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                        validation_data=(valid_inputs, valid_labels),
                        batch_size=batch_size,
                        callbacks=[early_stopping])

    path = os.path.dirname(os.path.realpath(__file__))
    path_models = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'saved_models_secir_groups_paper')
    if not os.path.isdir(path_models):
        os.mkdir(path_models)
    # path_models = '/localdata1/gnn_paper_2024/data_Iteration2/results/saved_models/with_age_groups'
    path_models = '/hpc_data/schm_a45/data_paper/'
    model.save(os.path.join(path_models, modelname))

    if (plot):
        plot_losses(history)
        df = get_test_statistic(test_inputs, test_labels, model, filename)
        print(df)
    return history


def plot_losses(history):
    """! Plots the losses of the model training.  

    @param history model training history. 

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


def get_test_statistic(test_inputs, test_labels, model, filename_df):
    """! Calculates the mean absolute percentage error based on the test dataset.   

    @param test_inputs inputs from test data.
    @param test_labels labels (output) from test data.
    @param model trained model. 

    """

    pred = model(test_inputs)
    pred = pred.numpy()
    test_labels = np.array(test_labels)

    diff = pred - test_labels
    relative_err = (abs(diff))/abs(test_labels)
    # reshape [batch, time, features] -> [features, time * batch]
    relative_err_transformed = relative_err.transpose(2, 0, 1).reshape(8, -1)
    relative_err_means_percentage = relative_err_transformed.mean(axis=1) * 100

    # same evaluation for rescaled data:
    pred = np.expm1(pred)
    test_labels = np.expm1(np.array(test_labels))

    diff_rescaled = pred - test_labels
    relative_err_rescaled = (abs(diff_rescaled))/abs(test_labels)
    # reshape [batch, time, features] -> [features, time * batch]
    relative_err_transformed_rescaled = relative_err_rescaled.transpose(
        2, 0, 1).reshape(8, -1)
    relative_err_means_percentage_rescaled = relative_err_transformed_rescaled.mean(
        axis=1) * 100

    # mean_percentage = pd.DataFrame(
    #    data=relative_err_means_percentage,
    #    index=[x for i, x in enumerate([str(compartment).split(
    #        '.')[1] for compartment in InfectionState.values()]) if i not in (3, 5)],
    #    columns=['MAPE_Scaled'])
    # mean_percentage['MAPE_rescaled'] = relative_err_means_percentage_rescaled

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    mean_percentage = pd.DataFrame(
        data=relative_err_means_percentage,
        index=infectionstates,
        columns=['MAPE_Scaled'])
    mean_percentage['MAPE_rescaled'] = relative_err_means_percentage_rescaled

    print('MAPE scaled data: ', mean_percentage['MAPE_Scaled'].mean(), "%")
    print('MAPE rescaled data: ', mean_percentage['MAPE_rescaled'].mean(), "%")

    # save the results as csv

    path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'testing_paper')
    # file_path = '/localdata1/gnn_paper_2024/data_Iteration2/results/testing'
    file_path = '/hpc_data/schm_a45/data_paper/'
    if not os.path.isdir(file_path):
        os.mkdir(file_path)
    file_path = os.path.join(file_path, filename_df)
    mean_percentage.to_csv(file_path)
    return mean_percentage


def split_data(inputs, labels, split_train=0.64,
               split_valid=0.16, split_test=0.2):
    """! Split data set in training, validation and testing data sets.

   @param inputs input dataset
   @param labels label dataset
   @param split_train Share of training data sets.
   @param split_valid Share of validation data sets.
   @param split_test Share of testing data sets.
   """

    if split_train + split_valid + split_test > 1 + 1e-10:
        raise ValueError(
            "Summed data set shares are greater than 1. Please adjust the values.")
    elif inputs.shape[0] != labels.shape[0] or inputs.shape[2] != labels.shape[2]:
        raise ValueError(
            "Number of batches or features different for input and labels")

    n = inputs.shape[0]
    n_train = int(n * split_train)
    n_valid = int(n * split_valid)
    n_test = n - n_train - n_valid

    inputs_train, inputs_valid, inputs_test = tf.split(
        inputs, [n_train, n_valid, n_test], 0)
    labels_train, labels_valid, labels_test = tf.split(
        labels, [n_train, n_valid, n_test], 0)

    data = {
        'train_inputs': inputs_train,
        'train_labels': labels_train,
        'valid_inputs': inputs_valid,
        'valid_labels': labels_valid,
        'test_inputs': inputs_test,
        'test_labels': labels_test
    }

    return data


def flat_input(input):
    """! Flatten input dimension

   @param input input array

   """
    dim = tf.reduce_prod(tf.shape(input)[1:])
    return tf.reshape(input, [-1, dim])


def split_contact_matrices(contact_matrices, split_train=0.64,
                           split_valid=0.16, split_test=0.2):
    """! Split dampings in train, valid and test

   @param contact_matrices contact matrices
   @param labels label dataset
   @param split_train ratio of train datasets
   @param split_valid ratio of validation datasets
   @param split_test ratio of test datasets
   """

    if split_train + split_valid + split_test != 1:
        ValueError("summed Split ratios not equal 1! Please adjust the values")

    n = contact_matrices.shape[0]
    n_train = int(n * split_train)
    n_valid = int(n * split_valid)
    n_test = n - n_train - n_valid

    contact_matrices_train, contact_matrices_valid, contact_matrices_test = tf.split(
        contact_matrices, [n_train, n_valid, n_test], 0)
    data = {
        "train": contact_matrices_train,
        "valid": contact_matrices_valid,
        "test": contact_matrices_test
    }

    return data


def split_damping_days(damping_days, split_train=0.64,
                       split_valid=0.16, split_test=0.2):
    """! Split damping days in train, valid and test

   @param damping_days damping days
   @param split_train ratio of train datasets
   @param split_valid ratio of validation datasets
   @param split_test ratio of test datasets
   """

    if split_train + split_valid + split_test != 1:
        ValueError("summed Split ratios not equal 1! Please adjust the values")
    damping_days = np.asarray(damping_days)
    n = damping_days.shape[0]
    n_train = int(n * split_train)
    n_valid = int(n * split_valid)
    n_test = n - n_train - n_valid

    damping_days_train, damping_days_valid, damping_days_test = tf.split(
        damping_days, [n_train, n_valid, n_test], 0)
    data = {
        "train": tf.reshape(damping_days_train, [n_train, 1]),
        "valid": tf.reshape(damping_days_valid, [n_valid, 1]),
        "test": tf.reshape(damping_days_test, [n_test, 1])
    }

    return data


def get_input_dim_lstm(path):
    """! Extract the dimensiond of the input data

   @param path path to the data 

   """
    file = open(path, 'rb')

    data = pickle.load(file)
    input_dim = data['inputs'].shape[2] + np.asarray(
        data['contact_matrix']).shape[1] * np.asarray(data['contact_matrix']).shape[2]+1

    return input_dim


if __name__ == "__main__":
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data_paper')

    dataset_name = 'data_secir_groups_30days_Germany_I_based_10k_1damp.pickle'
    filename = 'data_secir_groups_30days_Germany_I_based_10k_onedamp.csv'
    modelname = 'LSTM_groups_30days_onedamp_I_based_10k.h5'
    path = os.path.join(path_data, dataset_name)
    path = '/localdata1/gnn_paper_2024/data_Iteration2/one_population/with_agegroups_Germany/1damp/data_secir_groups_30days_I_based_Germany_10k_onedamp.pickle'

    max_epochs = 1500
    label_width = 30

    input_dim = get_input_dim_lstm(path)

    def lstm_multi_input_multi_output(label_width, num_age_groups=6):
        """! LSTM Network which uses multiple time steps as input and returns the 8 compartments for
        multiple time steps in the future.

        Input and output have shape [number of expert model simulations, time points in simulation,
        number of individuals in infection states].

        @param label_width Number of time steps in the output.
        """
        model = tf.keras.Sequential([
            tf.keras.layers.LSTM(128, return_sequences=False),
            tf.keras.layers.Dense(label_width * 8 * num_age_groups,
                                  kernel_initializer=tf.initializers.zeros()),
            tf.keras.layers.Reshape([label_width, 8 * num_age_groups])
        ])
        return model

    def lstm_multi_input_multi_output_Ibased(label_width, num_age_groups=6):
        """! LSTM Network which uses multiple time steps as input and returns the 8 compartments for
        multiple time steps in the future.

        Input and output have shape [number of expert model simulations, time points in simulation,
        number of individuals in infection states].

        @param label_width Number of time steps in the output.
        """
        model = tf.keras.Sequential([
            tf.keras.layers.LSTM(1024, return_sequences=False),
            tf.keras.layers.Dense(units=1024, activation='elu'),
            tf.keras.layers.Dense(label_width * 8 * num_age_groups,
                                  kernel_initializer=tf.initializers.zeros()),
            tf.keras.layers.Reshape([label_width, 8 * num_age_groups])
        ])
        return model

    model = lstm_multi_input_multi_output_Ibased(
        label_width)
    modeltype = 'timeseries'

    model_output = network_fit(
        path, model=model, modeltype=modeltype, modelname=modelname, filename=filename,
        max_epochs=max_epochs)
