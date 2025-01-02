#############################################################################
# Copyright (C) 2020-2024 MEmilio
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

# from memilio.simulation.osecir import InfectionState
# from memilio.surrogatemodel.ode_secir_groups import network_architectures


def network_fit(path_data, model, savename, filename_df,  max_epochs=30, early_stop=100, plot=True):
    """! Training and evaluation of a given model with mean squared error loss and Adam optimizer using the mean absolute error as a metric.

    @param path path of the dataset. 
    @param model Keras sequential model.
    @param max_epochs int maximum number of epochs in training. 
    @param early_stop Integer that forces an early stop of training if the given number of epochs does not give a significant reduction of validation loss. 

    """
    file = open(path_data, 'rb')
    data = pickle.load(file)

    # split data: 80% for training, 20% for testing
    # IMPORTANT: we need to have the same split as sor the grid search
    # --> the 20% test data have to be the 20% that we withhold from the grid
    # search wth cross validation
    train_inputs = data['inputs'][:int((0.8 * len(data['inputs'])))]
    train_labels = data['labels'][:int((0.8 * len(data['labels'])))]
    test_inputs = data['inputs'][int((0.8 * len(data['inputs']))):]
    test_labels = data['labels'][int((0.8 * len(data['labels']))):]

    # we need to split the training dataset into training and validation
    # the validation dataset is needed for model evalutaion during training
    # e.g. it is needed for early stopping
    # we take 20% from the training data as validation data

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')

    model.compile(
        loss=tf.keras.losses.MeanAbsolutePercentageError(),
        optimizer=tf.keras.optimizers.Adam(),
        metrics=[tf.keras.metrics.MeanAbsolutePercentageError()])

    history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                        validation_split=0.2,
                        callbacks=[early_stopping])

    # save the best model so we can use the weights afterward without needing
    # to train the model again

    path = os.path.dirname(os.path.realpath(__file__))
    path_models = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'saved_models_secir_groups_paper')
    # path_models = '/localdata1/gnn_paper_2024/data_Iteration2/results/saved_models/with_age_groups'
    path_models = 'hpc_data/schm_a45/data_paper/'
    if not os.path.isdir(path_models):
        os.mkdir(path_models)

    model.save(os.path.join(path_models, savename))

    if (plot):
        # plot_losses(history)
        df = get_test_statistic(test_inputs, test_labels, model, filename_df)
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

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']
    # mean_percentage = pd.DataFrame(
    #    data=relative_err_means_percentage,
    #    index=[x for i, x in enumerate([str(compartment).split(
    #        '.')[1] for compartment in InfectionState.values()]) if i not in (3, 5)],
    #    columns=['MAPE_Scaled'])
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
    path_models = '/hpc_data/schm_a45/data_paper/'

    if not os.path.isdir(file_path):
        os.mkdir(file_path)
    file_path = os.path.join(file_path, filename_df)
    mean_percentage.to_csv(file_path)
    return mean_percentage


if __name__ == "__main__":
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data_paper')

    # ataset_name = 'data_secir_groups_60days_Germany_10k_nodamp.pickle'
    # path_data = os.path.join(path_data, dataset_name)
    path_data = '/hpc_data/schm_a45/data_paper/data_with_agegroups/data_secir_groups_90days_I_based_Germany_100k_nodamp_2.pickle'
    max_epochs = 1500

    # model = network_architectures.lstm_multi_input_multi_output_Ibased(90)

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

    model = lstm_multi_input_multi_output_Ibased(90)
    savename = 'LSTM_NODAMP_90days_I_based_secirgroups_100k_2.h5'
    filename_df = 'LSTM_NODAMP_90days_I_based_secirgroups_100k_2.csv'
    model_output = network_fit(
        path_data, model=model, savename=savename, filename_df=filename_df,
        max_epochs=max_epochs)
