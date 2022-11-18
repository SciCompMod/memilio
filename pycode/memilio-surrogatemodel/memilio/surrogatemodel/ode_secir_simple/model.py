#############################################################################
# Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
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
import numpy as np
import pandas as pd
import os
import pickle
import tensorflow as tf
import matplotlib.pyplot as plt
from memilio.surrogatemodel.ode_secir_simple import data_generation
from memilio.surrogatemodel.ode_secir_simple import different_networks


def plotCol(
        inputs, labels, model=None, plot_col='InfectedSymptoms',
        max_subplots=8):
    """! Plot prediction of the model and label for one compartment.

    @param inputs test inputs for model prediction. 
    @param labels test labels. 
    @param model trained model. 
    @param plot_col string name of compartment to be plotted. 
    @param max_subplots number of plot to be plotted. 
    """

    input_width = inputs.shape[1]
    label_width = labels.shape[1]

    plt.figure(figsize=(12, 8))
    cols = np.array([
        'Susceptible', 'Exposed', 'InfectedNoSymptoms', 'InfectedSymptoms',
        'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead'])
    plot_col_index = np.where(cols == plot_col)[0][0]
    max_n = min(max_subplots, inputs.shape[0])

    for n in range(max_n):
        plt.subplot(max_n, 1, n+1)
        plt.ylabel(f'{plot_col}')

        input_array = inputs[n].numpy()
        label_array = labels[n].numpy()
        plt.plot(np.arange(0, input_width), input_array[:, plot_col_index],
                 label='Inputs', marker='.', zorder=-10)
        plt.scatter(
            np.arange(input_width, input_width + label_width),
            label_array[:, plot_col_index],
            edgecolors='k', label='Labels', c='#2ca02c', s=64)

        if model is not None:
            input_series = tf.expand_dims(inputs[n], axis=0)
            pred = model(input_series)
            pred = pred.numpy()
            plt.scatter(np.arange(input_width, input_width+pred.shape[-2]),
                        pred[0, :, plot_col_index],
                        marker='X', edgecolors='k', label='Predictions',
                        c='#ff7f0e', s=64)

    plt.xlabel('days')
    if os.path.isdir("plots") == False:
        os.mkdir("plots")
    plt.savefig('plots/evaluation_secir_simple_' + plot_col + '.pdf')


def network_fit(path, model, max_epochs=30, early_stop=500, plot=True):
    """! Training and evaluation of the model. 

    @param path path of the dataset. 
    @param model name of the model to be loaded from the file where the model architectures are saved. 
    @param max_epochs int maximum number of epochs in training. 
    @param early_stop int defines the number of epochs without significant change tolerated before forcing an early stop of training. 

    """

    if not os.path.isfile(os.path.join(path, 'data_secir_simple.pickle')):
        ValueError("no dataset found in path: " + path)

    file = open(os.path.join(path, 'data_secir_simple.pickle'), 'rb')

    data = pickle.load(file)
    data_splitted = data_generation.splitdata(data["inputs"], data["labels"])

    train_inputs = data_splitted["train_inputs"]
    train_labels = data_splitted["train_labels"]
    valid_inputs = data_splitted["valid_inputs"]
    valid_labels = data_splitted["valid_labels"]
    test_inputs = data_splitted["test_inputs"]
    test_labels = data_splitted["test_labels"]

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')

    model.compile(
        loss=tf.keras.losses.MeanSquaredError(),
        optimizer=tf.keras.optimizers.Adam(),
        metrics=[tf.keras.metrics.MeanAbsoluteError()])

    history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                        validation_data=(valid_inputs, valid_labels),
                        callbacks=[early_stopping])

    if(plot):
        plot_losses(history)
        plotCol(test_inputs, test_labels, model=model,
                plot_col='InfectedSymptoms', max_subplots=6)
        df = get_test_statistic(test_inputs, test_labels, model)
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
    plt.savefig('plots/losses_plot.pdf')
    plt.show()


def get_test_statistic(test_inputs, test_labels, model):
    """! Calculates the mean absolute percentage error based on the test dataset.   

    @param test_inputs inputs from test data.
    @param testlabels labels/outout from test data.
    @aram model trained model. 

    """

    pred = model(test_inputs)
    pred = pred.numpy()
    test_labels = np.array(test_labels)

    diff = pred - test_labels
    anteil = (abs(diff))/abs(test_labels)

    anteil_6 = anteil.transpose(2, 0, 1).reshape(6, -1)
    df = pd.DataFrame(data=anteil_6)
    df = df.transpose()

    mean_percentage = pd.DataFrame(
        data=(df.mean().values) * 100,
        index=['Exposed', 'Carrier', 'Infected', 'Hospitalized', 'ICU',
               'Dead'],
        columns=['Percentage Error'])

    return mean_percentage


if __name__ == "__main__":
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')
    max_epochs = 400

    model = "LSTM"
    if model == "Dense":
        model = different_networks.multilayer_multi_input()
    elif model == "LSTM":
        model = different_networks.lstm_multi_output(30)
    elif model == "CNN":
        model = different_networks.cnn_multi_output(30)

    model_output = network_fit(
        path_data, model=model,
        max_epochs=max_epochs)
