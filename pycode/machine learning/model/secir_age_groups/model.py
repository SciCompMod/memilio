import numpy as np
import pandas as pd
import os
import pickle
import tensorflow as tf
import matplotlib.pyplot as plt
from data_generator import splitdata, split_contact_matrices, flat_input
from different_networks import *


def plotCol(inputs, labels, model=None, plot_col='Infected', max_subplots=3):
    """! Plot prediction of the model and label for one compartment.

    @param inputs test inputs for model prediction. 
    @param labels test labels. 
    @param model trained model. 
    @param plot_col string name of compartment to be plotted. 
    @param max_subplots number of plot to be plotted. 
    """

    # 36 damping entries  and 48 age dependent compartments
    input_width = int((inputs.shape[1] - 36) / 48)
    #label_width = int(labels.shape[1] / 48)
    label_width = int(labels.shape[1] / 36)

    plt.figure(figsize=(12, 8))
    # cols = np.array(['Susceptible', 'Exposed', 'Carrier',
    #                'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead'])
    cols = np.array(['Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU', 'Dead'])

    plot_col_index = np.where(cols == plot_col)[0][0]
    max_n = min(max_subplots, inputs.shape[0])

    # predictions = model(inputs) # for just one input: input_series = tf.expand_dims(inputs[n], axis=0) -> model(input_series)
    for n in range(max_n):
        plt.subplot(max_n, 1, n+1)
        plt.ylabel(f'{plot_col}')

        input_array = inputs[n].numpy()
        label_array = labels[n].numpy()
        plt.plot(
            np.arange(0, input_width),
            input_array[plot_col_index: inputs.shape[1] - 36: 48],
            label='Inputs', marker='.', zorder=-10)
        plt.scatter(
            np.arange(input_width, input_width + label_width),
            #label_array[plot_col_index: -1: 48],
            label_array[plot_col_index: -1: 36],
            edgecolors='k', label='Labels', c='#2ca02c', s=64)

        if model is not None:
            input_series = tf.expand_dims(inputs[n], axis=0)
            pred = model(input_series)
            pred = pred[0].numpy()
            plt.scatter(np.arange(input_width, input_width+label_width),
                        # pred[plot_col_index:-1:48],
                        pred[plot_col_index:-1:36],
                        marker='X', edgecolors='k', label='Predictions',
                        c='#ff7f0e', s=64)

    plt.xlabel('days')
    plt.show()
    plt.savefig('evaluation_secir_simple_' + plot_col + '.pdf')


def network_fit(path, model, max_epochs=30, early_stop=3000):
    """! Training and evaluation of the model. 

    @param path path of the dataset. 
    @param model name of the model to be loaded from the file where the model architectures are saved. 
    @param max_epochs int maximum number of epochs in training. 
    @param early_stop int defines the number of epochs without significant change tolerated before forcing an early stop of training. 

    """

    if not os.path.isfile(os.path.join(path, 'data_secir_age_groups.pickle')):
        ValueError("no dataset found in path: " + path)

    file = open(os.path.join(path, 'data_secir_age_groups.pickle'), 'rb')

    data = pickle.load(file)
    data_splitted = splitdata(data["inputs"], data["labels"])

    train_inputs_compartments = flat_input(data_splitted["train_inputs"])
    train_labels = flat_input(data_splitted["train_labels"])
    valid_inputs_compartments = flat_input(data_splitted["valid_inputs"])
    valid_labels = flat_input(data_splitted["valid_labels"])
    test_inputs_compartments = flat_input(data_splitted["test_inputs"])
    test_labels = flat_input(data_splitted["test_labels"])

    contact_matrices = split_contact_matrices(tf.stack(data["contact_matrix"]))
    contact_matrices_train = flat_input(contact_matrices['train'])
    contact_matrices_valid = flat_input(contact_matrices['valid'])
    contact_matrices_test = flat_input(contact_matrices['test'])

    train_inputs = tf.concat(
        [tf.cast(train_inputs_compartments, tf.float32),
         tf.cast(contact_matrices_train, tf.float32)],
        axis=1, name='concat')
    valid_inputs = tf.concat(
        [tf.cast(valid_inputs_compartments, tf.float32),
         tf.cast(contact_matrices_valid, tf.float32)],
        axis=1, name='concat')
    test_inputs = tf.concat(
        [tf.cast(test_inputs_compartments, tf.float32),
         tf.cast(contact_matrices_test, tf.float32)],
        axis=1, name='concat')

    # train_inputs = tf.concat(
    #     [train_inputs_compartments, contact_matrices_train],
    #     axis=1, name='concat')
    # valid_inputs = tf.concat(
    #     [valid_inputs_compartments, contact_matrices_valid],
    #     axis=1, name='concat')
    # test_inputs = tf.concat(
    #     [test_inputs_compartments, contact_matrices_test],
    #     axis=1, name='concat')

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')

    model.compile(  # loss=tf.keras.losses.MeanSquaredError(),
        loss=tf.keras.losses.MeanAbsolutePercentageError(),
        optimizer=tf.keras.optimizers.Adam(),
        metrics=['mse', 'mae'])

    history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                        validation_data=(
                            valid_inputs, valid_labels),
                        callbacks=[early_stopping])

    plot_losses(history)
    # plotCol(test_inputs, test_labels, model=model,
    # lot_col='Susceptible', max_subplots=3)
    plotCol(test_inputs, test_labels, model=model,
            plot_col='Dead', max_subplots=6)
    plotCol(test_inputs, test_labels, model=model,
            plot_col='Hospitalized', max_subplots=6)
    plotCol(test_inputs, test_labels, model=model,
            plot_col='Infected', max_subplots=6)
    test_statistic(model, test_inputs, test_labels)

    return history

# simple benchmarking


def test_statistic(model, test_inputs, test_labels):
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

    anteil_6 = anteil.transpose(1, 0)
    df = pd.DataFrame(data=anteil_6)
    df = df.transpose()

    mean_percentage = pd.DataFrame(data=(df.mean().values)*100)
    print('Mean percentage error: ', mean_percentage.mean())


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
    plt.show()
    plt.savefig('losses plot.pdf')


if __name__ == "__main__":

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data_groups')

    max_epochs = 600

    network_fit(path_data, mlp_model(), max_epochs=max_epochs)
   # network_fit(path_data, cnn_model(), max_epochs=max_epochs)
