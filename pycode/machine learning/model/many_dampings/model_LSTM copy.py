import numpy as np
import pandas as pd
import os
import pickle
import tensorflow as tf
import matplotlib.pyplot as plt
from data_generator import getBaselineMatrix, splitdata, split_contact_matrices, flat_input
from different_models import *


def plotCol(inputs, labels, model=None, plot_col='Infected', max_subplots=3):

    # 72 damping entries  and 48 age dependent compartments
    input_width = int((inputs.shape[1] - 72) / 48)
    #label_width = int(labels.shape[1] / 48)
    label_width = int(labels.shape[1] / 48)

    plt.figure(figsize=(12, 8))
    # cols = np.array(['Susceptible', 'Exposed', 'Carrier',
    #                'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead'])
    cols = np.array(['Susceptible', 'Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead'])

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
            input_array[plot_col_index: inputs.shape[1] - 72: 48],
            label='Inputs', marker='.', zorder=-10)
        plt.scatter(
            np.arange(input_width, input_width + label_width),
            #label_array[plot_col_index: -1: 48],
            label_array[plot_col_index: -1: 48],
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

    if not os.path.isfile(os.path.join(path, 'data_secir_age_groups.pickle')):
        ValueError("no dataset found in path: " + path)

    # get data and split inputs from labels, contact matices and damping days
    file = open(os.path.join(path, 'data_secir_age_groups.pickle'), 'rb')

    data = pickle.load(file)
    data_splitted = splitdata(data["inputs"], data["labels"])

    train_inputs_compartments = (data_splitted["train_inputs"])
    train_labels = flat_input(data_splitted["train_labels"])
    valid_inputs_compartments = (data_splitted["valid_inputs"])
    valid_labels = flat_input(data_splitted["valid_labels"])
    test_inputs_compartments = (data_splitted["test_inputs"])
    test_labels = flat_input(data_splitted["test_labels"])

    n = np.array(data['damping_days']).shape[0]
    train_days = data['damping_days'][:int(n*0.7)]
    valid_days = data['damping_days'][int(n*0.7):int(n*0.9)]
    test_days = data['damping_days'][int(n*0.9):]

    label_days = data['labels'].shape[1]
    input_days = data['inputs'].shape[1]

    matrices = data['contact_matrix']
    days = data['damping_days']
    matrices_full = []
    for (i, j) in zip(matrices, days):
        matrices_run = i
        days_run = j
        days_diff = []
        k = 1
        days_diff.append(days_run[0])
        while k < len(days_run):
            days_diff.append(days_run[k]-days_run[k-1])
            k += 1
        days_diff.append((label_days+input_days)-days_run[k-1])
        baseline = getBaselineMatrix()
        matrix_array = np.array([baseline] * days_run[0])
        for (n, m) in zip(matrices_run, days_diff[1:]):
            array = np.array([n]*m)
            matrix_array = np.append(matrix_array, array)
        matrices_full.append(matrix_array)

    contact_matrices = split_contact_matrices(tf.stack(matrices_full))
    contact_matrices_train = flat_input(contact_matrices['train'])
    contact_matrices_valid = flat_input(contact_matrices['valid'])
    contact_matrices_test = flat_input(contact_matrices['test'])

    # concatenate the compartment data with contact matrices and damping days
    # to receive complete input data
    new_contact_train = []
    for i in contact_matrices_train:
        new_contact_train.extend([i for j in range(5)])

    new_contact_train = tf.reshape(
        tf.stack(new_contact_train),
        [train_inputs_compartments.shape[0],
         5, contact_matrices_train.shape[1]])

    new_damping_days_train = []
    for i in train_days:
        new_damping_days_train.extend([i for j in range(5)])
    new_damping_days_train = tf.reshape(
        tf.stack(new_damping_days_train),
        [train_inputs_compartments.shape[0],
         5, np.asarray(train_days).shape[1]])

    train_inputs = tf.concat(
        (tf.cast(train_inputs_compartments, tf.float16),
         tf.cast(new_contact_train, tf.float16),
         tf.cast(new_damping_days_train, tf.float16)),
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
         5, np.asarray(test_days).shape[1]])

    test_inputs = tf.concat(
        (tf.cast(test_inputs_compartments, tf.float16),
         tf.cast(new_contact_test, tf.float16),
         tf.cast(new_damping_days_test, tf.float16)),
        axis=2)

    new_contact_val = []
    for i in contact_matrices_valid:
        new_contact_val.extend([i for j in range(5)])

    new_contact_val = tf.reshape(tf.stack(new_contact_val), [
        contact_matrices_valid.shape[0], 5, contact_matrices_valid.shape[1]])

    new_damping_days_valid = []
    for i in valid_days:
        new_damping_days_valid.extend([i for j in range(5)])
    new_damping_days_valid = tf.reshape(
        tf.stack(new_damping_days_valid),
        [valid_inputs_compartments.shape[0],
         5, np.asarray(valid_days).shape[1]])

    valid_inputs = tf.concat(
        (tf.cast(valid_inputs_compartments, tf.float16),
         tf.cast(new_contact_val, tf.float16),
         tf.cast(new_damping_days_valid, tf.float16)),
        axis=2)

    # run the model

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')

    decayed_lr = tf.keras.optimizers.schedules.ExponentialDecay(
        initial_learning_rate=0.01,
        decay_steps=200,
        decay_rate=0.95,
        staircase=True)

    model.compile(  # loss=tf.keras.losses.MeanSquaredError(),
        loss=tf.keras.losses.MeanAbsolutePercentageError(),
        optimizer=tf.keras.optimizers.Adam(
            learning_rate=decayed_lr),
        metrics=['mse', 'mae'])

    history = model.fit(
        train_inputs, train_labels, epochs=max_epochs,
        validation_data=(valid_inputs, valid_labels),
        callbacks=[early_stopping], batch_size=32)

    plot_losses(history)
    test_statistic(model, test_inputs, test_labels)

    return history


def get_dimensions(path):
    # get data and split inputs from labels, contact matices and damping days
    file = open(os.path.join(path, 'data_secir_age_groups.pickle'), 'rb')

    data = pickle.load(file)

    number_of_dampings = np.asarray(data['damping_days'][0]).shape[0]

    input_dim = data['inputs'][0].shape[0]
    output_dim = (data['labels'][0].shape[0])*6*8

    return input_dim, output_dim


# simple benchmarking


def test_statistic(model, test_inputs, test_labels):
    pred = np.float16(model(test_inputs))

    test_labels = np.array(test_labels)

    diff = pred - test_labels
    anteil = (abs(diff))/abs(test_labels)

    anteil_6 = anteil.transpose(1, 0)
    df = pd.DataFrame(data=anteil_6)
    df = df.transpose()

    mean_percentage = pd.DataFrame(data=(df.mean().values)*100)
    print('Mean percentage error: ', mean_percentage.mean())


# Plot Performance


def plot_losses(history):
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.show()
    plt.savefig('losses plot.pdf')


print('x')
if __name__ == "__main__":

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'data_groups_with_damping_5damp')

    max_epochs = 800

    input_dim, output_dim = get_dimensions(path_data)
    # network_fit(path_data, lstm_multi_output(
    #    input_dim, output_dim), max_epochs=max_epochs)

    network_fit(path_data, cnn_lstm_hybrid_5damp(
        input_dim, output_dim), max_epochs=max_epochs)
