import numpy as np
import pandas as pd
import os
import pickle
import math
import tensorflow as tf
import matplotlib.pyplot as plt
from data_generator import getBaselineMatrix, splitdata, split_contact_matrices, flat_input
from different_models import *

# plots number of people in the chosen compartment (mean of all age groups)


def plotCol(inputs, labels, model=None, plot_col='Infected', max_subplots=3):

    label_width = int(labels.shape[1] / 48)

    plt.figure(figsize=(12, 8))
    # cols = np.array(['Susceptible', 'Exposed', 'Carrier',
    #                'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead'])
    cols = np.array(['Susceptible', 'Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead'])

    plot_col_index = np.where(cols == plot_col)[0][0]
    max_n = min(max_subplots, inputs.shape[0])

    for n in range(max_n):
        plt.subplot(max_n, 1, n+1)
        plt.ylabel(f'{plot_col}')

        label_array = labels[n].numpy()
        label_array = label_array.reshape(-1, 48)

        mean_comp = np.array([])
        for i in (label_array):
            mean_comp = np.append(mean_comp, np.mean(i[plot_col_index::8]))

        plt.scatter(
            np.arange(label_width), mean_comp,
            edgecolors='k', label='Labels', c='#2ca02c', s=64)

        if model is not None:
            input_series = tf.expand_dims(inputs[n], axis=0)
            pred = model(input_series)

            pred = np.asarray(pred).flatten().reshape(-1, 48)
            mean_comp_pred = np.array([])
            for i in (pred):
                mean_comp_pred = np.append(
                    mean_comp_pred, np.mean(i[plot_col_index::8]))

            plt.scatter(np.arange(label_width), mean_comp_pred,
                        marker='X', edgecolors='k', label='Predictions',
                        c='#ff7f0e', s=64)

    plt.xlabel('days')
    plt.show()
    plt.savefig('evaluation_secir_simple_' + plot_col + '.pdf')


def plotCol2(inputs, labels, model=None, plot_col='Infected', max_subplots=3):

    label_width = int(labels.shape[1] / 36)

    plt.figure(figsize=(12, 8))
    # cols = np.array(['Susceptible', 'Exposed', 'Carrier',
    #                'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead'])
    cols = np.array(['Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU',  'Dead'])

    plot_col_index = np.where(cols == plot_col)[0][0]
    max_n = min(max_subplots, inputs.shape[0])

    for n in range(max_n):
        plt.subplot(max_n, 1, n+1)
        plt.ylabel(f'{plot_col}')

        label_array = labels[n].numpy()
        label_array = label_array.reshape(-1, 36)

        mean_comp = np.array([])
        for i in (label_array):
            mean_comp = np.append(mean_comp, np.mean(i[plot_col_index::6]))

        plt.scatter(
            np.arange(label_width), mean_comp,
            edgecolors='k', label='Labels', c='#2ca02c', s=64)

        if model is not None:
            input_series = tf.expand_dims(inputs[n], axis=0)
            pred = model(input_series)

            pred = np.asarray(pred).flatten().reshape(-1, 36)
            mean_comp_pred = np.array([])
            for i in (pred):
                mean_comp_pred = np.append(
                    mean_comp_pred, np.mean(i[plot_col_index::6]))

            plt.scatter(np.arange(label_width), mean_comp_pred,
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

    contact_matrices = split_contact_matrices(tf.stack(data["contact_matrix"]))
    contact_matrices_train = flat_input(contact_matrices['train'])
    contact_matrices_valid = flat_input(contact_matrices['valid'])
    contact_matrices_test = flat_input(contact_matrices['test'])

    n = np.array(data['damping_days']).shape[0]
    train_days = data['damping_days'][:int(n*0.7)]
    valid_days = data['damping_days'][int(n*0.7):int(n*0.9)]
    test_days = data['damping_days'][int(n*0.9):]

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

    def step_decay(epoch):
        initial_lrate = 0.01
        drop = 0.8
        epochs_drop = 5.0
        lrate = initial_lrate * math.pow(drop,
                                         math.floor((1+epoch)/epochs_drop))
        return lrate
    lrate = tf.keras.callbacks.LearningRateScheduler(step_decay)
    callback_list = [lrate]

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
    plotCol2(test_inputs, test_labels, model=model,
             plot_col='Infected', max_subplots=6)
    plotCol2(test_inputs, test_labels, model=model,
             plot_col='Dead', max_subplots=6)
    plotCol2(test_inputs, test_labels, model=model,
             plot_col='Hospitalized', max_subplots=6)
    # plotCol2(test_inputs, test_labels, model=model,
    #        plot_col='Susceptible', max_subplots=6)
    plotCol2(test_inputs, test_labels, model=model,
             plot_col='ICU', max_subplots=6),
    plotCol2(test_inputs, test_labels, model=model,
             plot_col='Exposed', max_subplots=6),
    plotCol2(test_inputs, test_labels, model=model,
             plot_col='Carrier', max_subplots=6),
    # plotCol2(test_inputs, test_labels, model=model,
    #        plot_col='Recovered', max_subplots=6)

    return history


def get_dimensions(path):
    # get data and split inputs from labels, contact matices and damping days
    file = open(os.path.join(path, 'data_secir_age_groups.pickle'), 'rb')

    data = pickle.load(file)

    number_of_dampings = np.asarray(data['damping_days'][0]).shape[0]
    input_dim = (6*8+6*6*number_of_dampings+number_of_dampings)

    output_dim = (data['labels'][0].shape[0])*6*8

    return input_dim, output_dim


def get_dimensions_noSR(path):
    # get data and split inputs from labels, contact matices and damping days
    file = open(os.path.join(path, 'data_secir_age_groups.pickle'), 'rb')

    data = pickle.load(file)

    number_of_dampings = np.asarray(data['damping_days'][0]).shape[0]
    input_dim = (6*8+6*6*number_of_dampings+number_of_dampings)

    output_dim = (data['labels'][0].shape[0])*6*6

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


def get_testdata(path):
    file = open(os.path.join(path, 'data_secir_age_groups.pickle'), 'rb')

    data = pickle.load(file)
    data_splitted = splitdata(data["inputs"], data["labels"])
    test_inputs_compartments = (data_splitted["test_inputs"])
    test_labels = flat_input(data_splitted["test_labels"])
    contact_matrices = split_contact_matrices(tf.stack(data["contact_matrix"]))
    contact_matrices_test = flat_input(contact_matrices['test'])

    n = np.array(data['damping_days']).shape[0]
    test_days = data['damping_days'][int(n*0.9):]

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

    return test_inputs, test_labels


print('x')
if __name__ == "__main__":

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'data_groups_with_damping_6damp_10k_120days')

    max_epochs = 300

    input_dim, output_dim = get_dimensions_noSR(path_data)
    # network_fit(path_data, lstm_multi_output(
    #    input_dim, output_dim), max_epochs=max_epochs)

    network_fit(path_data, cnn_lstm_hybrid_3damp(
        input_dim, output_dim), max_epochs=max_epochs)
