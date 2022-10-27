from statistics import mean
import numpy as np
import pandas as pd
from datetime import date
from math import ceil
import random
import os
import pickle
from progress.bar import Bar  # pip install progess
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf
import matplotlib.pyplot as plt
from tensorflow import keras
from keras import Sequential
from keras.layers import Dense
from keras import backend as K
import seaborn as sns  # plot after normalization
from data_secir_age_groups import getBaselineMatrix, splitdata, split_contact_matrices, flat_input
from different_models import *
from tensorflow import keras
from data_secir_age_groups import run_secir_groups_simulation, get_population
from memilio.simulation import UncertainContactMatrix, ContactMatrix, Damping
from memilio.simulation.secir import SecirModel, simulate, AgeGroup, Index_InfectionState, SecirSimulation
from memilio.simulation.secir import InfectionState as State


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

    file = open(os.path.join(path, 'data_secir_age_groups.pickle'), 'rb')

    data = pickle.load(file)
    data_splitted = splitdata(data["inputs"], data["labels"])

    # train_inputs_compartments = flat_input(data_splitted["train_inputs"])
    # train_labels = flat_input(data_splitted["train_labels"])
    # valid_inputs_compartments = flat_input(data_splitted["valid_inputs"])
    # valid_labels = flat_input(data_splitted["valid_labels"])
    # test_inputs_compartments = flat_input(data_splitted["test_inputs"])
    # test_labels = flat_input(data_splitted["test_labels"])

# for lstm :

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


# for mlp and cnn (?)
    # train_inputs = tf.concat(
    #     [tf.cast(train_inputs_compartments, tf.float16),
    #      tf.cast(contact_matrices_train, tf.float16)],
    #     axis=1, name='concat')
    # train_inputs = np.expand_dims(train_inputs, -1)
    # valid_inputs = tf.concat(
    #     [tf.cast(valid_inputs_compartments, tf.float16),
    #      tf.cast(contact_matrices_valid, tf.float16)],
    #     axis=1, name='concat')
    # test_inputs = tf.concat(
    #     [tf.cast(test_inputs_compartments, tf.float16),
    #      tf.cast(contact_matrices_test, tf.float16)],
    #     axis=1, name='concat')

  # for lstm :

    new_contact_train = []
    for i in contact_matrices_train:
        new_contact_train.extend([i for j in range(5)])

    new_contact_train = tf.reshape(tf.stack(new_contact_train), [
                                   contact_matrices_train.shape[0], 5, 72])

    train_inputs = tf.concat(
        (tf.cast(train_inputs_compartments, tf.float16),
         tf.cast(new_contact_train, tf.float16)),
        axis=2)
    new_contact_test = []
    for i in contact_matrices_test:
        new_contact_test.extend([i for j in range(5)])

    new_contact_test = tf.reshape(tf.stack(new_contact_test), [
        contact_matrices_test.shape[0], 5, 72])

    test_inputs = tf.concat(
        (tf.cast(test_inputs_compartments, tf.float16),
         tf.cast(new_contact_test, tf.float16)),
        axis=2)
    new_contact_val = []
    for i in contact_matrices_valid:
        new_contact_val.extend([i for j in range(5)])

    new_contact_val = tf.reshape(tf.stack(new_contact_val), [
        contact_matrices_valid.shape[0], 5, 72])

    valid_inputs = tf.concat(
        (tf.cast(valid_inputs_compartments, tf.float16),
         tf.cast(new_contact_val, tf.float16)),
        axis=2)

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

    # mc = ModelCheckpoint("model_secir_age_group" + '.h5', monitor=['mse'],
    #                      mode='min', save_best_only=True)

    model.compile(  # loss=tf.keras.losses.MeanSquaredError(),
        loss=tf.keras.losses.MeanAbsolutePercentageError(),
        optimizer=tf.keras.optimizers.Nadam(learning_rate=0.0001),
        metrics=['mse', 'mae'])

    history = model.fit(
        train_inputs, train_labels, epochs=max_epochs,
        validation_data=(valid_inputs, valid_labels),
        callbacks=[early_stopping])

    # plot_losses(history)
    # plotCol(test_inputs, test_labels, model=model,
    # # lot_col='Susceptible', max_subplots=3)
    # plotCol(test_inputs, test_labels, model=model,
    #         plot_col='Dead', max_subplots=6)
    # plotCol(test_inputs, test_labels, model=model,
    #         plot_col='Hospitalized', max_subplots=6)
    # plotCol(test_inputs, test_labels, model=model,
    #         plot_col='Infected', max_subplots=6)
    test_statistic(model, test_inputs, test_labels)

    simulation(model)
    return history

# simple benchmarking


def test_statistic(model, test_inputs, test_labels):
    pred = np.float16(model(test_inputs))
    #pred = pred.numpy()
    test_labels = np.array(test_labels)

    diff = pred - test_labels
    anteil = (abs(diff))/abs(test_labels)

    anteil_6 = anteil.transpose(1, 0)
    df = pd.DataFrame(data=anteil_6)
    df = df.transpose()

    mean_percentage = pd.DataFrame(data=(df.mean().values)*100)
    print('Mean percentage error: ', mean_percentage.mean())


def timereps(reps, model, input):
    from time import time
    start = time()
    for i in range(0, reps):
        _ = model(input)
    end = time()
    time_passed = end - start
    print(time_passed)
    return (end - start) / reps

# Plot Performance


def plot_histories(histories):
    model_names = ["LSTM", "Dense"]
    count_names = 0
    interval = np.arange(len(histories))
    for x in interval:
        history = histories[x]
        width = 0.3

        train_mae = history.history['mean_absolute_error'][-1]
        valid_mae = history.history['val_mean_absolute_error'][-1]

        plt.bar(x - 0.17, train_mae, width, label='Train')
        plt.bar(x + 0.17, valid_mae, width, label='Validation')
        name = model_names[x]
        # plt.xticks(ticks=x, labels=[name,name],
        #         rotation=45)
        plt.ylabel(f'MAE')
        _ = plt.legend()
        count_names = count_names + 1

    plt.show()
    plt.savefig('evaluation_single_shot.pdf')


def plot_losses(history):
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.show()
    plt.savefig('losses plot.pdf')


def get_simulation_data():
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data_simulation')
    file = open(os.path.join(path_data, 'data_secir_age_groups.pickle'), 'rb')

    data = pickle.load(file)

    return data


def simulation(model):
    input_width = 5
    data = get_simulation_data()
    sim_inputs_compartments = flat_input(data['inputs'])
    sim_labels = (data['labels'])
    sim_contact_matrix = (data['contact_matrix'])
    matrices = np.asarray(sim_contact_matrix)
    baseline = getBaselineMatrix()

    i = 0

    mean_percentage_error = []
    mean_percentage_error_splitted = []
    while i < len(sim_inputs_compartments):

        t1 = data['days'][i][0]
        t2 = data['days'][i][1]
        t3 = data['days'][i][2]

        matrix1 = np.asarray((flat_input(matrices[i][0])))
        matrix2 = np.asarray((flat_input(matrices[i][1])))
        matrix3 = np.asarray((flat_input(matrices[i][2])))

        first_matrix = []
        first_matrix.extend(
            [np.concatenate((baseline, matrix1)).flatten(order='C')
             for j in range(5)])
        first_input = tf.concat(
            [tf.cast(tf.reshape(sim_inputs_compartments[i], [5, 48]),
                     tf.float16),
             tf.cast(
                tf.stack(first_matrix),
                tf.float16)],
            axis=1, name='concat')
        results = []
        first_output = model.predict(tf.reshape(first_input, [1, 5, 120]))
        first_days = t2-t1
        results = np.append(
            results, (first_output.reshape(960)[:(first_days*6*8)]))
        second_matrix = []
        second_matrix.extend(
            [np.concatenate((matrix1, matrix2)).flatten(order='C')
             for j in range(5)])

        second_input = tf.concat(
            [tf.cast(tf.reshape(results[-240:], [5, 48]),
                     tf.float16),
             tf.cast(
                second_matrix,
                tf.float16)],
            axis=1, name='concat')
        second_output = model.predict(tf.reshape(second_input, [1, 5, 120]))
        second_days = t3 - t2
        results = np.append(
            results, (second_output.reshape(960)[:(second_days*6*8)]))

        third_matrix = []
        third_matrix.extend(
            [np.concatenate((matrix2, matrix3)).flatten(order='C')
             for j in range(5)])

        third_input = tf.concat(
            [tf.cast(tf.reshape(results[-240:], [5, 48]),
                     tf.float16),
             tf.cast(
                third_matrix,
                tf.float16)],
            axis=1, name='concat')
        third_output = model.predict(tf.reshape(third_input, [1, 5, 120]))

        third_days = 20
        results = np.append(
            results, (third_output.reshape(960)[:(third_days*6*8)]))

        # fourth_matrix = []
        # fourth_matrix.extend(
        #     [np.concatenate((matrix2, matrix3)).flatten(order='C')
        #      for j in range(5)])

        # fourth_input = tf.concat(
        #     [tf.cast(tf.reshape(results[-240:], [5, 48]),
        #              tf.float16),
        #      tf.cast(
        #         fourth_matrix,
        #         tf.float16)],
        #     axis=1, name='concat')
        # fourth_output = model.predict(tf.reshape(fourth_input, [1, 5, 120]))
        # fourth_days = 20
        # results = np.append(
        #     results, (fourth_output.reshape(960)[:(fourth_days*6*8)]))
        labels_run = np.asarray(sim_labels[i].flatten())
        #fourth_days = int((len(labels_run)/48) - t3)

        diff = results-labels_run
        anteil = (abs(diff))/abs(labels_run)
        mean_percentage_error = np.append(mean_percentage_error, anteil.mean())
        mean_percentage_error = np.append(mean_percentage_error, anteil.mean())
        array = [
            #anteil[: t1 * 6 * 8].mean(),
            anteil[t1 * 6 * 8: t2 * 6 * 8].mean(),
            anteil[t2 * 6 * 8: t3 * 6 * 8].mean(),
            anteil[t3 * 6 * 8:].mean()]
        mean_percentage_error_splitted.append(array)

        i += 1

    MAPE_array = [

        np.asarray(mean_percentage_error_splitted).transpose()[0].mean(),
        np.asarray(mean_percentage_error_splitted).transpose()[1].mean(),
        np.asarray(mean_percentage_error_splitted).transpose()[2].mean()]
    # np.asarray(mean_percentage_error_splitted).transpose()[3].mean(),]
    plt.plot(MAPE_array)
    plt.savefig('MAPE_array')
    print('Array with the average loss after each damping:', MAPE_array)
    print('Simulation MAPE: ', mean_percentage_error.mean()*100)


# for cnn model :

    #     first_input = tf.concat(
    #         [tf.cast(sim_inputs_compartments[i],
    #                  tf.float16),
    #          tf.cast(
    #             np.concatenate((baseline, baseline)).flatten(
    #                 order='C'),
    #             tf.float16)],
    #         axis=0, name='concat')

    #     results = []
    #     first_output = model.predict(tf.reshape(first_input, [1, 312]))

    #     first_days = t1-input_width
    #     results = np.append(
    #         results, (first_output.reshape(960)[:(first_days*6*8)]))

    #     second_input = tf.concat(
    #         [tf.cast(results[-240:],
    #                  tf.float16),
    #          tf.cast(
    #             np.concatenate((baseline, matrix1)).flatten(
    #                 order='C'),
    #             tf.float16)],
    #         axis=0, name='concat')

    #     second_output = model.predict(tf.reshape(second_input, [1, 312]))
    #     second_days = t2-t1
    #     results = np.append(
    #         results, (second_output.reshape(960)[:(second_days*6*8)]))

    #     third_input = tf.concat(
    #         [tf.cast(results[-240:],
    #                  tf.float16),
    #          tf.cast(
    #             np.concatenate((matrix1, matrix2)).flatten(
    #                 order='C'),
    #             tf.float16)],
    #         axis=0, name='concat')

    #     third_output = model.predict(tf.reshape(third_input, [1, 312]))

    #     third_days = t3-t2
    #     results = np.append(
    #         results, (third_output.reshape(960)[:(third_days*6*8)]))

    #     ####
    #     fourth_input = tf.concat(
    #         [tf.cast(results[-240:],
    #                  tf.float16),
    #          tf.cast(
    #             np.concatenate((matrix2, matrix3)).flatten(
    #                 order='C'),
    #             tf.float16)],
    #         axis=0, name='concat')

    #     fourth_output = model.predict(tf.reshape(fourth_input, [1, 312]))
    #     labels_run = np.asarray(sim_labels[i].flatten())
    #     #fourth_days = int((len(labels_run)/48) - t3)
    #     fourth_days = 20
    #     results = np.append(
    #         results, (fourth_output.reshape(960)[:(fourth_days*6*8)]))

    # first_matrix = []
    # first_matrix.extend(
    #     [np.concatenate((baseline, baseline)).flatten(order='C')
    #      for j in range(5)])
    # first_input = tf.concat(
    #     [tf.cast(tf.reshape(sim_inputs_compartments[i], [5, 48]),
    #              tf.float16),
    #      tf.cast(
    #         tf.stack(first_matrix),
    #         tf.float16)],
    #     axis=1, name='concat')
    # results = []
    # first_output = model.predict(tf.reshape(first_input, [1, 5, 120]))
    # first_days = t1-input_width
    # results = np.append(
    #     results, (first_output.reshape(960)[:(first_days*6*8)]))
    # second_matrix = []
    # second_matrix.extend(
    #     [np.concatenate((baseline, matrix1)).flatten(order='C')
    #      for j in range(5)])

    # second_input = tf.concat(
    #     [tf.cast(tf.reshape(results[-240:], [5, 48]),
    #              tf.float16),
    #      tf.cast(
    #         second_matrix,
    #         tf.float16)],
    #     axis=1, name='concat')
    # second_output = model.predict(tf.reshape(second_input, [1, 5, 120]))
    # second_days = t2 - t1
    # results = np.append(
    #     results, (first_output.reshape(960)[:(second_days*6*8)]))

    # third_matrix = []
    # third_matrix.extend(
    #     [np.concatenate((matrix1, matrix2)).flatten(order='C')
    #      for j in range(5)])

    # third_input = tf.concat(
    #     [tf.cast(tf.reshape(results[-240:], [5, 48]),
    #              tf.float16),
    #      tf.cast(
    #         third_matrix,
    #         tf.float16)],
    #     axis=1, name='concat')
    # third_output = model.predict(tf.reshape(second_input, [1, 5, 120]))

    # third_days = t3-t2
    # results = np.append(
    #     results, (third_output.reshape(960)[:(third_days*6*8)]))

    # fourth_matrix = []
    # fourth_matrix.extend(
    #     [np.concatenate((matrix2, matrix3)).flatten(order='C')
    #      for j in range(5)])

    # fourth_input = tf.concat(
    #     [tf.cast(tf.reshape(results[-240:], [5, 48]),
    #              tf.float16),
    #      tf.cast(
    #         fourth_matrix,
    #         tf.float16)],
    #     axis=1, name='concat')
    # fourth_output = model.predict(tf.reshape(third_input, [1, 5, 120]))

    # labels_run = np.asarray(sim_labels[i].flatten())
    # #fourth_days = int((len(labels_run)/48) - t3)
    # fourth_days = 20
    # results = np.append(
    #     results, (fourth_output.reshape(960)[:(fourth_days*6*8)]))
    # MAPE_array = [
    #     np.asarray(mean_percentage_error_splitted).transpose()[0].mean(),
    #     np.asarray(mean_percentage_error_splitted).transpose()[1].mean(),
    #     np.asarray(mean_percentage_error_splitted).transpose()[2].mean(),
    #     np.asarray(mean_percentage_error_splitted).transpose()[3].mean()]
    # plt.plot(MAPE_array)
    # plt.savefig('MAPE_array')


# with t1 always= 5
    # input_width = 5
    # data = get_simulation_data()
    # sim_inputs_compartments = flat_input(data['inputs'])
    # sim_labels = (data['labels'])
    # sim_contact_matrix = (data['contact_matrix'])
    # matrices = np.asarray(sim_contact_matrix)
    # baseline = getBaselineMatrix()

    # i = 0

    # mean_percentage_error = []
    # mean_percentage_error_splitted = []
    # while i < len(sim_inputs_compartments):

    #     t1 = data['days'][i][0]
    #     t2 = data['days'][i][1]
    #     t3 = data['days'][i][2]

    #     matrix1 = np.asarray((flat_input(matrices[i][0])))
    #     matrix2 = np.asarray((flat_input(matrices[i][1])))
    #     matrix3 = np.asarray((flat_input(matrices[i][2])))

    #     first_matrix = []
    #     first_matrix.extend(
    #         [np.concatenate((baseline, matrix1)).flatten(order='C')
    #          for j in range(5)])
    #     first_input = tf.concat(
    #         [tf.cast(tf.reshape(sim_inputs_compartments[i], [5, 48]),
    #                  tf.float16),
    #          tf.cast(
    #             tf.stack(first_matrix),
    #             tf.float16)],
    #         axis=1, name='concat')
    #     results = []
    #     first_output = model.predict(tf.reshape(first_input, [1, 5, 120]))
    #     first_days = t2-t1
    #     results = np.append(
    #         results, (first_output.reshape(960)[:(first_days*6*8)]))
    #     second_matrix = []
    #     second_matrix.extend(
    #         [np.concatenate((matrix1, matrix2)).flatten(order='C')
    #          for j in range(5)])

    #     second_input = tf.concat(
    #         [tf.cast(tf.reshape(results[-240:], [5, 48]),
    #                  tf.float16),
    #          tf.cast(
    #             second_matrix,
    #             tf.float16)],
    #         axis=1, name='concat')
    #     second_output = model.predict(tf.reshape(second_input, [1, 5, 120]))
    #     second_days = t3 - t2
    #     results = np.append(
    #         results, (second_output.reshape(960)[:(second_days*6*8)]))

    #     third_matrix = []
    #     third_matrix.extend(
    #         [np.concatenate((matrix2, matrix3)).flatten(order='C')
    #          for j in range(5)])

    #     third_input = tf.concat(
    #         [tf.cast(tf.reshape(results[-240:], [5, 48]),
    #                  tf.float16),
    #          tf.cast(
    #             third_matrix,
    #             tf.float16)],
    #         axis=1, name='concat')
    #     third_output = model.predict(tf.reshape(third_input, [1, 5, 120]))

    #     third_days = 20
    #     results = np.append(
    #         results, (third_output.reshape(960)[:(third_days*6*8)]))

    #     labels_run = np.asarray(sim_labels[i].flatten())
    #     #fourth_days = int((len(labels_run)/48) - t3)

    #     diff = results-labels_run
    #     anteil = (abs(diff))/abs(labels_run)
    #     mean_percentage_error = np.append(mean_percentage_error, anteil.mean())
    #     mean_percentage_error = np.append(mean_percentage_error, anteil.mean())
    #     array = [

    #         anteil[t1 * 6 * 8: t2 * 6 * 8].mean(),
    #         anteil[t2 * 6 * 8: t3 * 6 * 8].mean(),
    #         anteil[t3 * 6 * 8:].mean()]
    #     mean_percentage_error_splitted.append(array)

    #     i += 1

    # MAPE_array = [
    #     np.asarray(mean_percentage_error_splitted).transpose()[0].mean(),
    #     np.asarray(mean_percentage_error_splitted).transpose()[1].mean(),
    #     np.asarray(mean_percentage_error_splitted).transpose()[2].mean()
    #    ]
    # plt.plot(MAPE_array)
    # plt.savefig('MAPE_array')

    # print('Simulation MAPE: ', mean_percentage_error.mean()*100)


print('x')

if __name__ == "__main__":
    # TODO: Save contact matrix depending on the damping.
    # In the actual state it might be enough to save the regular one and the damping

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data_groups_withSR')

    max_epochs = 10

    network_fit(path_data, lstm_multi_output(), max_epochs=max_epochs)
    #network_fit(path_data, cnn_model(), max_epochs=max_epochs)
