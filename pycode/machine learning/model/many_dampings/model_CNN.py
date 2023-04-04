import numpy as np
import pandas as pd
import os
import pickle
import tensorflow as tf
import matplotlib.pyplot as plt
from data_generator import getBaselineMatrix, splitdata, split_contact_matrices, flat_input
from different_models import *


def network_fit(path, model, max_epochs=30, early_stop=3000):

    if not os.path.isfile(
            os.path.join(path, 'data_secir_age_groups.pickle')):
        ValueError("no dataset found in path: " + path)

    file = open(
        os.path.join(path, 'data_secir_age_groups.pickle'),
        'rb')

    data = pickle.load(file)
    data_splitted = splitdata(data["inputs"], data["labels"])

    train_inputs_compartments = flat_input(data_splitted["train_inputs"])
    train_labels = flat_input(data_splitted["train_labels"])
    valid_inputs_compartments = flat_input(data_splitted["valid_inputs"])
    valid_labels = flat_input(data_splitted["valid_labels"])
    test_inputs_compartments = flat_input(data_splitted["test_inputs"])
    test_labels = flat_input(data_splitted["test_labels"])

    #contact_matrices = split_contact_matrices(tf.stack(data["contact_matrix"]))

    contact_matrices = split_contact_matrices(
        tf.stack(np.asarray(data['contact_matrix']).reshape(15000, 1, 6, 6)))
    contact_matrices_train = flat_input(contact_matrices['train'])
    contact_matrices_valid = flat_input(contact_matrices['valid'])
    contact_matrices_test = flat_input(contact_matrices['test'])

    # n = np.array(data['damping_days']).shape[0]
    # train_days = data['damping_days'][:int(n*0.7)]
    # valid_days = data['damping_days'][int(n*0.7):int(n*0.9)]
    # test_days = data['damping_days'][int(n*0.9):]

    #####
    # num_of_dampings = 3

    # new_contact_train = []
    # for i in contact_matrices_train:
    #     new_contact_train.extend([i for j in range(num_of_dampings)])

    # new_damping_days_train = []
    # for i in train_days:
    #     new_damping_days_train.extend([i for j in range(num_of_dampings)])

    # new_contact_test = []
    # for i in contact_matrices_test:
    #     new_contact_test.extend([i for j in range(num_of_dampings)])

    # new_damping_days_test = []
    # for i in test_days:
    #     new_damping_days_test.extend([i for j in range(num_of_dampings)])

    # new_contact_val = []
    # for i in contact_matrices_valid:
    #     new_contact_val.extend([i for j in range(num_of_dampings)])

    # new_damping_days_valid = []
    # for i in valid_days:
    #     new_damping_days_valid.extend([i for j in range(num_of_dampings)])

    #######

    train_inputs = tf.concat(
        [tf.cast(train_inputs_compartments, tf.float16),
         tf.cast(contact_matrices_train, tf.float16)],
        axis=1, name='concat')
    train_inputs = np.expand_dims(train_inputs, -1)
    valid_inputs = tf.concat(
        [tf.cast(valid_inputs_compartments, tf.float16),
         tf.cast(contact_matrices_valid, tf.float16)],
        axis=1, name='concat')
    test_inputs = tf.concat(
        [tf.cast(test_inputs_compartments, tf.float16),
         tf.cast(contact_matrices_test, tf.float16)],
        axis=1, name='concat')

    # train_inputs = tf.concat(
    #     [tf.cast(train_inputs_compartments, tf.float16),
    #      tf.cast(new_contact_train, tf.float16)],
    #     axis=1, name='concat')
    # train_inputs = np.expand_dims(train_inputs, -1)
    # valid_inputs = tf.concat(
    #     [tf.cast(valid_inputs_compartments, tf.float16),
    #      tf.cast(new_contact_val, tf.float16)],
    #     axis=1, name='concat')
    # test_inputs = tf.concat(
    #     [tf.cast(test_inputs_compartments, tf.float16),
    #      tf.cast(new_contact_test, tf.float16)],
    #     axis=1, name='concat')

    # run the model

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')

    model.compile(  # loss=tf.keras.losses.MeanSquaredError(),
        loss=tf.keras.losses.MeanAbsolutePercentageError(),
        optimizer=tf.keras.optimizers.Adam(
            learning_rate=0.0001),
        metrics=['mse', 'mae'])

    history = model.fit(
        train_inputs, train_labels, epochs=max_epochs,
        validation_data=(valid_inputs, valid_labels),
        callbacks=[early_stopping], batch_size=32)

    plot_losses(history)
    test_statistic(model, test_inputs, test_labels)

    return history


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


def get_dimensions(path):
    # get data and split inputs from labels, contact matices and damping days
    file = open(
        os.path.join(path, 'data_secir_age_groups.pickle'),
        'rb')

    data = pickle.load(file)

    number_of_dampings = np.asarray(data['damping_days'][0]).shape[0]
    input_dim = (6*8*5+6*6)

    output_dim = (data['labels'][0].shape[0])*6*8

    return input_dim, output_dim


print('x')
if __name__ == "__main__":

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'data_groups_with_damping_split')

    max_epochs = 200

    input_dim, output_dim = get_dimensions(path_data)
    network_fit(path_data, cnn_model(
        input_dim, output_dim), max_epochs=max_epochs)


# 1440 ouput
# 348 input
