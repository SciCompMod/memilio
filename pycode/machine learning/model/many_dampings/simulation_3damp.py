import os
import pickle
from model_LSTM import *
from different_models import *


def get_simulation_data():
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'data_groups_with_damping_9damp_sim')
    # path_data = os.path.join(os.path.dirname(os.path.realpath(
    #   os.path.dirname(os.path.realpath(path)))), 'data_simulation_fix')
    file = open(os.path.join(path_data, 'data_secir_age_groups.pickle'), 'rb')

    data = pickle.load(file)

    return data


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
    simulation(model)
    test_statistic(model, test_inputs, test_labels)
    return history


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

        damping_days = []
        for j in range(len(data['damping_days'][0])):
            damping_days.append('t' + str(j))

        for j in range(len(damping_days)):
            damping_days[j] = data['damping_days'][i][j]

        matrices_array = []
        for j in range(len(data['damping_days'][0])):
            matrices_array.append('matrix' + str(j))

        for j in range(len(matrices_array)):
            matrices_array[j] = np.asarray((flat_input(matrices[i][j])))

        dampings_in_input = 3

        first_matrix = []
        first_matrix.extend(
            [np.asarray(matrices_array[:dampings_in_input]).flatten(order='C')
             for j in range(5)])

        first_days = []
        first_days.extend(
            [np.asarray(damping_days[: dampings_in_input]).flatten(
                order='C') for j in range(5)])

        first_input = tf.concat(
            [tf.cast(tf.reshape(sim_inputs_compartments[i], [5, 48]),
                     tf.float16),
             tf.cast(
                tf.stack(first_matrix),
                tf.float16),
             tf.cast(first_days,  tf.float16)],
            axis=1, name='concat')

        results = []

        first_output = model.predict(tf.reshape(
            first_input, [1, first_input.shape[0], first_input.shape[1]]))

        first_days = damping_days[dampings_in_input]-input_width
        results = np.append(
            results, (first_output.reshape(output_dim)[:(first_days*6*8)]))

        second_matrix = []
        second_matrix.extend(
            [np.asarray(matrices_array[dampings_in_input:dampings_in_input*2]).flatten(order='C')
             for j in range(5)])

        second_days = []
        second_days.extend([np.asarray(
            damping_days[dampings_in_input:dampings_in_input*2]).flatten(order='C') for j in range(5)])

        second_input = tf.concat(
            [tf.cast(tf.reshape(results[(-5*6*8):], [5, 48]),
                     tf.float16),
             tf.cast(
                tf.stack(second_matrix),
                tf.float16),
             tf.cast(second_days,  tf.float16)],
            axis=1, name='concat')

        second_output = model.predict(
            tf.reshape(
                second_input,
                [1, second_input.shape[0],
                 second_input.shape[1]]))

        second_days = damping_days[dampings_in_input*2] - \
            damping_days[dampings_in_input]
        results = np.append(
            results, (second_output.reshape(output_dim)[:(second_days*6*8)]))

        third_matrix = []
        third_matrix.extend(
            [np.asarray(matrices_array[dampings_in_input*2:dampings_in_input*3]).flatten(order='C')
             for j in range(5)])

        third_days = []
        third_days.extend(
            [np.asarray(
                damping_days
                [dampings_in_input * 2: dampings_in_input * 3]).flatten(
                order='C') for j in range(5)])

        third_input = tf.concat(
            [tf.cast(tf.reshape(results[(-5*6*8):], [5, 48]),
                     tf.float16),
             tf.cast(
                tf.stack(third_matrix),
                tf.float16),
             tf.cast(third_days,  tf.float16)],
            axis=1, name='concat')

        third_output = model.predict(
            tf.reshape(
                third_input,
                [1, third_input.shape[0],
                 third_input.shape[1]]))

        third_days = damping_days[-1]-damping_days[dampings_in_input*2]
        # results = np.append(
        #     results, (third_output.reshape(output_dim)[:(third_days*6*8)]))
        results = np.append(
            results, (third_output.reshape(output_dim)))
        labels_run = np.asarray(sim_labels[i]).flatten()
        results = results[:labels_run.shape[0]]
        mean_percentage_error_splitted = []
        diff = results-labels_run
        anteil = (abs(diff))/abs(labels_run)
        mean_percentage_error = np.append(mean_percentage_error, anteil.mean())
        mean_percentage_error = np.append(mean_percentage_error, anteil.mean())

        i += 1
    print('Simulation MAPE: ', mean_percentage_error.mean()*100)

    # array = [

    #     anteil[t1 * 6 * 8: t2 * 6 * 8].mean(),
    #     anteil[t2 * 6 * 8: t3 * 6 * 8].mean(),
    #     anteil[t3 * 6 * 8:].mean()]
    # mean_percentage_error_splitted.append(array)

    # MAPE_array = [
    #     np.asarray(mean_percentage_error_splitted).transpose()[0].mean(),
    #     np.asarray(mean_percentage_error_splitted).transpose()[1].mean(),
    #     np.asarray(mean_percentage_error_splitted).transpose()[2].mean(),
    #     np.asarray(mean_percentage_error_splitted).transpose()[3].mean()]
    # plt.plot(MAPE_array)
    plt.savefig('MAPE_array')


print('x')
if __name__ == "__main__":
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'data_groups_with_damping_9damp_split_1mio')

    max_epochs = 50

    input_dim, output_dim = get_dimensions(path_data)
    # network_fit(path_data, lstm_multi_output(
    #    input_dim, output_dim), max_epochs=max_epochs)

    network_fit(path_data, cnn_lstm_hybrid_3damp(
        input_dim, output_dim), max_epochs=max_epochs)
