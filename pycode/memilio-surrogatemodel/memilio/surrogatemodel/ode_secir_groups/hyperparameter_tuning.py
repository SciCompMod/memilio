# from memilio.surrogatemodel.ode_secir_simple import network_architectures
# from memilio.surrogatemodel.ode_secir_simple.model import split_data, get_test_statistic
import os
import tensorflow as tf
import pickle
import pandas as pd
import time
from sklearn.model_selection import KFold
import numpy as np
from keras.optimizers import Adam, Nadam, SGD, Adagrad, RMSprop

# grid serach of activation and optimizer for our best model: the LSTM with 0 hidden layers and 128 neurons per layer
optimizers = ['Adam', 'Nadam', 'SGD', 'Adagrad', 'RMSProp']

parameters = []

for o in optimizers:
    parameters.append((o))

modelname = 'LSTM'

df_results = pd.DataFrame(
    columns=['model', 'optimizer',
             'mean_test_MAPE', 'kfold_train',
             'kfold_val', 'kfold_test', 'training_time',
             'train_losses', 'val_losses'])


def train_and_evaluate_model(param, max_epochs, filename, filename_df):

    optimizer = param

    if not os.path.isfile(os.path.join(path_data, 'data_secir_simple.pickle')):
        ValueError("no dataset found in path: " + path_data)

    file = open(os.path.join(path_data, filename), 'rb')

    data = pickle.load(file)
    inputs_grid_search = data['inputs'][:int((0.8 * len(data['inputs'])))]
    labels_grid_search = data['labels'][:int((0.8 * len(data['labels'])))]

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')

    kf = KFold(n_splits=5)
    train_losses = []
    val_losses = []

    losses_history_all = []
    val_losses_history_all = []

    start = time.perf_counter()
    for train_idx, val_idx in kf.split(inputs_grid_search):

        num_outputs = 8
        model = tf.keras.Sequential([
            tf.keras.layers.LSTM(128, return_sequences=False),
            tf.keras.layers.Dense(label_width * 8 * 6,
                                  kernel_initializer=tf.initializers.zeros()),
            tf.keras.layers.Reshape([label_width, 8 * 6])
        ])

        # Gather training and validation data based on the fold
        train_inputs = tf.gather(inputs_grid_search, indices=train_idx)
        train_labels = tf.gather(labels_grid_search, indices=train_idx)
        valid_inputs = tf.gather(inputs_grid_search, indices=val_idx)
        valid_labels = tf.gather(labels_grid_search, indices=val_idx)

        # start = time.perf_counter()

        model.compile(
            loss=tf.keras.losses.MeanAbsolutePercentageError(),
            optimizer=optimizer,
            metrics=[tf.keras.metrics.MeanAbsoluteError(), tf.keras.metrics.MeanAbsolutePercentageError()])

        history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                            validation_data=(valid_inputs, valid_labels),
                            callbacks=[early_stopping])

        train_losses.append(np.min(history.history['loss']))
        val_losses.append(np.min(history.history['val_loss']))
        losses_history_all.append(history.history['loss'])
        val_losses_history_all.append(history.history['val_loss'])

    elapsed = time.perf_counter() - start

    # Print out the results
    elapsed = time.perf_counter() - start
    print(f"Best train losses: {train_losses}")
    print(f"Best validation losses: {val_losses}")
    print("--------------------------------------------")
    print(f"K-Fold Train Score: {np.mean(train_losses)}")
    print(f"K-Fold Validation Score: {np.mean(val_losses)}")
    print(f"Time for training: {elapsed:.4f} seconds")
    print(f"Time for training: {elapsed / 60:.4f} minutes")

    # After cross-validation, we can test on the withhold dataset (outside of the loop)
    df_results.loc[len(df_results.index)] = [
        modelname, optimizer,  np.nan,  # Placeholder for test score
        np.mean(train_losses),
        np.mean(val_losses),
        np.nan,  # Placeholder for final test score
        elapsed / 60,
        [losses_history_all],
        [val_losses_history_all]
    ]

    path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'secir_groups_grid_search_paper')
    if not os.path.isdir(file_path):
        os.mkdir(file_path)
    file_path = os.path.join(file_path, filename_df)
    df_results.to_csv(file_path)


start_hyper = time.perf_counter()
max_epochs = 1500


path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'data_paper')

filename = "data_secir_groups_30days_Germany_10k_nodamp.pickle"

label_width = 30
early_stop = 100
filename_df = "groups_dataframe_optimizer_paper.csv"

for param in parameters:
    train_and_evaluate_model(param, max_epochs, filename, filename_df)

elapsed_hyper = time.perf_counter() - start_hyper
print(
    "Time for hyperparameter testing: {:.4f} minutes".format(
        elapsed_hyper / 60))
print(
    "Time for hyperparameter testing: {:.4f} hours".format(
        elapsed_hyper / 60 / 60))
