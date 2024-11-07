# from memilio.surrogatemodel.ode_secir_simple import network_architectures
# from memilio.surrogatemodel.ode_secir_simple.model import split_data, get_test_statistic
import os
import tensorflow as tf
import pickle
import pandas as pd
import time
from sklearn.model_selection import KFold
import numpy as np


path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'data_paper')

filename = "data_secir_groups_30days_I_based_Germany_10k_nodamp.pickle"
filename_df = "gridserach_secir_groups_30days_I_based_Germany_10k_nodamp.csv"

label_width = 30
early_stop = 100

# Define grid search parameters
hidden_layers = [0, 1, 2, 3, 4]
neurons_in_hidden_layer = [32, 64, 128, 512, 1024]
models = ["Dense", "CNN", "LSTM"]

parameters = [(layer, neuron_number, modelname)
              for layer in hidden_layers for neuron_number in neurons_in_hidden_layer for modelname in models]

# Create a DataFrame to store the results
df_results = pd.DataFrame(columns=['model', 'number_of_hidden_layers', 'number_of_neurons',
                                   'mean_test_MAPE', 'kfold_train', 'kfold_val',
                                   'kfold_test', 'training_time', 'train_losses', 'val_losses'])


if not os.path.isfile(os.path.join(path_data, filename)):
    raise ValueError(f"No dataset found in path: {path_data}")

with open(os.path.join(path_data, filename), 'rb') as file:
    data = pickle.load(file)

# Split the data: 80% will be used for grid search with cross-validation, and the remaining 20% is withheld for testing
inputs_grid_search = data['inputs'][:int((0.8 * len(data['inputs'])))]
labels_grid_search = data['labels'][:int((0.8 * len(data['labels'])))]
inputs_withhold = data['inputs'][int((0.8 * len(data['inputs']))):]
labels_withhold = data['labels'][int((0.8 * len(data['labels']))):]

# Function to train and evaluate the model using cross-validation


def train_and_evaluate_model(param, max_epochs):
    layer, neuron_number, modelname = param

    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                                      patience=early_stop,
                                                      mode='min')
    num_groups = 6
    kf = KFold(n_splits=5)
    train_losses = []
    val_losses = []
    rescaled_mape_scores = []  # New list to store rescaled MAPE for each fold

    losses_history_all = []
    val_losses_history_all = []

    start = time.perf_counter()

    # Custom metric for MAPE on rescaled data
    def rescaled_mape(y_true, y_pred):
        y_true_rescaled = np.expm1(y_true)
        y_pred_rescaled = np.expm1(y_pred)
        return 100 * np.mean(np.abs((y_true_rescaled - y_pred_rescaled) / y_true_rescaled))

    for train_idx, val_idx in kf.split(inputs_grid_search):
        # Model creation based on `modelname`
        if modelname == 'Dense':
            model = tf.keras.Sequential([tf.keras.layers.Flatten(),
                                         tf.keras.layers.Dense(units=neuron_number, activation='relu')])
            for _ in range(layer):
                model.add(tf.keras.layers.Dense(
                    units=neuron_number, activation='relu'))
            model.add(tf.keras.layers.Dense(
                units=label_width * 8 * num_groups))
            model.add(tf.keras.layers.Reshape([label_width, 8 * num_groups]))

        elif modelname == 'CNN':
            conv_size = 3
            num_outputs = 8 * num_groups
            model = tf.keras.Sequential([tf.keras.layers.Lambda(lambda x: x[:, -conv_size:, :]),
                                         tf.keras.layers.Conv1D(neuron_number, activation='relu', kernel_size=(conv_size))])
            for _ in range(layer):
                model.add(tf.keras.layers.Dense(
                    units=neuron_number, activation='relu'))
            model.add(tf.keras.layers.Dense(label_width * num_outputs,
                      kernel_initializer=tf.initializers.zeros()))
            model.add(tf.keras.layers.Reshape([label_width, num_outputs]))

        elif modelname == "LSTM":
            num_outputs = 8 * num_groups
            model = tf.keras.Sequential(
                [tf.keras.layers.LSTM(neuron_number, return_sequences=False)])
            for _ in range(layer):
                model.add(tf.keras.layers.Dense(
                    units=neuron_number, activation='relu'))
            model.add(tf.keras.layers.Dense(label_width * num_outputs,
                      kernel_initializer=tf.initializers.zeros()))
            model.add(tf.keras.layers.Reshape([label_width, num_outputs]))

        # Gather training and validation data based on the fold
        train_inputs = tf.gather(inputs_grid_search, indices=train_idx)
        train_labels = tf.gather(labels_grid_search, indices=train_idx)
        valid_inputs = tf.gather(inputs_grid_search, indices=val_idx)
        valid_labels = tf.gather(labels_grid_search, indices=val_idx)

        # Compile the model
        model.compile(loss=tf.keras.losses.MeanAbsolutePercentageError(),
                      optimizer=tf.keras.optimizers.Adam(),
                      metrics=[tf.keras.metrics.MeanAbsoluteError(), tf.keras.metrics.MeanAbsolutePercentageError()])

        # Train the model
        history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                            validation_data=(valid_inputs, valid_labels),
                            callbacks=[early_stopping])

        # Append minimum loss values
        train_losses.append(np.min(history.history['loss']))
        val_losses.append(np.min(history.history['val_loss']))
        losses_history_all.append(history.history['loss'])
        val_losses_history_all.append(history.history['val_loss'])

        # Calculate and store the rescaled MAPE for validation set
        val_pred = model.predict(valid_inputs)
        fold_rescaled_mape = rescaled_mape(valid_labels, val_pred)
        # Save each fold's rescaled MAPE
        rescaled_mape_scores.append(fold_rescaled_mape)

    elapsed = time.perf_counter() - start

    # Print out the results
    print(f"Best train losses: {train_losses}")
    print(f"Best validation losses: {val_losses}")
    print(f"Rescaled MAPE for each fold: {rescaled_mape_scores}")
    print("--------------------------------------------")
    print(f"K-Fold Train Score: {np.mean(train_losses)}")
    print(f"K-Fold Validation Score: {np.mean(val_losses)}")
    print(f"K-Fold Rescaled MAPE: {np.mean(rescaled_mape_scores)}")
    print(f"Time for training: {elapsed:.4f} seconds")
    print(f"Time for training: {elapsed / 60:.4f} minutes")

    # After cross-validation, save results to df_results
    df_results.loc[len(df_results.index)] = [
        modelname, layer, neuron_number, np.nan,  # Placeholder for test score
        np.mean(train_losses),
        np.mean(val_losses),
        np.mean(rescaled_mape_scores),  # Store the mean rescaled MAPE
        elapsed / 60,
        [losses_history_all],
        [val_losses_history_all]
    ]

    # Save df_results to CSV
    path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'secir_groups_grid_search_paper')
    if not os.path.isdir(file_path):
        os.mkdir(file_path)
    file_path = os.path.join(file_path, filename_df)
    #file_path = os.path.join('/localdata1/schm_a45/data_paper_November/')
    df_results.to_csv(file_path)


start_hyper = time.perf_counter()
max_epochs = 1500

for param in parameters:
    train_and_evaluate_model(param, max_epochs)

elapsed_hyper = time.perf_counter() - start_hyper
print(
    "Time for hyperparameter testing: {:.4f} minutes".format(
        elapsed_hyper / 60))
print(
    "Time for hyperparameter testing: {:.4f} hours".format(
        elapsed_hyper / 60 / 60))
