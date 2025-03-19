import tensorflow as tf
import numpy as np
import pandas as pd
import pickle
import time
from functools import partial
import concurrent.futures
from keras.activations import relu, elu, softmax, sigmoid, linear, tanh
from sklearn.model_selection import KFold
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
tf.config.threading.set_intra_op_parallelism_threads(1)
tf.config.threading.set_inter_op_parallelism_threads(1)

# Annahme: Die Funktion train_and_evaluate_model ist wie folgt definiert:


def train_and_evaluate_model(param, max_epochs, filename, filename_df):
    act = param[0]
    optimizer = param[1]
    activation_s = param[2]

    with open('/localdata1/gnn_paper_2024/data/one_population/with_agegroups_Germany/nodamp/data_secir_groups_30days_I_based_Germany_10k_nodamp.pickle', 'rb') as file:
        data = pickle.load(file)
    inputs_grid_search = data['inputs'][:int(0.8 * len(data['inputs']))]
    labels_grid_search = data['labels'][:int(0.8 * len(data['labels']))]

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
        model = tf.keras.Sequential([
            tf.keras.layers.LSTM(1024, return_sequences=False),
            tf.keras.layers.Dense(units=1024, activation=act),
            tf.keras.layers.Dense(label_width * 8 * 6,
                                  kernel_initializer=tf.initializers.zeros()),
            tf.keras.layers.Reshape([label_width, 8 * 6])
        ])

        train_inputs = tf.gather(inputs_grid_search, indices=train_idx)
        train_labels = tf.gather(labels_grid_search, indices=train_idx)
        valid_inputs = tf.gather(inputs_grid_search, indices=val_idx)
        valid_labels = tf.gather(labels_grid_search, indices=val_idx)

        model.compile(
            loss=tf.keras.losses.MeanAbsolutePercentageError(),
            optimizer=optimizer,
            metrics=[tf.keras.metrics.MeanAbsoluteError(),
                     tf.keras.metrics.MeanAbsolutePercentageError()]
        )

        history = model.fit(train_inputs, train_labels, epochs=max_epochs,
                            validation_data=(valid_inputs, valid_labels),
                            callbacks=[early_stopping],
                            verbose=0)  # Verbose kann angepasst werden

        train_losses.append(np.min(history.history['loss']))
        val_losses.append(np.min(history.history['val_loss']))
        losses_history_all.append(history.history['loss'])
        val_losses_history_all.append(history.history['val_loss'])

    elapsed = time.perf_counter() - start

    print(f"Parameter {activation_s} und {optimizer}: Train loss {np.mean(train_losses):.4f}, Val loss {np.mean(val_losses):.4f}, Zeit {elapsed/60:.2f} min")

    result = {
        'model': modelname,
        'optimizer': optimizer,
        'activation': activation_s,
        'mean_test_MAPE': np.nan,
        'kfold_train': np.mean(train_losses),
        'kfold_val': np.mean(val_losses),
        'kfold_test': np.nan,
        'training_time': elapsed / 60,
        'train_losses': losses_history_all,
        'val_losses': val_losses_history_all
    }
    return result


def worker(param, max_epochs, filename, filename_df):
    return train_and_evaluate_model(param, max_epochs, filename, filename_df)


if __name__ == '__main__':
    max_epochs = 1500
    filename = "data_secir_groups_30days_I_based_Germany_10k_nodamp.pickle"
    filename_df = "groups_I_based_dataframe_optimizer_paper.csv"
    modelname = 'LSTM'
    label_width = 30
    early_stop = 100

    activations_str = ['relu', 'elu', 'softmax', 'sigmoid', 'linear', 'tanh']
    activations = [relu, elu, softmax, sigmoid, linear, tanh]
    optimizers = ['Adam', 'Nadam', 'SGD']
    parameters = []
    for a in activations:
        for o in optimizers:
            for astr in activations_str:
                parameters.append((a, o, astr))

    worker_partial = partial(
        worker, max_epochs=max_epochs, filename=filename, filename_df=filename_df)

    start_hyper = time.perf_counter()
    with concurrent.futures.ProcessPoolExecutor(max_workers=18) as executor:
        results = list(executor.map(worker_partial, parameters))

    df_results = pd.DataFrame(results)
    file_path = os.path.join("/localdata1/zunk_he/memilio/saves", filename_df)
    df_results.to_csv(file_path, index=False)

    elapsed_hyper = time.perf_counter() - start_hyper
    print("Time for hyperparameter testing: {:.4f} minutes".format(
        elapsed_hyper / 60))
    print("Time for hyperparameter testing: {:.4f} hours".format(
        elapsed_hyper / 3600))
