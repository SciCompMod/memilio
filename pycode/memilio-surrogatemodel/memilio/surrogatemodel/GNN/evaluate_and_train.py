#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Agatha Schmidt, Henrik Zunker, Manuel Heger
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

import os
import pickle
import spektral
import time
import pandas as pd
import numpy as np


from tensorflow.keras.optimizers import Adam
from tensorflow.keras.losses import MeanAbsolutePercentageError
import tensorflow.keras.initializers as initializers

import tensorflow as tf
import memilio.surrogatemodel.GNN.network_architectures as network_architectures
from memilio.surrogatemodel.utils.helper_functions import (calc_split_index)

from spektral.data import MixedLoader
from spektral.layers import (GCSConv, GlobalAvgPool, GlobalAttentionPool, ARMAConv, AGNNConv,
                             APPNPConv, CrystalConv, GATConv, GINConv, XENetConv, GCNConv, GCSConv)
from spektral.utils.convolution import normalized_laplacian, rescale_laplacian


def create_dataset(path_cases, path_mobility, number_of_nodes=400):
    """
    Create a dataset from the given file path.
    :param path_cases: Path to the pickle file containing case data.
    :param path_mobility: Path to the directory containing mobility data.
    .param number_of_nodes: Number of nodes in the graph.
    :return: A Spektral Dataset object containing the processed data.
    """

    # Load data from pickle file
    file = open(
        path_cases, 'rb')
    data = pickle.load(file)
    # Extract inputs and labels from the loaded data
    inputs = data['inputs']
    labels = data['labels']

    len_dataset = len(inputs)

    if len_dataset == 0:  # check if dataset is empty
        raise ValueError(
            "Dataset is empty. Please provide a valid dataset with at least one sample.")
    # Calculate the flattened shape of inputs and labels
    shape_input_flat = np.asarray(
        inputs).shape[1]*np.asarray(inputs).shape[2]
    shape_labels_flat = np.asarray(
        labels).shape[1]*np.asarray(labels).shape[2]
    # Reshape inputs and labels to the required format
    new_inputs = np.asarray(
        inputs.transpose(0, 3, 1, 2)).reshape(
        len_dataset, number_of_nodes, shape_input_flat)
    new_labels = np.asarray(labels.transpose(0, 3, 1, 2)).reshape(
        len_dataset, number_of_nodes, shape_labels_flat)

    # Load mobility data and create adjacency matrix
    commuter_file = open(os.path.join(
        path_mobility, 'commuter_mobility_2022.txt'), 'rb')
    commuter_data = pd.read_csv(commuter_file, sep=" ", header=None)
    sub_matrix = commuter_data.iloc[:number_of_nodes, 0:number_of_nodes]

    adjacency_matrix = np.asarray(sub_matrix)

    adjacency_matrix[adjacency_matrix > 0] = 1  # make the adjacency binary

    node_features = new_inputs

    node_labels = new_labels

    #  Define a custom Dataset class
    class MyDataset(spektral.data.dataset.Dataset):
        def read(self):
            self.a = adjacency_matrix
            return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features, node_labels)]

            super().__init__(**kwargs)
    # Instantiate the custom dataset
    data = MyDataset()

    return data


def train_step(inputs, target, loss_fn, model, optimizer):
    '''Perform a single training step.
    :param inputs: Tuple (x, a) where x is the node features and a is the adjacency matrix.
    :param target: Ground truth labels.
    :param loss_fn: Loss function to use.
    :param model: The GNN model to train.
    :param optimizer: Optimizer to use for training.
    :return: Loss and accuracy for the training step.'''
    # Record operations for automatic differentiation
    with tf.GradientTape() as tape:
        predictions = model(inputs, training=True)
        loss = loss_fn(target, predictions) + \
            sum(model.losses)  # Add regularization losses
    # Compute gradients and update model weights
    gradients = tape.gradient(loss, model.trainable_variables)
    optimizer.apply_gradients(zip(gradients, model.trainable_variables))
    acc = tf.reduce_mean(loss_fn(target, predictions))
    return loss, acc


def evaluate(loader, model, loss_fn, retransform=False):
    '''Evaluate the model on a validation or test set.
    :param loader: Data loader for the evaluation set.
    :param model: The GNN model to evaluate.
    :param loss_fn: Loss function to use.
    :param retransform: Whether to apply inverse transformation to the outputs.
    :return: Losses and mean_loss for the evaluation set.'''
    output = []
    step = 0
    while step < loader.steps_per_epoch:
        step += 1
        inputs, target = loader.__next__()
        pred = model(inputs, training=False)
        if retransform:
            target = np.expm1(target)
            pred = np.expm1(pred)
        # Calculate loss and metrics
        outs = (
            loss_fn(target, pred),
            tf.reduce_mean(loss_fn(target, pred)),
            len(target),  # Keep track of batch size
        )
        output.append(outs)
        # Aggregate results at the end of the epoch
        if step == loader.steps_per_epoch:
            output = np.array(output)
            return np.average(output[:, :-1], 0, weights=output[:, -1])


def train_and_evaluate(data, batch_size, epochs,  model, loss_fn, optimizer, es_patience, save_results=False, save_name="model"):
    '''Train and evaluate the GNN model.
    :param data: The dataset to use for training and evaluation.
    :param batch_size: Batch size for training.
    :param epochs: Maximum number of epochs to train.
    :param model: The GNN model to train and evaluate.
    :param loss_fn: Loss function to use.
    :param optimizer: Optimizer to use for training.
    :param es_patience: Patience for early stopping.
    :param save_results: Whether to save the results to a file.
    :param save_name: Name to use when saving the results.
    :return: A dictionary containing training and evaluation results if save_results is False.'''
    n = len(data)
    # Determine split indices for training, validation, and test sets
    n_train, n_valid, n_test = calc_split_index(
        n, split_train=0.7, split_valid=0.2, split_test=0.1)

    # Split data into train, validation, and test sets
    train_data, valid_data, test_data = data[:n_train], data[n_train:n_train +
                                                             n_valid], data[n_train+n_valid:]

    # Define Data Loaders
    loader_tr = MixedLoader(
        train_data, batch_size=batch_size, epochs=epochs, shuffle=False)
    loader_val = MixedLoader(
        valid_data, batch_size=batch_size, shuffle=False)
    loader_test = MixedLoader(
        test_data, batch_size=n_test, shuffle=False)

    df = pd.DataFrame(columns=[
        "train_loss", "val_loss", "test_loss",
        "test_loss_orig", "training_time",
        "loss_history", "val_loss_history"])

    # Initialize variables to track best scores and histories
    test_scores = []
    test_scores_r = []
    train_losses = []
    val_losses = []

    losses_history_all = []
    val_losses_history_all = []

    best_val_loss = np.inf
    best_weights = None
    patience = es_patience
    results = []
    losses_history = []
    val_losses_history = []

    ################################################################################
    # Train model
    ################################################################################
    start = time.perf_counter()
    step = 0
    epoch = 0
    for batch in loader_tr:
        step += 1
        loss, acc = train_step(*batch, loss_fn, model, optimizer)
        results.append((loss, acc))

        # Compute validation loss and accuracy  at the end of each epoch
        if step == loader_tr.steps_per_epoch:
            step = 0
            epoch += 1
            val_loss, val_acc = evaluate(loader_val, model, loss_fn)
            print(
                "Ep. {} - Loss: {:.3f} - Acc: {:.3f} - Val loss: {:.3f} - Val acc: {:.3f}".format(
                    epoch, *np.mean(results, 0), val_loss, val_acc
                )
            )
            # Check if loss improved for early stopping
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience = es_patience
                print(f"New best val_loss {val_loss:.3f}")
                best_weights = model.get_weights()
            else:
                patience -= 1
                if patience == 0:
                    print(
                        "Early stopping (best val_loss: {})".format(
                            best_val_loss))
                    break
            results = []
            losses_history.append(loss)
            val_losses_history.append(val_loss)
    ################################################################################
    # Evaluate model
    ################################################################################
    # Load best model weights before evaluation on test set
    model.set_weights(best_weights)
    test_loss, test_acc = evaluate(loader_test, model, loss_fn)
    test_loss_r, test_acc_r = evaluate(
        loader_test, model, loss_fn, retransform=True)

    print(
        "Done. Test loss: {:.4f}. Test acc: {:.2f}".format(
            test_loss, test_acc))
    test_scores.append(test_loss)
    test_scores_r.append(test_loss_r)
    train_losses.append(np.asarray(losses_history).min())
    val_losses.append(np.asarray(val_losses_history).min())
    losses_history_all.append(np.asarray(losses_history))
    val_losses_history_all.append(val_losses_history)

    elapsed = time.perf_counter() - start

    # print out stats
    print(f"Best train losses: {train_losses} ")
    print(f"Best validation losses: {val_losses}")
    print(f"Test values: {test_scores}")
    print("--------------------------------------------")
    print(f"Train Score:{np.mean(train_losses)}")
    print(f"Validation Score:{np.mean(val_losses)}")
    print(f"Test Score (log): {np.mean(test_scores)}")
    print(f"Test Score (orig.): {np.mean(test_scores_r)}")

    print(f"Time for training: {elapsed:.4f} seconds")
    print(f"Time for training: {elapsed/60:.4f} minutes")

    ################################################################################
    # Save results
    ################################################################################
    if save_results:
        # save df
        df.loc[len(df.index)] = [  # layer_name, number_of_layer, channels,
            np.mean(train_losses),
            np.mean(val_losses),
            np.mean(test_scores),
            np.mean(test_scores_r),
            (elapsed / 60),
            [losses_history_all],
            [val_losses_history_all]]
        print(df)

        # Ensure that save_name has the .pickle extension
        if not save_name.endswith('.pickle'):
            save_name += '.pickle'

        # Save best weights as pickle
        path = os.path.dirname(os.path.realpath(__file__))
        file_path_w = os.path.join(path, 'saved_weights')
        # Ensure the directory exists
        if not os.path.isdir(file_path_w):
            os.mkdir(file_path_w)

        # Construct the full file path by joining the directory with save_name
        file_path_w = os.path.join(file_path_w, save_name)

        # Save the weights to the file
        with open(file_path_w, 'wb') as f:
            pickle.dump(best_weights, f)

        file_path_df = os.path.join(
            os.path.dirname(
                os.path.realpath(os.path.dirname(os.path.realpath(path)))),
            'model_evaluations_paper')
        if not os.path.isdir(file_path_df):
            os.mkdir(file_path_df)
        file_path_df = file_path_df+"/"+save_name.replace('.pickle', '.csv')
        df.to_csv(file_path_df)

    else:
        return {
            "model": save_name,
            "mean_train_loss": np.mean(train_losses),
            "mean_val_loss": np.mean(val_losses),
            "mean_test_loss": np.mean(test_scores),
            "mean_test_loss_orig": np.mean(test_scores_r),
            "training_time": elapsed/60,
            "train_losses": losses_history_all,
            "val_losses": val_losses_history_all
        }


if __name__ == "__main__":
    start_hyper = time.perf_counter()
    epochs = 10
    batch_size = 2
    es_patience = 10
    optimizer = Adam(learning_rate=0.001)
    loss_fn = MeanAbsolutePercentageError()

    # Generate the Dataset
    path_cases = "/localdata1/hege_mn/memilio/saves/GNN_data_30days_3dampings_classic5.pickle"
    path_mobility = '/localdata1/hege_mn/memilio/data/Germany/mobility'
    data = create_dataset(path_cases, path_mobility)

    # Define the model architecture
    def transform_a(adjacency_matrix):
        a = adjacency_matrix.numpy()
        a = rescale_laplacian(normalized_laplacian(a))
        return tf.convert_to_tensor(a, dtype=tf.float32)

    layer_types = [
        # Dense layer (only uses x)
        lambda: ARMAConv(512, activation='elu',
                         kernel_initializer=initializers.GlorotUniform(seed=None))
    ]
    num_repeat = [7]

    model_class = network_architectures.generate_model_class(
        "ARMA", layer_types, num_repeat, num_output=1440, transform=transform_a)

    model = model_class()

    # early stopping patience
    # name for df
    save_name = 'GNN_30days'  # name for model
    # for param in parameters:

    train_and_evaluate(data, batch_size, epochs, model,
                       loss_fn, optimizer, es_patience, save_name)

    elapsed_hyper = time.perf_counter() - start_hyper
    print(
        "Time for hyperparameter testing: {:.4f} minutes".format(
            elapsed_hyper / 60))
    print(
        "Time for hyperparameter testing: {:.4f} hours".format(
            elapsed_hyper / 60 / 60))
