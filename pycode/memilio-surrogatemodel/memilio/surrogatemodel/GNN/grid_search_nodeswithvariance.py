import os
import pickle
import spektral
import pandas as pd
import numpy as np
from sklearn.preprocessing import FunctionTransformer

# from first_attempt_withSpectral_withCV import preprocess, train_step, evaluate, test_evaluation
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
# import torch
# import torch.nn as nn
# import torch.nn.functional as F
import pickle
import spektral
import time

import numpy as np
import scipy.sparse as sp
import tensorflow as tf
from keras.layers import Dense
from keras.losses import MeanAbsolutePercentageError
from keras.metrics import mean_absolute_percentage_error
from keras.models import Model
# from keras.optimizers.legacy import Adam, Nadam, RMSprop, SGD, Adagrad
from keras.optimizers import Adam, Nadam, RMSprop, SGD, Adagrad

from sklearn.model_selection import KFold

from spektral.data import Dataset, DisjointLoader, Graph, Loader, BatchLoader, MixedLoader
from spektral.layers import GCSConv, GlobalAvgPool, GlobalAttentionPool, ARMAConv, AGNNConv, APPNPConv, CrystalConv, GATConv, GINConv, XENetConv, GCNConv, GCSConv
from spektral.transforms.normalize_adj import NormalizeAdj
from spektral.utils.convolution import gcn_filter, normalized_laplacian, rescale_laplacian, normalized_adjacency

tf.config.threading.set_intra_op_parallelism_threads(16)
tf.config.threading.set_inter_op_parallelism_threads(16)


# from memilio.simulation.secir import InfectionState


# load and prepare data
# path = os.path.dirname(os.path.realpath(__file__))
# path_data = os.path.join(
#    os.path.dirname(
#        os.path.realpath(os.path.dirname(os.path.realpath(path)))),
#   'data/data_GNN_nodamp_400pop_30days_2024_oldrange')

# file = open(os.path.join(path_data, 'data_secir_age_groups.pickle'), 'rb')
file = open(
    '/localdata1/gnn_paper_2024/data/GNNs/GNN_data_90days_nodeswithvariance_3damp_1k.pickle', 'rb')

data = pickle.load(file)


# len_dataset = data_secir['inputs'][0].shape[0]
# len_dataset = 800*5
# numer_of_nodes = np.asarray(data_secir['inputs']).shape[0]
numer_of_nodes = 400


# Split the data: 80% will be used for grid search with cross-validation, and the remaining 20% is withheld for testing
inputs_grid_search = data['inputs'][:int(0.8 * len(data['inputs']))]
labels_grid_search = data['labels'][:int(0.8 * len(data['labels']))]
inputs_withhold = data['inputs'][int(0.8 * len(data['inputs'])):]
labels_withhold = data['labels'][int(0.8 * len(data['labels'])):]

len_dataset = len(inputs_grid_search)

shape_input_flat = np.asarray(
    inputs_grid_search).shape[1]*np.asarray(inputs_grid_search).shape[2]
shape_labels_flat = np.asarray(
    labels_grid_search).shape[1]*np.asarray(labels_grid_search).shape[2]


new_inputs = np.asarray(
    inputs_grid_search.transpose(0, 3, 1, 2)).reshape(
    len_dataset, numer_of_nodes, shape_input_flat)
new_labels = np.asarray(labels_grid_search.transpose(0, 3, 1, 2)).reshape(
    len_dataset, numer_of_nodes, shape_labels_flat)

n_days = int(new_labels.shape[2]/48)

path = os.path.dirname(os.path.realpath(__file__))
path_data = '/localdata1/zunk_he/memilio/data/mobility'
commuter_file = open(os.path.join(
    path_data, 'commuter_mobility.txt'), 'rb')
commuter_data = pd.read_csv(commuter_file, sep=" ", header=None)
sub_matrix = commuter_data.iloc[:numer_of_nodes, 0:numer_of_nodes]


adjacency_matrix = np.asarray(sub_matrix)

adjacency_matrix[adjacency_matrix > 0] = 1  # make the adjacency binary
node_features = new_inputs

node_labels = new_labels


layers = [GCNConv, APPNPConv, GATConv]
# layers = [ARMAConv]
number_of_layers = [1, 2, 3]
# number_of_layers = [2,3,4]
number_of_neurons = [32, 128, 512, 1024]


parameters = []
for l in layers:
    for nl in number_of_layers:
        for nn in number_of_neurons:
            parameters.append((l, nl, nn))

df = pd.DataFrame(
    columns=['layer', 'number_of_layers', 'number_of_neurons',
             'learning_rate', 'kfold_train',
             'kfold_val', 'kfold_test', 'training_time', 'all_testscores', 'test_rescaled', 'all_val_scores', 'all_trainscores'])


def train_and_evaluate_model(
        epochs, learning_rate, param):

    layer, number_of_layer, number_of_n = param

    class MyDataset(spektral.data.dataset.Dataset):
        def read(self):
            if layer == GCNConv:
                self.a = gcn_filter(adjacency_matrix)
            elif layer == ARMAConv:
                self.a = rescale_laplacian(
                    normalized_laplacian(adjacency_matrix))
            elif layer == APPNPConv:
                self.a = gcn_filter(adjacency_matrix)
            elif layer == GATConv:
                self. a = adjacency_matrix

            # self.a = normalized_adjacency(adjacency_matrix)
            # self.a = rescale_laplacian(normalized_laplacian(adjacency_matrix))

            return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features, node_labels)]

            super().__init__(**kwargs)

    # data = MyDataset()
    data = MyDataset(transforms=NormalizeAdj())
    batch_size = 32
    epochs = epochs
    es_patience = 100  # Patience for early stopping

    if number_of_layer == 1:
        class Net(Model):
            def __init__(self):
                super().__init__()
                self.conv1 = layer(number_of_n, activation='relu')
                self.dense = Dense(data.n_labels, activation="linear")

            def call(self, inputs, mask=None):
                x, a = inputs
                a = tf.convert_to_tensor(a, dtype=tf.float32)

                x = self.conv1([x, a])

                output = self.dense(x)

                return output

    elif number_of_layer == 2:
        class Net(Model):
            def __init__(self):
                super().__init__()
                self.conv1 = layer(number_of_n, activation='relu')
                self.conv2 = layer(number_of_n, activation='relu')
                self.dense = Dense(data.n_labels, activation="linear")

            def call(self, inputs, mask=None):
                x, a = inputs
                a = np.asarray(a)

                x = self.conv1([x, a])
                x = self.conv2([x, a])

                output = self.dense(x)

                return output

    elif number_of_layer == 3:
        class Net(Model):
            def __init__(self):
                super().__init__()

                self.conv1 = layer(number_of_n, activation='relu')
                self.conv2 = layer(number_of_n, activation='relu')
                self.conv3 = layer(number_of_n, activation='relu')
                self.dense = Dense(data.n_labels, activation="linear")

            def call(self, inputs, mask=None):
                x, a = inputs
                a = np.asarray(a)

                x = self.conv1([x, a])
                x = self.conv2([x, a])
                x = self.conv3([x, a])

                output = self.dense(x)
                return output

    elif number_of_layer == 4:
        class Net(Model):
            def __init__(self):
                super().__init__()

                self.conv1 = layer(number_of_n, activation='relu')
                self.conv2 = layer(number_of_n, activation='relu')
                self.conv3 = layer(number_of_n, activation='relu')
                self.conv4 = layer(number_of_n, activation='relu')
                self.dense = Dense(data.n_labels, activation="linear")

            def call(self, inputs, mask=None):
                x, a = inputs
                a = np.asarray(a)

                x = self.conv1([x, a])
                x = self.conv2([x, a])
                x = self.conv3([x, a])
                x = self.conv4([x, a])

                output = self.dense(x)
                return output

    elif number_of_layer == 5:
        class Net(Model):
            def __init__(self):
                super().__init__()

                self.conv1 = layer(number_of_n, activation='relu')
                self.conv2 = layer(number_of_n, activation='relu')
                self.conv3 = layer(number_of_n, activation='relu')
                self.conv4 = layer(number_of_n, activation='relu')
                self.conv5 = layer(number_of_n, activation='relu')
                self.dense = Dense(data.n_labels, activation="linear")

            def call(self, inputs, mask=None):
                x, a = inputs
                a = np.asarray(a)

                x = self.conv1([x, a])
                x = self.conv2([x, a])
                x = self.conv3([x, a])
                x = self.conv4([x, a])
                x = self.conv5([x, a])

                output = self.dense(x)

                return output

    # optimizer = Adam(learning_rate=learning_rate)
    optimizer = Adam(learning_rate=learning_rate)
    loss_fn = MeanAbsolutePercentageError()

    def train_step(inputs, target):
        optimizer = Adam(learning_rate=learning_rate)
        with tf.GradientTape() as tape:
            predictions = model(inputs, training=True)
            loss = loss_fn(target, predictions) + sum(model.losses)
        gradients = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(gradients, model.trainable_variables))
        acc = tf.reduce_mean(
            mean_absolute_percentage_error(
                target, predictions))
        return loss, acc

    def evaluate(loader):
        output = []
        step = 0
        while step < loader.steps_per_epoch:
            step += 1
            inputs, target = loader.__next__()
            pred = model(inputs, training=False)
            outs = (
                loss_fn(target, pred),
                tf.reduce_mean(mean_absolute_percentage_error(target, pred)),
                len(target),  # Keep track of batch size
            )
            output.append(outs)
            if step == loader.steps_per_epoch:
                output = np.array(output)
                return np.average(output[:, :-1], 0, weights=output[:, -1])

    def evaluate_r(loader):
        output = []
        step = 0
        while step < loader.steps_per_epoch:
            step += 1
            inputs, target = loader.__next__()
            pred = model(inputs, training=False)
            outs = (
                loss_fn(np.expm1(target), np.expm1(pred)),
                tf.reduce_mean(mean_absolute_percentage_error(
                    np.expm1(target), np.expm1(pred))),
                len(target),  # Keep track of batch size
            )
            output.append(outs)
            if step == loader.steps_per_epoch:
                output = np.array(output)
                return np.average(output[:, :-1], 0, weights=output[:, -1])

    n_days = int(new_labels.shape[2]/48)

    kf = KFold(n_splits=5)
    train_idxs = []
    test_idxs = []
    for i, (train_index, test_index) in enumerate(kf.split(data)):

        train_idxs.append(train_index)
        test_idxs.append(test_index)

    test_scores = []
    train_losses = []
    test_rescaled = []
    val_losses = []

    losses_history_all = []
    val_losses_history_all = []

    start = time.perf_counter()
    for train_idx, test_idx in zip(train_idxs, test_idxs):
        model = Net()

        data_tr = data[train_idx[:(int(0.8*len(train_idx)))]]
        data_va = data[train_idx[(int(0.8*len(train_idx))):]]
        data_te = data[test_idx]

        # Data loaders
        loader_tr = MixedLoader(
            data_tr, batch_size=batch_size, epochs=epochs)
        loader_va = MixedLoader(data_va,  batch_size=batch_size)
        loader_te = MixedLoader(data_te, batch_size=data_te.n_graphs)

        epoch = step = 0
        best_val_loss = np.inf
        best_weights = None
        patience = es_patience
        results = []
        losses_history = []
        val_losses_history = []

        start = time.perf_counter()

        for batch in loader_tr:
            step += 1
            loss, acc = train_step(*batch)
            results.append((loss, acc))
            if step == loader_tr.steps_per_epoch:
                step = 0
                epoch += 1

                # Compute validation loss and accuracy
                val_loss, val_acc = evaluate(loader_va)
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
        model.set_weights(best_weights)  # Load best model
        test_loss, test_acc = evaluate(loader_te)
        test_mape_rescaled = evaluate_r(loader_te)
        # test_MAPE = test_evaluation(loader_te)
        # print(test_MAPE)

        print(
            "Done. Test loss: {:.4f}. Test acc: {:.2f}".format(
                test_loss, test_acc))
        test_scores.append(test_loss)
        test_rescaled.append(test_mape_rescaled)
        train_losses.append(np.asarray(losses_history).min())
        val_losses.append(np.asarray(val_losses_history).min())
        losses_history_all.append(losses_history)
        val_losses_history_all.append(val_losses_history)

        elapsed = time.perf_counter() - start

    # plot the losses
    # plt.figure()
    # plt.plot(np.asarray(losses_history_all).mean(axis=0), label='train loss')
    # plt.plot(np.asarray(val_losses_history_all).mean(axis=0), label='val loss')
    # plt.xlabel('Epoch')
    # plt.ylabel('Loss ( MAPE)')
    # plt.title('Loss for' + str(layer))
    # plt.legend()
    # plt.savefig('losses'+str(layer)+'.png')

    # print out stats
    print(f"Best train losses: {train_losses} ")
    print(f"Best validation losses: {val_losses}")
    print(f"Test values: {test_scores}")
    print("--------------------------------------------")
    print(f"K-Fold Train Score:{np.mean(train_losses)}")
    print(f"K-Fold Validation Score:{np.mean(val_losses)}")
    print(f"K-Fold Test Score: {np.mean(test_scores)}")
    print(f"K-Fold ORIGINAL DATA Test Score: {np.mean(test_rescaled)}")

    print(f"Time for training: {elapsed:.4f} seconds")
    print(f"Time for training: {elapsed/60:.4f} minutes")

    df.loc[len(df.index)] = [layer, number_of_layer, number_of_n,
                             learning_rate, np.mean(train_losses),
                             np.mean(val_losses),
                             np.mean(test_scores),
                             (elapsed / 60), test_scores, test_rescaled, val_losses, train_losses]
    # [np.asarray(losses_history_all).mean(axis=0)],
    # [np.asarray(val_losses_history_all).mean(axis=0)]]

    path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'model_evaluations_paper')
    if not os.path.isdir(file_path):
        os.mkdir(file_path)
    file_path = file_path+filename
    df.to_csv(file_path)


start_hyper = time.perf_counter()
epochs = 1500
filename = '/GNN_gridsearch_withCV_90days_nodeswithvariance_3damp.csv'
for param in parameters:
    train_and_evaluate_model(epochs, 0.001, param)

elapsed_hyper = time.perf_counter() - start_hyper
print(
    "Time for hyperparameter testing: {:.4f} minutes".format(
        elapsed_hyper / 60))
print(
    "Time for hyperparameter testing: {:.4f} hours".format(
        elapsed_hyper / 60 / 60))
