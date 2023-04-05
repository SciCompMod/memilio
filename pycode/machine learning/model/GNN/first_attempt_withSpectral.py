import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import pickle
import spektral

import numpy as np
import scipy.sparse as sp
import tensorflow as tf
from keras.layers import Dense
from keras.losses import MeanAbsolutePercentageError
from keras.metrics import mean_absolute_percentage_error
from keras.models import Model
from keras.optimizers import Adam

from spektral.data import Dataset, DisjointLoader, Graph, Loader, BatchLoader, MixedLoader
from spektral.layers import GCSConv, GlobalAvgPool, GlobalAttentionPool, ARMAConv, AGNNConv, APPNPConv, CrystalConv, GATConv, GINConv, XENetConv
from spektral.transforms.normalize_adj import NormalizeAdj


from memilio.simulation.secir import InfectionState

# add to setup: torch, spektral


path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.dirname(os.path.realpath(path)))),
    'data/data_GNN_nodamp_401pop_1k')

file = open(os.path.join(path_data, 'data_secir_age_groups.pickle'), 'rb')
data_secir = pickle.load(file)

len_dataset = data_secir['inputs'][0].shape[0]
numer_of_nodes = np.asarray(data_secir['inputs']).shape[0]
shape_input_flat = np.asarray(
    data_secir['inputs']).shape[2]*np.asarray(data_secir['inputs']).shape[3]
shape_labels_flat = np.asarray(
    data_secir['labels']).shape[2]*np.asarray(data_secir['labels']).shape[3]


new_inputs = np.asarray(
    data_secir['inputs']).reshape(
    len_dataset, numer_of_nodes, shape_input_flat)
new_labels = np.asarray(data_secir['labels']).reshape(
    len_dataset, numer_of_nodes, shape_labels_flat)

######## open commuter data #########
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(path,
                         'data')
commuter_file = open(os.path.join(
    path_data, 'commuter_migration_scaled.txt'), 'rb')
commuter_data = pd.read_csv(commuter_file, sep=" ", header=None)
sub_matrix = commuter_data.iloc[:numer_of_nodes, 0:numer_of_nodes]


adjacency_matrix = np.asarray(sub_matrix)

node_features = new_inputs

node_labels = new_labels


class MyDataset(spektral.data.dataset.Dataset):
    def read(self):
        # self.a = sp.csr_matrix(adjacency_matrix)
        self.a = adjacency_matrix

        return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features, node_labels)]

        super().__init__(**kwargs)


# data = MyDataset()
data = MyDataset(transforms=NormalizeAdj())
batch_size = 32
epochs = 200
es_patience = 20  # Patience for early stopping

# Train/valid/test split
idxs = np.random.permutation(len(data))
split_va, split_te = int(0.8 * len(data)), int(0.9 * len(data))
idx_tr, idx_va, idx_te = np.split(idxs, [split_va, split_te])
data_tr = data[idx_tr]
data_va = data[idx_va]
data_te = data[idx_te]

# # Data loaders
# loader_tr = DisjointLoader(data_tr, batch_size=batch_size, epochs=epochs)
# loader_va = DisjointLoader(data_va, batch_size=batch_size)
# loader_te = DisjointLoader(data_te, batch_size=batch_size)

# Data loaders
loader_tr = MixedLoader(
    data_tr, batch_size=batch_size, epochs=epochs)
loader_va = MixedLoader(data_va,  batch_size=batch_size)
loader_te = MixedLoader(data_te, batch_size=data_te.n_graphs)

# self.conv1 = GCSConv(32, activation="relu")
# self.conv2 = GCSConv(32, activation="relu")
# self.conv3 = GCSConv(32, activation="relu")
# #self.global_pool = GlobalAttentionPool(3)
# self.dense = Dense(data.n_labels, activation="linear")

# self.global_pool = GlobalAvgPool()
# x = self.global_pool([x])

################################################################################
# Build model
################################################################################


class Net(Model):
    def __init__(self):
        super().__init__()

        ##### for mixed mode ####
        self.conv1 = GATConv(32,     activation="relu")
        self.conv2 = GATConv(32,   activation="relu")
        self.conv3 = GATConv(32,   activation="relu")

        self.dense = Dense(data.n_labels, activation="linear")

    def call(self, inputs):
        x, a = inputs
        a = np.asarray(a)

        x = self.conv1([x, a])
        x = self.conv2([x, a])
        x = self.conv3([x, a])

        output = self.dense(x)

        return output


# decayed_lr = tf.keras.optimizers.schedules.ExponentialDecay(
#         initial_learning_rate=0.01,
#         decay_steps=200,
#         decay_rate=0.95,
#         staircase=True)
learning_rate = 0.001
model = Net()
optimizer = Adam(learning_rate=learning_rate)
loss_fn = MeanAbsolutePercentageError()


# @tf.function(input_signature=loader_tr.tf_signature(),
# experimental_relax_shapes=True)
def train_step(inputs, target):
    with tf.GradientTape() as tape:
        predictions = model(inputs, training=True)
        loss = loss_fn(target, predictions) + sum(model.losses)
    gradients = tape.gradient(loss, model.trainable_variables)
    optimizer.apply_gradients(zip(gradients, model.trainable_variables))
    acc = tf.reduce_mean(mean_absolute_percentage_error(target, predictions))
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


def test_evaluation(loader):

    inputs, target = loader.__next__()
    pred = model(inputs, training=False)

    mean_per_batch = []

    for batch_p, batch_t in zip(pred, target):
        MAPE_v = []
        for v_p, v_t in zip(batch_p, batch_t):

            pred_ = tf.reshape(v_p, (26, 48))
            target_ = tf.reshape(v_t, (26, 48))

            diff = pred_ - target_
            relative_err = (abs(diff))/abs(target_)
            relative_err_transformed = np.asarray(
                relative_err).transpose().reshape(8, -1)
            relative_err_means_percentage = relative_err_transformed.mean(
                axis=1) * 100

            MAPE_v.append(relative_err_means_percentage)

        mean_per_batch.append(np.asarray(MAPE_v).transpose().mean(axis=1))

    mean_percentage = pd.DataFrame(
        data=np.asarray(mean_per_batch).transpose().mean(axis=1),
        index=[str(compartment).split('.')[1]
               for compartment in InfectionState.values()],
        columns=['Percentage Error'])

    return mean_percentage


epoch = step = 0
best_val_loss = np.inf
best_weights = None
patience = es_patience
results = []
losses_history = []
val_losses_history = []
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
            print("New best val_loss {:.3f}".format(val_loss))
            best_weights = model.get_weights()
        else:
            patience -= 1
            if patience == 0:
                print("Early stopping (best val_loss: {})".format(best_val_loss))
                break
        results = []
        losses_history.append(loss)
        val_losses_history.append(val_loss)

################################################################################
# Evaluate model
################################################################################
model.set_weights(best_weights)  # Load best model
test_loss, test_acc = evaluate(loader_te)
test_MAPE = test_evaluation(loader_te)
print(test_MAPE)

print("Done. Test loss: {:.4f}. Test acc: {:.2f}".format(test_loss, test_acc))
