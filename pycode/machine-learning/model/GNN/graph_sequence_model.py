# Panagopolus Transfer GNN for Pandemic Forecasting
# Given a sequence of Graphs G1, G2,... Gt that corresponds to a sequence of dates they utilize a MPNN ar each time step to obtain a
# sequence of representations H1, H2, ...Ht. These representations are then fed into a LSTM network  (two LSTM layers)
# then H1,...Ht and the G1,...Gt are passed to an output layer similar to MPNN


# "Following Kipf and Welling (2017), we normalize the adjacency matrix A such that the
# sum of the weights of the incoming edges of each node is equal to 1."


import os
import pickle
import spektral
import torch

import torch.nn as nn
import numpy as np
import pandas as pd

from spektral.layers.convolutional import GCNConv
# from torch_geometric.nn import GCNConv
import time

import tensorflow as tf
from keras.layers import Dense
from keras.losses import MeanAbsolutePercentageError, MeanSquaredError
from keras.metrics import mean_absolute_percentage_error
from keras.models import Model
from keras.optimizers import Adam

from sklearn.model_selection import KFold

from spektral.data import BatchLoader, MixedLoader

from data_generator_withoutdampings_all401 import get_population
from memilio.simulation.secir import InfectionState


# GAR
from keras.layers import TimeDistributed, Input, LSTM, Lambda, Dense, Reshape, Concatenate, Dropout, Flatten, Reshape
from keras.models import Model
from keras.regularizers import l2
from sklearn.model_selection import train_test_split
from spektral.layers import GlobalAttentionPool


from sklearn.preprocessing import FunctionTransformer

######### Load data #############################################
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.dirname(os.path.realpath(path)))),
    'data/data_GNN_nodamp_400pop_1k_30days')

file = open(os.path.join(path_data, 'data_secir_age_groups.pickle'), 'rb')
data_secir = pickle.load(file)

len_dataset = data_secir['inputs'][0].shape[0]
number_of_nodes = np.asarray(data_secir['inputs']).shape[0]
num_days = np.asarray(data_secir['inputs']).shape[2]
num_features = np.asarray(data_secir['inputs']).shape[3]

shape_labels_flat = np.asarray(
    data_secir['labels']).shape[2]*np.asarray(data_secir['labels']).shape[3]


new_inputs = np.asarray(
    data_secir['inputs']).reshape(
    len_dataset, number_of_nodes, num_days, num_features)
new_labels = np.asarray(data_secir['labels']).reshape(
    len_dataset, number_of_nodes, shape_labels_flat)

# open commuter data
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(path,
                         'data')
commuter_file = open(os.path.join(
    path_data, 'commuter_migration_scaled.txt'), 'rb')
commuter_data = pd.read_csv(commuter_file, sep=" ", header=None)
sub_matrix = commuter_data.iloc[:number_of_nodes, 0:number_of_nodes]


# adjacency_matrix_01 = adjacency_matrix[adjacency_matrix > 0] = 1
# adjacency_matrix_01 = sub_matrix[sub_matrix > 0] = 1
adjacency_matrix = np.asarray(sub_matrix.copy())
# adjacency_matrix[adjacency_matrix > 0] = 1
node_labels = new_labels

# edge features
df = pd.DataFrame(data=np.asarray(sub_matrix))

# df = pd.DataFrame(data=adjacency_list)
edge_features = df.rename_axis('Source')\
    .reset_index()\
    .melt('Source', value_name='Weight', var_name='Target')\
    .reset_index(drop=True)

edge_features = edge_features[
    edge_features.Weight != 0]

transformer = FunctionTransformer(np.log1p, validate=True)
scaled_edges = transformer.transform(
    edge_features['Weight'].values.reshape(-1, 1))


########### create graph sequence #########################################
# input sequence
datasets_input = []
input_array = []
for i in range(len_dataset):
    input_i = new_inputs[i].reshape(
        new_inputs[i].shape[1],
        new_inputs[i].shape[0],
        new_inputs[i].shape[2])
    node_features_i = input_i
    input_array.append(node_features_i)

    class MyDataset(spektral.data.dataset.Dataset):
        def read(self):
            self.a = spektral.utils.convolution.gcn_filter(adjacency_matrix)
            # self.a = np.asarray(adjacency_list)
            # self.e = scaled_edges.reshape(np.asarray(edge_features['Weight']).shape[0], 1)
            # self.e = scaled_edges
            # return [spektral.data.Graph(x=x, e=self.e, a=self.a, y=y) for x, y in zip(node_features_i, node_labels)]
            return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features_i, node_labels)]

            super().__init__(**kwargs)

    data = MyDataset()

    datasets_input.append(data)


############## define the model ###########################
# for spektral gcn layer

class MPNN(Model):
    def __init__(self, nfeat, nhid, nout, n_nodes, dropout):
        super().__init__()
        # self.n_nodes = n_nodes

        # self.batch_size = batch_size
        self.nhid = nhid

        self.conv1 = GCNConv(nfeat, activation='relu')
        self.conv2 = GCNConv(nfeat, activation='relu')
        self.dense1 = Dense(nfeat, activation='linear')
        self.dense2 = Dense(nfeat, activation='linear')

        self.lstm = LSTM(nout)

        self.dense3 = Dense(nout*n_nodes, activation='linear')
        self.reshape = Reshape([n_nodes, nout])

    def call(self, inputs):
        lst = list()
        # print(x.shape)
        # print(adj.shape)
        graphs, a, e, = inputs
        for graph in graphs:
            x = tf.convert_to_tensor(graph)

            # x, a, e = inputs
            # a = spektral.utils.convolution.gcn_filter(a)
            # lst.append(ident)

            # x = x[:,mob_feats]x

            # x = xt.index_select(1, mob_feats)
            # lst.append(tf.convert_to_tensor(x))

            x = self.conv1([x, a])
            # lst.append(x)

            x = self.conv2([x, a])
            # lst.append(x)

            # x = tf.concat(lst, axis=1)

            x = (self.dense1(x))

            x = (self.dense2(x))
            lst.append(x)

        x = self.lstm(np.asarray(lst))

        x = self.dense3(np.expand_dims((np.asarray(x).flatten()), axis=0))
        x = self.reshape(x)

        # np als outout ein problem ?? -> sorgt das für das gradient problem ?
        return np.squeeze(x)

#######################################################################################


learning_rate = 0.001

optimizer = Adam(learning_rate=learning_rate)
loss_fn = MeanAbsolutePercentageError()


# @tf.function(input_signature=loader_tr.tf_signature(),
# experimental_relax_shapes=True)
def train_step(inputs, target):
    with tf.GradientTape() as tape:
        predictions = model(inputs)
        loss = loss_fn(target, predictions) + sum(model.losses)
    gradients = tape.gradient(loss, model.trainable_variables)
    optimizer.apply_gradients(zip(gradients, model.trainable_variables))
    acc = tf.reduce_mean(mean_absolute_percentage_error(target, predictions))
    return loss, acc


def evaluate(inputs, target):
    output = []

    pred = model(inputs, training=False)
    outs = (
        loss_fn(target, pred),
        tf.reduce_mean(mean_absolute_percentage_error(target, pred)),
        len(target),  # Keep track of batch size
    )
    output.append(outs)

    output = np.array(output)
    return np.average(output[:, :-1], 0, weights=output[:, -1])


def test_evaluation(inputs, target, n_days):

    pred = model(inputs, training=False)

    mean_per_batch = []

    for batch_p, batch_t in zip(pred, target):
        MAPE_v = []
        for v_p, v_t in zip(batch_p, batch_t):

            pred_ = tf.reshape(v_p, (n_days, 48))
            target_ = tf.reshape(v_t, (n_days, 48))

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


########## instead if dataloader: transform graph data to necessary form ############
node_features = np.asarray(input_array)
adjacency_matrix_array = np.repeat(
    np.dstack([adjacency_matrix]).transpose(),
    len_dataset, axis=0)
edge_features = np.expand_dims(
    np.repeat(np.dstack([adjacency_matrix]).transpose(),
              len_dataset, axis=0),
    axis=0).reshape(
    len_dataset, number_of_nodes, number_of_nodes, 1)
input_tuple = (node_features, adjacency_matrix_array, edge_features)


# Train/valid/test split
idxs = np.random.permutation(len(node_labels))
split_va, split_te = int(0.8 * len(node_labels)), int(0.9 * len(node_labels))
idx_tr, idx_va, idx_te = np.split(idxs, [split_va, split_te])

data_tr_labels = node_labels[idx_tr]
data_va_labels = node_labels[idx_va]
data_te_labels = node_labels[idx_te]

data_tr_features = node_features[idx_tr]
data_va_features = node_features[idx_va]
data_te_features = node_features[idx_te]

data_tr_input_tuple = (
    data_tr_features, adjacency_matrix_array[: split_va],
    edge_features[: split_va])
data_va_input_tuple = (
    data_va_features, adjacency_matrix_array[split_va: split_te],
    edge_features[split_va: split_te])
data_te_input_tuple = (
    data_te_features, adjacency_matrix_array[split_te:],
    edge_features[split_te:])


n_nodes = node_features.shape[2]
# n_samples =4
# n_features =5
dropout = 0.3
nhid = 64
nfeat = node_features.shape[3]
nout = data.n_labels
n_days = node_labels.shape[1]

epochs = 10
es_patience = 10
lr = 0.001
batch_size = 32
steps_per_epoch = 23
loss_fn = MeanAbsolutePercentageError()
optimizer = Adam(lr=lr)


# input_tuple_i = [input_tuple[0][0], input_tuple[1][0], input_tuple[2][0]]


model = MPNN(nfeat, nhid, nout, n_nodes, dropout)

# pred = model(input_tuple_i)


start = time.perf_counter()


epoch = step = 0
best_val_loss = np.inf
best_weights = None
patience = es_patience
results = []
losses_history = []
val_losses_history = []
# start = time.perf_counter()
for x_t, a_t, e_t, labels_t, x_v, a_v, e_v, labels_v in zip(
        data_tr_input_tuple[0],
        data_tr_input_tuple[1],
        data_tr_input_tuple[2],
        data_tr_labels,
        data_va_input_tuple[0],
        data_va_input_tuple[1],
        data_va_input_tuple[2],
        data_va_labels):

    step += 1
    inputs_tr = (x_t, a_t, e_t)
    loss, acc = train_step(inputs_tr, labels_t)
    results.append((loss, acc))
    if step == steps_per_epoch:
        step = 0
        epoch += 1

        # Compute validation loss and accuracy
        inputs_va = (x_v, a_v, e_v)
        val_loss, val_acc = evaluate(inputs_va, labels_v)
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
                print(
                    "Early stopping (best val_loss: {})".format(
                        best_val_loss))
                break
        results = []
        losses_history.append(loss)
        val_losses_history.append(val_loss)
elapsed = time.perf_counter() - start
################################################################################
# Evaluate model
################################################################################

start = time.perf_counter()
model.set_weights(best_weights)  # Load best model


test_loss, test_acc = evaluate(data_te_input_tuple)

test_MAPE = test_evaluation(data_te_input_tuple, n_days)
print(
    "Done. Test loss: {:.4f}. Test acc: {:.2f}".format(
        test_loss, test_acc))

elapsed_test = time.perf_counter() - start


print("Time for training: {:.4f} seconds".format(elapsed))
print("Time for training: {:.4f} minutes".format(elapsed/60))

print("Time for testing: {:.4f} seconds".format(elapsed_test))
print("Time for testing: {:.4f} minutes".format(elapsed_test/60))


# # split the data

# kf = KFold(n_splits=4)
# train_idxs = []
# test_idxs = []
# for i, (train_index, test_index) in enumerate(kf.split(input_tuple[0])):

#     train_idxs.append(train_index)
#     test_idxs.append(test_index)


# for train_idx, test_idx in zip(train_idxs, test_idxs):

#     data_tr_input = (input_tuple[0][train_idx[:(int(0.8*len(train_idx)))]], input_tuple[1][train_idx[:(int(0.8*len(train_idx)))]], input_tuple[2][train_idx[:(int(0.8*len(train_idx)))]])
#     data_va_input = (input_tuple[0][train_idx[(int(0.8*len(train_idx))):]], input_tuple[1][train_idx[:(int(0.8*len(train_idx)))]], input_tuple[2][train_idx[:(int(0.8*len(train_idx)))]])
#     data_te_input = (input_tuple[0][test_idx], input_tuple[1][test_idx], input_tuple[2][test_idx])


#     data_tr_label = node_labels[train_idx[:(int(0.8*len(train_idx)))]]
#     data_va_label = node_labels[train_idx[(int(0.8*len(train_idx))):]]
#     data_te_label = node_labels[test_idx]


######################################
# x, a, e, = input_tuple
# x_t = torch.tensor(x[0])
# a_t = torch.tensor(a[0])
# e_t = torch.tensor(e[0])
# pred = model(x_t, a_t, e_t)


################# NGAR                      #############################

# N = 20                      # Number of nodes
# F = 48                       # Dimension of the node attributes
# dropout_rate = 0.0
# l2_reg = 5e-4              # Weight for l2 regularization
# X_in = data_tr_input[0][0]
# filter_in = data_tr_input[0][1]
# batch_size = 32
# mape = MeanAbsolutePercentageError()

# # Convolutional block
# conc1 = Concatenate()([X_in, filter_in])
# gc1 = TimeDistributed(
#     Lambda(lambda x_:
#                GraphConv(128, activation='relu', kernel_regularizer=l2(l2_reg), use_bias=True)([x_[..., :-N], x_[..., -N:]])
#                )
# )(conc1)
# gc1 = Dropout(dropout_rate)(gc1)
# conc2 = Concatenate()([gc1, filter_in])
# gc2 = TimeDistributed(
#     Lambda(lambda x_:
#                GraphConv(128, activation='relu', kernel_regularizer=l2(l2_reg), use_bias=True)([x_[..., :-N], x_[..., -N:]])
#                )
# )(conc2)

# # pool = TimeDistributed(NodeAttentionPool())(gc2)
# pool = TimeDistributed(GlobalAttentionPool(128))(gc2)
#     # pool = Lambda(lambda x_: K.reshape(x_, (-1, ts, N * 128)))(gc2)

# # Recurrent block
# lstm = LSTM(256, return_sequences=True)(pool)
# lstm = LSTM(256)(lstm)

# # Dense block
# # dense1 = BatchNormalization()(lstm)
# # dense1 = Dropout(dropout_rate)(dense1)
# dense1 = Dense(256, activation='relu', kernel_regularizer=l2(l2_reg))(lstm)
# # dense2 = BatchNormalization()(dense1)
# # dense2 = Dropout(dropout_rate)(dense2)
# dense2 = Dense(512, activation='relu', kernel_regularizer=l2(l2_reg))(dense1)

# adj_out = Dense(N * N, activation='sigmoid')(dense2)
# adj_out = Reshape((N, N), name='ADJ')(adj_out)
# nf_out = Dense(N * F, activation='linear')(dense2)
# nf_out = Reshape((N, F), name='NF')(nf_out)

# # Build model
# model = Model(inputs=[X_in, filter_in],
#                   outputs=[nf_out, adj_out])
# model.compile('adam',
#                   [mape],
#                   metrics=[mape])

# x, a, e = data_tr_input[0]
# x_label = data_tr_label[0]
# a_label = a

# x_v, a_v, e_v = data_va_input[0]
# va_x_label = data_va_label[0]
# va_a_label = a_v

# # Train model
# validation_data = [[x_v, a_v], [va_x_label, va_a_label]]
# model.fit([x, a],
#               [x_label, a_label],
#               epochs=epochs,
#               batch_size=batch_size,
#               validation_data=validation_data)


#     epoch = step = 0
#     best_val_loss = np.inf
#     best_weights = None
#     patience = es_patience
#     results = []
#     losses_history = []
#     val_losses_history = []
#     #start = time.perf_counter()


#     loss, acc = train_step(data_tr_input, data_tr_label)
#     results.append((loss, acc))

#     # Compute validation loss and accuracy
#     val_loss, val_acc = evaluate(data_va_input, data_va_label)
#     print(
#                 "Ep. {} - Loss: {:.3f} - Acc: {:.3f} - Val loss: {:.3f} - Val acc: {:.3f}".format(
#                     epoch, *np.mean(results, 0), val_loss, val_acc
#                 )
#             )

#     # Check if loss improved for early stopping
#     if val_loss < best_val_loss:
#                 best_val_loss = val_loss
#                 patience = es_patience
#                 print("New best val_loss {:.3f}".format(val_loss))
#                 best_weights = model.get_weights()

#     losses_history.append(loss)
#     val_losses_history.append(val_loss)
#     #elapsed = time.perf_counter() - start
#     ################################################################################
#     # Evaluate model
#     ################################################################################
#     model.set_weights(best_weights)  # Load best model
#     test_loss, test_acc = evaluate(data_te_input, data_te_label)
#     test_MAPE = test_evaluation(data_te_input, data_te_label)
#     print(test_MAPE)

#     print(
#         "Done. Test loss: {:.4f}. Test acc: {:.2f}".format(
#             test_loss, test_acc))
# #     test_scores.append(test_loss)
# #     train_losses.append(np.asarray(losses_history).min())
# #     val_losses.append(np.asarray(val_losses_history).min())

# # elapsed = time.perf_counter() - start

# # print("Best train losses: {} ".format(train_losses))
# # print("Best validation losses: {}".format(val_losses))
# # print("Test values: {}".format(test_scores))
# # print("--------------------------------------------")
# # print("K-Fold Train Score:{}".format(np.mean(train_losses)))
# # print("K-Fold Validation Score:{}".format(np.mean(val_losses)))
# # print("K-Fold Test Score: {}".format(np.mean(test_scores)))

# # print("Time for training: {:.4f} seconds".format(elapsed))
# # print("Time for training: {:.4f} minutes".format(elapsed/60))
