import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import pickle
import spektral
import scipy
import time
from functools import partial

import numpy as np
import scipy.sparse as sp
import tensorflow as tf
from keras.layers import Dense
from keras.losses import MeanAbsolutePercentageError
from keras.metrics import mean_absolute_percentage_error
from keras.models import Model
from keras.optimizers import Adam

from sklearn.model_selection import KFold

from spektral.data import Dataset, DisjointLoader, Graph, Loader, BatchLoader, MixedLoader
from spektral.layers import GATConv, GINConv, XENetConv, XENetConvBatch, ECCConv, CrystalConv, CensNetConv, SRCPool, MinCutPool

from spektral.transforms.normalize_adj import NormalizeAdj
from spektral.utils.convolution import gcn_filter, line_graph, incidence_matrix
from spektral.utils.sparse import sp_matrix_to_sp_tensor
from utils.gnn_utils import incidence_matrix

from data_generator_withoutdampings_all401 import get_population
from memilio.simulation.secir import InfectionState

from sklearn.preprocessing import FunctionTransformer

# add to setup: torch, spektral


def getBaselineMatrix():
    """! loads the baselinematrix
    """

    baseline_contact_matrix0 = os.path.join(
        "./data/contacts/baseline_home.txt")
    baseline_contact_matrix1 = os.path.join(
        "./data/contacts/baseline_school_pf_eig.txt")
    baseline_contact_matrix2 = os.path.join(
        "./data/contacts/baseline_work.txt")
    baseline_contact_matrix3 = os.path.join(
        "./data/contacts/baseline_other.txt")

    baseline = np.loadtxt(baseline_contact_matrix0) \
        + np.loadtxt(baseline_contact_matrix1) + \
        np.loadtxt(baseline_contact_matrix2) + \
        np.loadtxt(baseline_contact_matrix3)

    return baseline


path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.dirname(os.path.realpath(path)))),
    'data/data_GNN_nodamp_20pop_1k')

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


adjacency_matrix = np.asarray(sub_matrix.copy())
adjacency_matrix[adjacency_matrix > 0] = 1

node_features = new_inputs

node_labels = new_labels

# # my calculation
# population = get_population()
# baseline = getBaselineMatrix()
# pop_sum_per_agegroup = pd.DataFrame(data=population).sum()

# contacts_per_agegroup = baseline.sum(axis = 1)

# mean_contact_pP= sum(pop_sum_per_agegroup*contacts_per_agegroup)/sum(pop_sum_per_agegroup)


# mean_contacts_per_pop=(np.asarray(population).mean(axis=1))*mean_contact_pP

# mean_contacts_per_pop=(np.asarray(population).mean(axis=1))*10

# # add mean contacts per pop to adjacency matrix
# adjacency_list = adjacency_matrix.tolist()
# for i in range(0, len(adjacency_matrix)):
#     adjacency_list[i][i] = mean_contacts_per_pop[i]


# mean_contacts_per_pop = []
# for pop in population:
#     sum_pop = np.asarray(pop).sum()*10
#     mean_contacts_per_pop.append(sum_pop)

# # add baseline contact matrix to adjacency matrix
# adjacency_list = adjacency_matrix.tolist()
# for i in range(0, len(adjacency_matrix)):
#     adjacency_list[i][i] = baseline.mean()

# add mean contacts per pop to adjacency matrix
# adjacency_list = adjacency_matrix.tolist()
# for i in range(0, len(adjacency_matrix)):
#     adjacency_list[i][i] = mean_contacts_per_pop[i]

# # add edge with weight 1 to adjacency matrix (??? better solution?)
# adjacency_matrix_withselfedges = adjacency_matrix
# for i in range(0, len(adjacency_matrix)):
#     adjacency_matrix_withselfedges[i][i] = 1

# create edge features from adjacency matrix
df = pd.DataFrame(data=np.asarray(sub_matrix))

#df = pd.DataFrame(data=adjacency_list)
edge_features = df.rename_axis('Source')\
    .reset_index()\
    .melt('Source', value_name='Weight', var_name='Target')\
    .reset_index(drop=True)

edge_features = edge_features[
    edge_features.Weight != 0]


transformer = FunctionTransformer(np.log1p, validate=True)
scaled_edges = transformer.transform(
    edge_features['Weight'].values.reshape(-1, 1))


# df = pd.DataFrame(data=adjacency_list)
# edge_features_withoutselfedges = df.rename_axis('Source')\
#     .reset_index()\
#     .melt('Source', value_name='Weight', var_name='Target')\
#     .query('Source != Target')\
#     .reset_index(drop=True)

# edge_features_withoutselfedges = edge_features_withoutselfedges[
#     edge_features_withoutselfedges.Weight != 0]


# df = pd.DataFrame(data=adjacency_list)
# edge_features = df.rename_axis('Source')\
#     .reset_index()\
#     .melt('Source', value_name='Weight', var_name='Target')\
#     .reset_index(drop=True)

# edge_features_correct = pd.concat(
#     [edge_features.loc[edge_features['Source'] == edge_features['Target']],
#      edge_features_withoutselfedges]).sort_index()


# class MyDataset(spektral.data.dataset.Dataset):
#     def read(self):
#         # self.a = sp.csr_matrix(adjacency_matrix)
#         self.a = adjacency_matrix
#         self.e = np.asarray(edge_features_correct['Weight']).reshape(
#             np.asarray(edge_features_correct['Weight']).shape[0], 1)
#         return [spektral.data.Graph(x=x, y=y, a=self.a, e=self.e, ) for x, y in zip(node_features, node_labels)]

#         super().__init__(**kwargs)


# for mixed mode


class MyDataset(spektral.data.dataset.Dataset):
    def read(self):
        self.a = adjacency_matrix
        #self.a = np.asarray(adjacency_list)
        # #self.e = np.asarray(edge_features['Weight']).reshape(
        # np.asarray(edge_features['Weight']).shape[0], 1)
        self.e = scaled_edges.reshape(np.asarray(
            edge_features['Weight']).shape[0], 1)

        return [spektral.data.Graph(x=x, y=y, e=self.e, a=self.a) for x, y in zip(node_features, node_labels)]

        super().__init__(**kwargs)


# class MyDataset(spektral.data.dataset.Dataset):
#     def read(self):
#         # self.a = sp.csr_matrix(adjacency_matrix)
#         self.a = adjacency_matrix
#         self.e = np.asarray(edge_features_withoutselfedges['Weight']).reshape(
#             np.asarray(edge_features_withoutselfedges['Weight']).shape[0], 1)
#         return [spektral.data.Graph(x=x, y=y, e=self.e) for x, y in zip(node_features, node_labels)]
#         super().__init__(**kwargs)
# data = MyDataset()
#data = MyDataset(transforms=NormalizeAdj())
data = MyDataset()
batch_size = 32
epochs = 500
es_patience = 50  # Patience for early stopping

# Train/valid/test split
# idxs = np.random.permutation(len(data))
# split_va, split_te = int(0.8 * len(data)), int(0.9 * len(data))
# idx_tr, idx_va, idx_te = np.split(idxs, [split_va, split_te])
# data_tr = data[idx_tr]
# data_va = data[idx_va]
# data_te = data[idx_te]

# # Data loaders
# loader_tr = DisjointLoader(data_tr, batch_size=batch_size, epochs=epochs)
# loader_va = DisjointLoader(data_va, batch_size=batch_size)
# loader_te = DisjointLoader(data_te, batch_size=batch_size)

# Data loaders
# loader_tr = BatchLoader(
#     data_tr, batch_size=batch_size, epochs=epochs)
# loader_va = BatchLoader(data_va,  batch_size=batch_size)
# loader_te = BatchLoader(data_te, batch_size=data_te.n_graphs)

# loader_tr = MixedLoader(
#     data_tr, batch_size=batch_size, epochs=epochs)
# loader_va = MixedLoader(data_va,  batch_size=batch_size)
# loader_te = MixedLoader(data_te, batch_size=data_te.n_graphs)


def preprocess(adjacency):
    laplacian = gcn_filter(adjacency)
    incidence = incidence_matrix(adjacency)
    edge_laplacian = gcn_filter(line_graph(incidence).numpy())

    return laplacian, edge_laplacian, incidence


#laplacian, edge_laplacian, incidence = preprocess(adjacency_matrix)
#matrix_tuple = [laplacian, edge_laplacian, incidence]


################################################################################
# Build model
################################################################################
class Net(Model):
    def __init__(self):
        super().__init__()

        #self.conv1 = ECCConv(24,   activation="relu")
        initializer = tf.keras.initializers.GlorotUniform(seed = 42)
        self.conv1 = XENetConvBatch(
            32, 240, 1,  activation="relu", kernel_initializer=initializer)

        #self.srcpool = SRCPool()
        #self.sagpool = MinCutPool(k = 50)

        # self.conv2 = XENetConv(
        #     32, 240, 1,  activation="relu")

        # self.conv1 = XENetConvBatch(
        #     32, 240, 1,    activation="relu")
        # self.conv2 = XENetConvBatch(
        #     32, 240, 1,  activation="relu")
        # self.conv3 = XENetConv(32, 64,64,  activation="relu")

        # self.global_pool = GlobalAvgPool()
        self.dense = Dense(data.n_labels, activation="linear")

    def call(self, inputs):
        # x, a, e = inputs
        x, a, e = inputs
        # e = np.asarray(edge_features['Weight']).reshape(
        #     np.asarray(edge_features['Weight']).shape[0], 1)

        # a = np.asarray(a)
        # e = np.asarray(e)
        #x,e= self.conv1([x, a, e])

        # e = np.asarray(edge_features['Weight']).reshape(
        #     np.asarray(edge_features['Weight']).shape[0], 1)

        x, e = self.conv1([x, a, e])
        #x_pool, a_pool = self.srcpool([x, a])
        #x_pool, a_pool = self.sagpool([x, a])


        #x,e= self.conv1([x, matrix_tuple,e])

        # x, e = self.conv2([x, a, e])
        # x = self.conv3([x, a])
        # x = self.global_pool([x])
        #output = self.dense(x_pool)
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


kf = KFold(n_splits=5)
train_idxs = []
test_idxs = []
for i, (train_index, test_index) in enumerate(kf.split(data)):

    train_idxs.append(train_index)
    test_idxs.append(test_index)

test_scores = []
train_losses = []
val_losses = []

start = time.perf_counter()


for train_idx, test_idx in zip(train_idxs, test_idxs):

    learning_rate = 0.001
    model = Net()
    optimizer = Adam(learning_rate=learning_rate)
    loss_fn = MeanAbsolutePercentageError()


    data_tr = data[train_idx[:(int(0.8*len(train_idx)))]]
    data_va = data[train_idx[(int(0.8*len(train_idx))):]]
    data_te = data[test_idx]

    # Data loaders
    loader_tr = BatchLoader(
        data_tr, batch_size=batch_size, epochs=epochs)
    loader_va = BatchLoader(data_va,  batch_size=batch_size)
    loader_te = BatchLoader(data_te, batch_size=data_te.n_graphs)

    epoch = step = 0
    best_val_loss = np.inf
    best_weights = None
    patience = es_patience
    results = []
    losses_history = []
    val_losses_history = []
    #start = time.perf_counter()
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
                    print(
                        "Early stopping (best val_loss: {})".format(
                            best_val_loss))
                    break
            results = []
            losses_history.append(loss)
            val_losses_history.append(val_loss)
    #elapsed = time.perf_counter() - start
    ################################################################################
    # Evaluate model
    ################################################################################
    model.set_weights(best_weights)  # Load best model
    test_loss, test_acc = evaluate(loader_te)
    test_MAPE = test_evaluation(loader_te)
    print(test_MAPE)

    print(
        "Done. Test loss: {:.4f}. Test acc: {:.2f}".format(
            test_loss, test_acc))
    test_scores.append(test_loss)
    train_losses.append(np.asarray(losses_history).min())
    val_losses.append(np.asarray(val_losses_history).min())


    ##### set random weight for next run 
    # model.set_weights 
     

elapsed = time.perf_counter() - start

print("Best train losses: {} ".format(train_losses))
print("Best validation losses: {}".format(val_losses))
print("Test values: {}".format(test_scores))
print("--------------------------------------------")
print("K-Fold Train Score:{}".format(np.mean(train_losses)))
print("K-Fold Validation Score:{}".format(np.mean(val_losses)))
print("K-Fold Test Score: {}".format(np.mean(test_scores)))

print("Time for training: {:.4f} seconds".format(elapsed))
print("Time for training: {:.4f} minutes".format(elapsed/60))
