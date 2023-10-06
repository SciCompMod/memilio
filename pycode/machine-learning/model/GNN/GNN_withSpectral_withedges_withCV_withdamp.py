import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
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
from keras.optimizers.legacy import Adam, Nadam, RMSprop, SGD, Adagrad

from sklearn.preprocessing import FunctionTransformer
from sklearn.model_selection import KFold

from spektral.data import Dataset, DisjointLoader, Graph, Loader, BatchLoader, MixedLoader
from spektral.layers import GCSConv, GlobalAvgPool, GlobalAttentionPool, ARMAConv, AGNNConv, APPNPConv, CrystalConv, GATConv, GINConv, XENetConv, GCNConv
from spektral.transforms.normalize_adj import NormalizeAdj
from spektral.layers import GATConv, GINConv, XENetConv, XENetConvBatch
from spektral.utils.convolution import gcn_filter, normalized_laplacian, rescale_laplacian, normalized_adjacency

from memilio.simulation.secir import InfectionState





path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.dirname(os.path.realpath(path)))),
    'data_GNN_400pop_one_var_damp_30days_1k_withmatrix')

file = open(os.path.join(path_data, 'data_secir_age_groups.pickle'), 'rb')
data_secir = pickle.load(file)


len_dataset = data_secir['inputs'][0].shape[0]
numer_of_nodes = np.asarray(data_secir['inputs']).shape[0]
damping_days = data_secir['damping_day'][0]
damping_coeffs = data_secir['damping_coeff'][0]
matrices = data_secir['damped_matrix'][0]
num_runs = np.asarray(data_secir['inputs']).shape[1]
number_of_dampings = 1 ##change to flexible 

# shape_input_flat = np.asarray(
#     data_secir['inputs']).shape[2]*np.asarray(data_secir['inputs']).shape[3]+ number_of_dampings+ number_of_dampings + number_of_dampings*36
shape_input_flat = np.asarray(
     data_secir['inputs']).shape[2]*np.asarray(data_secir['inputs']).shape[3]
shape_labels_flat = np.asarray(
    data_secir['labels']).shape[2]*np.asarray(data_secir['labels']).shape[3]


new_inputs = np.asarray(
    data_secir['inputs']).reshape(
    len_dataset, numer_of_nodes, shape_input_flat)
new_labels = np.asarray(data_secir['labels']).reshape(
    len_dataset, numer_of_nodes, shape_labels_flat)



## adding damping informations to inputs 
new_inputs_with_damp=[]
for run,damping_day, damping_coeff, matrix in zip(new_inputs, damping_days, damping_coeffs,matrices): 
    for node in run  :
        new_inputs_with_damp.append(np.append(node,[np.append([damping_day, damping_coeff], matrix.flatten())]))
new_inputs_with_damp = np.asarray(new_inputs_with_damp).reshape(num_runs, numer_of_nodes,np.asarray(new_inputs_with_damp).shape[1])

######## open commuter data #########
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(path,
                         'data')
commuter_file = open(os.path.join(
    path_data, 'commuter_migration_scaled.txt'), 'rb')
commuter_data = pd.read_csv(commuter_file, sep=" ", header=None)
sub_matrix = commuter_data.iloc[:numer_of_nodes, 0:numer_of_nodes]


adjacency_matrix = np.asarray(sub_matrix)

# adjacency_matrix[adjacency_matrix > 0] = 1
node_features = new_inputs_with_damp

node_labels = new_labels



# create edge features from adjacency matrix
df = pd.DataFrame(data=np.asarray(sub_matrix))

# df = pd.DataFrame(data=adjacency_list)
edge_features = df.rename_axis('Source')\
    .reset_index()\
    .melt('Source', value_name='Weight', var_name='Target')\
    .reset_index(drop=True)


edge_features = edge_features[
    edge_features.Weight != 0]  # 157864 edges



transformer = FunctionTransformer(np.log1p, validate=True)
scaled_edges = transformer.transform(
    edge_features['Weight'].values.reshape(-1, 1))
edge_features['Weight']=scaled_edges

class MyDataset(spektral.data.dataset.Dataset):
    def read(self):
        self.a = adjacency_matrix
        # self.a = normalized_adjacency(adjacency_matrix)
        # self.a = rescale_laplacian(normalized_laplacian(adjacency_matrix))
        self.e = scaled_edges.reshape(np.asarray(
            edge_features['Weight']).shape[0], 1)
        return [spektral.data.Graph(x=x, y=y, e = self.e, a = self.a) for x, y in zip(node_features, node_labels)]

        super().__init__(**kwargs)


# data = MyDataset()
data = MyDataset(transforms=NormalizeAdj())
batch_size = 32
epochs = 2
es_patience = 100  # Patience for early stopping

shape_node_features = node_features.shape[2]
shape_node_labels = node_labels.shape[2]
################################################################################
# Build model
################################################################################
class Net(Model):
    def __init__(self):
        super().__init__()

        # self.conv1 = ECCConv(24,   activation="relu")
        # initializer = tf.keras.initializers.GlorotUniform(seed=42)
        self.conv1 = XENetConvBatch(
            32, shape_node_features, 1,  activation="relu")

        # self.srcpool = SRCPool()
        # self.sagpool = MinCutPool(k = 50)

        # self.conv2 = XENetConv(
        #    32, 240, 1,  activation="relu")

        # self.conv1 = XENetConvBatch(
        #     32, 240, 1,    activation="relu")
        # self.conv2 = XENetConvBatch(
        #     32, 240, 1,  activation="relu")
        # self.conv3 = XENetConv(32, 64,64,  activation="relu")

        # self.global_pool = GlobalAvgPool()
        self.dense = Dense(shape_node_labels, activation="linear")

    def call(self, inputs):
        # x, a, e = inputs
        x, a, e = inputs
        # e = np.asarray(edge_features['Weight']).reshape(
        #     np.asarray(edge_features['Weight']).shape[0], 1)

        # a = np.asarray(a)
        # e = np.asarray(e)
        # x,e= self.conv1([x, a, e])

        # e = np.asarray(edge_features['Weight']).reshape(
        #     np.asarray(edge_features['Weight']).shape[0], 1)

        x, e = self.conv1([x, a, e])
        # x_pool, a_pool = self.srcpool([x, a])
        # x_pool, a_pool = self.sagpool([x, a])

        # x,e= self.conv1([x, matrix_tuple,e])

        # x, e = self.conv2([x, a, e])
        # x = self.conv3([x, a])
        # x = self.global_pool([x])
        # output = self.dense(x_pool)
        output = self.dense(x)

        return output
    



learning_rate = 0.001

optimizer = Adam(learning_rate=learning_rate)
#optimizer = Adagrad(learning_rate=learning_rate)
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
        pred = model.predict(inputs)
        outs = (
            loss_fn(target, pred),
            tf.reduce_mean(mean_absolute_percentage_error(target, pred)),
            len(target),  # Keep track of batch size
        )
        output.append(outs)
        if step == loader.steps_per_epoch:
            output = np.array(output)
            return np.average(output[:, :-1], 0, weights=output[:, -1])

n_days = int(new_labels.shape[2]/48)

def test_evaluation(loader):

    inputs, target = loader.__next__()
    pred = model.predict(inputs)

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
            
                # delete the two confirmed compartments from InfectionStates
            compartment_array = []
            for compartment in InfectionState.values():
                compartment_array.append(compartment) 
            index = [3,5]
            compartments_cleaned= np.delete(compartment_array, index)

            MAPE_v.append(relative_err_means_percentage)

        mean_per_batch.append(np.asarray(MAPE_v).transpose().mean(axis=1))

    mean_percentage = pd.DataFrame(
        data=np.asarray(mean_per_batch).transpose().mean(axis=1),
        index=[str(compartment).split('.')[1]
               for compartment in compartments_cleaned],
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

losses_history_all = []
val_losses_history_all = []

start = time.perf_counter()
for train_idx, test_idx in zip(train_idxs, test_idxs):
    model = Net()

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
    losses_history_all.append(losses_history)
    val_losses_history_all.append(val_losses_history)

elapsed = time.perf_counter() - start


# #plot the losses
# plt.figure()
# plt.plot(np.asarray(losses_history_all).mean(axis=0), label='train loss')
# plt.plot(np.asarray(val_losses_history_all).mean(axis=0), label='val loss')
# plt.xlabel('Epoch')
# plt.ylabel('Loss ( MAPE)')
# plt.title('Loss for GCN Adagrad ')
# plt.legend()
# plt.savefig('losses_GCN_Adagrad.png')


# print out stats
print("Best train losses: {} ".format(train_losses))
print("Best validation losses: {}".format(val_losses))
print("Test values: {}".format(test_scores))
print("--------------------------------------------")
print("K-Fold Train Score:{}".format(np.mean(train_losses)))
print("K-Fold Validation Score:{}".format(np.mean(val_losses)))
print("K-Fold Test Score: {}".format(np.mean(test_scores)))

print("Time for training: {:.4f} seconds".format(elapsed))
print("Time for training: {:.4f} minutes".format(elapsed/60))
