
import os
import numpy as np
import pandas as pd

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

#from memilio.simulation.secir import InfectionState

layer_name='ARMAConv'
number_of_layers = 3
number_of_channels = 128
layer = ARMAConv


file = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/data_GNN_400pop_one_var_damp_100days_1k_withmatrix/data_secir_age_groups.pickle','rb')
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

n_days = int(new_labels.shape[2]/48)

######## open commuter data #########
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(path,
                         'data')
commuter_file = open(os.path.join(
    path_data, 'commuter_migration_scaled.txt'), 'rb')
commuter_data = pd.read_csv(commuter_file, sep=" ", header=None)
sub_matrix = commuter_data.iloc[:numer_of_nodes, 0:numer_of_nodes]


adjacency_matrix = np.asarray(sub_matrix)

#adjacency_matrix[adjacency_matrix > 0] = 1
node_features = new_inputs

node_labels = new_labels

class MyDataset(spektral.data.dataset.Dataset):
        def read(self):              

            self.a = rescale_laplacian(
                normalized_laplacian(adjacency_matrix))

            return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features, node_labels)]

            super().__init__(**kwargs)

    # data = MyDataset()
data = MyDataset(transforms=NormalizeAdj())
batch_size = 32


class Net(Model):
            def __init__(self):
                super().__init__()

                self.conv1 = layer(number_of_channels, activation='elu')
                self.conv2 = layer(number_of_channels, activation='elu')
                self.conv3 = layer(number_of_channels, activation='elu')
                self.dense = Dense(data.n_labels, activation="linear")

            def call(self, inputs):
                x, a = inputs
                a = np.asarray(a)

                x = self.conv1([x, a])
                x = self.conv2([x, a])
                x = self.conv2([x, a])

                output = self.dense(x)

                return output


idxs = np.arange(len(data))#same for each run 
split_va, split_te = int(0.8 * len(data)), int(0.9 * len(data))
idx_tr, idx_va, idx_te = np.split(idxs, [split_va, split_te])

data_tr = data[idx_tr]
data_va = data[idx_va]
data_te = data[idx_te]
epochs = 2


loader_tr = MixedLoader(
            data_tr, batch_size=batch_size, epochs=epochs)
loader_va = MixedLoader(data_va,  batch_size=batch_size)
loader_te = MixedLoader(data_te, batch_size=data_te.n_graphs) 

inputs, target = loader_tr.__next__()

model = Net()
model(inputs, training=False)
model.load_weights('/home/schm_a45/Documents/code3/memilio/pycode/machine-learning/ARMAConv_1damp_saved_model_test')