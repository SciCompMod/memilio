import pickle
import os
import numpy as np
import pandas as pd
import spektral
from spektral.utils.convolution import normalized_laplacian, rescale_laplacian
from spektral.transforms.normalize_adj import NormalizeAdj
from spektral.data import MixedLoader
from spektral.layers import ARMAConv
from keras.layers import Dense
from keras.models import Model


file = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/data_GNN_nodamp_90days/data_secir_age_groups.pickle', 'rb')
data_secir = pickle.load(file)

len_dataset = data_secir['inputs'].shape[0]
number_of_nodes = data_secir['inputs'].shape[1]
######## open commuter data #########
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(path,
                         'data')
commuter_file = open(os.path.join(
    path_data, 'commuter_migration_scaled.txt'), 'rb')
commuter_data = pd.read_csv(commuter_file, sep=" ", header=None)
sub_matrix = commuter_data.iloc[:number_of_nodes, 0:number_of_nodes]

adjacency_matrix = np.asarray(sub_matrix)
adjacency_matrix[adjacency_matrix > 0] = 1


class Net(Model):
    def __init__(self):
        super().__init__()
        self.conv1 = ARMAConv(1024, order=2, iterations=2,
                              share_weights=True,  activation='elu')
        self.dense = Dense(data.n_labels, activation="linear")

    def call(self, inputs):
        x, a = inputs
        a = np.asarray(a)
        x = self.conv1([x, a])
        output = self.dense(x)

        return output


shape_input_flat = np.asarray(
    data_secir['inputs']).shape[2]*np.asarray(data_secir['inputs']).shape[3]
shape_labels_flat = np.asarray(
    data_secir['labels']).shape[2]*np.asarray(data_secir['labels']).shape[3]

n_input_days = data_secir['inputs'].shape[2]
n_label_days = data_secir['labels'].shape[2]

n_age_groups = 6
n_compartments = 8

new_inputs = np.asarray(
    data_secir['inputs']).reshape(
    len_dataset, number_of_nodes, n_input_days*n_age_groups*n_compartments)
new_labels = np.asarray(data_secir['labels']).reshape(
    len_dataset, number_of_nodes, n_label_days*n_age_groups*n_compartments)
n_days = int(new_labels.shape[2]/48)

node_features = new_inputs
node_labels = new_labels


class MyDataset(spektral.data.dataset.Dataset):
    def read(self):

        self.a = rescale_laplacian(
            normalized_laplacian(adjacency_matrix))

        return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features, node_labels)]

        super().__init__(**kwargs)


data = MyDataset(transforms=NormalizeAdj())

idxs = np.arange(len(data))  # same for each run
split_va, split_te = int(0.8 * len(data)), int(0.9 * len(data))
idx_tr, idx_va, idx_te = np.split(idxs, [split_va, split_te])

data_tr = data[idx_tr]
data_va = data[idx_va]
data_te = data[idx_te]
epochs = 2
batch_size = data_te.n_graphs


loader_tr = MixedLoader(
    data_tr, batch_size=batch_size, epochs=epochs)
loader_va = MixedLoader(data_va,  batch_size=batch_size)
loader_te = MixedLoader(data_te, batch_size=data_te.n_graphs)

inputs, target = loader_tr.__next__()

model = Net()

inputs, target = loader_te.__next__()

# saving input and labels
with open("inputs_90days_paper_new.pickle", 'wb') as f:
    pickle.dump(inputs, f)
with open("labels_90days_paper_new.pickle", 'wb') as f:
    pickle.dump(target, f)
