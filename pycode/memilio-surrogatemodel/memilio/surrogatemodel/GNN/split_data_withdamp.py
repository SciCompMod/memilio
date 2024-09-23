
import pickle
import os
import numpy as np
import pandas as pd
import spektral
from spektral.utils.convolution import normalized_laplacian, rescale_laplacian
from spektral.transforms.normalize_adj import NormalizeAdj
from spektral.data import MixedLoader
from spektral.data import MixedLoader
from spektral.layers import ARMAConv
from spektral.transforms.normalize_adj import NormalizeAdj
from spektral.utils.convolution import normalized_laplacian, rescale_laplacian
from keras.layers import Dense
from keras.losses import MeanAbsolutePercentageError
from keras.metrics import mean_absolute_percentage_error
from keras.models import Model

# load data and model with damping

# file = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/data_GNN_400pop_4var_damp_100days_1k_w/GNN_400pop_damp_w.pickle','rb')
file = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/data_GNN_with_5_dampings/data_secir_age_groups.pickle', 'rb')
data_secir = pickle.load(file)
len_dataset = 1000

number_of_nodes = 400

shape_input_flat = np.asarray(
    data_secir['inputs']).shape[1]*np.asarray(data_secir['inputs']).shape[2]
shape_labels_flat = np.asarray(
    data_secir['labels']).shape[1]*np.asarray(data_secir['labels']).shape[2]

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

damping_factors = data_secir['damping_coeff']  # one for every run
damping_days = data_secir['damping_day']  # one for every rum
contact_matrices = data_secir['damped_matrix']

n_runs = new_inputs.shape[0]
n_pop = new_inputs.shape[1]
n_dampings = np.asarray(data_secir['damping_day']).shape[1]
# n_dampings = 1

# inputs_with_damp = np.dstack((new_inputs,(np.asarray(damping_factors).reshape([n_runs,n_pop,n_dampings])),
#                                (np.asarray(damping_days).reshape([n_runs,n_pop,n_dampings])),
#                                (np.asarray(contact_matrices).reshape([n_runs,n_pop,36*n_dampings]))))

inputs_with_damp = np.dstack((new_inputs, np.repeat(np.asarray(damping_factors)[:, np.newaxis, :], 400, axis=1),
                              (np.repeat(np.asarray(damping_days)
                               [:, np.newaxis, :], 400, axis=1)),
                              (np.repeat(np.asarray(contact_matrices)[:, np.newaxis, :], 400, axis=1).reshape([n_runs, n_pop, 36*n_dampings]))))

# for one damp:
# damping_factors_reshaped =  np.tile(np.asarray(damping_factors)[:, np.newaxis], (1, 400))
# damping_days_reshaped = np.tile(np.asarray(damping_days)[:, np.newaxis], (1, 400))
# contact_matrices_reshaped = np.repeat(np.asarray(contact_matrices)[:, np.newaxis, :, :], 400, axis=1).reshape([n_runs,n_pop,-1])
# inputs_with_damp = np.concatenate((new_inputs,
#                                   np.expand_dims(damping_factors_reshaped,axis = 2),
#                                   np.expand_dims(damping_days_reshaped, axis = 2),
#                                   contact_matrices_reshaped), axis = 2  )


node_features = inputs_with_damp

node_labels = new_labels
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


class MyDataset(spektral.data.dataset.Dataset):
    def read(self):

        self.a = rescale_laplacian(
            normalized_laplacian(adjacency_matrix))

        return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features, node_labels)]

        super().__init__(**kwargs)


# class Net(Model):
#     def __init__(self):
#         super().__init__()
#         self.conv1 = ARMAConv(1024, order=2, iterations=2,
#                               share_weights=True,  activation='elu')
#         self.dense = Dense(data.n_labels, activation="linear")

#     def call(self, inputs):
#         x, a = inputs
#         a = np.asarray(a)
#         x = self.conv1([x, a])
#         output = self.dense(x)

#         return output


# data = MyDataset()
data = MyDataset(transforms=NormalizeAdj())
batch_size = 32


idxs = np.arange(len(data))  # same for each run
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

#model = Net()

# saving input and labels
with open("/hpc_data/schm_a45/GNN_data/paper/inputs_100days_5damp_paper.pickle", 'wb') as f:
    pickle.dump(inputs, f)
with open("/hpc_data/schm_a45/GNN_data/paper/labels_100days_5damp_paper.pickle", 'wb') as f:
    pickle.dump(target, f)

