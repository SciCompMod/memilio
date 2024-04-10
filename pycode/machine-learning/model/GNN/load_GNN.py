
import os
import pandas as pd
import json 
import pickle
import spektral
import time
import numpy as np
import scipy.sparse as sp
import tensorflow as tf
import matplotlib.pyplot as plt

from keras.layers import Dense
from keras.metrics import mean_absolute_percentage_error
from keras.models import Model

from spektral.data import  MixedLoader
from spektral.layers import ARMAConv
from spektral.transforms.normalize_adj import NormalizeAdj
from spektral.utils.convolution import gcn_filter, normalized_laplacian, rescale_laplacian, normalized_adjacency

from memilio.surrogatemodel.ode_secir_groups.data_generation_nodamp import get_population
from memilio.simulation.secir import InfectionState


#from memilio.simulation.secir import InfectionState

layer_name='ARMAConv'
number_of_layers = 3
number_of_channels = 128
layer = ARMAConv


######## open commuter data #########
numer_of_nodes = 400 
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(path,
                         'data')
commuter_file = open(os.path.join(
    path_data, 'commuter_migration_scaled.txt'), 'rb')
commuter_data = pd.read_csv(commuter_file, sep=" ", header=None)
sub_matrix = commuter_data.iloc[:numer_of_nodes, 0:numer_of_nodes]


adjacency_matrix = np.asarray(sub_matrix)

#adjacency_matrix[adjacency_matrix > 0] = 1



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




## load data and model with damping 

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

damping_factors = data_secir['damping_coeff'] # one for every run
damping_days = data_secir['damping_day'] # one for every rum 
contact_matrices = data_secir['damped_matrix']

n_runs = new_inputs.shape[0]
n_pop = new_inputs.shape[1]
#n_dampings = np.asarray(data_secir['damping_day']).shape[2]
n_dampings = 1


inputs_with_damp = np.dstack((new_inputs,(np.asarray(damping_factors).reshape([n_runs,n_pop,n_dampings])),
                               (np.asarray(damping_days).reshape([n_runs,n_pop,n_dampings])),
                               (np.asarray(contact_matrices).reshape([n_runs,n_pop,36*n_dampings]))))


node_features = inputs_with_damp

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
#model.load_weights('/home/schm_a45/Documents/code3/memilio/pycode/machine-learning/ARMAConv_1damp_saved_model_test')

with open("/home/schm_a45/Documents/Code/memilio/memilio/best_weights_ARMAConv_test.pickle", "rb") as fp:
     b = pickle.load(fp)
best_weights = b
model.set_weights(best_weights)

inputs, target = loader_te.__next__()
   
pred = model(inputs, training=False)


############################ for ARMAConv model without dampings #########################################################

## generate and save input and labels file
file = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/data_GNN_nodamp_400pop_1k_60days_w/data_secir_age_groups.pickle', 'rb')
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

node_features = new_inputs
node_labels = new_labels

class MyDataset(spektral.data.dataset.Dataset):
        def read(self):              

            self.a = rescale_laplacian(
                normalized_laplacian(adjacency_matrix))

            return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features, node_labels)]

            super().__init__(**kwargs)

data = MyDataset(transforms=NormalizeAdj())

idxs = np.arange(len(data))#same for each run 
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

# saving input and labels 
with open("inputs_60days_w.pickle", 'wb') as f:
      pickle.dump(inputs, f) 
with open("labels_60days_w.pickle", 'wb') as f:
      pickle.dump(target, f) 

###################################################################################################

# load the saved input and label file 
file_i = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_inputs/inputs_60days_w.pickle', 'rb')
test_inputs  = pickle.load(file_i)

file_l = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_labels/labels_60days_w.pickle', 'rb')
test_labels  = pickle.load(file_l)

node_features = test_inputs[0]
class MyDataset(spektral.data.dataset.Dataset):
        def read(self):              

            self.a = rescale_laplacian(
                normalized_laplacian(adjacency_matrix))

            return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features, test_labels)]

            super().__init__(**kwargs)

data = MyDataset(transforms=NormalizeAdj())
batch_size = node_features.shape[0]

model = Net()
model(test_inputs, training=False)
#model.load_weights('/home/schm_a45/Documents/code3/memilio/pycode/machine-learning/ARMAConv_1damp_saved_model_test')

with open("/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_model_weights/best_weights_ARMAConv_60days_nodamp_w.pickle", "rb") as fp:
     b = pickle.load(fp)
best_weights = b
model.set_weights(best_weights)
   
pred = model(test_inputs, training=False)

pred_reversed = np.expm1(pred)
labels_reversed = np.expm1(test_labels)

mae =  np.mean((abs(labels_reversed - pred_reversed)), axis = 0).transpose()
mape = 100*np.mean(abs((test_labels - pred)/test_labels), axis = 0).transpose()

mae_400 = np.mean(mae, axis = 0)
mape_400 = np.mean(mape, axis = 0)

n_days = 60 

def reshape_for_plot(array):
    #reshaped = []
    #for sample in array:
        sample = array.transpose().reshape(400,n_days,48)
        per_node =  []
        for node in sample:
            compartment_sum = []
            for day in node:
                    compartment_sum.append(day.reshape(-1, 8).transpose().mean(axis = 1))
                    
            per_node.append(np.asarray(compartment_sum).transpose().mean(axis = 1))
        #reshaped.append(per_node)
        return per_node

mape_reshaped = reshape_for_plot(mape)
mae_reshaped = reshape_for_plot(mae)


# pred_r = pred_reversed[0]
# p_48 = pred_r.reshape(400,60,48)[0][0]
# p_48.reshape(-1, 8).transpose().sum(axis = 1)

################################

path_population = os.path.abspath(
                r"data//pydata//Germany//county_population.json")
populations = get_population(path_population )
populations = np.asarray(populations).sum(axis = 1)

df = pd.DataFrame({'population':populations, 'mape':mape_400})

# from scipy.stats import pearsonr
# df_pearson = pd.DataFrame(pearsonr(df['population'], df['mape']), 
#                             index = ['pearson-coeff', 'p-value'], 
#                             columns = ['resultat_test'])


node1 = 324     # Berlin --> biggest city
node2 = 74       # Mettmann --> medium county 
node3 = 334   # Havelland --> small county (

def compartments_two_nodescomparison(node_id1, node_id2, node_id3):

    plt.figure().clf()

    compartment_array = []
    for compartment in InfectionState.values():
                compartment_array.append(compartment) 
    index = [3,5]
    compartments_cleaned= np.delete(compartment_array, index)

    compartmentnames=[str(compartment).split('.')[1]
               for compartment in compartments_cleaned]

    compartment_errors = [mape_reshaped[node_id1], mape_reshaped[node_id2], mape_reshaped[node_id3],
                             mae_reshaped[node_id1], mae_reshaped[node_id2],  mae_reshaped[node_id3]] 


    width = 0.25  # the width of the bars
    multiplier = 0
    
    fig, ((ax1, ax2,ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, sharey=True, figsize =(6,8), layout='constrained')
    axs = [ax1, ax2, ax3, ax4, ax5, ax6]
    xlabels=['MAPE', 'MAPE', 'MAPE', 'MAE', 'MAE', 'MAE']
    titles = ['Berlin MAPE', 'Mettmann MAPE', 'Havelland MAPE', 'Berlin MAE', 'Mettmann MAE', 'Havelland MAE']
    p_values = [2, 2, 2, 1000, 1000, 1000]

    for ax, xlabel, title , pvalue, measurement in zip(axs, xlabels, titles, p_values, compartment_errors):
                x = np.arange(len(compartmentnames))  
                rects = ax.barh(x, measurement.round(4), width, label=title)

                large_bars = [p if p > pvalue else '' for p in measurement.round(4)] # for MAPE
                small_bars = [p if p <= pvalue else '' for p in measurement.round(4)] # for MAPE
                
                large_bars = [p if p > pvalue else '' for p in measurement.round(2)] # for MAE
                small_bars = [p if p <= pvalue else '' for p in measurement.round(2)] # for MAE

                ax.bar_label(rects, small_bars,
                  padding=3, color='black', fontsize = 8 )
                ax.bar_label(rects, large_bars,
                                        padding=-40, color='white', fontsize = 8)

                multiplier += 1

                # Add some text for labels, title and custom x-axis tick labels, etc.
                ax.set_xlabel(xlabel)
                ax.set_title(title)                  
                ax.set_yticks(x + width * 0.2) 
                ax.set_yticklabels(compartmentnames)

    plt.savefig("GNN_compartmentdistribution.png")