
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
from sklearn.preprocessing import FunctionTransformer

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
number_of_layers = 4
number_of_channels = 512
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

adjacency_matrix[adjacency_matrix > 0] = 1



class Net(Model):
            def __init__(self):
                super().__init__()

                self.conv1 = layer(number_of_channels, order = 3, activation='relu')
                self.conv2 = layer(number_of_channels, order = 3, activation='relu')
                self.conv3 = layer(number_of_channels, order = 3, activation='relu')
                self.conv4 = layer(number_of_channels, order = 3, activation='relu')
                self.dense = Dense(data.n_labels, activation="linear")

            def call(self, inputs):
                x, a = inputs
                a = np.asarray(a)

                x = self.conv1([x, a])
                x = self.conv2([x, a])
                x = self.conv3([x, a])
                x = self.conv4([x,a])

                output = self.dense(x)

                return output




################### load data and model with damping 

#file = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/data_GNN_400pop_4var_damp_100days_1k_w/GNN_400pop_damp_w.pickle','rb')
file = open('/localdata1/schm_a45/GNN_data/dampings/data_GNN_100days_400pop_4damp_2024/GNN_400pop_damp_w.pickle', 'rb')
data_secir = pickle.load(file)


len_dataset = 1000

numer_of_nodes = 400
shape_input_flat = np.asarray(
    data_secir['inputs'][1]).shape[1]*np.asarray(data_secir['inputs'][1]).shape[2]
shape_labels_flat = np.asarray(
    data_secir['labels'][1]).shape[1]*np.asarray(data_secir['labels'][1]).shape[2]


### forgot to scale data in data generation 
transformer = FunctionTransformer(np.log1p, validate=True)
# inputs = np.asarray(
#     data['inputs']).transpose(
            #     2, 0, 1)
inputs = np.asarray(
data_secir['inputs'][1]).transpose(2,0,1,3).reshape(48,-1)
scaled_inputs = transformer.transform(inputs)
scaled_inputs = scaled_inputs.transpose().reshape(len_dataset, 400,5, 48)

labels = np.asarray(
            data_secir['labels'][1]).transpose(2,0,1,3).reshape(48,-1)
scaled_labels = transformer.transform(labels)
scaled_labels = scaled_labels.transpose().reshape(len_dataset, 400,100, 48)

new_inputs = np.asarray(
    scaled_inputs).reshape(
    len_dataset, numer_of_nodes, shape_input_flat)
new_labels = np.asarray(scaled_labels).reshape(
    len_dataset, numer_of_nodes, 100*48)


n_days = int(new_labels.shape[2]/48)

damping_factors = data_secir['damping_coeff'][0] # one for every run
damping_days = data_secir['damping_day'][0]# one for every rum 
contact_matrices = data_secir['damped_matrix'][0]

n_runs = new_inputs.shape[0]
n_pop = new_inputs.shape[1]
n_dampings = np.asarray(data_secir['damping_day']).shape[2]
#n_dampings = 1

# inputs_with_damp = np.dstack((new_inputs,(np.asarray(damping_factors).reshape([n_runs,n_pop,n_dampings])),
#                                (np.asarray(damping_days).reshape([n_runs,n_pop,n_dampings])),
#                                (np.asarray(contact_matrices).reshape([n_runs,n_pop,36*n_dampings]))))

inputs_with_damp = np.dstack((new_inputs,np.repeat(np.asarray(damping_factors)[:, np.newaxis, :], 400, axis=1),
                               (np.repeat(np.asarray(damping_days)[:, np.newaxis, :], 400, axis=1)),
                              (np.repeat(np.asarray(contact_matrices)[:, np.newaxis, :], 400, axis=1).reshape([n_runs,n_pop,36*n_dampings]))))

### for one damp: 
# damping_factors_reshaped =  np.tile(np.asarray(damping_factors)[:, np.newaxis], (1, 400))
# damping_days_reshaped = np.tile(np.asarray(damping_days)[:, np.newaxis], (1, 400))
# contact_matrices_reshaped = np.repeat(np.asarray(contact_matrices)[:, np.newaxis, :, :], 400, axis=1).reshape([n_runs,n_pop,-1])
# inputs_with_damp = np.concatenate((new_inputs,
#                                   np.expand_dims(damping_factors_reshaped,axis = 2),
#                                   np.expand_dims(damping_days_reshaped, axis = 2),
#                                   contact_matrices_reshaped), axis = 2  )


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

inputs, target = loader_te.__next__() 



# saving input and labels 
with open("inputs_100days_4damp_graphdata.pickle", 'wb') as f:
      pickle.dump(inputs, f) 
with open("labels_100days_4damp_graphdata.pickle", 'wb') as f:
      pickle.dump(target, f) 


#model.load_weights('/home/schm_a45/Documents/code3/memilio/pycode/machine-learning/ARMAConv_1damp_saved_model_test')

# with open("/home/schm_a45/Documents/Code/memilio/memilio/best_weights_ARMAConv_test.pickle", "rb") as fp:
#      b = pickle.load(fp)
# best_weights = b
# model.set_weights(best_weights)

#pred = model(inputs, training=False)

############################ for ARMAConv model without dampings #########################################################

## generate and save input and labels file
#file = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/data_GNN_nodamp_400pop_1k_120days_w/data_secir_age_groups.pickle', 'rb')
file = open('/localdata1/schm_a45/GNN_data/days/data_GNN_nodamp_400pop_60days_2024/data_secir_age_groups.pickle', 'rb')
data_secir = pickle.load(file)

len_dataset = data_secir['inputs'][0].shape[0]
#numer_of_nodes = np.asarray(data_secir['inputs']).shape[0]
numer_of_nodes = 400
shape_input_flat = np.asarray(
    data_secir['inputs'][0]).shape[2]*np.asarray(data_secir['inputs'][0]).shape[3]
shape_labels_flat = np.asarray(
    data_secir['labels'][0]).shape[2]*np.asarray(data_secir['labels'][0]).shape[3]

new_inputs = np.asarray(
    data_secir['inputs'][0]).reshape(
    len_dataset, numer_of_nodes, shape_input_flat)
new_labels = np.asarray(data_secir['labels'][0]).reshape(
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
with open("inputs_120days_graphdata.pickle", 'wb') as f:
      pickle.dump(inputs, f) 
with open("labels_120days_graphdata.pickle", 'wb') as f:
      pickle.dump(target, f) 

###################################################################################################

# load the saved input and label file 
file_i = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_inputs_graphdata/inputs_60days_graphdata.pickle', 'rb')
test_inputs  = pickle.load(file_i)

file_l = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_labels_graphdata/labels_60days_graphdata.pickle', 'rb')
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

with open("/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_model_weights_graphdata/best_weights_ARMAConv_60days_nodamp_graphdata.pickle", "rb") as fp:
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

def compartments_nodescomparison(node_id1, node_id2, node_id3):

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
    
    fig, ((ax1, ax2,ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, sharey=True, figsize =(8,8), layout='constrained')
    axs = [ax1, ax2, ax3, ax4, ax5, ax6]
    xlabels=['MAPE', 'MAPE', 'MAPE', 'MAE', 'MAE', 'MAE']
    titles = ['Berlin MAPE', 'Mettmann MAPE', 'Havelland MAPE', 'Berlin MAE', 'Mettmann MAE', 'Havelland MAE']
    p_values = [2, 10, 2, 200, 1000, 800]

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

    plt.savefig("GNN_compartmentdistribution_MAEtraining_graphdata.png")


def populations_metrics_sctter(populations, mape_400, mae_400):
    plt.figure().clf()
    fig, axs = plt.subplots(1,2 ,figsize=(8,5))
    axs[0].scatter(populations, mape_400)
    axs[0].set_xlabel('Population size')
    axs[0].set_ylabel('MAPE')
    axs[0].set_xlim(0,500000)
    axs[0].set_title('Scatterplot of population vs MAPE ')

    axs[1].scatter(populations, mae_400)
    axs[1].set_xlabel('Population size')
    axs[1].set_ylabel('MAE')
    axs[1].set_xlim(0,500000)
    axs[1].set_title('Scatterplot of population vs MAE')
    plt.savefig("GNN_pop_mae_mape_60days.png")

def MAPE_MAE_distribution(mape_400, mae_400):
    plt.figure().clf()
    fig, axs = plt.subplots(1,2 ,figsize=(8,5))
    axs[0].hist(mape_400, bins = 40)
    axs[0].set_title('Hist MAPE ')
    axs[0].set_xlabel('MAPE')
    axs[0].set_ylabel('Count')

    axs[1].hist(mae_400, bins = 40 )
    axs[1].set_xlabel('MAE')
    axs[1].set_ylabel('Count')
    axs[1].set_title('Hist MAE ')
    plt.savefig("hist_MAPE_MAE_60w.png")


def plot_mobility():

    adjacency_matrix = np.asarray(sub_matrix)
    
    df_commuter = pd.DataFrame(data = adjacency_matrix)
    # Calculate sum of each row
    row_sums = df_commuter.sum(axis=1)
    # Calculate sum of each column
    column_sums = df_commuter.sum(axis=0)
    
    # Annahme: rows sind die herausfahrenden und columns sind die hereinfahrenden individuem (getestet an berlin und Havelland)
    df = pd.DataFrame({'incomming':column_sums, 'outgoing': row_sums, 'mape': mape_400, 'mae':mae_400})
    df['sum_commuters'] = df['incomming']+ df['outgoing']
    df['pop_size']  = populations
    df['relative_traffic'] = df['sum_commuters']/ df['pop_size']


    plt.figure().clf()
    fig, ax = plt.subplots()
    ax.scatter(df['relative_traffic'], df['mape'])
    ax.set_title('realive traffic vs. MAPE ')
    #ax.set_xlabel('number of commuters')
    ax.set_ylabel('MAPE')
    plt.savefig("relative_traffic_vs_mape.png")

    # plt.figure().clf()
    # fig, axs = plt.subplots(1,2 ,figsize=(8,5))
    # axs[0].scatter(df['incomming'], df['mape'])
    # axs[0].set_title('Incomming commuters vs. MAPE ')
    # axs[0].set_xlabel('number of commuters')
    # axs[0].set_ylabel('MAPE')

    # axs[1].scatter(df['incomming'], df['mae'])
    # axs[1].set_title('Incomming commuters vs. MAE ')
    # axs[1].set_xlabel('number of commuters')
    # axs[1].set_ylabel('MAE')
    # plt.savefig("incomming_commuters_vs_error.png")

    
def absolute_relative_traffic():
    adjacency_matrix = np.asarray(sub_matrix)

    df_commuter = pd.DataFrame(data=adjacency_matrix)
    row_sums = df_commuter.sum(axis=1)
    column_sums = df_commuter.sum(axis=0)

    df = pd.DataFrame({'incoming': column_sums, 'outgoing': row_sums, 'mape': mape_400, 'mae': mae_400})
    df['sum_commuters'] = df['incoming'] + df['outgoing']
    df['pop_size'] = populations
    df['relative_traffic'] = df['sum_commuters'] / df['pop_size']

    # Plotting
    plt.figure().clf()
    fig, axs = plt.subplots(2, 1, figsize=(10, 8))  # Create a 2x1 grid layout

    axs[0].scatter(df['sum_commuters'], df['mape'])
    axs[0].set_title('Traffic vs. MAPE')
    axs[0].set_xlabel('Traffic')
    axs[0].set_ylabel('MAPE')

    axs[1].scatter(df['relative_traffic'], df['mape'])
    axs[1].set_title('Relative Traffic vs. MAPE')
    axs[1].set_xlabel('Relative Traffic')
    axs[1].set_ylabel('MAPE')

    plt.subplots_adjust(hspace=0.5)  # Adjust the vertical spacing between subplots

    plt.savefig("absolute_and_relative_traffic_vs_mape.png")

    plt.savefig("absolute_and_relative_traffic_vs_mape.png")

def normalrange_vs_adjustedrange_plotdays():
    w_30 = pd.DataFrame(data = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframes_w/GNNtype1_ARMA_30days_w_mittel.csv'))
    #w_30= w_30[w_60.columns]
    w_60 = pd.DataFrame(data = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframes_w/GNNtype1_ARMA_60days_w.csv'))
    w_90 = pd.DataFrame(data = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframes_w/GNNtype1_ARMA_90days_w.csv'))
    w_120 = pd.DataFrame(data = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframes_w/GNNtype1_ARMA_120days_w.csv'))

    df_w = pd.concat([w_30, w_60, w_90, w_120])
    df_w['days'] = [30, 60, 90, 120]

    b_30 = 2.14
    b_60 = pd.DataFrame(data = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_gridsearch_2024/GNNtype1_ARMA_60days.csv'))
    b_90 = pd.DataFrame(data = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_gridsearch_2024/GNNtype1_ARMA_90days.csv'))
    b_120 = pd.DataFrame(data = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_gridsearch_2024/GNNtype1_ARMA_120days.csv'))
    df_baseline = pd.concat([b_60, b_90, b_120])
    

    MAPE = {'baseline': np.insert(df_baseline['kfold_test'],0, b_30), 
            'adjusted':df_w['kfold_test']}

    days = df_w['days'].values
    #df_bar=df_opt[['optimizer',  'kfold_test']]


    x = np.arange(1,len(days)+1)  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained')

    for attribute, measurement in MAPE.items():
        offset = width * multiplier
        rects = ax.bar(x +offset, measurement.round(2), width, label=attribute)
        ax.bar_label(rects, padding=3, fontsize = 10)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.

    #ax.set_yticks(x+width, layers)
    ax.legend(loc='upper right', ncols=3)
    ax.set_xticks(x + width, days)
    #ax.set_xticks(x, labels=x, fontsize = 12)
    
    ax.set_ylabel('test MAPE')
    ax.set_xlabel('number of days in prediction')
    #ax.set_title('')
    plt.savefig("GNN_days_baseline_vs_adjusted.png")



# compartment error plots 
# for secir simple
def compartment_error_simple():
                    
    layer_name='ARMAConv'
    number_of_layers = 4
    number_of_channels = 512
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

    class Net(Model):
                def __init__(self):
                    super().__init__()

                    self.conv1 = layer(number_of_channels,order = 3,  activation='relu')
                    self.conv2 = layer(number_of_channels, order = 3, activation='relu')
                    self.conv3 = layer(number_of_channels, order = 3,  activation='relu')
                    self.conv4 = layer(number_of_channels, order = 3,  activation='relu')
                    self.dense = Dense(data.n_labels, activation="linear")

                def call(self, inputs):
                    x, a = inputs
                    a = np.asarray(a)

                    x = self.conv1([x, a])
                    x = self.conv2([x, a])
                    x = self.conv3([x, a])
                    x = self.conv4([x,a])

                    output = self.dense(x)

                    return output

    
    df_plot = pd.DataFrame(columns = ['dampings', 'test_MAPE'])
    dampings =[1,2,3,4,5]
    
    #dampings = [1,2,3,4]

    path_inputs ='/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_inputs_graphdata'
    path_labels = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_labels_graphdata'
    path_weights = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_model_weights_graphdata'


    filenames_i = ['inputs_30days_graphdata.pickle', 'inputs_60days_graphdata.pickle', 'inputs_90days_graphdata.pickle']
        

    filenames_l = ['labels_30days_graphdata.pickle', 'labels_60days_graphdata.pickle', 'labels_90days_graphdata.pickle']
          
    saved_weights = ['best_weights_ARMAConv_30days_nodamp_graphdata.pickle', 'best_weights_ARMAConv_60days_nodamp_graphdata.pickle', 
                     'best_weights_ARMAConv_90days_nodamp_graphdata.pickle']
        
    plotnames = ['150_GNN_MAE_and_MAPE_30', '150_GNN_MAE_and_MAPE_60', '150_GNN_MAE_and_MAPE_90']
    days = [30,60,90]
    for name_i, name_l, weights, d , plotname in zip(filenames_i, filenames_l, saved_weights, days, plotnames):
    # load the saved input and label file 
        file_i = open(os.path.join(path_inputs, name_i), 'rb')
        test_inputs  = pickle.load(file_i)

        file_l = open(os.path.join(path_labels, name_l), 'rb')
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

        with open(os.path.join(path_weights, weights), "rb") as fp:
            b = pickle.load(fp)
        best_weights = b
        model.set_weights(best_weights)
        
        pred = model(test_inputs, training=False)

        pred_reversed = np.expm1(pred)
        labels_reversed = np.expm1(test_labels)

        mae =  np.mean((abs(labels_reversed - pred_reversed)), axis = 0).transpose()
        mape = 100*np.mean(abs((test_labels - pred)/test_labels), axis = 0).transpose()

        
        mape_mean_nodes = np.mean(mape, axis = 1)
        mae_mean_nodes = np.mean(mae, axis = 1)

        #def reshape_for_plot(array):
        #reshaped = []
    #for sample in array:
        sample_mape = mape_mean_nodes.transpose().reshape(d,48)
        compartment_sum_mape = []
        for day in sample_mape:                  
            compartment_sum_mape.append(day.reshape(-1, 8).transpose().mean(axis = 1))
                            
        sample_mae = mae_mean_nodes.transpose().reshape(d,48)
        compartment_sum_mae = []
        for day in sample_mae:                  
            compartment_sum_mae.append(day.reshape(-1, 8).transpose().mean(axis = 1))


        mae_plot = np.asarray(compartment_sum_mae).transpose()
        mape_plot = np.asarray(compartment_sum_mape).transpose()

        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(10,13))
        
        #fig, axes = plt.subplots(nrows=2, ncols=4, sharey=False)
        axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
        infectionstates = ['Susceptible','Exposed', 'InfectedNoSymptoms', 'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']
        for ax, c, ms, ma in zip(axes, infectionstates, mae_plot, mape_plot):
                
                color = 'tab:blue'
                ax.plot(ms, label ='MAE', color = color)
                ax.set_xlabel('Number of days')
                ax.set_ylabel('MAE', color = color)
                ax.tick_params(axis='y', labelcolor=color)
                ax.set_title(c, fontsize = 10)

                ax2 = ax.twinx()
                color = 'tab:green'
                ax2.set_ylabel('MAPE', color = color)
                ax2.plot(ma, color = color, label = 'MAPE' , linestyle = '--')
                ax2.tick_params(axis='y', labelcolor=color)
                fig.tight_layout()

        
        ax7.set_xlabel('Days')
        ax8.set_xlabel('Days')
            
        plt.savefig(plotname)
