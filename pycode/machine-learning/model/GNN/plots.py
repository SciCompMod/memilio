
import os
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import pickle
import json 
import tensorflow as tf 
from spektral.layers import ARMAConv
from keras.models import Model
from keras.layers import Dense


# def get_population(path="data\pydata\Germany\county_population_dim401.json"):
def get_population(path="/home/schm_a45/Documents/Code/memilio/memilio/data/pydata/Germany/county_population.json"):

    with open(path) as f:
        data = json.load(f)
    population = []
    for data_entry in data:
        population_county = []
        population_county.append(
            data_entry['<3 years'] + data_entry['3-5 years'] / 2)
        population_county.append(data_entry['6-14 years'])
        population_county.append(
            data_entry['15-17 years'] + data_entry['18-24 years'] +
            data_entry['25-29 years'] + data_entry['30-39 years'] / 2)
        population_county.append(
            data_entry['30-39 years'] / 2 + data_entry['40-49 years'] +
            data_entry['50-64 years'] * 2 / 3)
        population_county.append(
            data_entry['65-74 years'] + data_entry['>74 years'] * 0.2 +
            data_entry['50-64 years'] * 1 / 3)
        population_county.append(
            data_entry['>74 years'] * 0.8)

        population.append(population_county)
    return population 


def plot_mutliple_damp():

    # path = os.path.dirname(os.path.realpath(__file__))
    # path_data = os.path.join(
    #     os.path.dirname(
    #         os.path.realpath(os.path.dirname(os.path.realpath(path)))),
    #     'dataframes_dampingexperiments')
    path_data = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframes_dampingexperiments'

    filenames = ['GNNtype1_ARMA_100days_1damp.csv', 'GNNtype1_ARMA_100days_2damp.csv', 'GNNtype1_ARMA_100days_3damp.csv', 
                 'GNNtype1_ARMA_100days_4damp.csv',
                 'GNNtype1_ARMA_100days_5damp.csv']
    
    df_plot = pd.DataFrame(columns = ['dampings', 'test_MAPE'])
    dampings =[1,2,3,4,5]
    #dampings = [1,2,3,4]
    for damp, filename in zip(dampings, filenames):

        file = open(os.path.join(path_data, filename),'rb')
        data_secir = pd.read_csv(file)
        df_plot.loc[len(df_plot)] = [damp, data_secir['kfold_test'][0]]
        df_plot['dampings']= df_plot['dampings'].astype('int')

    plt.figure().clf() 
    fig, ax = plt.subplots()
    rects = ax.bar(df_plot['dampings'].round(0), df_plot['test_MAPE'].round(4))

    ax.set_ylabel('Test MAPE')
    ax.set_xlabel('Number of dampings')
    ax.set_xticks(dampings)
    ax.set_title('Test MAPE for different number of dampings')

    #large_bars = [p if p > 2 else '' for p in df_bar['kfold_test'].round(4)]
    small_bars = [p if p <= 10 else '' for p in df_plot['test_MAPE'].round(4)]

    ax.bar_label(rects, small_bars,
                  padding=5, color='black')
    #ax.bar_label(rects, large_bars,
    #              padding=-40, color='white')
    ax.margins(y=0.1)

    plt.show()
    plt.savefig("damping_ARMAConv.png")

def  populationsize_mape():
        foldername = 'data_GNN_400pop_one_var_damp_100days_1k_withmatrix'
        filename = 'data_secir_age_groups.pickle'
        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                        os.path.dirname(os.path.realpath(path)))), foldername)
        
        file = open(os.path.join(path_data,filename), 'rb')
        data = pickle.load(file)

        data['inputs'] = (data['inputs'])
        data['labels'] = ((data['labels']))
        data['damping_coeff'] =(data['damping_coeff'])
        data['damping_day'] = ((data['damping_day']))
        data['damped_matrix'] = ((data['damped_matrix']))

        len_dataset = data['inputs'][0].shape[0]
        numer_of_nodes = np.asarray(data['inputs']).shape[0]
        shape_input_flat = np.asarray(
                data['inputs']).shape[2]*np.asarray(data['inputs']).shape[3]
        shape_labels_flat = np.asarray(
                data['labels']).shape[2]*np.asarray(data['labels']).shape[3]


        new_inputs = np.asarray(
                data['inputs']).reshape(
                len_dataset, numer_of_nodes, shape_input_flat)
        new_labels = np.asarray(data['labels']).reshape(
                len_dataset, numer_of_nodes, shape_labels_flat)

        n_days = int(new_labels.shape[2]/48)

        damping_factors = data['damping_coeff'] # one for every run
        damping_days = data['damping_day'] # one for every rum 
        contact_matrices = data['damped_matrix']

        n_runs = new_inputs.shape[0]
        n_pop = new_inputs.shape[1]

        n_dampings = 1  # for one damping
        #n_dampings = np.asarray(data['damping_day']).shape[2] # for more than one damping


        inputs_with_damp = np.dstack((new_inputs,(np.asarray(damping_factors).reshape([n_runs,n_pop,n_dampings])),
                                    (np.asarray(damping_days).reshape([n_runs,n_pop,n_dampings])),
                                    (np.asarray(contact_matrices).reshape([n_runs,n_pop,36*n_dampings]))))
        


       # n_labels = 

        

        path_population = os.path.abspath(
                r"data//pydata//Germany//county_population.json")
        populations = get_population(path_population )
        populations = np.asarray(populations).sum(axis = 1)

        # load saved model 


        foldername = 'data_GNN_nodamp_4000pop_1k_90days_24'
        filename = 'data_secir_age_groups.pickle'
        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                        os.path.dirname(os.path.realpath(path)))), foldername)
        
        file = open(os.path.join(path_data,filename), 'rb')
        data = pickle.load(file)

        data['inputs'] = (data['inputs'])
        data['labels'] = ((data['labels']))

        #new_model = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/ARMAConv_100days_one_damp_saved_model')
        #pred = new_model.predict(test_inputs)
        



        layer = ARMAConv
        number_of_layer = 3
        channels = 128

        class Net(Model):
                def __init__(self):
                    super().__init__()

                    self.conv1 = layer(channels, activation='elu')
                    self.conv2 = layer(channels, activation='elu')
                    self.conv3 = layer(channels, activation='elu')
                    self.dense = Dense(1440, activation="linear")

                def call(self, inputs):
                    x, a = inputs
                    #a = np.asarray(a)

                    x = self.conv1([x, a])
                    x = self.conv2([x, a])
                    x = self.conv2([x, a])

                    output = self.dense(x)

                    return output

        model = Net()



        f = open('/home/schm_a45/Documents/Code/memilio/memilio/GNN_90days.json')
    
        model_json = json.load(f)
        from keras.models import model_from_json
        model = model_from_json(model_json) 
        
        #model.load_weights("/home/schm_a45/Documents/code3/memilio/pycode/machine-learning/ARMAConv_90days_saved_model_test.h5")
        model.build(input_shape = (100,400,240))


def plot_commuter_data():
    commuterdata = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/data/commuter_migration_scaled.txt', sep = " ", header = None)

    fig, ax = plt.subplots()
    
    im = ax.imshow(commuterdata.values)

    # Show all ticks and label them with the respective list entries
    #ax.set_xticks(np.arange(len(commuterdata.columns)), labels=commuterdata.columns)
    #ax.set_yticks(np.arange(len(commuterdata.index)), labels=commuterdata.index)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")


    cbar_kw = {}        
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel('commuters', rotation=-90, va="bottom")        


    ax.set_title("Visualization of commuter matrix")
    fig.tight_layout()
    plt.show()
    plt.savefig("commuter_matrix.png")
    import networkx as nx 


    commuterdata[commuterdata <1] = 0
    G = nx.from_numpy_array(np.asarray(commuterdata))

    pos = nx.spring_layout(G)
    plt.figure(figsize=(8,8))
    nx.draw(G, node_size = 5 ,  node_color='blue', edge_color='gray', alpha=0.5, with_labels=False)
    plt.savefig('networkx_graph_network.png')