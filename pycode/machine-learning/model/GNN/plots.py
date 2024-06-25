
import os
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib.colors import ListedColormap
import pickle
import json 
import tensorflow as tf 
import ast
from spektral.layers import ARMAConv
from keras.models import Model
from keras.layers import Dense
import seaborn as sns
from memilio.simulation.secir import InfectionState
from keras.layers import Dense
from keras.metrics import mean_absolute_percentage_error
from keras.models import Model
import spektral
from spektral.data import  MixedLoader
from spektral.layers import ARMAConv
from spektral.transforms.normalize_adj import NormalizeAdj
from spektral.utils.convolution import gcn_filter, normalized_laplacian, rescale_laplacian, normalized_adjacency


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



def lineplot_number_of_days():

        #ilenames =[ "GNNtype1_ARMA_30days_w_mittel.csv", "GNNtype1_ARMA_60days_w.csv", "GNNtype1_ARMA_90days_w.csv", "GNNtype1_ARMA_120days_w.csv"]
        filenames = ['ARMA_4_512_nodamp_30days_graphdata_withCV_newrange.csv', 'ARMA_4_512_nodamp_60days_graphdata_withCV_newrange.csv', 
                     'ARMA_4_512_nodamp_90days_graphdata_withCV_newrange.csv']
        #days = [30,60,90,120]
        days = [30,60,90]
        MAPE = []
        
        for file in filenames: 
                
                path = os.path.dirname(os.path.realpath(__file__))
                path_data = os.path.join(os.path.dirname(os.path.realpath(
                        os.path.dirname(os.path.realpath(path)))), 'model_evaluations_graphdata')
                
                filename = os.path.join(path_data, file)

                df = pd.read_csv(filename)

                MAPE.append(np.mean(df['kfold_test']))
     
                
        plt.figure().clf()
        fig, ax = plt.subplots(figsize =(8,5))
        ax.plot(days, MAPE,  marker = 'o' )
        ax.set_xticks(days)
        ax.set_xlabel('Number of days')
        ax.set_ylabel('MAPE')
        ax.set_title('MAPE for number of days to be predicted')
        plt.savefig("plot_days_GNN_w.png")



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


def boxplot_inputs():

        file = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/data_GNN_nodamp_400pop_1k_60days_1_24/data_secir_age_groups.pickle', 'rb')
        #file_w2 = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/data_GNN_nodamp_400pop_1k_30days_w1/data_secir_age_groups.pickle', 'rb')
        file_w = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/data_GNN_nodamp_400pop_1k_60days_w/data_secir_age_groups.pickle', 'rb') 
       
        data = pickle.load(file)
        data = np.expm1(data['inputs'])
        data_w = pickle.load(file_w)
        data_w = np.expm1(data_w['inputs'])
        #data_w2 = pickle.load(file_w2)
        #data_w2 = np.expm1(data_w2['inputs'])


        node1 = 324     # Berlin --> biggest city
        node2 = 334     # Havelland --> small county (?)
        node3 = 74      # Mettmann --> medium county (?)


        compartment_array = []
        for compartment in InfectionState.values():
                compartment_array.append(compartment)
        index = [3,5]
        compartments_cleaned= np.delete(compartment_array, index)                

        index = [str(compartment).split('.')[1] for compartment in compartments_cleaned]

        node_id = 334
        #summarize age group information 
        
        indices = [0,1,2,3,4,5,6,7]             
            
        sum_all = []
        for run in data[node_id]:
                sum_run = []
                for day in run:  
                    sum_day = []                    
                    for i in indices:
                        x = day[i::8]
                        sum_day.append(x.sum())
                    sum_run.append(sum_day)
                sum_all.append(sum_run)


        sum_all_w = []
        for run in data_w[node_id]:
                sum_run = []
                for day in run:  
                    sum_day = []                    
                    for i in indices:
                        x = day[i::8]
                        sum_day.append(x.sum())
                    sum_run.append(sum_day)
                sum_all_w.append(sum_run)
            	    
        # sum_all_w2 = []
        # for run in data_w2[node_id]:
        #         sum_run = []
        #         for day in run:  
        #             sum_day = []                    
        #             for i in indices:
        #                 x = day[i::8]
        #                 sum_day.append(x.sum())
        #             sum_run.append(sum_day)
        #         sum_all_w2.append(sum_run)
            	    

        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(7, 9))
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

        #for ax, compartment, d, dw1, dw2 in zip(axes, index, np.asarray(sum_all).transpose(), np.asarray(sum_all_w).transpose(), np.asarray(sum_all_w2).transpose()):
        for ax, compartment, d, dw1 in zip(axes, index, np.asarray(sum_all).transpose(), np.asarray(sum_all_w).transpose()):
          
                d_df = pd.DataFrame(data = d)
                d_df = pd.melt(d_df.transpose(),var_name="Day") 
                d_df['type'] = 'dataset 1'  

                dw1_df = pd.DataFrame(data = dw1)
                dw1_df = pd.melt(dw1_df.transpose(), var_name="Day") 
                dw1_df['type'] = 'dataset 2'  
        
                # dw2_df = pd.DataFrame(data = dw2)
                # dw2_df = pd.melt(dw2_df.transpose(), var_name="Day") 
                # dw2_df['type'] = 'w2' 
                #df_all = pd.concat([d_df, dw1_df, dw2_df], ignore_index=True)
                df_all = pd.concat([d_df, dw1_df], ignore_index=True)

                sns.boxplot(ax = ax, x='Day', y='value', data = df_all, hue = 'type', palette = 'Set1', width = 0.8, legend = 'auto')
                ax.set_title(compartment, fontsize = 10)
                ax.set_xticklabels(np.arange(1,6))
                ax.legend().set_visible(False)


                handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.55, 0.001), frameon=False)
        plt.tight_layout()
                 
        plt.savefig("GNN_boxplot_input_compartments_b_w_Havelland_60.png")


def plot_mutliple_damp_w():
    
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


    filenames_i = ['inputs_100days_1damp_graphdata.pickle', 'inputs_100days_2damp_graphdata.pickle', 'inputs_100days_3damp_graphdata.pickle', 
                   'inputs_100days_4damp_graphdata.pickle', 'inputs_100days_5damp_graphdata.pickle']
    

    filenames_l = ['labels_100days_1damp_graphdata.pickle', 'labels_100days_2damp_graphdata.pickle', 'labels_100days_3damp_graphdata.pickle', 
                   'labels_100days_4damp_graphdata.pickle', 'labels_100days_5damp_graphdata.pickle']
    
    saved_weights = ['best_weights_ARMAConv_100days_1damp_graphdata.pickle', 'best_weights_ARMAConv_100days_2damp_graphdata.pickle', 
                     'best_weights_ARMAConv_100days_3damp_graphdata.pickle','best_weights_ARMAConv_100days_4damp_graphdata.pickle', 
                     'best_weights_ARMAConv_100days_5damp_graphdata.pickle']
        
    for name_i, name_l, weights, d in zip(filenames_i, filenames_l, saved_weights, dampings):
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
        df_plot.loc[len(df_plot)] = [d, mape.mean()]



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
    plt.savefig("damping_ARMAConv_graphdata_noCV.png")


def plot_mutliple_damp():

    # path = os.path.dirname(os.path.realpath(__file__))
    # path_data = os.path.join(
    #     os.path.dirname(
    #         os.path.realpath(os.path.dirname(os.path.realpath(path)))),
    #     'dataframes_dampingexperiments')
    path_data = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model_evaluation_graphdata'

    filenames = ['GNNt_ARMA_100days_1damp_graphdata_withCV.csv', 'GNNt_ARMA_100days_2damp_graphdata_withCV.csv', 
                 'GNNt_ARMA_100days_3damp_graphdata_withCV.csv',
                 'GNNt_ARMA_100days_4damp_graphdata_withCV.csv',
                 'GNNt_ARMA_100days_5damp_graphdata_withCV.csv']
    
    df_plot = pd.DataFrame(columns = ['dampings', 'test_MAPE'])
    dampings =[1,2,3,4,5]
    #dampings = [1,2,3]
    for damp, filename in zip(dampings, filenames):

        file = open(os.path.join(path_data, filename),'rb')
        data_secir = pd.read_csv(file)
        df_plot.loc[len(df_plot)] = [damp, data_secir['kfold_test'][0]]
    df_plot['dampings']= df_plot['dampings'].astype('int')

    plt.figure().clf() 
    fig, ax = plt.subplots()
    rects = ax.bar(df_plot['dampings'].round(0), df_plot['test_MAPE'].round(4))

    ax.set_ylabel('MAPE')
    ax.set_xlabel('Number of dampings')
    ax.set_xticks(dampings)
    #ax.set_title('Mean cros MAPE for different number of dampings')

    #large_bars = [p if p > 2 else '' for p in df_bar['kfold_test'].round(4)]
    small_bars = [p if p <= 10 else '' for p in df_plot['test_MAPE'].round(4)]

    ax.bar_label(rects, small_bars,
                  padding=5, color='black')
    #ax.bar_label(rects, large_bars,
    #              padding=-40, color='white')
    ax.margins(y=0.1)

    plt.show()
    plt.savefig("damping_ARMAConv_graphdata_withCV.png")


def ARMA_dampings_scatter():


    path_data = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model_evaluation_graphdata'

    filenames = [
        'GNNt_ARMA_100days_1damp_graphdata_withCV.csv',
        'GNNt_ARMA_100days_2damp_graphdata_withCV.csv', 
        'GNNt_ARMA_100days_3damp_graphdata_withCV.csv',
        'GNNt_ARMA_100days_4damp_graphdata_withCV.csv',
        'GNNt_ARMA_100days_5damp_graphdata_withCV.csv'
    ]

    df_plot = pd.DataFrame(columns=['dampings', 'test_MAPE'])
    dampings = [1, 2, 3, 4, 5]

    for damp, filename in zip(dampings, filenames):
        file = open(os.path.join(path_data, filename), 'rb')
        data_secir = pd.read_csv(file)
        df_plot.loc[len(df_plot)] = [damp, data_secir['all_test_scores'].values]
        df_plot['dampings'] = df_plot['dampings'].astype('int')

    # Convert string representations of arrays to actual arrays
    df_plot['test_MAPE'] = df_plot['test_MAPE'].apply(lambda x: np.array(ast.literal_eval(x[0])))

    # Explode the dataframe based on the 'all_test_scores' column
    df_expanded = df_plot.explode('test_MAPE')
    df_expanded.reset_index(drop=True, inplace=True)

    # Define positions
    array = ['1st split', '2nd split', '3rd split', '4th split', '5th split']
    df_expanded['position'] = np.tile(array, len(df_plot))

    # Define your own custom colors
    custom_colors = ['red', 'blue', 'green', 'orange', 'purple']

    # Create a ListedColormap using the custom colors
    custom_cmap = ListedColormap(custom_colors)

    # Plot the scatter plot with color coding based on 'position'
    plt.figure(figsize=(6, 5))  # Adjust figure size as needed
    for i, position in enumerate(array):
        plt.scatter(df_expanded[df_expanded['position'] == position]['dampings'], df_expanded[df_expanded['position'] == position]['test_MAPE'], c=custom_colors[i], label=position, cmap=custom_cmap)

    plt.legend()
    plt.ylabel('MAPE')
    plt.xlabel('Number of Dampings')
    plt.xticks(ticks=dampings, labels=dampings)  # Set the x-ticks and labels explicitly
    plt.tight_layout()  # Adjust layout to prevent label overlap
    plt.savefig("ARMAConv_CV_testscores_damp.png")



def MAPE_MAE_ARMA():
        
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

    
    df_plot = pd.DataFrame(columns = ['days', 'test_MAPE', 'test_MAE'])

    path_inputs ='/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_inputs_graphdata'
    path_labels = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_labels_graphdata'
    path_weights = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_model_weights_graphdata'


    filenames_i = ['inputs_30days_graphdata.pickle', 'inputs_60days_graphdata.pickle', 'inputs_90days_graphdata.pickle']             

    filenames_l = ['labels_30days_graphdata.pickle', 'labels_60days_graphdata.pickle', 'labels_90days_graphdata.pickle']
           
    saved_weights = ['best_weights_ARMAConv_30days_nodamp_graphdata.pickle', 'best_weights_ARMAConv_60days_nodamp_graphdata.pickle', 
                     'best_weights_ARMAConv_90days_nodamp_graphdata.pickle']
    
    days = [30,60,90]
        
    for name_i, name_l, weights, d in zip(filenames_i, filenames_l, saved_weights, days):
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
        df_plot.loc[len(df_plot)] = [d, mape.mean(), mae.mean()]
       



    plt.figure().clf()
    fig, ax1 = plt.subplots(figsize =(8,4))

    color = 'tab:red'
    ax1.plot(days, df_plot['test_MAPE'], label ='mape', color = color, marker = 'o' )
    ax1.set_xticks(days)
    ax1.set_xlabel('Number of days')
    ax1.set_ylabel('MAPE', color = color)
    ax1.tick_params(axis='y', labelcolor=color)
    #ax1.set_title('MAPE for number of days to be predicted')
        
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('MAE', color = color)
    ax2.plot(days, df_plot['test_MAE'], color = color, label = 'groups', marker = 'x')
    ax2.tick_params(axis='y', labelcolor=color)
    fig.legend(loc= 'lower left', fontsize = 9)
    fig.tight_layout()
        
    plt.savefig("ARMA_days_MAE_MAPE.png")



def plot_MAE_MAPE_compartment_ARMA():
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

    
    df_plot = pd.DataFrame(columns = ['days', 'test_MAPE', 'test_MAE'])

    path_inputs ='/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_inputs_graphdata'
    path_labels = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_labels_graphdata'
    path_weights = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_model_weights_graphdata'


    filenames_i = ['inputs_30days_graphdata.pickle', 'inputs_60days_graphdata.pickle', 'inputs_90days_graphdata.pickle']             

    filenames_l = ['labels_30days_graphdata.pickle', 'labels_60days_graphdata.pickle', 'labels_90days_graphdata.pickle']
           
    saved_weights = ['best_weights_ARMAConv_30days_nodamp_graphdata.pickle', 'best_weights_ARMAConv_60days_nodamp_graphdata.pickle', 
                     'best_weights_ARMAConv_90days_nodamp_graphdata.pickle']
    
    days = [30,60,90]
        
    plotnames = ['30_GNN_MAE_and_MAPE_graphdata.png','60_GNN_MAE_and_MAPE_graphdata.png', '90_GNN_MAE_and_MAPE_graphdata.png' ]
    for name_i, name_l, weights, d ,plotname in zip(filenames_i, filenames_l, saved_weights, days, plotnames):
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
        df_plot.loc[len(df_plot)] = [d, mape.mean(), mae.mean()]
       

                
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(10,13))
        
        #fig, axes = plt.subplots(nrows=2, ncols=4, sharey=False)
        axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
        infectionstates = ['Susceptible','Exposed', 'InfectedNoSymptoms', 'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']
        for ax, c, ms, ma in zip(axes, infectionstates, mae, mape):
                
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


                #ax.plot(ms, label ='MSE')
                #ax.plot(ma, label='MAPE', linestyle = '--')#

                #ax.set_title(c, fontsize = 10)
                
                #ax.legend(loc='upper right', ncols=3)
        
        ax7.set_xlabel('Days')
        ax8.set_xlabel('Days')
            
        plt.savefig(plotname)

     


def ARMA_days_scatter():


    path_data = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model_evaluations_graphdata'

    filenames = [
        'ARMA_4_512_nodamp_30days_graphdata_withCV_newrange.csv',
        'ARMA_4_512_nodamp_60days_graphdata_withCV_newrange.csv', 
         'ARMA_4_512_nodamp_90days_graphdata_withCV_newrange.csv']

    df_plot = pd.DataFrame(columns=['days', 'test_MAPE'])
    days = [30,60,90]

    for damp, filename in zip(days, filenames):
        file = open(os.path.join(path_data, filename), 'rb')
        data_secir = pd.read_csv(file)
        df_plot.loc[len(df_plot)] = [damp, data_secir['all_test_scores'].values]
        df_plot['days'] = df_plot['days'].astype('int')

    # Convert string representations of arrays to actual arrays
    df_plot['test_MAPE'] = df_plot['test_MAPE'].apply(lambda x: np.array(ast.literal_eval(x[0])))

    # Explode the dataframe based on the 'all_test_scores' column
    df_expanded = df_plot.explode('test_MAPE')
    df_expanded.reset_index(drop=True, inplace=True)

    # Define positions
    array = ['1st split', '2nd split', '3rd split', '4th split', '5th split']
    df_expanded['position'] = np.tile(array, len(df_plot))

    # Define your own custom colors
    custom_colors = ['red', 'blue', 'green', 'orange', 'purple']

    # Create a ListedColormap using the custom colors
    custom_cmap = ListedColormap(custom_colors)

    # Plot the scatter plot with color coding based on 'position'
    plt.figure(figsize=(6, 5))  # Adjust figure size as needed
    for i, position in enumerate(array):
        plt.scatter(df_expanded[df_expanded['position'] == position]['days'], df_expanded[df_expanded['position'] == position]['test_MAPE'], c=custom_colors[i], label=position, cmap=custom_cmap)

    plt.legend()
    plt.ylabel('MAPE')
    plt.xlabel('Number of Days')
    plt.xticks(ticks=days, labels=days)  # Set the x-ticks and labels explicitly
    plt.tight_layout()  # Adjust layout to prevent label overlap
    plt.savefig("ARMAConv_CV_testscores_days.png")
