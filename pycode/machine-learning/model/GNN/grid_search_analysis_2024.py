import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import pickle
import re
from matplotlib.colors import ListedColormap
from sklearn.model_selection import KFold
import seaborn as sns
from matplotlib.lines import Line2D


df_2024 = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_hyperparameters_2024/dataframe_hyperparameter_tuning_30days_oldrange_noCV.csv')
df_ARMA = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_hyperparameters_2024/dataframe_ARMAConv_30days_oldrange_noCV.csv')
#df_XENET = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_hyperparameters_we/dataframe_XENEtConv_Baseline_noCV.csv')
df_part1 = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_gridsearch_2024/dataframe_grid_search_networkarchitecture_GNN_type1.csv')
df_part2 = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_gridsearch_2024/dataframe_grid_searc_networkarchitecture_GNN_type1_part2.csv')
df_combined = pd.concat([df_part1 ,df_part2], ignore_index = True )

df_type1_best_binary = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_gridsearch_2024/GNNtype1_bestmodels_baseline_binaryadjacency.csv')
df_type1_best_nonbinary = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_gridsearch_2024/GNNtype1_bestmodels_baseline_nonbinaryadjacency.csv')

def heatmap(df):
    df_heatmap1 = pd.DataFrame(data = df.loc[(df['layer'] == 'GCNConv') & (df['number_of_layers']==3)][['activation', 'optimizer', 'kfold_test']])
    df_heatmap1= df_heatmap1.pivot(index='activation', columns='optimizer', values='kfold_test')

    fig, ax = plt.subplots()
    im = ax.imshow(df_heatmap1.values)

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(df_heatmap1.columns)), labels=df_heatmap1.columns)
    ax.set_yticks(np.arange(len(df_heatmap1.index)), labels=df_heatmap1.index)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(df_heatmap1.index)):
        for j in range(len(df_heatmap1.columns)):
            text = ax.text(j, i, np.around(df_heatmap1.values, decimals=2)[i, j],
                        ha="center", va="center", color="w")


    cbar_kw = {}        
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel('MAPE', rotation=-90, va="bottom")        


    ax.set_title("Activation function and optimizer for 3 GCNConv layers")
    fig.tight_layout()
    plt.show()
    plt.savefig("heatmap_activation_optimizer.png")


def heatmaps(df):


    df_1 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GCNConv')][['number_of_layers', 'channels', 'kfold_test']])
    df_1= df_1.pivot(index='number_of_layers', columns='channels', values='kfold_test')
    
    df_2 = pd.DataFrame(data =  df.loc[(df['layer'] == 'ARMAConv')][['number_of_layers', 'channels', 'kfold_test']])
    df_2= df_2.pivot(index='number_of_layers', columns='channels', values='kfold_test')
    
    df_3 = pd.DataFrame(data =  df.loc[(df['layer'] == 'APPNPConv')][['number_of_layers', 'channels', 'kfold_test']])
    df_3= df_3.pivot(index='number_of_layers', columns='channels', values='kfold_test')
    
    df_4 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GATConv')][['number_of_layers', 'channels', 'kfold_test']])
    df_4= df_4.pivot(index='number_of_layers', columns='channels', values='kfold_test')

        
    df_5 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GCSConv')][['number_of_layers', 'channels', 'kfold_test']])
    df_5= df_5.pivot(index='number_of_layers', columns='channels', values='kfold_test')


    plt.figure().clf() 
    layers = df['layer'].unique()



    fig, axs = plt.subplots(nrows = 3, ncols = 2, sharex=False, figsize = (7,9), constrained_layout = True)



    for ax, df_heatmap, name  in zip(axs.flat, [df_1 ,df_2, df_3, df_4, df_5], layers):
        
        im = ax.imshow(df_heatmap.values, cmap ='summer' )

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(df_heatmap.columns)), labels=df_heatmap.columns)
        ax.set_yticks(np.arange(len(df_heatmap.index)), labels=df_heatmap.index)

        ax.set_ylabel('number of layers')
        ax.set_xlabel('channels')

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(df_heatmap.index)):
            for j in range(len(df_heatmap.columns)):
                text = ax.text(j, i, np.around(df_heatmap.values, decimals=2)[i, j],
                            ha="center", va="center", color="k")
                
        ax.set_title('Model = '+name) 



    fig.colorbar(im, ax = axs, shrink=0.75, label = 'MAPE')
    
    fig.delaxes(axs[2][1])
    
    #plt.rcParams.update({'font.size': 25})
    plt.show()
    plt.savefig("GNN1_layers_channels_summer.png")


def barplot_baseline(df, df_GSC):

 
    # get performance for basic one layer model with relu and Adam 
    df_baseline = df.loc[(df['number_of_layers']=='1' )& (df['activation']=='relu') & (df['optimizer']=='Adam') ][['layer', 'kfold_test', 'kfold_train']]
    df_GSC = df_GSC[['layer', 'kfold_test', 'kfold_train']]
    df_baseline = pd.concat([df_baseline, df_GSC])
    df_baseline[['kfold_train', 'kfold_test']]=df_baseline[['kfold_train', 'kfold_test']].astype('float').round(2)
       

    MAPE = {'train MAPE': df_baseline['kfold_train'], 
            'test MAPE':df_baseline['kfold_test']}




    layers = df_baseline['layer'].values
    #df_bar=df_opt[['optimizer',  'kfold_test']]


    x = np.arange(1,len(layers)+1)  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained')

    for attribute, measurement in MAPE.items():
        offset = width * multiplier
        rects = ax.bar(x +offset, measurement, width, label=attribute)
        ax.bar_label(rects, padding=3, fontsize = 10)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.

    #ax.set_yticks(x+width, layers)
    ax.legend(loc='upper right', ncols=3)
    ax.set_xticks(x + width, layers)
    #ax.set_xticks(x, labels=x, fontsize = 12)
    


    ax.set_ylabel('test MAPE')
    ax.set_xlabel('Layers')
    ax.set_title('Mean Test MAPE for types of layers')
    plt.savefig("bar_GNN1_baselineLayers.png")

def plot_losses(df_XENET):
    train = df_XENET['train_losses']
    train =  [i.replace('[', '') for i in list(train[0].split())]
    train = train[:-1]
    train = list(np.float_(train))

    val = df_XENET['val_losses']
    val =  [i.replace('[', '') for i in list(val[0].split())]
    val = val[:-1]
    val = list(np.float_(val))

    plt.figure().clf() 
    plt.plot(train, label = 'train losses')
    plt.plot(val, label = 'validation losses')
    plt.legend()
    plt.ylabel('MAPE')
    plt.xlabel('Epochs')
    plt.title('Losses for XENEt baseline')
    plt.savefig("losses_XENEt.png")



def plot_binary_vs_nonbinary_adjacency(df_binary, df_nonbinary):
    df_binary = df_binary[['layer', 'number_of_layers', 'channels', 'kfold_train', 'kfold_test']]
    df_nonbinary = df_nonbinary[['layer', 'number_of_layers', 'channels', 'kfold_train', 'kfold_test']]

    MAPE = {'binary A': df_binary['kfold_test'], 
            'non-binary A':df_nonbinary['kfold_test']}


    layers = df_binary['layer'].values
    #df_bar=df_opt[['optimizer',  'kfold_test']]


    x = np.arange(1,len(layers)+1)  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0.5

    fig, ax = plt.subplots(layout='constrained')

    for attribute, measurement in MAPE.items():
        offset = width * multiplier
        rects = ax.bar(x +offset, round(measurement,2), width, label=attribute)
        ax.bar_label(rects, padding=3, fontsize = 10)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.

    #ax.set_yticks(x+width, layers)
    ax.legend(loc='upper left', ncols=3)
    ax.set_xticks(x + width, layers)
    #ax.set_xticks(x, labels=x, fontsize = 12)
    

    ax.set_ylabel('test MAPE')
    ax.set_xlabel('Layers')
    ax.set_title('Comparison of test MAPE for binary vs non-binary adjacency matrix')
    plt.savefig("GNN1_adjacency.png")


#plot_max_min_losses(min_loss_train, min_loss_val, max_loss_train, max_loss_val)
#heatmap()
#plot_layer()

#for i in df.columns[1:5]:
#    print('Mean Test MAPE by ' + i ,':',  pd.DataFrame(data = df[[i, 'kfold_test']]).groupby([i]).mean())


def boxplot_layers(df_combined):
    df = df_combined
    plt.figure().clf()   

    df_1 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GCNConv')][['optimizer', 'activation', 'kfold_test']])
    df_1= df_1.pivot(index='number_of_layers', columns='channels', values='kfold_test')
    
    df_2 = pd.DataFrame(data =  df.loc[(df['layer'] == 'ARMAConv')][['number_of_layers', 'channels', 'kfold_test']])
    df_2= df_2.pivot(index='number_of_layers', columns='channels', values='kfold_test')
    
    df_3 = pd.DataFrame(data =  df.loc[(df['layer'] == 'APPNPConv')][['number_of_layers', 'channels', 'kfold_test']])
    df_3= df_3.pivot(index='number_of_layers', columns='channels', values='kfold_test')
    
    df_4 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GATConv')][['number_of_layers', 'channels', 'kfold_test']])
    df_4= df_4.pivot(index='number_of_layers', columns='channels', values='kfold_test')
        
    df_5 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GCSConv')][['number_of_layers', 'channels', 'kfold_test']])
    df_5= df_5.pivot(index='number_of_layers', columns='channels', values='kfold_test')

    dataframes = [df_1, df_2, df_3, df_4, df_5]
    models = ['GCNConv', 'ARMAConv', 'APPNPConv', 'GATConv' ,'GCSConv']
    
    fig, ((ax1, ax2, ax3),(ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(9, 4), sharex=False, sharey=False,)
    axs = [ax1, ax2, ax3, ax4, ax5, ax6]

    for ax, model, df in zip(axs, models, dataframes):

        ax.boxplot(df)
        ax.set_title(model)
        ax.yaxis.grid(True)
        ax.set_xticks([y + 1 for y in range(len(df))],
                        labels=['1', '2', '3'])
        ax.set_xlabel('number of layers')
        ax.set_ylabel('test MAPE') 
    
    fig.tight_layout()
    fig.delaxes(ax6)

    plt.savefig("boxplot_GNN_layers.png")


def heatmaps_act_opt(df):
    df_1 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GCNConv') & (df['number_of_layers']==2)& (df['channels']==512)][['optimizer', 'activation', 'kfold_test']])
    df_1= df_1.pivot(index='optimizer', columns='activation', values='kfold_test')
    
    df_2 = pd.DataFrame(data =  df.loc[(df['layer'] == 'ARMAConv')&(df['number_of_layers']==3)& (df['channels']==128)][['optimizer', 'activation', 'kfold_test']])
    df_2= df_2.pivot(index='optimizer', columns='activation', values='kfold_test')
    
    df_3 = pd.DataFrame(data =  df.loc[(df['layer'] == 'APPNPConv')&(df['number_of_layers']==3)& (df['channels']==512)][['optimizer', 'activation', 'kfold_test']])
    df_3= df_3.pivot(index='optimizer', columns='activation', values='kfold_test')
    
    df_4 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GATConv')&(df['number_of_layers']==3)& (df['channels']==1024)][['optimizer', 'activation', 'kfold_test']])
    df_4= df_4.pivot(index='optimizer', columns='activation', values='kfold_test')
        
    df_5 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GCSConv')&(df['number_of_layers']==3)& (df['channels']==512)][['optimizer', 'activation', 'kfold_test']])
    df_5= df_5.pivot(index='optimizer', columns='activation', values='kfold_test')


    plt.figure().clf() 
    layers = df['layer'].unique()



    fig, axs = plt.subplots(nrows = 3, ncols = 2, sharex=False, figsize = (7,9), constrained_layout = True)



    for ax, df_heatmap, name  in zip(axs.flat, [df_1 ,df_2, df_3, df_4, df_5], layers):
        
        im = ax.imshow(df_heatmap.values, cmap ='summer_r' )

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(df_heatmap.columns)), labels=df_heatmap.columns)
        ax.set_yticks(np.arange(len(df_heatmap.index)), labels=df_heatmap.index)

        ax.set_ylabel('optimizer')
        ax.set_xlabel('activation')

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(df_heatmap.index)):
            for j in range(len(df_heatmap.columns)):
                text = ax.text(j, i, np.around(df_heatmap.values, decimals=2)[i, j],
                            ha="center", va="center", color="k")
                
        ax.set_title('Model = '+name) 



    fig.colorbar(im, ax = axs, shrink=0.75, label = 'Test MAPE')
    
    fig.delaxes(axs[2][1])
    
    #plt.rcParams.update({'font.size': 25})
    plt.show()
    plt.savefig("GNN1_optimizer_activation.png")


def ARMA_heatmap():
    df = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/machine-learning/dataframe_gridsearch_2024/ARMAConv_act_opt_run2.csv')        
        
    optimizer = []
    for i in df['optimizer']:
        optimizer.append(re.sub('[^a-zA-Z]+', '', i.split(' ')[0].split('.')[-1]))

    df['optimizer'] = optimizer
    df_heatmap1 = pd.DataFrame(data = df[['activation', 'optimizer', 'kfold_test']])
    df_heatmap1= df_heatmap1.pivot(index='activation', columns='optimizer', values='kfold_test').transpose()
    #df_heatmap1 = df_heatmap1[['elu', 'linear', 'relu']]

    fig, ax = plt.subplots()
    im = ax.imshow(df_heatmap1.values , cmap = 'summer_r')

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(df_heatmap1.columns)), labels=df_heatmap1.columns)
    ax.set_yticks(np.arange(len(df_heatmap1.index)), labels=df_heatmap1.index)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(df_heatmap1.index)):
        for j in range(len(df_heatmap1.columns)):
            text = ax.text(j, i, np.around(df_heatmap1.values, decimals=2)[i, j],
                        ha="center", va="center", color="black")


    cbar_kw = {}        
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel('MAPE', rotation=-90, va="bottom")        


    #ax.set_title("Activation function and optimizer for best ARMAConv architechture")
    fig.tight_layout()
    plt.show()
    plt.savefig("ARMA_heatmap_activation_optimizer_run2.png")


def grid_search_2014(df_2024):
    df = df_2024
    layers = []
    for i in df['layer']:
        layers.append(re.sub('[^a-zA-Z]+', '', i.split('.')[-1]))

    df['layer'] = layers
    df['kfold_test'] = df['kfold_test'].astype('float')

    df_1 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GCNConv')][['number_of_layers', 'number_of_neurons', 'kfold_test']])
    df_1= df_1.pivot(index='number_of_layers', columns='number_of_neurons', values='kfold_test')
    
    df_2 = pd.DataFrame(data =  df.loc[(df['layer'] == 'ARMAConv')][['number_of_layers', 'number_of_neurons', 'kfold_test']])
    df_2= df_2.pivot(index='number_of_layers', columns='number_of_neurons', values='kfold_test')
    
    df_3 = pd.DataFrame(data =  df.loc[(df['layer'] == 'APPNPConv')][['number_of_layers', 'number_of_neurons', 'kfold_test']])
    df_3= df_3.pivot(index='number_of_layers', columns='number_of_neurons', values='kfold_test')
    
    df_4 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GATConv')][['number_of_layers', 'number_of_neurons', 'kfold_test']])
    df_4= df_4.pivot(index='number_of_layers', columns='number_of_neurons', values='kfold_test')


    plt.figure().clf() 
    layers = df['layer'].unique()


    fig, axs = plt.subplots(nrows = 2, ncols = 2, sharex=False, figsize = (7,7), constrained_layout = True)

    for ax, df_heatmap, name  in zip(axs.flat, [df_1 ,df_2, df_3, df_4], layers):
        
        im = ax.imshow(df_heatmap.values, cmap ='summer_r' )

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(df_heatmap.columns)), labels=df_heatmap.columns)
        ax.set_yticks(np.arange(len(df_heatmap.index)), labels=df_heatmap.index)

        ax.set_ylabel('number of layers')
        ax.set_xlabel('channels')

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(df_heatmap.index)):
            for j in range(len(df_heatmap.columns)):
                text = ax.text(j, i, np.around(df_heatmap.values, decimals=2)[i, j],
                            ha="center", va="center", color="k")
                
        ax.set_title('Model = '+name) 



    fig.colorbar(im, ax = axs, shrink=0.75, label = 'MAPE')
    
    #plt.rcParams.update({'font.size': 25})
    plt.show()
    plt.savefig("GNN1_layers_channels_2024.png")

def heatmap_ARMAConv():
    df_ARMA_1 = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_hyperparameters_2024/dataframe_ARMAConv_30days_oldrange_noCV.csv')
    df_ARMA_2 = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_hyperparameters_2024/dataframe_ARMAConv_30days_oldrange_noCV_part2.csv')
    df = pd.concat([df_ARMA_1, df_ARMA_2])

    df_plot = pd.DataFrame(data =  df[['number_of_layers', 'number_of_neurons', 'kfold_test']])
    df_plot = df_plot.pivot(index='number_of_layers', columns='number_of_neurons', values='kfold_test')
    
    fig, ax = plt.subplots()
    im = ax.imshow(df_plot.values, cmap ='summer_r')

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(df_plot.columns)), labels=df_plot.columns)
    ax.set_yticks(np.arange(len(df_plot.index)), labels=df_plot.index)
    ax.set_ylabel('number of layers')
    ax.set_xlabel('channels')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(df_plot.index)):
        for j in range(len(df_plot.columns)):
            text = ax.text(j, i, np.around(df_plot.values, decimals=3)[i, j],
                        ha="center", va="center", color="black", fontsize = 14 )

    cbar_kw = {}        
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel('MAPE', rotation=-90, va="bottom")    

    ax.set_title("")
    fig.tight_layout()
    plt.show()
    plt.savefig("ARNMAConv_enhanced_gridsearch.png")


def plot_data_differentsplits():
    
        file = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/data_GNN_nodamp_400pop_30days_2024_oldrange/data_secir_age_groups.pickle', 'rb')
        data_secir = pickle.load(file)

        kf = KFold(n_splits=5)
        train_idxs = []
        test_idxs = []
        for i, (train_index, test_index) in enumerate(kf.split(data_secir['labels'][0])):

            train_idxs.append(train_index)
            test_idxs.append(test_index)

        first_test_split = data_secir['labels'][0][test_idxs[0]]
        last_test_split = data_secir['labels'][0][test_idxs[4]]

        first_array = np.asarray(first_test_split.transpose().reshape(48,-1))
        all_comp_first = []
        index = [0,1,2,3,4,5,6,7]
        for i in index :
            comparray = first_array[i::8]
            all_comp_first.append(comparray.flatten())

        last_array = np.asarray(last_test_split.transpose().reshape(48,-1))
        all_comp_last = []
        index = [0,1,2,3,4,5,6,7]
        for i in index :
            comparray = last_array[i::8]
            all_comp_last.append(comparray.flatten())

        infectionstates = ['Susceptible','Exposed', 'InfectedNoSymptoms', 'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Receovered', 'Dead']

        #index = [str(compartment).split('.')[1] for compartment in  infectionstates ]
       
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(7, 9))
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

        #for ax, compartment, d, dw1, dw2 in zip(axes, index, np.asarray(sum_all).transpose(), np.asarray(sum_all_w).transpose(), np.asarray(sum_all_w2).transpose()):
        for ax, compartment, data in zip(axes, infectionstates, all_comp_last):
          
                ax.hist(data, bins = 30)
                ax.set_title(compartment, fontsize = 10)
                handles, labels = ax.get_legend_handles_labels()

        fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.55, 0.001), frameon=False)
        plt.tight_layout()
                 
        plt.savefig("GNN_hist_compartments_las.png")


def heatmap_ARMAConv_withCV():
    df_ARMA = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_hyperparameters_2024/dataframe_ARMAConv_30days_oldrange_withCV.csv')
  
    df_plot = pd.DataFrame(data =  df_ARMA[['number_of_layers', 'number_of_neurons', 'kfold_test']])
    df_plot = df_plot.pivot(index='number_of_layers', columns='number_of_neurons', values='kfold_test')
    
    fig, ax = plt.subplots()
    im = ax.imshow(df_plot.values, cmap ='summer_r')

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(df_plot.columns)), labels=df_plot.columns)
    ax.set_yticks(np.arange(len(df_plot.index)), labels=df_plot.index)
    ax.set_ylabel('number of layers')
    ax.set_xlabel('channels')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(df_plot.index)):
        for j in range(len(df_plot.columns)):
            text = ax.text(j, i, np.around(df_plot.values, decimals=3)[i, j],
                        ha="center", va="center", color="black", fontsize = 14 )

    cbar_kw = {}        
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel('MAPE', rotation=-90, va="bottom")    

    ax.set_title("")
    fig.tight_layout()
    plt.show()
    plt.savefig("ARNMAConv_withCV.png")


def scatterplot_ARMA():
    # Read dataframe from CSV
    df_ARMA = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_hyperparameters_2024/dataframe_ARMAConv_30days_oldrange_withCV.csv')

    # Extract relevant columns
    df_plot = pd.DataFrame(data=df_ARMA[['number_of_layers', 'number_of_neurons', 'all_testscores']])

    # Convert string representations of arrays to actual arrays
    df_plot['all_testscores'] = df_plot['all_testscores'].apply(lambda x: np.array(eval(x)))

    # Explode the dataframe based on the 'all_testscores' column
    df_expanded = df_plot.explode('all_testscores')
    df_expanded.reset_index(drop=True, inplace=True)

    # Create a new column 'summary' by concatenating 'number_of_layers' and 'number_of_neurons'
    df_expanded['summary'] = df_expanded.apply(lambda row: f"{row['number_of_layers']}, {row['number_of_neurons']}", axis=1)

    # Define positions
    array = ['1st split', '2nd split', '3d split', '4th split', '5th split']
    df_expanded['position'] = np.tile(array, 9)

    # Define your own custom colors
    custom_colors = ['red', 'blue', 'green', 'orange', 'purple']

    # Create a ListedColormap using the custom colors
    custom_cmap = ListedColormap(custom_colors)

    # Plot the scatter plot with color coding based on 'position'
    plt.figure(figsize=(6, 5))  # Adjust figure size as needed
    for i, position in enumerate(array):
        plt.scatter(df_expanded[df_expanded['position'] == position]['summary'], df_expanded[df_expanded['position'] == position]['all_testscores'], c=custom_colors[i], label=position, cmap=custom_cmap)

    plt.legend()
    plt.ylabel('MAPE')
    plt.xlabel('model architecture')  # Fixed x-label
    plt.xticks(rotation=45)
    plt.tight_layout()  # Adjust layout to prevent label overlap
    plt.savefig("ARMAConv_CV_testscores.png")


def ARMA_op_act_heatmap_cv_3_512():
    df_1 = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/machine-learning/dataframe_gridsearch_2024/ARMAConv_act_opt_2024_withCV.csv')        
    df_2 = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/machine-learning/dataframe_gridsearch_2024/ARMAConv_act_opt_2024_withCV_part2.csv')

    df = pd.concat([df_1, df_2], axis = 0)
    optimizer = []
    for i in df['optimizer']:
        optimizer.append(re.sub('[^a-zA-Z]+', '', i.split(' ')[0].split('.')[-1]))

    df['optimizer'] = optimizer

    df_heatmap1 = pd.DataFrame(data = df[['activation', 'optimizer', 'kfold_test']])
    df_heatmap1= df_heatmap1.pivot(index='activation', columns='optimizer', values='kfold_test').transpose()
    #df_heatmap1 = df_heatmap1[['elu', 'linear', 'relu']]

    fig, ax = plt.subplots()
    im = ax.imshow(df_heatmap1.values , cmap = 'summer_r')

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(df_heatmap1.columns)), labels=df_heatmap1.columns)
    ax.set_yticks(np.arange(len(df_heatmap1.index)), labels=df_heatmap1.index)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(df_heatmap1.index)):
        for j in range(len(df_heatmap1.columns)):
            text = ax.text(j, i, np.around(df_heatmap1.values, decimals=2)[i, j],
                        ha="center", va="center", color="black", fontsize = 15)


    cbar_kw = {}        
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel('MAPE', rotation=-90, va="bottom")        


    #ax.set_title("Activation function and optimizer for best ARMAConv architechture")
    fig.tight_layout()
    plt.show()
    plt.savefig("ARMA_heatmap_activation_optimizer_CV_3_512.png")


def ARMA_op_act_heatmap_nocv_4_512():
    df = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/machine-learning/dataframe_gridsearch_2024/ARMAConv_act_opt_2024_nocv_4_512.csv')        

    optimizer = []
    for i in df['optimizer']:
        optimizer.append(re.sub('[^a-zA-Z]+', '', i.split(' ')[0].split('.')[-1]))

    df['optimizer'] = optimizer

    df_heatmap1 = pd.DataFrame(data = df[['activation', 'optimizer', 'kfold_test']])
    df_heatmap1= df_heatmap1.pivot(index='activation', columns='optimizer', values='kfold_test').transpose()
    #df_heatmap1 = df_heatmap1[['elu', 'linear', 'relu']]

    fig, ax = plt.subplots()
    im = ax.imshow(df_heatmap1.values , cmap = 'summer_r')

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(df_heatmap1.columns)), labels=df_heatmap1.columns)
    ax.set_yticks(np.arange(len(df_heatmap1.index)), labels=df_heatmap1.index)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(df_heatmap1.index)):
        for j in range(len(df_heatmap1.columns)):
            text = ax.text(j, i, np.around(df_heatmap1.values, decimals=2)[i, j],
                        ha="center", va="center", color="black", fontsize = 15)


    cbar_kw = {}        
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel('MAPE', rotation=-90, va="bottom")        


    #ax.set_title("Activation function and optimizer for best ARMAConv architechture")
    fig.tight_layout()
    plt.show()
    plt.savefig("ARMA_heatmap_activation_optimizer_noCV_4_512.png")

def plot_ARMA_modifications_scatter():
    df = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_hyperparameters_2024/dataframe_ARMAConv_4_512_30days_oldrange_withCV_hyperparameters.csv')
        # Extract relevant columns
    df_plot = pd.DataFrame(data=df[['configuration', 'all_testscores']])

    # Convert string representations of arrays to actual arrays
    df_plot['all_testscores'] = df_plot['all_testscores'].apply(lambda x: np.array(eval(x)))
    df_plot = df_plot[:5]

    # Explode the dataframe based on the 'all_testscores' column
    df_expanded = df_plot.explode('all_testscores')
    df_expanded.reset_index(drop=True, inplace=True)

    # Define positions
    array = ['1st split', '2nd split', '3d split', '4th split', '5th split']
    df_expanded['position'] = np.tile(array, len(df_plot))

    # Define your own custom colors
    custom_colors = ['red', 'blue', 'green', 'orange', 'purple']

    # Create a ListedColormap using the custom colors
    custom_cmap = ListedColormap(custom_colors)

    # Plot the scatter plot with color coding based on 'position'
    plt.figure(figsize=(6, 5))  # Adjust figure size as needed
    for i, position in enumerate(array):
        plt.scatter(df_expanded[df_expanded['position'] == position]['configuration'], df_expanded[df_expanded['position'] == position]['all_testscores'], c=custom_colors[i], label=position, cmap=custom_cmap)


    plt.legend()
    plt.ylabel('MAPE')
    plt.xlabel('model architecture')  # Fixed x-label
    plt.xticks(rotation=45)
    plt.tight_layout()  # Adjust layout to prevent label overlap
    plt.savefig("ARMAConv_CV_testscores_configurations.png")


def plot_ARMA_modifications_barplot():


# Read the data from CSV
    df = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_hyperparameters_2024/dataframe_ARMAConv_4_512_30days_oldrange_withCV_hyperparameters.csv')

    # Convert kfold_test values to float
    df['kfold_test'] = df['kfold_test'].astype(float)

    # Select relevant columns and first 5 rows
    df_plot = df[['kfold_test', 'configuration']][:5]

    # Append a row with entries ['4.24', 'baseline']
    df_plot.loc[len(df_plot)] = ['4.24', 'baseline']

    # Extract values for the bars and convert to float
    values = df_plot['kfold_test'].astype(float)

    # Create the bar plot
    fig, ax = plt.subplots(layout='constrained')
    bars = ax.barh(df_plot['configuration'], values)

    # Add annotations to the bars
    for i, bar in enumerate(bars):
        ax.text(bar.get_width() - 0.6, bar.get_y() + bar.get_height()/2, f'{values[i]:.2f}',
                va='center', ha='left', fontsize=10, color='white', fontweight='bold')

    # Customize plot
    #ax.legend(['test MAPE'], loc='upper right', fontsize=8, bbox_to_anchor=(1.2, 1))
    ax.set_ylabel('Model Configuration')
    ax.set_xlabel('Test MAPE')

    # Save plot
    plt.savefig("ARMAConv_configurations_bar.png", bbox_inches='tight')
    plt.show()
