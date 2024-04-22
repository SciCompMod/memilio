import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import re




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