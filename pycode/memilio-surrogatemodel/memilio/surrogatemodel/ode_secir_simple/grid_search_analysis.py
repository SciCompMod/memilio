import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl
import re

df = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_simple_grid_searchdataframe')



# we do k forld cross validation in the grid search ( 5 folds) with early stop
# each fold can have a different number of epochs 
# we only plot the losses from the first fold
################# best model ###################
# get validation losses 
min_row = df.iloc[df['kfold_test'].idxmin()]

first_fold_minval = min_row.val_losses.split(']')[0]

min_loss_val  =  [i.replace('[', '') for i in first_fold_minval.split(',')]
min_loss_val  =  [i.replace(']', '') for i in min_loss_val]
min_loss_val = np.asarray(list(map(float, min_loss_val)))

# get training losses 
first_fold_mintrain = min_row.train_losses.split(']')[0]
min_loss_train  =  [i.replace('tf.Tensor: shape=(), dtype=float32, numpy=', '') for i in first_fold_mintrain.split('<')]
symbols = ['>', ',', '[', ']']
for s in symbols: 
    min_loss_train  =  [i.replace(s, '') for i in min_loss_train]
min_loss_train = np.asarray(list(map(float, min_loss_train[1:])))

################# worst model  ####################
max_row = df.iloc[df['kfold_test'].idxmax()]
first_fold_maxval = max_row.val_losses.split(']')[0]
max_loss_val  =  [i.replace('[', '') for i in first_fold_maxval.split(',')]
max_loss_val  =  [i.replace(']', '') for i in max_loss_val]
max_loss_val = np.asarray(list(map(float, max_loss_val)))

# get training losses 
first_fold_maxtrain = max_row.train_losses.split(']')[0]
max_loss_train  =  [i.replace('tf.Tensor: shape=(), dtype=float32, numpy=', '') for i in first_fold_maxtrain.split('<')]
symbols = ['>', ',', '[', ']']
for s in symbols: 
    max_loss_train  =  [i.replace(s, '') for i in max_loss_train]
max_loss_train = np.asarray(list(map(float, max_loss_train[1:])))

def plot_max_min_losses(min_loss_train, min_loss_val, max_loss_train, max_loss_val):
    plt.figure().clf() 
    plt.plot(min_loss_train, label = 'Train Best',  linestyle='--')
    plt.plot(min_loss_val, label = 'Validation Best ', linestyle='-')
    plt.plot(max_loss_train, label = 'Train Worst', linestyle='-.')
    plt.plot(max_loss_val, label = 'Validation Worst', linestyle=':')
    plt.legend()

    plt.title('Comparison of Losses for Best and Worst Model')
    plt.xlabel('Number of epochs')
    plt.ylabel('MAPE Loss')
    plt.savefig("min_max_losses_secir_simple.png")





def heatmap():
    df_heatmap1 = pd.DataFrame(data =  df.loc[(df['model'] == 'LSTM')][['number_of_hidden_layers', 'number_of_neurons', 'mean_test_MAPE']])
    df_heatmap1= df_heatmap1.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='mean_test_MAPE')

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


    ax.set_title("Number of layers and neurons for LSTM")
    fig.tight_layout()
    plt.show()
    plt.savefig("heatmap_layers_neurons_LSTM.png")







def heatmaps():
    df_heatmap1 = pd.DataFrame(data =  df.loc[(df['model'] == 'Dense')][['number_of_hidden_layers', 'number_of_neurons', 'mean_test_MAPE']])
    df_heatmap1= df_heatmap1.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='mean_test_MAPE')

    df_heatmap2 = pd.DataFrame(data =  df.loc[(df['model'] == 'CNN')][['number_of_hidden_layers', 'number_of_neurons', 'mean_test_MAPE']])
    df_heatmap2= df_heatmap2.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='mean_test_MAPE')

    df_heatmap3 = pd.DataFrame(data =  df.loc[(df['model'] == 'LSTM')][['number_of_hidden_layers', 'number_of_neurons', 'mean_test_MAPE']])
    df_heatmap3= df_heatmap3.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='mean_test_MAPE')

    fig, (ax0, ax1, ax2) = plt.subplots(nrows = 1, ncols = 3, sharex=True, figsize = (25,10))

    for ax, df_heatmap, name  in zip([ax0, ax1, ax2], [df_heatmap1 ,df_heatmap2, df_heatmap3], ['MLP', 'CNN', 'LSTM']):
        im = ax.imshow(df_heatmap.values)

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(df_heatmap.columns)), labels=df_heatmap.columns)
        ax.set_yticks(np.arange(len(df_heatmap.index)), labels=df_heatmap.index)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(df_heatmap.index)):
            for j in range(len(df_heatmap.columns)):
                text = ax.text(j, i, np.around(df_heatmap.values, decimals=2)[i, j],
                            ha="center", va="center", color="w")


        cbar_kw = {}        
        cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        cbar.ax.set_ylabel('MAPE', rotation=-90, va="bottom")        


        ax.set_title("Number of layers and neurons for"+name)
    fig.tight_layout()
    plt.show()
    plt.savefig("heatmap_layers_neurons_all.png")



    df_heatmap1 = pd.DataFrame(data =  df.loc[(df['model'] == 'Dense')][['number_of_hidden_layers', 'number_of_neurons', 'mean_test_MAPE']])
    df_heatmap1= df_heatmap1.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='mean_test_MAPE')

    df_heatmap2 = pd.DataFrame(data =  df.loc[(df['model'] == 'CNN')][['number_of_hidden_layers', 'number_of_neurons', 'mean_test_MAPE']])
    df_heatmap2= df_heatmap2.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='mean_test_MAPE')

    df_heatmap3 = pd.DataFrame(data =  df.loc[(df['model'] == 'LSTM')][['number_of_hidden_layers', 'number_of_neurons', 'mean_test_MAPE']])
    df_heatmap3= df_heatmap3.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='mean_test_MAPE')

    fig, axs = plt.subplots(nrows = 2, ncols = 2, sharex=True, figsize = (20,20))

    for ax, df_heatmap, name  in zip(axs.flat, [df_heatmap1 ,df_heatmap2, df_heatmap3], ['MLP', 'CNN', 'LSTM']):
        
        im = ax.imshow(df_heatmap.values)

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(df_heatmap.columns)), labels=df_heatmap.columns)
        ax.set_yticks(np.arange(len(df_heatmap.index)), labels=df_heatmap.index)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(df_heatmap.index)):
            for j in range(len(df_heatmap.columns)):
                text = ax.text(j, i, np.around(df_heatmap.values, decimals=2)[i, j],
                            ha="center", va="center", color="w")


        cbar_kw = {}        
        cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        cbar.ax.set_ylabel('MAPE', rotation=-90, va="bottom")        


        ax.set_title("Number of layers and neurons for "+name)
    fig.tight_layout()
    plt.rcParams.update({'font.size': 15})
    plt.show()
    plt.savefig("heatmap_layers_neurons_all.png")






def plot_layer():
    plt.figure().clf() 
    fig, axs = plt.subplots(nrows = 5, ncols = 1, sharex=True, layout='constrained', figsize = (15, 20))


    for neurons, ax in zip(df['number_of_neurons'].unique(),axs.flat):
        df_layer_plot = pd.DataFrame(data = df.loc[(df['number_of_neurons'] == neurons)][['number_of_hidden_layers','model',  'mean_test_MAPE']])
        df_layer_plot= df_layer_plot.pivot(index='number_of_hidden_layers', columns='model', values='mean_test_MAPE')
    
        linestyles = ['--', '-', ':']
        markers = []
        for m in Line2D.markers:
                try:
                    if len(m) == 1 and m != " ":
                        markers.append(m)
                except TypeError:
                    pass



        for i, ls, m in zip(df_layer_plot.columns, linestyles, markers): 
            ax.plot(df_layer_plot[i], label = i,  linestyle=ls, marker = m )
        ax.set_title('Influence of number of hidden layers with '+str(neurons)+ ' neurons')
        ax.legend()
        ax.set_xlabel('Number of layers')
        ax.set_ylabel('MAPE loss')
        


        #plt.title('Influence of number of hidden layers with '+str(neurons)+ ' neurons')
    #plt.legend()
    #plt.xlabel('Number of layers')
    #plt.ylabel('MAPE loss')
    plt.rcParams.update({'font.size': 10})
    plt.savefig("layers_simple.png")


def plot_neurons():
    plt.figure().clf() 
    fig, axs = plt.subplots(nrows = 5, ncols = 1, sharex=True, layout='constrained', figsize = (15, 30))


    for layers, ax in zip(df['number_of_hidden_layers'].unique(),axs.flat):
        df_layer_plot = pd.DataFrame(data = df.loc[(df['number_of_hidden_layers'] == layers)][['number_of_neurons','model',  'mean_test_MAPE']])
        df_layer_plot= df_layer_plot.pivot(index='number_of_neurons', columns='model', values='mean_test_MAPE')
    
        linestyles = ['--', '-', ':']
        markers = []
        for m in Line2D.markers:
                try:
                    if len(m) == 1 and m != " ":
                        markers.append(m)
                except TypeError:
                    pass



        for i, ls, m in zip(df_layer_plot.columns, linestyles, markers): 
            ax.plot(df_layer_plot[i], label = i,  linestyle=ls, marker = m )
        ax.set_title(str(layers)+ ' hidden layers')
        ax.legend()
        ax.set_xlabel('Number of layers')
        ax.set_ylabel('MAPE loss')
        

    plt.rcParams.update({'font.size': 8})
    plt.savefig("neurons_simple.png")



plot_max_min_losses(min_loss_train, min_loss_val, max_loss_train, max_loss_val)
heatmap()
plot_layer()

for i in df.columns[1:5]:
    print('Mean Test MAPE by ' + i ,':',  pd.DataFrame(data = df[[i, 'kfold_test']]).groupby([i]).mean())
