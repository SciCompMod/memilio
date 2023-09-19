import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl
import re

df = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/machine-learning/model/GNN/backup_dataframe_hyperparameter_tuning_full.csv')


# saving the array of losses in the df was a little bit problematic, and I made some mistakes
# this is why I have to go through some preparation steps in order to get the values in form of an array

##### replace layer and optimizer description
layers = []
for i in df['layer']:
    layers.append(re.sub('[^a-zA-Z]+', '', i.split('.')[-1]))


optimizer = []
for i in df['optimizer']:
    optimizer.append(re.sub('[^a-zA-Z]+', '', i.split(' ')[0].split('.')[-1]))

df['layer'] = layers
df['optimizer'] = optimizer


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
    plt.savefig("min_max_losses_.png")





def heatmap():
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




def plot_layer():
    plt.figure().clf() 
    df_layer_plot = pd.DataFrame(data = df.loc[(df['activation'] == 'elu') & (df['optimizer']=='Adam')][['layer', 'number_of_layers', 'kfold_test']])
    df_layer_plot= df_layer_plot.pivot(index='number_of_layers', columns='layer', values='kfold_test')
    linestyles = ['--', '-', '.', ':']

    markers = []
    for m in Line2D.markers:
        try:
            if len(m) == 1 and m != " ":
                markers.append(m)
        except TypeError:
            pass

    for i, ls, m in zip(df_layer_plot.columns, linestyles, markers): 
        plt.plot(df_layer_plot[i], label = i,  linestyle=ls, marker = m )


    plt.title('Layers and number of layers')
    plt.legend()
    plt.xlabel('Number of layers')
    plt.ylabel('MAPE loss')
    plt.savefig("layers.png")


plot_max_min_losses(min_loss_train, min_loss_val, max_loss_train, max_loss_val)
heatmap()
plot_layer()

for i in df.columns[1:5]:
    print('Mean Test MAPE by ' + i ,':',  pd.DataFrame(data = df[[i, 'kfold_test']]).groupby([i]).mean())
