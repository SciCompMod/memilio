import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import re

df = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_hyperparameters/dataframe_frid_search_GNN_noedges_full_concat.csv')
df_GSC = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_gridsearch_2024/baseline_GSCConv.csv')
#df_new = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/dataframe_hyperparameters/dataframe_hyperparameter_tuning_all_concatenated_new.csv')
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
df.drop([80], inplace = True)

df['kfold_test'] = df['kfold_test'].astype('float')


df_ARMA =  df.loc[df['layer'] == 'ARMAConv']
df_ARMA = df_ARMA[['number_of_layers', 'activation', 'optimizer', 'kfold_test']].loc[df['number_of_layers'] == '2']


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


def barplot(df):
       

    df_1 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GCNConv')][['number_of_layers', 'kfold_test']])
    df_1 = df_1.groupby('number_of_layers').mean()
    
    df_2 = pd.DataFrame(data =  df.loc[(df['layer'] == 'ARMAConv')][['number_of_layers', 'kfold_test']])
    df_2 = df_2.groupby('number_of_layers').mean()
    
    df_3 = pd.DataFrame(data =  df.loc[(df['layer'] == 'APPNPConv')][['number_of_layers',  'kfold_test']])
    df_3 = df_3.groupby('number_of_layers').mean()
    
    df_4 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GATConv')][['number_of_layers', 'kfold_test']])
    df_4 = df_4.groupby('number_of_layers').mean()

    MAPE = {
        'GCNConv': np.squeeze(df_1.values).round(2).tolist(), 
        'ARMAConv':np.squeeze(df_2.values).round(2).tolist(), 
        'APPNPConv':np.squeeze(df_3.values).round(2).tolist(), 
        'GATConv': np.squeeze(df_4.values).round(2).tolist(), 
    }

    layers = df_1.index.values
    #df_bar=df_opt[['optimizer',  'kfold_test']]


    x = np.arange(1,len(layers)+1)  # the label locations
    width = 0.15  # the width of the bars
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

    ax.set_xticks(x, labels=x, fontsize = 12)
    


    ax.set_ylabel('test MAPE')
    ax.set_xlabel('Number of layers')
    ax.set_title('Mean Test MAPE for different number of layers')
    plt.savefig("bar_GNN1_layer.png")

def heatmap(df):
    #df_heatmap1 = pd.DataFrame(data = df.loc[(df['layer'] == 'GCNConv') & (df['number_of_layers']==3)][['activation', 'optimizer', 'kfold_test']])
    df['kfold_test'] = df['kfold_test'].astype(float)
    df_heatmap1= df.pivot(index='activation', columns='optimizer', values='kfold_test')

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


# def plot_layer():
#     plt.figure().clf() 
#     df_layer_plot = pd.DataFrame(data = df.loc[(df['activation'] == 'elu') & (df['optimizer']=='Adam')][['layer', 'number_of_layers', 'kfold_test']])
#     df_layer_plot= df_layer_plot.pivot(index='number_of_layers', columns='layer', values='kfold_test')
#     linestyles = ['--', '-', '.', ':']

#     markers = []
#     for m in Line2D.markers:
#         try:
#             if len(m) == 1 and m != " ":
#                 markers.append(m)
#         except TypeError:
#             pass

#     for i, ls, m in zip(df_layer_plot.columns, linestyles, markers): 
#         plt.plot(df_layer_plot[i], label = i,  linestyle=ls, marker = m )


#     plt.title('Layers and number of layers')
#     plt.legend()
#     plt.xlabel('Number of layers')
#     plt.ylabel('MAPE loss')
#     plt.savefig("layers.png")


def heatmaps(df):

    df_1 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GCNConv')& (df['number_of_layers']=='2')][['activation', 'optimizer', 'kfold_test']])
    df_1= df_1.pivot(index='activation', columns='optimizer', values='kfold_test')
    
    df_2 = pd.DataFrame(data =  df.loc[(df['layer'] == 'ARMAConv')& (df['number_of_layers']=='2')][['activation', 'optimizer', 'kfold_test']])
    df_2= df_2.pivot(index='activation', columns='optimizer', values='kfold_test')
    
    df_3 = pd.DataFrame(data =  df.loc[(df['layer'] == 'APPNPConv')& (df['number_of_layers']=='2')][['activation', 'optimizer', 'kfold_test']])
    df_3= df_3.pivot(index='activation', columns='optimizer', values='kfold_test')
    
    df_4 = pd.DataFrame(data =  df.loc[(df['layer'] == 'GATConv')& (df['number_of_layers']=='2')][['activation', 'optimizer', 'kfold_test']])
    df_4= df_4.pivot(index='activation', columns='optimizer', values='kfold_test')

    plt.figure().clf() 
    layers = df['layer'].unique()

    fig, axs = plt.subplots(nrows = 2, ncols = 2, sharex=False,figsize = (8,8),  constrained_layout = True)

    for ax, df_heatmap, name  in zip(axs.flat, [df_1 ,df_2, df_3, df_4], layers):
        
        im = ax.imshow(df_heatmap.values, cmap ='summer_r' )

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(df_heatmap.columns)), labels=df_heatmap.columns)
        ax.set_yticks(np.arange(len(df_heatmap.index)), labels=df_heatmap.index)

        ax.set_ylabel('activation')
        ax.set_xlabel('optimizer')

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
    
    #fig.delaxes(axs[1][1])
    
    #plt.rcParams.update({'font.size': 25})
    plt.title('GNN activation and optimizer for models with 2 layers')
    plt.show()
    plt.savefig("GNN_heatmap_activ_opt_2layer.png")



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

plot_max_min_losses(min_loss_train, min_loss_val, max_loss_train, max_loss_val)
heatmap()
#plot_layer()

for i in df.columns[1:5]:
    print('Mean Test MAPE by ' + i ,':',  pd.DataFrame(data = df[[i, 'kfold_test']]).groupby([i]).mean())
