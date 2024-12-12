import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import ast
import pickle
import re
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib.colors import Normalize
from matplotlib.colorbar import Colorbar

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from matplotlib.colors import Normalize


def heatmap_gridsearch_results(df_gridsearch, dimension1, dimension2, dimension3, savename, vmin=None, vmax=None):
    """! Creates a heatmap for all four GNN layer types. 

    @param df_gridsearch
    @param dimension1 
    @param dimension2 
    @param dimension3 
    @param savename
    @param vmin
    @return vmax
    """
    df = df_gridsearch
    layers = []
    for i in df['layer']:
        layers.append(re.sub('[^a-zA-Z]+', '', i.split('.')[-1]))
    df['layer'] = layers

    # Create pivot tables for each layer type
    df_heatmap1 = pd.DataFrame(data=df.loc[(df['layer'] == 'ARMAConv')][[
                               dimension1, dimension2, dimension3]])
    df_heatmap1 = df_heatmap1.pivot(
        index=dimension1,  columns=dimension2, values=dimension3)

    df_heatmap2 = pd.DataFrame(data=df.loc[(df['layer'] == 'GCNConv')][[
        dimension1, dimension2, dimension3]])
    df_heatmap2 = df_heatmap2.pivot(
        index=dimension1, columns=dimension2, values=dimension3)

    df_heatmap3 = pd.DataFrame(data=df.loc[(df['layer'] == 'APPNPConv')][[
                               dimension1, dimension2, dimension3]])
    df_heatmap3 = df_heatmap3.pivot(
        index=dimension1, columns=dimension2, values=dimension3)

    df_heatmap4 = pd.DataFrame(data=df.loc[(df['layer'] == 'GATConv')][[
                               dimension1, dimension2, dimension3]])
    df_heatmap4 = df_heatmap4.pivot(
        index=dimension1, columns=dimension2, values=dimension3)

    # Combine all data to find global min and max
    all_values = np.concatenate([
        df_heatmap1.values.flatten(),
        df_heatmap2.values.flatten(),
        df_heatmap3.values.flatten(),
        df_heatmap4.values.flatten()
    ])

    # Set clipping range (if not provided, use data range)
    if vmin is None:
        vmin = np.nanmin(all_values)
    if vmax is None:
        vmax = np.nanmax(all_values)

    # Clip data to the specified range for color normalization
    norm = Normalize(vmin=vmin, vmax=vmax)

    # Create subplots
    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=False,
                            figsize=(20, 20), constrained_layout=True)

    # Generate heatmaps
    for ax, df_heatmap, name in zip(
        axs.flat, [df_heatmap1, df_heatmap2, df_heatmap3, df_heatmap4],
        ['ARMAConv', 'GCNConv', 'GATConv', 'APPNPConv']
    ):
        # Clip values for color normalization
        clipped_data = np.clip(df_heatmap.values, vmin, vmax)

        # Create heatmap with clipped data for colormap
        im = ax.imshow(clipped_data, cmap='RdYlGn_r', norm=norm)

        plt.rcParams.update({'font.size': 30})
        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(df_heatmap.columns)),
                      labels=df_heatmap.columns, fontsize=25)
        ax.set_yticks(np.arange(len(df_heatmap.index)),
                      labels=df_heatmap.index, fontsize=25)

        ax.set_ylabel('Number of Hidden Layers', fontsize=25)
        ax.set_xlabel('Number of Neurons per Layer', fontsize=25)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(df_heatmap.index)):
            for j in range(len(df_heatmap.columns)):
                # Get the real value (unclipped) for annotation
                real_value = df_heatmap.values[i, j]

                # Annotate with the real value
                text = ax.text(j, i, np.around(real_value, decimals=2),
                               ha="center", va="center", color="k", fontsize=25)

        ax.set_title('Model = '+name, fontsize=30)

    # Add a single colorbar for all heatmaps
    cbar = fig.colorbar(im, ax=axs, location='right',
                        shrink=0.75, label='Validation MAPE')

    # Save and show the figure
    plt.show()
    plt.savefig(savename, bbox_inches='tight')


def ARMA_days_scatter(df):
    layers = []
    for i in df['layer']:
        layers.append(re.sub('[^a-zA-Z]+', '', i.split('.')[-1]))
    df['layer'] = layers
    df_plot = df.loc[(df['layer'] == 'ARMAConv')]

    df_plot = df[['number_of_layers',
                  'number_of_neurons', 'all_testscores']]

    df_plot['all_testscores'] = df_plot['all_testscores'].apply(
        lambda x: np.asarray(ast.literal_eval(x)))

    # Explode the dataframe based on the 'all_test_scores' column
    df_expanded = df_plot.explode('all_testscores')
    df_expanded.reset_index(drop=True, inplace=True)

    # Define positions
    array = ['1st split', '2nd split', '3rd split', '4th split', '5th split']
    df_expanded['position'] = np.tile(array, len(df_plot))

    # Create a new column by combining the first two columns
    df_expanded['summary'] = df_expanded['number_of_layers'].astype(
        str) + ',' + df_expanded['number_of_neurons'].astype(str)

    # Define your own custom colors
    custom_colors = ['red', 'blue', 'green', 'orange', 'purple']

    # Create a ListedColormap using the custom colors
    custom_cmap = ListedColormap(custom_colors)

    # Plot the scatter plot with color coding based on 'position'
    plt.figure(figsize=(20, 8))  # Adjust figure size as needed
    for i, position in enumerate(array):
        plt.scatter(df_expanded[df_expanded['position'] == position]['summary'], df_expanded[df_expanded['position']
                    == position]['all_testscores'], c=custom_colors[i], label=position, cmap=custom_cmap)

    plt.legend()
    plt.ylabel('MAPE')
    plt.xlabel('Number of Days')
    # Set the x-ticks and labels explicitly
    plt.xticks(ticks=df_expanded['summary'],
               labels=df_expanded['summary'], rotation=90)
    plt.tight_layout()  # Adjust layout to prevent label overlap
    plt.savefig("ARMAConv_CV_testscores_days_equalnodes.png")


def heatmap_ARMAConv(df_gridsearch, savename, dimension1, dimension2, dimension3, vmin=None, vmax=None):
    """! Creates a heatmap for all four GNN layer types. 

    @param df_gridsearch
    @param dimension1 
    @param dimension2 
    @param dimension3 
    @param savename
    @param vmin
    @return vmax
    """
    df = df_gridsearch

    # layers = []
    # for i in df['layer']:
    #     layers.append(re.sub('[^a-zA-Z]+', '', i.split('.')[-1]))
    # df['layer'] = layers
    # df = df.loc[df['layer'] == 'ARMAConv']

    # Create pivot tables for each layer type
    df_heatmap = pd.DataFrame(
        data=df[[dimension1, dimension2, dimension3]])
    df_heatmap = df_heatmap.pivot(
        index=dimension1,  columns=dimension2, values=dimension3)

    # Create subplots
    fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)

    # Create heatmap with clipped data for colormap
    im = ax.imshow(df_heatmap, cmap='RdYlGn_r')

    plt.rcParams.update({'font.size': 30})
    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(df_heatmap.columns)),
                  labels=df_heatmap.columns, fontsize=25)
    ax.set_yticks(np.arange(len(df_heatmap.index)),
                  labels=df_heatmap.index, fontsize=25)

    ax.set_ylabel('Number of parallel GCS stacks', fontsize=25)
    ax.set_xlabel('Number of iterations per stack ', fontsize=25)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(df_heatmap.index)):
        for j in range(len(df_heatmap.columns)):
            # Get the real value (unclipped) for annotation
            real_value = df_heatmap.values[i, j]

            # Annotate with the real value
            text = ax.text(j, i, np.around(real_value, decimals=2),
                           ha="center", va="center", color="k", fontsize=25)

    ax.set_title('Model = ARMAConv', fontsize=30)

    # Add a single colorbar for all heatmaps
    cbar = fig.colorbar(im, ax=ax, location='right',
                        shrink=0.75, label='Validation MAPE')

    # Save and show the figure
    plt.show()
    plt.savefig(savename, bbox_inches='tight')


# HEATMAP FOUR MODELS

path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'model_evaluations_paper')


# TEST MAPE
dimension1 = 'number_of_layers'
dimension2 = 'number_of_neurons'
dimension3 = 'kfold_test'
# dimension3 = 'training_time'
# equal nodes
filename = 'GNN_gridsearch_withCV_equalnodes.csv'  # filename 4x3 all models
df = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename)))
df_2 = pd.DataFrame(data=pd.read_csv(os.path.join(
    path_data, 'GNN_gridsearch_ARMA_samenodes_2.csv')))  # etxneded grid search ARMAConv
df_3 = pd.DataFrame(data=pd.read_csv(os.path.join(
    path_data, 'GNN_gridsearch_ARMA_equalnodes_fillNAN.csv')))  # fill NAN from ARMAConv
df_all = pd.concat([df, df_2, df_3])
savename = "GNN_gridsearch_equalnodes.png"

heatmap_gridsearch_results(
    df_all, dimension1, dimension2, dimension3, savename, vmin=None, vmax=70)

# nodes with variance
filename = 'GNN_gridsearch_withCV_nodeswithvariance.csv'  # 3x4 heatmaps
df = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename)))
df_2 = pd.DataFrame(data=pd.read_csv(os.path.join(
    path_data, 'GNN_gridsearch_ARMA_nodeswithvariance_2.csv')))  # extended grid search for ARMAConv
df_3 = pd.DataFrame(data=pd.read_csv(os.path.join(
    path_data, 'GNN_gridsearch_ARMA_nodeswithvariance_part3.csv')))
df_all = pd.concat([df, df_2, df_3])
savename = "GNN_gridsearch_nodeswithvariance.png"

heatmap_gridsearch_results(
    df_all, dimension1, dimension2, dimension3, savename, vmin=None, vmax=70)

# HEATMAP FOR ARMAConv HYPERPARAMETERS
dimension1 = 'order'
dimension2 = 'iterations'
dimension3 = 'kfold_test'

# equal nodes

filename1 = 'GNN_gridsearch_ARMAConv_hyperparameters_4thrun.csv'
filename2 = 'GNN_gridsearch_ARMAConv_hyperparameters_transposed.csv'
filename3 = 'GNN_gridsearch_ARMAConv_hyperparameters_5thrun.csv'
filename4 = 'GNN_gridsearch_ARMAConv_hyperparameters_transposed_part2.csv'

df1 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename1)))
df2 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename2)))
df3 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename3)))
df4 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename4)))
# df4 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename4)))
df_all = pd.concat([df1, df2, df3, df4])


savename = 'ARMAConv_hyper_equalnodes.png'

heatmap_ARMAConv(df_all, savename, dimension1, dimension2, dimension3,
                 vmin=None, vmax=None)


# nodes with variance

filename1 = 'GNN_gridsearch_ARMAConv_hyperparameters_nodeswithvariance.csv'
filename2 = 'GNN_gridsearch_ARMAConv_hyperparameters_transposed.csv'
filename3 = 'GNN_gridsearch_ARMAConv_hyperparameters_transposed_part2.csv'

df1 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename1)))
df2 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename2)))
df3 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename3)))
savename = 'ARMAConv_hyper_nodeswithvariance.png'
df_all = pd.concat([df1, df2, df_3])
heatmap_ARMAConv(df_all, dimension1, dimension2, dimension3, savename,
                 vmin=None, vmax=None)
