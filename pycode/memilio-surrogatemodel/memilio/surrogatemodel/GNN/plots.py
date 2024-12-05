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


def heatmap_gridsearch_results(df_gridsearch, savename, vmin=None, vmax=None):
    df = df_gridsearch
    layers = []
    for i in df['layer']:
        layers.append(re.sub('[^a-zA-Z]+', '', i.split('.')[-1]))
    df['layer'] = layers

    # Create pivot tables for each layer type
    df_heatmap1 = pd.DataFrame(data=df.loc[(df['layer'] == 'ARMAConv')][[
                               'number_of_layers', 'number_of_neurons', 'kfold_test']])
    df_heatmap1 = df_heatmap1.pivot(
        index='number_of_layers', columns='number_of_neurons', values='kfold_test')

    df_heatmap2 = pd.DataFrame(data=df.loc[(df['layer'] == 'GCNConv')][[
        'number_of_layers', 'number_of_neurons', 'kfold_test']])
    df_heatmap2 = df_heatmap2.pivot(
        index='number_of_layers', columns='number_of_neurons', values='kfold_test')

    df_heatmap3 = pd.DataFrame(data=df.loc[(df['layer'] == 'APPNPConv')][[
                               'number_of_layers', 'number_of_neurons', 'kfold_test']])
    df_heatmap3 = df_heatmap3.pivot(
        index='number_of_layers', columns='number_of_neurons', values='kfold_test')

    df_heatmap4 = pd.DataFrame(data=df.loc[(df['layer'] == 'GATConv')][[
                               'number_of_layers', 'number_of_neurons', 'kfold_test']])
    df_heatmap4 = df_heatmap4.pivot(
        index='number_of_layers', columns='number_of_neurons', values='kfold_test')

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


def ARMA_days_scatter(df_plot):
    df_plot = df[['number_of_layers',
                  'number_of_neurons', 'all_testscores']][:12]
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
    plt.figure(figsize=(6, 5))  # Adjust figure size as needed
    for i, position in enumerate(array):
        plt.scatter(df_expanded[df_expanded['position'] == position]['summary'], df_expanded[df_expanded['position']
                    == position]['all_testscores'], c=custom_colors[i], label=position, cmap=custom_cmap)

    plt.legend()
    plt.ylabel('MAPE')
    plt.xlabel('Number of Days')
    # Set the x-ticks and labels explicitly
    plt.xticks(ticks=df_plot['summary'], labels=df_plot['summary'])
    plt.tight_layout()  # Adjust layout to prevent label overlap
    plt.savefig("ARMAConv_CV_testscores_days_equalnodes.png")


# HEATMAP
filename = 'GNN_gridsearch_withCV_nodeswithvariance.csv'
# filename = 'GNN_gridsearch_withCV_equalnodes.csv'
# filename = 'gridserach_secir_groups_30days_I_based_Germany_10k_nodamp.csv'
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'model_evaluations_paper')
filepath = os.path.join(path_data, filename)
df = pd.DataFrame(data=pd.read_csv(filepath))
savename = "GNN_gridsearch_nodeswithvariance.png"
# savename = "GNN_gridsearch_equalnodes.png"
heatmap_gridsearch_results(df, savename, vmin=None, vmax=110)
