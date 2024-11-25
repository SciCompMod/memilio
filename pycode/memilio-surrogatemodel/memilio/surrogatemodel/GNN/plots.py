import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pickle
import re
import matplotlib.gridspec as gridspec


def heatmap_gridsearch_results(df_gridsearch, savename):
    df = df_gridsearch
    layers = []
    for i in df['layer']:
        layers.append(re.sub('[^a-zA-Z]+', '', i.split('.')[-1]))

    df['layer'] = layers
    plt.figure().clf()
    df_heatmap1 = pd.DataFrame(data=df.loc[(df['layer'] == 'ARMAConv')][[
                               'number_of_layers', 'number_of_neurons', 'kfold_val']])
    df_heatmap1 = df_heatmap1.pivot(
        index='number_of_layers', columns='number_of_neurons', values='kfold_val')

    # df_heatmap2 = pd.DataFrame(data=df.loc[(df['layer'] == 'GCNConv')][[
    #                            'number_of_layers', 'number_of_neurons', 'kfold_val']])
    # df_heatmap2 = df_heatmap2.pivot(
    #     index='number_of_layers', columns='number_of_neurons', values='kfold_val')

    # df_heatmap3 = pd.DataFrame(data=df.loc[(df['layer'] == 'GATConv')][[
    #                            'number_of_layers', 'number_of_neurons', 'kfold_val']])
    # df_heatmap3 = df_heatmap3.pivot(
    #     index='number_of_layers', columns='number_of_neurons', values='kfold_val')

    # df_heatmap4 = pd.DataFrame(data=df.loc[(df['layer'] == 'APPNPConv')][[
    #                            'number_of_layers', 'number_of_neurons', 'kfold_val']])
    # df_heatmap4 = df_heatmap3.pivot(
    #     index='number_of_layers', columns='number_of_neurons', values='kfold_val')

    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=False,
                            figsize=(20, 20), constrained_layout=True)

    # for ax, df_heatmap, name in zip(axs.flat, [df_heatmap1, df_heatmap2, df_heatmap3, df_heatmap4], ['ARMAConv', 'GCNConv', 'GATConv', 'APPNPConv']):
    for ax, df_heatmap, name in zip(axs.flat, [df_heatmap1], ['ARMAConv']):
        im = ax.imshow(df_heatmap.values, cmap='RdYlGn_r')
        plt.rcParams.update({'font.size': 30})
        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(df_heatmap.columns)),
                      labels=df_heatmap.columns, fontsize=25)
        ax.set_yticks(np.arange(len(df_heatmap.index)),
                      labels=df_heatmap.index, fontsize=25)

        ax.set_ylabel('number of hidden layers', fontsize=25)
        ax.set_xlabel('number of neurons per layer', fontsize=25)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(df_heatmap.index)):
            for j in range(len(df_heatmap.columns)):
                text = ax.text(j, i, np.around(df_heatmap.values, decimals=2)[i, j],
                               ha="center", va="center", color="k", fontsize=25)

        ax.set_title('Model = '+name, fontsize=30)

    fig.colorbar(im, ax=axs, shrink=0.75, label='Validation MAPE')
    # fig.delaxes(axs[1][1])
    plt.show()
    plt.savefig(savename)


# HEATMAP
filename = 'GNN_gridsearch_withCV_equalnodes.csv'
# filename = 'gridserach_secir_groups_30days_I_based_Germany_10k_nodamp.csv'
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'model_evaluations_paper')
filepath = os.path.join(path_data, filename)
df = pd.DataFrame(data=pd.read_csv(filepath))
savename = "GNN_gridsearch_nodeswithvariance.png"
