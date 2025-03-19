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
from matplotlib.cm import ScalarMappable
import matplotlib.ticker as ticker


def heatmap_gridsearch_results(df_gridsearch, dimension1, dimension2, dimension3, savename, vmin=None, vmax=None, size_tikz=40, size_labels=40, size_title=50, size_legend=20):
    """!
    Creates heatmaps for four GNN layer types without an integrated colorbar, and then
    saves separate figures containing a vertical and a horizontal colorbar.

    @param df_gridsearch: DataFrame with grid search results
    @param dimension1: DataFrame column to use for the heatmap's rows
    @param dimension2: DataFrame column to use for the heatmap's columns
    @param dimension3: DataFrame column to use for the heatmap's cell values
    @param savename: Filename to save the main heatmap figure (e.g. "heatmap.png").
                     Two additional files will be created by appending suffixes.
    @param vmin: Minimum value for colormap scaling (if None, computed from data)
    @param vmax: Maximum value for colormap scaling (if None, computed from data)
    @param size_tikz: Font size for tick labels.
    @param size_labels: Font size for axis labels.
    @param size_title: Font size for titles.
    @param size_legend: Font size for legend.
    """

    # Copy the DataFrame so the original is not modified
    df = df_gridsearch.copy()

    # Clean up the layer column (e.g. "models.layers.ARMAConv" -> "ARMAConv")
    layers = []
    for i in df['layer']:
        layers.append(re.sub('[^a-zA-Z]+', '', i.split('.')[-1]))
    df['layer'] = layers

    # --- Create pivot tables for each layer type ---
    # ARMAConv: only keep specific rows and columns
    df_heatmap1 = pd.DataFrame(data=df.loc[(df['layer'] == 'ARMAConv')][[
                               dimension1, dimension2, dimension3]])
    df_heatmap1 = df_heatmap1.pivot(
        index=dimension1, columns=dimension2, values=dimension3)
    df_heatmap1 = df_heatmap1.loc[[1, 2, 3], [32, 128, 512, 1024]]

    # GCNConv
    df_heatmap2 = pd.DataFrame(data=df.loc[(df['layer'] == 'GCNConv')][[
                               dimension1, dimension2, dimension3]])
    df_heatmap2 = df_heatmap2.pivot(
        index=dimension1, columns=dimension2, values=dimension3)

    # APPNPConv
    df_heatmap3 = pd.DataFrame(data=df.loc[(df['layer'] == 'APPNPConv')][[
                               dimension1, dimension2, dimension3]])
    df_heatmap3 = df_heatmap3.pivot(
        index=dimension1, columns=dimension2, values=dimension3)

    # GATConv
    df_heatmap4 = pd.DataFrame(data=df.loc[(df['layer'] == 'GATConv')][[
                               dimension1, dimension2, dimension3]])
    df_heatmap4 = df_heatmap4.pivot(
        index=dimension1, columns=dimension2, values=dimension3)

    # --- Determine global vmin and vmax across all heatmaps ---
    all_values = np.concatenate([
        df_heatmap1.values.flatten(),
        df_heatmap2.values.flatten(),
        df_heatmap3.values.flatten(),
        df_heatmap4.values.flatten()
    ])
    if vmin is None:
        vmin = np.nanmin(all_values)
    if vmax is None:
        vmax = np.nanmax(all_values)
    norm = Normalize(vmin=vmin, vmax=vmax)

    # --- Create heatmap subplots (without a colorbar) ---
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(
        20, 20), constrained_layout=True)

    # Note: To match the intended labels, we re-order the data so that:
    # ARMAConv -> df_heatmap1, GCNConv -> df_heatmap2,
    # GATConv -> df_heatmap4, APPNPConv -> df_heatmap3.
    for ax, df_heatmap, name in zip(
        axs.flat,
        [df_heatmap1, df_heatmap2, df_heatmap4, df_heatmap3],
        ['ARMAConv', 'GCNConv', 'GATConv', 'APPNPConv']
    ):
        # Clip the data to ensure that the colormap spans vmin to vmax
        clipped_data = np.clip(df_heatmap.values, vmin, vmax)
        im = ax.imshow(clipped_data, cmap='RdYlGn_r', norm=norm)

        # Set tick labels
        ax.set_xticks(np.arange(len(df_heatmap.columns)))
        ax.set_xticklabels(df_heatmap.columns, fontsize=size_tikz)
        ax.set_yticks(np.arange(len(df_heatmap.index)))
        ax.set_yticklabels(df_heatmap.index, fontsize=size_tikz)
        ax.set_xlabel('Number of Channels per Layer', fontsize=size_labels)
        ax.set_ylabel('Number of Hidden Layers', fontsize=size_labels)
        plt.setp(ax.get_xticklabels(), rotation=45,
                 ha="right", rotation_mode="anchor")

        # Annotate each cell with its numeric value
        for i in range(len(df_heatmap.index)):
            for j in range(len(df_heatmap.columns)):
                real_value = df_heatmap.values[i, j]
                ax.text(j, i, np.around(real_value, decimals=2),
                        ha="center", va="center", color="k", fontsize=size_tikz-10)
        ax.set_title('Model = ' + name, fontsize=size_title)

    # Save the main heatmap figure (without any colorbar)
    save_dir = '/localdata1/zunk_he/memilio/new_heatmaps/'
    plt.savefig(save_dir+savename, bbox_inches='tight', dpi=300)
    plt.close(fig)

    # --- Create a ScalarMappable for the colorbar using the same colormap and norm ---
    sm = ScalarMappable(norm=norm, cmap='RdYlGn_r')
    sm.set_array([])

    # --- Save a vertical colorbar in its own file ---
    vertical_fig, vertical_ax = plt.subplots(figsize=(1, 10))
    cbar = vertical_fig.colorbar(sm, cax=vertical_ax,
                                 orientation='vertical', label='Validation MAPE')
    cbar.set_label('Validation MAPE', size=size_labels)
    cbar.ax.tick_params(labelsize=size_tikz)
    for label in cbar.ax.get_yticklabels():
        label.set_fontsize(size_tikz)
    base, ext = os.path.splitext(savename)
    vertical_savename = base + '_vertical_colorbar' + ext
    plt.savefig(save_dir+vertical_savename, bbox_inches='tight', dpi=300)
    plt.close(vertical_fig)

    # --- Save a horizontal colorbar in its own file ---
    horizontal_fig, horizontal_ax = plt.subplots(figsize=(10, 1))
    cbar = horizontal_fig.colorbar(sm, cax=horizontal_ax,
                                   orientation='horizontal', label='Validation MAPE')
    cbar.set_label('Validation MAPE', size=size_labels)
    cbar.ax.tick_params(labelsize=size_tikz)
    for label in cbar.ax.get_xticklabels():
        label.set_fontsize(size_tikz)
    horizontal_savename = base + '_horizontal_colorbar' + ext
    plt.savefig(save_dir+horizontal_savename, bbox_inches='tight', dpi=300)
    plt.close(horizontal_fig)


def heatmap_gridsearch_ARMAConv_results(
    df_gridsearch, dimension1, dimension2, dimension3, savename,
    vmin=None, vmax=None, size_labels=40, size_ticks=40, size_text=32, size_title=50, time=False
):
    """
    Creates a heatmap for ARMAConv grid search results and saves the figure.

    The function filters the DataFrame for rows corresponding to the ARMAConv layer,
    pivots the data using the specified dimensions, determines global vmin/vmax if not provided,
    and then plots and annotates the heatmap.

    Parameters:
        df_gridsearch (pd.DataFrame): DataFrame with grid search results, which must include a 'layer' column.
        dimension1 (str): DataFrame column to use for the heatmap's rows (e.g., 'Number of Hidden Layers').
        dimension2 (str): DataFrame column to use for the heatmap's columns (e.g., 'Number of Channels per Layer').
        dimension3 (str): DataFrame column to use for the heatmap's cell values (e.g., 'Validation MAPE').
        savename (str): Filename (with extension) to save the heatmap figure.
        vmin (float, optional): Minimum value for colormap scaling. If None, it is computed from the data.
        vmax (float, optional): Maximum value for colormap scaling. If None, it is computed from the data.
        size_labels (int, optional): Font size for axis labels.
        size_ticks (int, optional): Font size for tick labels.
        size_text (int, optional): Font size for annotations inside heatmap cells.
        size_title (int, optional): Font size for subplot title.
        time (bool, optional): If True, the heatmap is for training time instead of MAPE.
    """

    # Copy and clean the DataFrame (strip layer strings like "models.layers.ARMAConv" -> "ARMAConv")
    df = df_gridsearch.copy()
    df['layer'] = df['layer'].apply(
        lambda x: re.sub('[^a-zA-Z]+', '', x.split('.')[-1]))

    # Define static heatmap values for ARMAConv
    if not time:
        df_heatmap = pd.DataFrame([
            [51.36, 51.30, 46.42, 46.74, 47.85],
            [23.08, 40.72, 28.71, 28.94, 47.28],
            [29.98, 28.65, 20.71, 39.50, 29.02],
            [30.47, 43.76, 27.56, 21.64, 35.71],
            [21.31, 37.22, 11.56, 30.96, 35.96],
            [21.98, 28.65,  8.98, 11.98, 29.47],
            [43.85, 35.91, 9.03, 30.66, 30.62]
        ], index=[1, 2, 3, 4, 5, 6, 7], columns=[32, 128, 512, 1024, 2048])
    else:
        df_heatmap = pd.DataFrame([
            [35.84, 38.58, 41.18, 32.01, 36.85],
            [50.57, 42.06, 28.16, 32.38, 47.25],
            [61.69, 30.02, 29.74, 61.21, 40.98],
            [37.73, 32.20, 36.38, 29.57, 61.68],
            [20.98, 44.03, 54.14, 39.28, 72.99],
            [81.15, 48.76, 58.38, 86.14, 138.81],
            [48.05, 46.29, 31.46, 82.58, 77.84]
        ], index=[1, 2, 3, 4, 5, 6, 7], columns=[32, 128, 512, 1024, 2048])

    # Determine global vmin/vmax if not provided
    all_values = df_heatmap.values.flatten()
    if vmin is None:
        vmin = np.nanmin(all_values)
    if vmax is None:
        vmax = np.nanmax(all_values)
    norm = Normalize(vmin=vmin, vmax=vmax)

    # Plot the heatmap
    fig, ax = plt.subplots(figsize=(15, 15), constrained_layout=True)

    # Clip data to ensure colormap spans exactly vmin to vmax
    clipped_data = np.clip(df_heatmap.values, vmin, vmax)
    im = ax.imshow(clipped_data, cmap='RdYlGn_r', norm=norm)

    # Set tick labels
    ax.set_xticks(np.arange(len(df_heatmap.columns)))
    ax.set_xticklabels(df_heatmap.columns, fontsize=size_ticks)
    ax.set_yticks(np.arange(len(df_heatmap.index)))
    ax.set_yticklabels(df_heatmap.index, fontsize=size_ticks)
    ax.set_xlabel('Number of Channels per Layer', fontsize=size_labels)
    ax.set_ylabel('Number of Hidden Layers', fontsize=size_labels)
    plt.setp(ax.get_xticklabels(), rotation=45,
             ha="right", rotation_mode="anchor")

    # Annotate each cell with its value
    for i in range(len(df_heatmap.index)):
        for j in range(len(df_heatmap.columns)):
            val = df_heatmap.values[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center",
                    va="center", color="k", fontsize=size_text)

    ax.set_title('Model = ARMAConv', fontsize=size_title)

    # Save and close the figure
    save_dir = '/localdata1/zunk_he/memilio/new_heatmaps/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    plt.savefig(os.path.join(save_dir, savename), bbox_inches='tight', dpi=300)
    plt.close(fig)


def heatmap_gridsearch_time_results(
    df_gridsearch, dimension1, dimension2, dimension3, savename,
    vmin=None, vmax=None,
    size_labels=40, size_ticks=40, size_title=50, size_text=32, size_colorbar=30
):
    """!
    Creates heatmaps for four GNN layer types without an integrated colorbar, and then
    saves separate figures containing a vertical and a horizontal colorbar.

    @param df_gridsearch: DataFrame with grid search results
    @param dimension1: DataFrame column to use for the heatmap's rows
    @param dimension2: DataFrame column to use for the heatmap's columns
    @param dimension3: DataFrame column to use for the heatmap's cell values
    @param savename: Filename to save the main heatmap figure (e.g. "heatmap.png").
                     Two additional files will be created by appending suffixes.
    @param vmin: Minimum value for colormap scaling (if None, computed from data)
    @param vmax: Maximum value for colormap scaling (if None, computed from data)
    @param size_labels: Font size for axis labels.
    @param size_ticks: Font size for tick labels.
    @param size_title: Font size for subplot titles.
    @param size_text: Font size for annotations inside heatmap cells.
    @param size_colorbar: Font size for colorbar labels and ticks.
    """

    # Copy the DataFrame so the original is not modified
    df = df_gridsearch.copy()

    # Clean up the layer column (e.g. "models.layers.ARMAConv" -> "ARMAConv")
    df['layer'] = df['layer'].apply(lambda x: ''.join(
        filter(str.isalpha, x.split('.')[-1])))

    # --- Create pivot tables for each layer type ---
    layer_names = ["ARMAConv", "GCNConv", "APPNPConv", "GATConv"]
    df_heatmaps = {layer: df[df["layer"] == layer].pivot(
        index=dimension1, columns=dimension2, values=dimension3) for layer in layer_names}

    # ARMAConv: Keep specific rows and columns
    if "ARMAConv" in df_heatmaps:
        df_heatmaps["ARMAConv"] = df_heatmaps["ARMAConv"].loc[[
            1, 2, 3], [32, 128, 512, 1024]]

    # --- Determine global vmin and vmax across all heatmaps ---
    all_values = np.concatenate([df_heatmap.values.flatten()
                                for df_heatmap in df_heatmaps.values()])
    if vmin is None:
        vmin = np.nanmin(all_values)
    if vmax is None:
        vmax = np.nanmax(all_values)
    norm = Normalize(vmin=vmin, vmax=vmax)

    # --- Create heatmap subplots (without a colorbar) ---
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(
        20, 20), constrained_layout=True)

    # Order layers for plotting
    layer_order = ["ARMAConv", "GCNConv", "GATConv", "APPNPConv"]

    for ax, layer in zip(axs.flat, layer_order):
        df_heatmap = df_heatmaps[layer]

        # Clip the data to ensure that the colormap spans vmin to vmax
        clipped_data = np.clip(df_heatmap.values, vmin, vmax)
        im = ax.imshow(clipped_data, cmap='RdYlGn_r', norm=norm)

        # Set tick labels
        ax.set_xticks(np.arange(len(df_heatmap.columns)))
        ax.set_xticklabels(df_heatmap.columns, fontsize=size_ticks)
        ax.set_yticks(np.arange(len(df_heatmap.index)))
        ax.set_yticklabels(df_heatmap.index, fontsize=size_ticks)
        ax.set_xlabel('Number of Channels per Layer', fontsize=size_labels)
        ax.set_ylabel('Number of Hidden Layers', fontsize=size_labels)
        plt.setp(ax.get_xticklabels(), rotation=45,
                 ha="right", rotation_mode="anchor")

        # Annotate each cell with its numeric value
        for i in range(len(df_heatmap.index)):
            for j in range(len(df_heatmap.columns)):
                ax.text(j, i, np.around(df_heatmap.values[i, j], decimals=2),
                        ha="center", va="center", color="k", fontsize=size_text)

        ax.set_title(f'Model = {layer}', fontsize=size_title)

    # Save the main heatmap figure (without any colorbar)
    save_dir = '/localdata1/zunk_he/memilio/new_heatmaps/'
    plt.savefig(os.path.join(save_dir, savename), bbox_inches='tight', dpi=300)
    plt.close(fig)

    # --- Create a ScalarMappable for the colorbar using the same colormap and norm ---
    sm = ScalarMappable(norm=norm, cmap='RdYlGn_r')
    sm.set_array([])

    # --- Save a vertical colorbar in its own file ---
    vertical_fig, vertical_ax = plt.subplots(figsize=(1, 10))
    label_name = 'Training Time per Fold (in minutes)'
    cbar = vertical_fig.colorbar(
        sm, cax=vertical_ax, orientation='vertical', label=label_name)
    cbar.set_label(label_name, size=size_colorbar)
    cbar.ax.tick_params(labelsize=size_colorbar)
    # Ensure all ticks are in the correct size
    for label in cbar.ax.get_yticklabels():
        label.set_fontsize(size_colorbar)
    # Reduce the number of ticks
    cbar.locator = ticker.MaxNLocator(nbins=5)
    cbar.update_ticks()
    base, ext = os.path.splitext(savename)
    vertical_savename = base + '_time_vertical_colorbar' + ext
    plt.savefig(os.path.join(save_dir, vertical_savename),
                bbox_inches='tight', dpi=300)
    plt.close(vertical_fig)

    # --- Save a horizontal colorbar in its own file ---
    horizontal_fig, horizontal_ax = plt.subplots(figsize=(10, 1))
    cbar = horizontal_fig.colorbar(
        sm, cax=horizontal_ax, orientation='horizontal', label=label_name)
    cbar.set_label(label_name, size=size_colorbar)
    cbar.ax.tick_params(labelsize=size_ticks)
    horizontal_savename = base + '_time_horizontal_colorbar' + ext
    plt.savefig(os.path.join(save_dir, horizontal_savename),
                bbox_inches='tight', dpi=300)
    plt.close(horizontal_fig)


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


def heatmap_ARMAConv_hyper(
    df_gridsearch, dimension1, dimension2, dimension3, savename,
    vmin=None, vmax=None,
    size_labels=36, size_ticks=40, size_text=35, size_title=50
):
    """!
    Creates a heatmap for ARMAConv hyperparameter tuning results.

    @param df_gridsearch: DataFrame with grid search results (not used as data is hardcoded).
    @param dimension1: Column name representing rows of the heatmap.
    @param dimension2: Column name representing columns of the heatmap.
    @param dimension3: Column name representing the cell values in the heatmap.
    @param savename: Filename to save the heatmap figure.
    @param vmin: Minimum value for colormap scaling (if None, computed from data).
    @param vmax: Maximum value for colormap scaling (if None, computed from data).
    @param size_labels: Font size for axis labels.
    @param size_ticks: Font size for tick labels.
    @param size_text: Font size for annotations inside heatmap cells.
    @param size_title: Font size for the subplot title.
    """

    # Hard-coded data dictionary
    data = {
        1: [9.03, 20.16, 26.07, 32.83],
        2: [27.13, 19.15, 36.6, 36.17],
        3: [43.74, 33.96, 43.43, 43.94],
        4: [43.56, 39.66, 44.18, 44.17]
    }

    # Convert dictionary to DataFrame
    df_wide = pd.DataFrame(data, index=[1, 2, 3, 4])
    df_wide.index.name = dimension1  # Row label
    df_wide.columns.name = dimension2  # Column label

    # Convert from wide to long format
    df_long = df_wide.reset_index().melt(
        id_vars=dimension1, var_name=dimension2, value_name=dimension3
    )

    # Pivot for heatmap structure
    df_heatmap = df_long.pivot(
        index=dimension1, columns=dimension2, values=dimension3)

    # Determine vmin and vmax if not provided
    all_values = df_heatmap.values.flatten()
    if vmin is None:
        vmin = np.nanmin(all_values)
    if vmax is None:
        vmax = np.nanmax(all_values)
    norm = Normalize(vmin=vmin, vmax=vmax)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)

    # Create heatmap
    im = ax.imshow(df_heatmap, cmap='RdYlGn_r', norm=norm)

    # Set tick labels
    ax.set_xticks(np.arange(len(df_heatmap.columns)))
    ax.set_xticklabels(df_heatmap.columns, fontsize=size_ticks)
    ax.set_yticks(np.arange(len(df_heatmap.index)))
    ax.set_yticklabels(df_heatmap.index, fontsize=size_ticks)

    ax.set_xlabel('Number of iterations per stack', fontsize=size_labels)
    ax.set_ylabel('Number of parallel GCS stacks', fontsize=size_labels)
    ax.set_title('Model = ARMAConv', fontsize=size_title)

    # Rotate the tick labels
    plt.setp(ax.get_xticklabels(), rotation=45,
             ha="right", rotation_mode="anchor")

    # Annotate each cell
    for i in range(len(df_heatmap.index)):
        for j in range(len(df_heatmap.columns)):
            ax.text(j, i, f"{df_heatmap.values[i, j]:.2f}",
                    ha="center", va="center", color="k", fontsize=size_text)

    # Save the figure
    save_dir = "/localdata1/zunk_he/memilio/new_heatmaps/"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    plt.savefig(os.path.join(save_dir, savename), bbox_inches='tight', dpi=300)
    plt.close(fig)


def heatmap_ARMAConv_time_hyper(df_gridsearch, dimension1, dimension2, dimension3, savename, vmin=None, vmax=None):
    """! Creates a heatmap for all four GNN layer types. 

    @param df_gridsearch
    @param dimension1 
    @param dimension2 
    @param dimension3 
    @param savename
    @param vmin
    @return vmax
    """
    # Hard-coded dictionary
    data = {
        32:   [35.84, 50.57, 61.69, 37.73, 20.98, 81.15, 48.05],
        128:  [38.58, 42.06, 30.02, 32.20, 44.03, 48.76, 46.29],
        512:  [41.18, 28.16, 29.74, 36.38, 54.14, 58.38, 31.46],
        1024: [32.01, 32.38, 61.21, 29.57, 39.28, 86.14, 82.58],
        2048: [36.85, 47.25, 40.98, 61.68, 72.99, 138.81, 77.84],
    }

    df_wide = pd.DataFrame(data, index=[1, 2, 3, 4, 5, 6, 7])
    df_wide.index.name = dimension1  # row label
    df_wide.columns.name = dimension2  # column label

    # Convert from wide to long so we have dimension1, dimension2, dimension3
    df_long = df_wide.reset_index().melt(
        id_vars=dimension1,
        var_name=dimension2,
        value_name=dimension3
    )

    # Pivot so that dimension1 is the row index, dimension2 is the column index,
    # and dimension3 is the cell value
    df_heatmap = df_long.pivot(
        index=dimension1,
        columns=dimension2,
        values=dimension3
    )

    # layers = []
    # for i in df['layer']:
    #     layers.append(re.sub('[^a-zA-Z]+', '', i.split('.')[-1]))
    # df['layer'] = layers
    # df = df.loc[df['layer'] == 'ARMAConv']

    # # Create pivot tables for each layer type
    # df_heatmap = pd.DataFrame(
    #     data=df[[dimension1, dimension2, dimension3]])
    # df_heatmap = df_heatmap.pivot(
    #     index=dimension1,  columns=dimension2, values=dimension3)

    # Create subplots
    fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)

    # Create heatmap with clipped data for colormap
    im = ax.imshow(df_heatmap, cmap='RdYlGn_r', vmin=vmin, vmax=vmax)

    plt.rcParams.update({'font.size': 30})
    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(df_heatmap.columns)),
                  labels=df_heatmap.columns, fontsize=25)
    ax.set_yticks(np.arange(len(df_heatmap.index)),
                  labels=df_heatmap.index, fontsize=25)

    ax.set_ylabel('Number of Layers', fontsize=25)
    ax.set_xlabel('Number of Channels per Layer ', fontsize=25)

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
    # cbar = fig.colorbar(im, ax=ax, location='right',
    #                     shrink=0.75, label='Validation MAPE')

    # Save and show the figure
    # plt.show()
    plt.savefig("/localdata1/zunk_he/memilio/new_heatmaps/" +
                savename, bbox_inches='tight')


def heatmap_ARMAConv(df_gridsearch, dimension1, dimension2, dimension3, savename, vmin=None, vmax=None):
    """! Creates a heatmap for all four GNN layer types. 

    @param df_gridsearch
    @param dimension1 
    @param dimension2 
    @param dimension3 
    @param savename
    @param vmin
    @return vmax
    """
    # Hard-coded dictionary
    data = {
        32:   [51.36, 23.08, 29.98, 30.47, 21.31, 21.98, 43.85],
        128:  [51.30, 40.72, 28.65, 43.76, 37.22, 28.65, 35.91],
        512:  [46.42, 28.71, 20.71, 27.56, 11.56,  8.98, 9.03],
        1024: [46.74, 28.94, 39.50, 21.64, 30.96, 11.98, 30.66],
        2048: [47.85, 47.28, 29.02, 35.71, 35.96, 29.47, 30.62]
    }

    # Convert dict to a "wide" DataFrame: rows = 1..7, columns = 32..2048
    # We'll label the index as dimension1, columns as dimension2
    df_wide = pd.DataFrame(data, index=range(1, 8))
    df_wide.index.name = dimension1
    df_wide.columns.name = dimension2

    # Convert from wide to long so we have (dimension1, dimension2, dimension3)
    df_long = df_wide.reset_index().melt(
        id_vars=dimension1,
        var_name=dimension2,
        value_name=dimension3
    )

    # Pivot so dimension1 becomes rows, dimension2 becomes columns
    df_heatmap = df_long.pivot(
        index=dimension1,
        columns=dimension2,
        values=dimension3
    )

    # Create a figure
    fig, ax = plt.subplots(figsize=(9, 9), constrained_layout=True)

    # Plot the heatmap using the passed-in vmin/vmax
    im = ax.imshow(df_heatmap, cmap='RdYlGn_r', vmin=vmin, vmax=vmax)

    # Increase font size
    plt.rcParams.update({'font.size': 14})

    # Show all ticks and label them
    ax.set_xticks(np.arange(len(df_heatmap.columns)),
                  labels=df_heatmap.columns)
    ax.set_yticks(np.arange(len(df_heatmap.index)),   labels=df_heatmap.index)

    ax.set_xlabel('Number of Layers', fontsize=16)
    ax.set_ylabel('Number of Channels per Layer', fontsize=16)
    ax.set_title('Model = ARMAConv', fontsize=18)

    # Rotate column labels
    plt.setp(ax.get_xticklabels(), rotation=45,
             ha="right", rotation_mode="anchor")

    # Annotate each cell
    for i in range(len(df_heatmap.index)):
        for j in range(len(df_heatmap.columns)):
            val = df_heatmap.values[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center",
                    va="center", color="k", fontsize=25)

    ax.set_title('Model = ARMAConv', fontsize=30)

    # Add a single colorbar for all heatmaps
    # cbar = fig.colorbar(im, ax=ax, location='right',
    #                     shrink=0.75, label='Validation MAPE')

    # Save and show the figure
    # plt.show()
    plt.savefig("/localdata1/zunk_he/memilio/new_heatmaps/" +
                savename, bbox_inches='tight')


# HEATMAP FOUR MODELS

path_data = '/localdata1/gnn_paper_2024/data/results/grid_search/GNN/'


# TEST MAPE
dimension1 = 'number_of_layers'
dimension2 = 'number_of_neurons'
dimension3 = 'kfold_test'
# dimension3 = 'training_time'
# equal nodes
# filename = 'GNN_gridsearch_withCV_equalnodes.csv'  # filename 4x3 all models
# df = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename)))
# df_2 = pd.DataFrame(data=pd.read_csv(os.path.join(
#     path_data, 'GNN_gridsearch_ARMA_samenodes_2.csv')))  # etxneded grid search ARMAConv
# df_3 = pd.DataFrame(data=pd.read_csv(os.path.join(
#     path_data, 'GNN_gridsearch_ARMA_equalnodes_fillNAN.csv')))  # fill NAN from ARMAConv
# df_all = pd.concat([df, df_2, df_3])
# savename = "GNN_gridsearch_equalnodes.png"

# heatmap_gridsearch_results(
#     df_all, dimension1, dimension2, dimension3, savename, vmin=None, vmax=70)

# nodes with variance
filename = 'GNN_gridsearch_withCV_nodeswithvariance.csv'  # 3x4 heatmaps
df = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename)))
df_2 = pd.DataFrame(data=pd.read_csv(os.path.join(
    path_data, 'GNN_gridsearch_ARMA_nodeswithvariance_2.csv')))  # extended grid search for ARMAConv
df_3 = pd.DataFrame(data=pd.read_csv(os.path.join(
    path_data, 'GNN_gridsearch_ARMA_nodeswithvariance_part3.csv')))
df_all = pd.concat([df, df_2, df_3])
savename = "GNN_gridsearch_nodeswithvariance.png"

# heatmap_gridsearch_results(
#     df_all, dimension1, dimension2, dimension3, savename, vmin=10, vmax=60)

dimension3 = 'training_time'
savename = "GNN_gridsearch_time_nodeswithvariance.png"
# heatmap_gridsearch_time_results(
#     df_all, dimension1, dimension2, dimension3, savename, vmin=20, vmax=65)

# heatmap for ARMAConv gridsearch
dimension3 = 'kfold_test'
savename = "GNN_gridsearch_nodeswithvariance_ARMAConv"
# heatmap_gridsearch_ARMAConv_results(
#     df_all, dimension1, dimension2, dimension3, savename, vmin=10, vmax=60)

# heatmap for ARMAConv time  gridsearch
dimension3 = 'training_time'
savename = "GNN_gridsearch_time_nodeswithvariance_ARMAConv"
# heatmap_gridsearch_ARMAConv_results(
#     df_all, dimension1, dimension2, dimension3, savename, vmin=20, vmax=65, time=True)

# HEATMAP FOR ARMAConv HYPERPARAMETERS
dimension1 = 'order'
dimension2 = 'iterations'
dimension3 = 'kfold_test'

# equal nodes

# filename1 = 'GNN_gridsearch_ARMAConv_hyperparameters_4thrun.csv'
# filename2 = 'GNN_gridsearch_ARMAConv_hyperparameters_transposed.csv'
# filename3 = 'GNN_gridsearch_ARMAConv_hyperparameters_5thrun.csv'
# filename4 = 'GNN_gridsearch_ARMAConv_hyperparameters_transposed_part2.csv'

# df1 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename1)))
# df2 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename2)))
# df3 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename3)))
# df4 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename4)))
# # df4 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename4)))
# df_all = pd.concat([df1, df2, df3, df4])


# savename = 'ARMAConv_hyper_equalnodes.png'

# heatmap_ARMAConv(df_all, savename, dimension1, dimension2, dimension3,
#                  vmin=None, vmax=None)


# nodes with variance

filename1 = 'GNN_gridsearch_ARMAConv_hyperparameters_nodeswithvariance.csv'
filename2 = 'GNN_gridsearch_ARMAConv_hyperparameters_transposed.csv'
filename3 = 'GNN_gridsearch_ARMAConv_hyperparameters_transposed_part2.csv'

df1 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename1)))
df2 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename2)))
df3 = pd.DataFrame(data=pd.read_csv(os.path.join(path_data, filename3)))
df_all = pd.concat([df1, df2, df3]).reset_index(drop=True)
savename = 'ARMAConv_hyper_nodeswithvariance.png'
heatmap_ARMAConv_hyper(df_all, dimension1, dimension2, dimension3, savename,
                       vmin=10, vmax=60)
