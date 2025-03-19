import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pickle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from memilio.simulation.osecir import InfectionState
import seaborn as sns

filename = 'data_secir_simple_90days_10k.pickle'

# import data
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'data_paper')
path_data = '/localdata1/gnn_paper_2024/data/one_population/without_agegroups/'

if not os.path.isfile(os.path.join(path_data, filename)):
    ValueError("no dataset found in path: " + path)

file = open(os.path.join(path_data, filename), 'rb')
data = pickle.load(file)

data['inputs'] = np.expm1(np.asarray(data['inputs']))
data['labels'] = np.expm1(np.asarray(data['labels']))

# create lineplot for n data samples  with 8 subplots


def lineplots_compartments(inputs_reversed, labels_reversed, num_plots):
    plt.clf()
    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    for i in range(num_plots):
        plt.clf()
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
            nrows=4, ncols=2, sharey=False, figsize=(10, 13), constrained_layout=True)

        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

        for ax, c, inp, lab in zip(axes, infectionstates, inputs_reversed[i].transpose(), labels_reversed[i].transpose()):
            # Plot the first 5 days from inp (in blue)
            ax.plot(np.arange(1, 6),
                    inp[:5], color='blue', label='inputs (first 5 days)')

            # Plot the remaining days from lab (in orange) starting from day 6
            ax.plot(np.arange(6, 96), lab, color='orange',
                    label='labels (days 6-95)')

            ax.set_xlabel('Number of days')
            ax.set_ylabel('Number of individuals')
            ax.set_title(c, fontsize=10)

        ax7.set_xlabel('Time')
        ax8.set_xlabel('Time')

        lines = []
        line_labels = []
        for ax in fig.axes:
            Line, Label = ax.get_legend_handles_labels()
            lines.extend(Line)
            line_labels.extend(Label)

        fig.legend(lines[:2], line_labels[:2], loc='lower center')
        fig.suptitle(
            'Predicted values and labels for compartments', fontsize=16)

        print('Plot No.'+str(i) + ' done.')
        plt.savefig("/localdata1/gnn_paper_2024/images/no_agegroups_90days_I_based_8plt/compartment_lines_noagegroups_90days_I_based_8plt_10k_paper_no" + str(i) + ".png")

# create lineplot with 8 lines in one plot


def SINGLE_lineplot_compartments(inputs_reversed, labels_reversed, num_plots):
    plt.clf()
    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    # Define different colors for each infection state
    colors = ['blue', 'orange', 'green', 'red',
              'purple', 'brown', 'pink', 'gray']

    for i in range(num_plots):
        plt.clf()
        # One plot for all infection states
        fig, ax = plt.subplots(figsize=(10, 8))

        # Loop through each infection state and plot both inputs and labels
        for idx, (c, inp, lab) in enumerate(zip(infectionstates, inputs_reversed[i].transpose(), labels_reversed[i].transpose())):
            # Plot the first 5 days from inp (in a specific color)
            ax.plot(np.arange(1, 6),
                    inp[:5], color=colors[idx], label=f'{c} (first 5 days)')

            # Plot the remaining days from lab (in the same color but dashed)
            ax.plot(np.arange(6, 96), lab,
                    color=colors[idx], linestyle='dashed', label=f'{c} (days 6-95)')

        # Set axis labels and title
        ax.set_xlabel('Number of days')
        ax.set_ylabel('Number of individuals')
        ax.set_title('Predicted values and labels for all compartments')

        # Add a legend showing the lines for each infection state
        ax.legend(loc='upper right')

        # Save the figure
        plt.savefig(
            f"/localdata1/gnn_paper_2024/images/no_agegroups_90days_I_based_1plt/compartment_lines_noagegroups_90days_I_based1plt_10k_paper_no{i}.png")
        print(f'Plot No. {i} done.')


def SINGLE_lineplot_compartments_log_and_nolog(inputs_reversed, labels_reversed, num_plots, num_days):

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    # Define different colors for each infection state
    colors = ['blue', 'orange', 'green', 'red',
              'purple', 'brown', 'pink', 'gray']

    for i in range(num_plots):
        plt.clf()

        # Create a figure with two subplots (1 row, 2 columns)
        # Adjust the figure size as needed
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

        # First subplot (original data)
        for idx, (c, inp, lab) in enumerate(zip(infectionstates, inputs_reversed[i].transpose(), labels_reversed[i].transpose())):
            # Plot the first 5 days from inp (in a specific color)
            ax1.plot(np.arange(1, 6),
                     inp[:5], color=colors[idx], label=f'{c} (first 5 days)')

            # Plot the remaining days from lab (in the same color but dashed)
            ax1.plot(np.arange(6, num_days+1), lab,
                     color=colors[idx], linestyle='dashed', label=f'{c} (days 6-' + str(num_days))

        ax1.set_xlabel('Number of days')
        ax1.set_ylabel('Number of individuals')
        ax1.set_title('Predicted values and labels (Original Data)')
        ax1.legend(loc='upper right')

        # Second subplot (data transformed using np.log1p)
        for idx, (c, inp, lab) in enumerate(zip(infectionstates, inputs_reversed[i].transpose(), labels_reversed[i].transpose())):
            # Apply np.log1p transformation before plotting
            ax2.plot(np.arange(1, 6), np.log1p(
                inp[:5]), color=colors[idx], label=f'{c} (first 5 days, log1p)')
            ax2.plot(np.arange(6, num_days+1), np.log1p(lab),
                     color=colors[idx], linestyle='dashed', label=f'{c} (days 6-'+str(num_days) + ', log1p)')

        ax2.set_xlabel('Number of days')
        ax2.set_ylabel('log(1 + Number of individuals)')
        ax2.set_title('Predicted values and labels (log1p Transformed)')
        ax2.legend(loc='upper right')

        # Adjust layout for better spacing between subplots
        plt.tight_layout()

        # Save the figure
        plt.savefig(
            f"/localdata1/gnn_paper_2024/images/no_agegroups_90days_2plt/compartment_lines_noagegroups_90days_2plt_10k_paper_no{i}.png")
        print(f'Plot No. {i} done.')


# create boxplot for n input datasamples
def boxplot_inputs(days=60, size_tikz=10, size_labels=12, size_title=14, size_legend=12):
    file = open(
        f'/localdata1/gnn_paper_2024/data/one_population/without_agegroups/data_secir_simple_{days}days_10k.pickle', 'rb')
    file_I = open(
        f'/localdata1/gnn_paper_2024/data/one_population/without_agegroups/data_secir_simple_{days}days_I_based_10k.pickle', 'rb')

    data = pickle.load(file)
    data = np.expm1(data['inputs'])
    data_w = pickle.load(file_I)
    data_w = np.expm1(data_w['inputs'])

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
        nrows=4, ncols=2, sharey=False, figsize=(8, 10))
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

    for ax, compartment, d, d_i in zip(axes, infectionstates, np.asarray(data).transpose(),
                                       np.asarray(data_w).transpose()):

        d_df = pd.DataFrame(data=d)
        d_df = pd.melt(d_df.transpose(), var_name="Day")
        d_df['type'] = 'Outbreak'

        d_i_df = pd.DataFrame(data=d_i)
        d_i_df = pd.melt(d_i_df.transpose(), var_name="Day")
        d_i_df['type'] = 'Persistent threat'

        df_all = pd.concat([d_df, d_i_df], ignore_index=True)

        sns.boxplot(ax=ax, x='Day', y='value', data=df_all,
                    hue='type', palette='Set1', width=0.8, legend='auto')
        ax.set_title(compartment, fontsize=size_title)
        ax.legend().set_visible(False)
        ax.tick_params(axis='both', which='major', labelsize=size_tikz)
        ax.set_xlabel(ax.get_xlabel(), fontsize=size_labels)
        ax.set_ylabel(ax.get_ylabel(), fontsize=size_labels)

        handles, labels = ax.get_legend_handles_labels()

    # save legend in separate file
    save_dir = "/localdata1/gnn_paper_2024/images/without_spatial_res/boxplots_inputdata/"
    fig_legend = plt.figure(figsize=(8, 2))
    fig_legend.legend(handles, labels, loc='center', ncol=3,
                      frameon=False, fontsize=size_legend)
    fig_legend.savefig(save_dir + "boxplot_legend.png")
    plt.close(fig_legend)

    plt.tight_layout()
    plt.savefig(
        save_dir + "boxplot_input_compartments_noagegroups_comparison.png")


# boxplot_inputs(days=60)


# create boxplot for a single input datasample
def boxplot_inputs_single(label='', days=60, size_tikz=10, size_labels=12, size_title=14):
    file = open(
        f'/localdata1/gnn_paper_2024/data/one_population/without_agegroups/data_secir_simple_{days}days_I_based_10k.pickle', 'rb')

    data = pickle.load(file)
    data = np.expm1(data['inputs'])

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
        nrows=4, ncols=2, sharey=False, figsize=(8, 10))
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

    for ax, compartment, d in zip(axes, infectionstates, np.asarray(data).transpose()):

        d_df = pd.DataFrame(data=d)
        d_df = pd.melt(d_df.transpose(), var_name="Day")
        d_df['type'] = label

        sns.boxplot(ax=ax, x='Day', y='value', data=d_df,
                    hue='type', palette='Set1', width=0.8, legend='auto')
        ax.set_title(compartment, fontsize=size_title)
        ax.legend().set_visible(False)
        ax.tick_params(axis='both', which='major', labelsize=size_tikz)
        ax.set_xlabel(ax.get_xlabel(), fontsize=size_labels)
        ax.set_ylabel(ax.get_ylabel(), fontsize=size_labels)

    save_dir = "/localdata1/gnn_paper_2024/images/without_spatial_res/boxplots_inputdata/"
    plt.tight_layout()
    plt.savefig(save_dir + "boxplot_input_compartments_noagegroups_I_based.png")


boxplot_inputs_single()


def heatmap_gridsearch_results(df_gridsearch, savename, size_labels=30):
    """
    Creates three heatmaps (MLP, CNN, LSTM) side-by-side without an integrated colorbar.
    Then saves two additional figures containing a vertical and a horizontal colorbar.

    :param df_gridsearch: DataFrame with grid search results. Must have columns:
                          ['model', 'number_of_hidden_layers', 'number_of_neurons', 'kfold_val'].
    :param savename: Filename to save the main heatmap figure (e.g. "heatmap.png").
                     Two additional files will be created by appending suffixes.
    :param size_labels: Font size for axis labels and tick labels.
    """

    # --- Pivot tables for each model ---
    df_heatmap1 = df_gridsearch.loc[df_gridsearch['model'] == 'Dense',
                                    ['number_of_hidden_layers', 'number_of_neurons', 'kfold_val']]
    df_heatmap1 = df_heatmap1.pivot(index='number_of_hidden_layers',
                                    columns='number_of_neurons',
                                    values='kfold_val')

    df_heatmap2 = df_gridsearch.loc[df_gridsearch['model'] == 'CNN',
                                    ['number_of_hidden_layers', 'number_of_neurons', 'kfold_val']]
    df_heatmap2 = df_heatmap2.pivot(index='number_of_hidden_layers',
                                    columns='number_of_neurons',
                                    values='kfold_val')

    df_heatmap3 = df_gridsearch.loc[df_gridsearch['model'] == 'LSTM',
                                    ['number_of_hidden_layers', 'number_of_neurons', 'kfold_val']]
    df_heatmap3 = df_heatmap3.pivot(index='number_of_hidden_layers',
                                    columns='number_of_neurons',
                                    values='kfold_val')

    # --- Compute global vmin and vmax across all 3 heatmaps ---
    all_values = np.concatenate([df_heatmap1.values.flatten(),
                                 df_heatmap2.values.flatten(),
                                 df_heatmap3.values.flatten()])
    vmin = np.nanmin(all_values)
    vmax = np.nanmax(all_values)
    norm = Normalize(vmin=vmin, vmax=vmax)

    # --- Create the main figure with 3 subplots (no colorbar) ---
    fig, axs = plt.subplots(nrows=1, ncols=3, sharex=False,
                            figsize=(25, 10), constrained_layout=True)

    for ax, df_heatmap, name in zip(
        axs.flat,
        [df_heatmap1, df_heatmap2, df_heatmap3],
        ['MLP', 'CNN', 'LSTM']
    ):
        # Plot each pivot table with the global norm
        im = ax.imshow(df_heatmap.values, cmap='RdYlGn_r', norm=norm)

        # Ticks/labels
        ax.set_xticks(np.arange(len(df_heatmap.columns)))
        ax.set_xticklabels(df_heatmap.columns, fontsize=25)
        ax.set_yticks(np.arange(len(df_heatmap.index)))
        ax.set_yticklabels(df_heatmap.index, fontsize=25)

        ax.set_ylabel('Number of Hidden Layers', fontsize=size_labels)
        ax.set_xlabel('Units per Layer', fontsize=size_labels)
        plt.setp(ax.get_xticklabels(), rotation=45,
                 ha="right", rotation_mode="anchor")

        # Annotate cells
        for i in range(len(df_heatmap.index)):
            for j in range(len(df_heatmap.columns)):
                val = df_heatmap.values[i, j]
                ax.text(j, i, f"{val:.2f}", ha="center",
                        va="center", color="k", fontsize=25)

        ax.set_title(f'Model = {name}', fontsize=30)

    # Save the main figure (NO colorbar in it)
    savedir = "/localdata1/zunk_he/memilio/new_heatmaps/simple"
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    main_path = os.path.join(savedir, savename)
    plt.savefig(main_path, dpi=300)
    plt.close(fig)

    # --- Now create separate figures for vertical & horizontal colorbars ---
    # Create a scalar mappable using the same cmap and norm
    sm = ScalarMappable(norm=norm, cmap='RdYlGn_r')
    sm.set_array([])  # We won't pass real data, just the norm/cmap

    # Build filenames for vertical/horizontal colorbars
    base, ext = os.path.splitext(savename)
    vertical_savename = f"{base}_vertical_colorbar{ext}"
    horizontal_savename = f"{base}_horizontal_colorbar{ext}"

    # 1) Vertical colorbar
    vert_fig, vert_ax = plt.subplots(figsize=(1, 10))
    cbar_v = vert_fig.colorbar(sm, cax=vert_ax, orientation='vertical')
    cbar_v.set_label('Validation MAPE', size=20)
    cbar_v.ax.tick_params(labelsize=15)
    vert_fig.savefig(os.path.join(savedir, vertical_savename),
                     bbox_inches='tight')
    plt.close(vert_fig)

    # 2) Horizontal colorbar
    horiz_fig, horiz_ax = plt.subplots(figsize=(10, 1))
    cbar_h = horiz_fig.colorbar(sm, cax=horiz_ax, orientation='horizontal')
    cbar_h.set_label('Validation MAPE', size=20)
    cbar_h.ax.tick_params(labelsize=15)
    horiz_fig.savefig(os.path.join(
        savedir, horizontal_savename), bbox_inches='tight')
    plt.close(horiz_fig)

    # lineplots_compartments(data['inputs'], data['labels'], 100)
    # SINGLE_lineplot_compartments(data['inputs'], data['labels'], 100)

    # boxplot_inputs()

    # boxplot_inputs_single(file)

    # SINGLE_lineplot_compartments_log_and_nolog(
    #    data['inputs'], data['labels'], 100, 95)


def plot_optimizer():
    df = pd.read_csv(
        '/home/schm_a45/Documents/Code/memilio_test/memilio/pycode/memilio-surrogatemodel/memilio/secir_simple_grid_search_paper/dataframe_optimizer_paper.csv')
    df_plot = df[['optimizer', 'kfold_val']]
    fig, ax = plt.subplots()
    ax.bar(df_plot['optimizer'], df_plot['kfold_val'])
    ax.set_ylabel('Validation MAPE')
    ax.set_xlabel('Optimizer')
    ax.set_title('Validation MAPE depending on optimizer ')
    ax.bar_label(ax.containers[0], label_type='edge')
    # pad the spacing between the number and the edge of the figure
    ax.margins(y=0.1)

    plt.show()
    plt.savefig("secir_simple_optimizer.png")


# secir simple
filepath = os.path.join(
    "/localdata1/gnn_paper_2024/data/results/grid_search/without_agegroups/dataframe_30days_I_based_10k_nodamp.csv")
df = pd.DataFrame(data=pd.read_csv(filepath))
savename = "heatmap_secir_simple_30days_10k_paper.png"
heatmap_gridsearch_results(df, savename)

# secir groups
filepath = os.path.join(
    "/localdata1/gnn_paper_2024/data/results/grid_search/with_agegroups/gridserach_secir_groups_30days_I_based_Germany_10k_nodamp.csv")
df = pd.DataFrame(data=pd.read_csv(filepath))
savename = "heatmap_secir_groups_30days_10k_paper.png"
heatmap_gridsearch_results(df, savename)


# plot_optimizer()
