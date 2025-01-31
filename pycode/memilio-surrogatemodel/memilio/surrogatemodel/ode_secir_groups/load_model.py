import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pickle
import matplotlib.gridspec as gridspec
from memilio.simulation.osecir import InfectionState
import seaborn as sns
import tensorflow as tf
from matplotlib.ticker import ScalarFormatter

days = 30
path_data = "/localdata1/gnn_paper_2024/data/one_population/with_agegroups_Germany/nodamp/"
filename = f'data_secir_groups_{days}days_I_based_Germany_10k_nodamp.pickle'

# import data
if not os.path.isfile(os.path.join(path_data, filename)):
    ValueError("no dataset found in path: " + path_data)

file = open(os.path.join(path_data, filename), 'rb')
data = pickle.load(file)

# data['inputs'] = np.expm1(np.asarray(data['inputs']))
# data['labels'] = np.expm1(np.asarray(data['labels']))

test_inputs = data['inputs'][int(0.8 * len(data['inputs'])):]
test_labels = data['labels'][int(0.8 * len(data['labels'])):]


# load trained model
# new_model = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple')
new_model = tf.keras.models.load_model(
    f'/localdata1/gnn_paper_2024/data/results/saved_models/saved_models_secir_groups_paper//LSTM_NODAMP_{days}days_I_based_secirgroups_10k.h5')

pred = new_model.predict(test_inputs)

pred_reversed = np.expm1(pred)
labels_reversed = np.expm1(np.asarray(test_labels))

mape_per_day = 100*np.mean(abs((test_labels - pred)/test_labels), axis=0)
mape_reversed_per_day = 100 * \
    np.mean(abs((labels_reversed - pred_reversed)/labels_reversed), axis=0)


def mape_per_day_compartments_per_agegroup(mape_per_day, mape_reversed_per_day, age_groups, savename):

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    num_age_groups = len(age_groups)
    num_compartments = len(infectionstates)

    # Reshape MAPE data into (days, num_age_groups, num_compartments)
    reshaped_mape = mape_per_day.reshape(
        days, num_age_groups, num_compartments)
    reshaped_mape_reversed = mape_reversed_per_day.reshape(
        days, num_age_groups, num_compartments)

    for group_idx, group in enumerate(age_groups):
        plt.clf()
        fig, axes = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(
            10, 13), constrained_layout=True)
        axes = axes.flatten()  # Flatten the 2D array

        # Extract MAPE values for the current age group
        mape_group = reshaped_mape[:, group_idx, :]
        mape_reversed_group = reshaped_mape_reversed[:, group_idx, :]

        for ax, compartment, m, mr in zip(axes, infectionstates, mape_group.T, mape_reversed_group.T):
            ax.plot(m, color='blue', label='MAPE (scaled)')
            ax.plot(mr, color='orange', label='MAPE (log scale)')
            ax.set_ylabel('MAPE (orig. scale)')
            ax.set_title(compartment, fontsize=10)
            # Use logarithmic scale
            ax.set_yscale('log')

        # Set x-axis labels for the last row only
        axes[-2].set_xlabel('Day')
        axes[-1].set_xlabel('Day')

        lines, labels = axes[0].get_legend_handles_labels()
        fig.legend(lines, labels, loc='lower center', ncol=2)

        fig.suptitle(
            f'MAPE for Age Group {group}: Log ({np.round(np.mean(mape_group), 4)}%) | Orig. ({np.round(np.mean(mape_reversed_group), 4)}%)', fontsize=14)

        # Save plot in a directory specific to the age group
        plot_dir = f"/localdata1/gnn_paper_2024/images/lineplots_MAPE_per_day/with_agegroups/secir_with_agegroup_{days}days_I_based/"
        os.makedirs(plot_dir, exist_ok=True)
        save_path = os.path.join(
            plot_dir, savename + "_group_" + group + ".png")
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Saved MAPE plot for age group {group} at: {save_path}")

        plt.close(fig)


def lineplots_pred_labels_per_agegroup(pred_reversed, labels_reversed, test_inputs, age_groups, num_plots, savename):

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    num_age_groups = len(age_groups)
    num_compartments = len(infectionstates)

    # Reshape data: (num_plots, days, num_age_groups, num_compartments)
    reshaped_pred = pred_reversed.reshape(
        2000, days, num_age_groups, num_compartments)
    reshaped_labels = labels_reversed.reshape(
        2000, days, num_age_groups, num_compartments)
    reshaped_inputs = np.expm1(np.asarray(test_inputs)).reshape(
        2000, -1, num_age_groups, num_compartments)

    for i in range(num_plots):
        for group_idx, group in enumerate(age_groups):
            fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(
                10, 13), constrained_layout=True, dpi=300)
            axes = axes.flatten()

            # Extract data for the current age group
            pred_group = reshaped_pred[i, :, group_idx, :]
            labels_group = reshaped_labels[i, :, group_idx, :]
            inputs_group = reshaped_inputs[i, :, group_idx, :]

            # scientific notation
            sci_formatter = ScalarFormatter(useMathText=True)
            sci_formatter.set_powerlimits((6, 6))

            for ax, state, pred, label_data, input_data in zip(axes, infectionstates, pred_group.T, labels_group.T, inputs_group.T):
                pred_plot = np.insert(pred, 0, input_data[-1])
                label_plot = np.insert(label_data, 0, input_data[-1])

                # Plot lines
                ax.plot(
                    np.arange(1, 6), input_data[:5], color='black', label='Inputs', linewidth=3)
                ax.plot(np.arange(5, days + 6), label_plot,
                        color='red', label='Labels', linewidth=3)
                ax.plot(np.arange(5, days + 6), pred_plot, color='deepskyblue',
                        label='Predictions', linestyle='--', linewidth=2)

                # Compute and add text (MAPE)
                not_log_mape = np.round(
                    100 * np.mean(abs((label_data - pred) / label_data)), 4)
                log_mape = np.round(
                    100 * np.mean(abs((np.log1p(label_data) - np.log1p(pred)) / np.log1p(label_data))), 4)

                textstr = '\n'.join((
                    f'MAPE (log scale): {log_mape}%',
                    f'MAPE (orig. scale): {not_log_mape}%'
                ))

                props = dict(boxstyle='round,pad=0.3',
                             facecolor='lightgray', edgecolor='black', alpha=0.8)
                ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                        verticalalignment='top', horizontalalignment='left', bbox=props)

                ax.set_ylabel('Number of individuals')
                ax.set_title(state)
                ax.yaxis.set_major_formatter(sci_formatter)

            # X-axis labels for last row only
            axes[-2].set_xlabel('Day')
            axes[-1].set_xlabel('Day')

            lines, labels = axes[0].get_legend_handles_labels()
            fig.legend(lines, labels, loc='lower center',
                       ncol=3, frameon=False, fontsize=10)

            # Save plo
            plot_dir = f"/localdata1/gnn_paper_2024/images/lineplots_compartments/with_agegroups/secir_with_agegroup_{days}days_I_based/"
            os.makedirs(plot_dir, exist_ok=True)

            save_path = os.path.join(
                plot_dir, f"pred_labels_group_{group}_plot_{i}.png")
            plt.savefig(save_path, bbox_inches='tight')
            print(
                f"Saved {group} Plot Nr {i}")

            plt.close(fig)  # Close figure to prevent memory issues


# def SINGLE_lineplot_per_group(inputs_reversed, labels_reversed, num_plots, input_width, label_width):
#     num_days = input_width + label_width
#     infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
#                        'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

#     age_groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']
#     colors = ['blue', 'orange', 'green', 'red',
#               'purple', 'brown', 'pink', 'gray']

#     # Reshape arrays to have age group as a separate dimension
#     # Reshape data to enable plotting by age group
#     reshaped_array_input = inputs_reversed.reshape(10000, input_width, 6, 8)
#     reshaped_array_labels = labels_reversed.reshape(10000, label_width, 6, 8)

#     for i in range(num_plots):
#         for group_idx, group in enumerate(age_groups):
#             fig = plt.figure(figsize=(15, 10), constrained_layout=True)
#             grid = gridspec.GridSpec(4, 2, wspace=0.3, hspace=0.4)

#             # Extract data for this group
#             input_per_group = reshaped_array_input[i, :, group_idx, :]
#             label_per_group = reshaped_array_labels[i, :, group_idx, :]

#             for j, compartment in enumerate(infectionstates):
#                 ax = fig.add_subplot(grid[j])

#                 # Extract relevant compartment data
#                 inp = input_per_group[:, j]
#                 lab = label_per_group[:, j]

#                 # Plot input and label data
#                 ax.plot(
#                     np.arange(1, 6), inp[:5], color=colors[j], linestyle='dashed', label="Input")
#                 ax.plot(np.arange(6, num_days + 1), lab,
#                         color=colors[j], label="Label")

#                 ax.set_title(compartment)
#                 if j % 2 == 0:
#                     ax.set_ylabel('Number of individuals')
#                 if j >= 6:
#                     ax.set_xlabel('Days')

#             # Save each figure per age group
#             plot_dir = "/localdata1/gnn_paper_2024/images/lineplots_compartments/with_agegroups/secir_with_agegroup_30days_I_based/"
#             if not os.path.exists(plot_dir):
#                 os.makedirs(plot_dir)

#             save_path = os.path.join(
#                 plot_dir, f"compartment_lines_group_{group}_plot_{i}.png")
#             plt.savefig(save_path, bbox_inches='tight')
#             print(f'Plot for group {group} - No. {i} saved.')

#             plt.close(fig)  # Close figure to prevent memory issues


def boxplot_inputs_single(file_path, savename):

    file = open(file_path, 'rb')
    data = pickle.load(file)
    data = np.expm1(data['inputs'])

    reshaped_array_input = data.reshape(10000, 5, 6, 8)
    reshaped_array_input = reshaped_array_input.transpose(2, 0, 1, 3)

    # compartment_array = []
    # for compartment in InfectionState.values():
    #        compartment_array.append(compartment)
    # infectionstates = [str(compartment).split('.')[1] for compartment in compartment_array]
    agegroups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']
    for group_input, agegroup in zip(reshaped_array_input, agegroups):

        infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                           'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
            nrows=4, ncols=2, sharey=False, figsize=(8, 10))
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

        for ax, compartment, d in zip(axes, infectionstates, np.asarray(group_input).transpose()):

            d_df = pd.DataFrame(data=d)
            d_df = pd.melt(d_df.transpose(), var_name="Day")
            d_df['type'] = 'b'

            sns.boxplot(ax=ax, x='Day', y='value', data=d_df,
                        hue='type', palette='Set1', width=0.8, legend='auto')
            ax.set_title(compartment, fontsize=10)
            ax.legend().set_visible(False)

            handles, labels = ax.get_legend_handles_labels()
            fig.legend(handles, labels, loc='upper right', ncol=3,
                       bbox_to_anchor=(0.7, 1), frameon=False)
        plt.tight_layout()

        plt.savefig(
            savename+agegroup+".png")


def heatmap_gridsearch_results(df_gridsearch, savename):
    df = df_gridsearch

    plt.figure().clf()
    df_heatmap1 = pd.DataFrame(data=df.loc[(df['model'] == 'Dense')][[
                               'number_of_hidden_layers', 'number_of_neurons', 'kfold_val']])
    df_heatmap1 = df_heatmap1.pivot(
        index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_val')

    df_heatmap2 = pd.DataFrame(data=df.loc[(df['model'] == 'CNN')][[
                               'number_of_hidden_layers', 'number_of_neurons', 'kfold_val']])
    df_heatmap2 = df_heatmap2.pivot(
        index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_val')

    df_heatmap3 = pd.DataFrame(data=df.loc[(df['model'] == 'LSTM')][[
                               'number_of_hidden_layers', 'number_of_neurons', 'kfold_val']])
    df_heatmap3 = df_heatmap3.pivot(
        index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_val')

    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=False,
                            figsize=(20, 20), constrained_layout=True)

    for ax, df_heatmap, name in zip(axs.flat, [df_heatmap1, df_heatmap2, df_heatmap3], ['MLP', 'CNN', 'LSTM']):

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
    fig.delaxes(axs[1][1])
    plt.show()
    plt.savefig(savename)


def barplot_hyperparameters(filename):
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'secir_groups_grid_search_paper')
    filepath = os.path.join(path_data, filename)
    df = pd.DataFrame(data=pd.read_csv(filepath))

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
    plt.savefig("secir_groups_optimizer.png")


def heatmap_activation_optimiizer(filename, filename_2):
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'secir_groups_grid_search_paper')
    filepath = os.path.join(path_data, filename)
    df = pd.DataFrame(data=pd.read_csv(filepath))
    # plot for part2
    filepath = os.path.join(path_data, filename_2)
    df_2 = pd.DataFrame(data=pd.read_csv(filepath))
    df_concat = pd.conct(df, df_2)
    df_plot = df_concat[['optimizer', 'activation', 'kfold_test']]

    plt.figure().clf()

    df_plot = df_plot.pivot(
        index='optimizer', columns='activation', values='kfold_val')

    fig, ax = plt.subplots()
    im = ax.imshow(df_plot.values, cmap='RdYlGn_r')

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(df_plot.columns)), labels=df_plot.columns)
    ax.set_yticks(np.arange(len(df_plot.index)), labels=df_plot.index)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(df_plot.index)):
        for j in range(len(df_plot.columns)):
            text = ax.text(j, i, np.around(df_plot.values, decimals=2)[i, j],
                           ha="center", va="center", color="black")

    cbar_kw = {}
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel('MAPE', rotation=-90, va="bottom")

    ax.set_title("Activation anf Optimizer for LSTM")
    fig.tight_layout()
    plt.show()
    plt.savefig("heatmap_activation_optimizer_LSTM_GroupsNN_I_based.png")


# line plots
title_add = f'_{days}d_Ibased_10k_noDamp'
savename = 'mape_per_day_GroupsNN_log_and_nonlog_mape' + title_add
age_groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']
# mape_per_day_compartments_per_agegroup(
#     mape_per_day, mape_reversed_per_day, age_groups, savename)


savename_2 = f'compartments_secir_noagegroup_{days}days_I_based' + title_add
lineplots_pred_labels_per_agegroup(pred_reversed, labels_reversed,
                                   test_inputs, age_groups, 50, savename_2)

# # BOXPLOT INPUTS
# file_path = '/home/schm_a45/Documents/Code/memilio_test/memilio/pycode/memilio-surrogatemodel/memilio/data_paper/data_secir_groups_30days_I_based_Germany_10k.pickle'
# savename = 'boxplot_withagegroups_I_based_10k_'
# boxplot_inputs_single(file_path, savename)

# # HEATMAP
# ibased = True
# age_groups = True
# if ibased:
#     filename = 'gridserach_secir_groups_30days_I_based_Germany_10k_nodamp.csv'
# else:
#     filename = 'gridsearch_data_secir_groups_30days_Germany_10k_nodamp.csv'
# path = os.path.dirname(os.path.realpath(__file__))
# path_data = os.path.join(os.path.dirname(os.path.realpath(
#     os.path.dirname(os.path.realpath(path)))), 'secir_groups_grid_search_paper')
# filepath = os.path.join(path_data, filename)
# df = pd.DataFrame(data=pd.read_csv(filepath))
# title_groups_add = "noagegroups"
# if age_groups:
#     title_groups_add = "withagegroups"

# if ibased:
#     title_add = title_groups_add + '_30d_Ibased_10k_noDamp.png'
# else:
#     title_add = title_groups_add + '_30d_initial_10k_noDamp.png'
# savename = "heatmap_secir_withagegroups" + title_add
# # NODAMP_heatmap_secir_withagegroups_30days_I_based_10k_paper

# heatmap_gridsearch_results(df, savename)

# # HYPERPARAMETERS
# filename = 'groups_dataframe_optimizer_paper.csv'
# barplot_hyperparameters(filename)

# # HYPERPARAMETERS I_BASED
# filename = 'groups_I_based_dataframe_optimizer_paper.csv'
# filename_2 = 'groups_I_based_dataframe_optimizer_paper_part2.csv'
# heatmap_activation_optimiizer(filename, filename_2)
