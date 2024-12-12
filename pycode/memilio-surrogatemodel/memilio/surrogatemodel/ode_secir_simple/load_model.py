import os
import tensorflow as tf
import pickle
from memilio.surrogatemodel.ode_secir_simple.model import split_data, get_test_statistic
# from memilio.surrogatemodel.ode_secir_groups.data_generation_nodamp import get_population
import numpy as np
from memilio.simulation.osecir import InfectionState
import matplotlib.pyplot as plt
import pandas as pd
# import seaborn as sns
from matplotlib.ticker import ScalarFormatter

# load data path = os.path.dirname(os.path.realpath(__file__))

# import data
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'data_paper')

path_data = "/localdata1/gnn_paper_2024/data_Iteration2/one_population/without_agegroups"
filename = "data_secir_simple_90days_I_based_10k.pickle"

# if not os.path.isfile(os.path.join(path_data, filename)):
#    ValueError("no dataset found in path: " + path_data)

file = open(os.path.join(path_data, filename), 'rb')

data = pickle.load(file)


test_inputs = data['inputs'][int((0.8 * len(data['inputs']))):]
test_labels = data['labels'][int((0.8 * len(data['labels']))):]


# load trained model
# new_model = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple')
new_model = tf.keras.models.load_model(
    '/localdata1/gnn_paper_2024/data_Iteration2/saved_models/saved_models_secir_simple_paper/LSTM_90days_secirsimple_I_based_10k.h5')

pred = new_model.predict(test_inputs)

pred_reversed = np.expm1(pred)
labels_reversed = np.expm1(np.asarray(test_labels))

mape_per_day = 100*np.mean(abs((test_labels - pred)/test_labels), axis=0)
mape_reversed_per_day = 100 * \
    np.mean(abs((labels_reversed - pred_reversed)/labels_reversed), axis=0)

n_runs = test_labels.shape[0]
n_days = test_labels.shape[1]
n_compartments = test_labels.shape[2]

mape_per_run = 100*np.mean(abs((np.asarray(test_labels).reshape(n_runs, n_days*n_compartments) - np.asarray(pred).reshape(
    n_runs, n_days*n_compartments))/np.asarray(test_labels).reshape(n_runs, n_days*n_compartments)), axis=1)

mape_reversed_per_run = 100 * \
    np.mean(abs((np.asarray(labels_reversed).reshape(n_runs, n_days*n_compartments) - np.asarray(pred_reversed).reshape(
        n_runs, n_days*n_compartments))/np.asarray(labels_reversed).reshape(n_runs, n_days*n_compartments)), axis=1)


def lineplots_compartments(mape_per_day, mape_reversed_per_day, savename):

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    plt.clf()
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
        nrows=4, ncols=2, sharey=False, figsize=(10, 13), constrained_layout=True)

    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

    for ax, c, m, mr in zip(axes, infectionstates, mape_per_day.transpose(), mape_reversed_per_day.transpose()):

        ax.plot(m, color='blue', label='Test MAPE scaled')

        ax.plot(mr, color='orange',
                label='Test MAPE not scaled')

        # ax.set_xlabel('Day')
        ax.set_ylabel('Test MAPE')
        ax.set_title(c, fontsize=10)

    ax7.set_xlabel('Day')
    ax8.set_xlabel('Day')

    lines = []
    line_labels = []
    for ax in fig.axes:
        Line, Label = ax.get_legend_handles_labels()
        lines.extend(Line)
        line_labels.extend(Label)

    fig.legend(lines[:2], line_labels[:2], loc='lower center')
    fig.suptitle(
        'Test MAPE for scaled ('+str(np.round(np.mean(mape_per_day), 4))+'%) and unscaled data ('+str(np.round(np.mean(mape_reversed_per_day), 4))+'%)', fontsize=16)

    plt.savefig(savename)


def lineplots_compartments_twoaxes(mape_per_day, mape_reversed_per_day, savename):
    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    plt.clf()
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
        nrows=4, ncols=2, sharey=False, figsize=(10, 13), constrained_layout=True)

    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

    for ax, c, m, mr in zip(axes, infectionstates, mape_per_day.transpose(), mape_reversed_per_day.transpose()):
        # Plot the first line (mape_per_day) on the primary y-axis (left)
        ax.plot(m, color='blue', label='Test MAPE scaled')
        ax.set_ylabel('Test MAPE (Scaled)', color='blue')
        ax.tick_params(axis='y', labelcolor='blue')

        # Create a secondary y-axis for the second line (mape_reversed_per_day)
        ax2 = ax.twinx()
        ax2.plot(mr, color='orange', label='Test MAPE not scaled')
        ax2.set_ylabel('Test MAPE (Not Scaled)', color='red')
        ax2.tick_params(axis='y', labelcolor='red')

        ax.set_title(c, fontsize=10)

    # Set x-axis labels for the bottom plots
    ax7.set_xlabel('Day')
    ax8.set_xlabel('Day')
    # Collect unique legend information for the combined legend
    lines, labels = [], []
    for ax in fig.axes:
        for line, label in zip(*ax.get_legend_handles_labels()):
            if label not in labels:  # Only add if label is unique
                lines.append(line)
                labels.append(label)

    fig.legend(lines, labels, loc='lower center')
    fig.suptitle(
        'Test MAPE for scaled ('+str(np.round(np.mean(mape_per_day), 4))+'%) and unscaled data ('+str(np.round(np.mean(mape_reversed_per_day), 4))+'%)', fontsize=16)

    plt.savefig(savename)


# def lineplots_pred_labels(pred_reversed, labels_reversed, num_plots):
#     from matplotlib.transforms import Bbox

#     infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
#                        'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

#     def get_best_text_position(ax, textstr, props, candidate_positions):
#         """
#         Determine the best position for the text box to minimize overlap with plot lines.
#         """
#         best_position = candidate_positions[0]
#         min_overlap = float('inf')  # Initialize with a very large number

#         for pos in candidate_positions:
#             x_pos, y_pos = pos
#             # Simulate the bounding box of the text at this position
#             renderer = ax.figure.canvas.get_renderer()
#             temp_text = ax.text(
#                 x_pos, y_pos, textstr, transform=ax.transAxes, fontsize=10,
#                 bbox=props, visible=False  # Don't draw it; just calculate position
#             )
#             text_bbox = temp_text.get_window_extent(renderer=renderer)
#             temp_text.remove()  # Remove the temporary text box

#             # Convert bbox to data coordinates
#             text_bbox_data = Bbox(ax.transData.inverted().transform(text_bbox))

#             # Check overlap with plotted lines
#             overlap = 0
#             for line in ax.get_lines():
#                 line_data = line.get_xydata()
#                 in_bbox = (
#                     (line_data[:, 0] >= text_bbox_data.xmin) &
#                     (line_data[:, 0] <= text_bbox_data.xmax) &
#                     (line_data[:, 1] >= text_bbox_data.ymin) &
#                     (line_data[:, 1] <= text_bbox_data.ymax)
#                 )
#                 overlap += in_bbox.sum()  # Count the points within the text box

#             if overlap < min_overlap:
#                 min_overlap = overlap
#                 best_position = pos  # Update the best position

#         return best_position

#     for i in range(num_plots):
#         plt.clf()
#         plt.tight_layout()
#         fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
#             nrows=4, ncols=2, sharey=False, figsize=(10, 13), constrained_layout=True, dpi=300)

#         axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
#         # Formatter for scientific notation (e.g., x10^6)
#         sci_formatter = ScalarFormatter(useMathText=True)
#         sci_formatter.set_powerlimits((6, 6))  # Force le6 notation (1e6)

#         # Candidate positions in Axes coordinates
#         candidate_positions = [
#             (0.05, 0.95),  # Top-left
#             (0.95, 0.95),  # Top-right
#             (0.05, 0.05),  # Bottom-left
#             (0.95, 0.05),  # Bottom-right
#         ]

#         for ax, c, p, l, input in zip(axes, infectionstates, pred_reversed[i].transpose(), labels_reversed[i].transpose(), np.expm1(np.asarray(test_inputs))[i].transpose()):
#             ax.plot(np.arange(1, 6), input, color='black',
#                     label='Inputs', linewidth=3)

#             ax.plot(np.arange(6, pred_reversed.shape[1] + 6), l, color='orange',
#                     label='Labels', linewidth=3)

#             ax.plot(np.arange(6, pred_reversed.shape[1]+6), p, color='blue', label='Predictions',
#                     linestyle='--', linewidth=2)

#             # Prepare text box properties and text
#             textstr = '\n'.join((
#                 'not log MAPE = ' +
#                 str(np.round(100 * np.mean(abs((l - p) / l)), 4)) + '%',
#                 'log MAPE: ' +
#                 str(np.round(
#                     100 * np.mean(abs((np.log1p(l) - np.log1p(p)) / np.log1p(l))), 4)) + '%'
#             ))
#             props = dict(boxstyle='round,pad=0.3',
#                          facecolor='lightgray', edgecolor='black', alpha=0.8)

#             # Dynamically determine the best position for the text box
#             best_position = get_best_text_position(
#                 ax, textstr, props, candidate_positions)
#             best_position = (0.95, 0.95)
#             ax.text(*best_position, textstr, transform=ax.transAxes, fontsize=10,
#                     verticalalignment='top', bbox=props)

#             ax.set_ylabel('Number of individuals')
#             ax.set_title(c)
#             ax.yaxis.set_major_formatter(sci_formatter)

#         ax7.set_xlabel('Day')
#         ax8.set_xlabel('Day')

#         lines = []
#         line_labels = []
#         for ax in fig.axes:
#             Line, Label = ax.get_legend_handles_labels()
#             lines.extend(Line)
#             line_labels.extend(Label)

#         fig.legend(lines[:3], line_labels[:3], loc='lower center')

#         print('Plot No.' + str(i) + ' done.')
#         plt.savefig("/localdata1/gnn_paper_2024/images/pred_labels_trajectories/no_agegroups_90days_I_based_pred_labels_withinput/compartment_lines_noagegroups_90days_I_based_pred_labels_10k_paper_no" + str(i) + ".png")


def lineplots_pred_labels(pred_reversed, labels_reversed, num_plots):
    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    for i in range(num_plots):
        # Reset figure
        plt.clf()
        fig, axes = plt.subplots(
            nrows=4, ncols=2, figsize=(10, 13), constrained_layout=True, dpi=300
        )
        axes = axes.flatten()  # Flatten the array of subplots for easier iteration

        # Formatter for scientific notation (e.g., x10^6)
        sci_formatter = ScalarFormatter(useMathText=True)
        # Force scientific notation (e.g., x10^6)
        sci_formatter.set_powerlimits((6, 6))

        for ax, state, pred, label, input_data in zip(
                axes, infectionstates, pred_reversed[i].transpose(),
                labels_reversed[i].transpose(), np.expm1(np.asarray(test_inputs))[i].transpose()):

            # Plot lines
            ax.plot(np.arange(1, 6), input_data,
                    color='black', label='Inputs', linewidth=3)
            ax.plot(np.arange(
                6, pred_reversed.shape[1] + 6), label, color='orange', label='Labels', linewidth=3)
            ax.plot(np.arange(6, pred_reversed.shape[1] + 6), pred, color='blue', label='Predictions',
                    linestyle='--', linewidth=2)

            # Compute and add annotation text
            not_log_mape = np.round(
                100 * np.mean(abs((label - pred) / label)), 4)
            log_mape = np.round(
                100 * np.mean(abs((np.log1p(label) - np.log1p(pred)) / np.log1p(label))), 4)
            textstr = f"not log MAPE = {not_log_mape}%\nlog MAPE = {log_mape}%"
            props = dict(boxstyle='round,pad=0.3',
                         facecolor='lightgray', edgecolor='black', alpha=0.8)

            # Dynamically position the text box in the upper-right corner
            ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                    verticalalignment='top', horizontalalignment='right', bbox=props)

            # Configure axes
            ax.set_ylabel('Number of individuals')
            ax.set_title(state)
            ax.yaxis.set_major_formatter(sci_formatter)

        # Add shared x-axis labels
        axes[-2].set_xlabel('Day')
        axes[-1].set_xlabel('Day')

        # Consolidate the legend
        lines, labels = axes[0].get_legend_handles_labels()
        fig.legend(lines, labels, loc='lower center',
                   ncol=3, frameon=False, fontsize=10)

        # Save the figure
        save_path = f"/localdata1/gnn_paper_2024/images/pred_labels_trajectories/no_agegroups_90days_I_based_pred_labels_withinput/compartment_lines_noagegroups_90days_I_based_pred_labels_10k_paper_no{i}.png"
        plt.savefig(save_path, bbox_inches='tight')
        print(f'Plot No. {i} saved')


def histogram_mape_per_run(mape_per_run, savename):

    plt.figure().clf()

    fig, ax = plt.subplots()
    # axs[0].ticklabel_format(style='plain')
    ax.hist(mape_per_run, bins=100)
    ax.set_xlabel('Test MAPE')
    ax.set_ylabel('Frequency')
    ax.set_title('Histogram for Test MAPE per run')
    # ax.set_xticks(np.arange(0,3500000,500000))

    plt.savefig(savename)


def mape_log_original(num_plots, pred_reversed, labels_reversed, test_inputs):

    log_mapes = []
    original_mape = []

    for i in range(num_plots):
        # Reset figure
        plt.clf()
        fig, axes = plt.subplots(
            nrows=4, ncols=2, figsize=(10, 13), constrained_layout=True, dpi=300
        )
        axes = axes.flatten()  # Flatten the array of subplots for easier iteration

        # Formatter for scientific notation (e.g., x10^6)
        sci_formatter = ScalarFormatter(useMathText=True)
        # Force scientific notation (e.g., x10^6)
        sci_formatter.set_powerlimits((6, 6))

        for ax,  pred, label, input_data in zip(
                axes,  pred_reversed[i].transpose(),
                labels_reversed[i].transpose(), np.expm1(np.asarray(test_inputs))[i].transpose()):

            # Plot lines
            ax.plot(np.arange(1, 6), input_data,
                    color='black', label='Inputs', linewidth=3)
            ax.plot(np.arange(
                6, pred_reversed.shape[1] + 6), label, color='orange', label='Labels', linewidth=3)
            ax.plot(np.arange(6, pred_reversed.shape[1] + 6), pred, color='blue', label='Predictions',
                    linestyle='--', linewidth=2)

            # Compute and add annotation text
            not_log_mape = np.round(
                100 * np.mean(abs((label - pred) / label)), 4)
            original_mape.append(not_log_mape)
            log_mape = np.round(
                100 * np.mean(abs((np.log1p(label) - np.log1p(pred)) / np.log1p(label))), 4)
            log_mapes.append(log_mape)
    df = pd.DataFrame(data={'log': log_mapes, 'original': original_mape})
    plt.clf()
    plt.figure(figsize=(5, 5))
    plt.scatter(df['log'], df['original'], s=3)
    plt.xlabel('MAPE(log)')
    plt.ylabel('MAPE (original)')
    plt.savefig('SimpleNN_log_original_90days.png')


def lineplots_pred_labels_selected_plot(pred_reversed, labels_reversed, plotID):
    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    plt.clf()
    fig, axes = plt.subplots(
        nrows=4, ncols=2, figsize=(10, 13), constrained_layout=True, dpi=300
    )
    axes = axes.flatten()  # Flatten the array of subplots for easier iteration

    # Formatter for scientific notation (e.g., x10^6)
    sci_formatter = ScalarFormatter(useMathText=True)
    # Force scientific notation (e.g., x10^6)
    sci_formatter.set_powerlimits((6, 6))

    for ax, state, pred, label, input_data in zip(
            axes, infectionstates, pred_reversed[plotID].transpose(),
            labels_reversed[plotID].transpose(), np.expm1(np.asarray(test_inputs))[plotID].transpose()):

        # Plot lines
        ax.plot(np.arange(1, 6), input_data,
                color='black', label='Inputs', linewidth=3)
        ax.plot(np.arange(
                6, pred_reversed.shape[1] + 6), label, color='orange', label='Labels', linewidth=3)
        ax.plot(np.arange(6, pred_reversed.shape[1] + 6), pred, color='blue', label='Predictions',
                linestyle='--', linewidth=2)

        # Compute and add annotation text
        not_log_mape = np.round(
            100 * np.mean(abs((label - pred) / label)), 4)
        log_mape = np.round(
            100 * np.mean(abs((np.log1p(label) - np.log1p(pred)) / np.log1p(label))), 4)
        textstr = f"not log MAPE = {not_log_mape}%\nlog MAPE = {log_mape}%"
        props = dict(boxstyle='round,pad=0.3',
                     facecolor='lightgray', edgecolor='black', alpha=0.8)

        # Dynamically position the text box in the upper-right corner
        ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right', bbox=props)

        # Configure axes
        ax.set_ylabel('Number of individuals')
        ax.set_title(state)
        ax.yaxis.set_major_formatter(sci_formatter)
        ax.set_yscale('log')

    # Add shared x-axis labels
    axes[-2].set_xlabel('Day')
    axes[-1].set_xlabel('Day')

    # Consolidate the legend
    lines, labels = axes[0].get_legend_handles_labels()
    fig.legend(lines, labels, loc='lower center',
               ncol=3, frameon=False, fontsize=10)

    # Save the figure
    save_path = f"/localdata1/gnn_paper_2024/images/pred_labels_trajectories/compartment_lines_noagegroups_90days_I_based_pred_labels_10k_paper_no{plotID}.png"
    plt.savefig(save_path, bbox_inches='tight')
    print(f'Plot No. {plotID} saved')


savename = 'secir_noagegroup_90days_I_based_MAPEs'
lineplots_compartments(mape_per_day, mape_reversed_per_day, savename)


savename_2 = 'secir_noagegroup_90days_I_based_MAPEs_two_axex'
lineplots_compartments_twoaxes(mape_per_day, mape_reversed_per_day, savename_2)

lineplots_pred_labels(pred_reversed, labels_reversed, num_plots=100)

savename = 'histogram_mape_per_run_simpleNN.png'
histogram_mape_per_run(mape_per_run, savename)
