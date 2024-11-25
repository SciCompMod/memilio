import os
import tensorflow as tf
import pickle
from memilio.surrogatemodel.ode_secir_simple.model import split_data, get_test_statistic
# from memilio.surrogatemodel.ode_secir_groups.data_generation_nodamp import get_population
import numpy as np
from memilio.simulation.osecir import InfectionState
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import ScalarFormatter

# load data path = os.path.dirname(os.path.realpath(__file__))

# import data
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'data_paper')

path_data = "/localdata1/gnn_paper_2024/data_Iteration2/one_population/without_agegroups/"
filename = "data_secir_simple_90days_I_based_10k.pickle"

if not os.path.isfile(os.path.join(path_data, filename)):
    ValueError("no dataset found in path: " + path_data)

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


def lineplots_pred_labels(pred_reversed, labels_reversed, num_plots):

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    for i in range(num_plots):

        plt.clf()
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
            nrows=4, ncols=2, sharey=False, figsize=(10, 13), constrained_layout=True)

        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
        # Formatter for scientific notation (e.g., x10^6)
        sci_formatter = ScalarFormatter(useMathText=True)
        sci_formatter.set_powerlimits((6, 6))  # Force le6 notation (1e6)
        for ax, c, p, l, input in zip(axes, infectionstates, pred_reversed[i].transpose(), labels_reversed[i].transpose(), np.expm1(np.asarray(test_inputs))[i].transpose()):
            ax.plot(np.arange(1, 6), input, color='black',
                    label='Inputs', linewidth=2)

            ax.plot(np.arange(6, 96), l, color='orange',
                    label='Labels', linewidth=2)

            ax.plot(np.arange(6, 96), p, color='blue', label='Predictions',
                    linestyle='--', linewidth=1)

            # ax.set_xlabel('Day')
            ax.set_ylabel('Number of individuals')
            ax.set_title(c+' / MAPE = ' + str(np.round(100 *
                         np.mean(abs((l - p)/l)), 4))+'%'+"//scaled MAPE: " + str(np.round(100 *
                                                                                           np.mean(abs((np.log1p(l) - np.log1p(p))/np.log1p(l))), 4)), fontsize=10)
            ax.yaxis.set_major_formatter(sci_formatter)
        ax7.set_xlabel('Day')
        ax8.set_xlabel('Day')

        lines = []
        line_labels = []
        for ax in fig.axes:
            Line, Label = ax.get_legend_handles_labels()
            lines.extend(Line)
            line_labels.extend(Label)

        fig.legend(lines[:3], line_labels[:3], loc='lower center')

        print('Plot No.'+str(i) + ' done.')
        plt.savefig("/localdata1/gnn_paper_2024/images/pred_labels_trajectories/no_agegroups_90days_I_based_pred_labels_withinput/compartment_lines_noagegroups_90days_I_based_pred_labels_10k_paper_no" + str(i) + ".png")


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


savename = 'secir_noagegroup_90days_I_based_MAPEs'
lineplots_compartments(mape_per_day, mape_reversed_per_day, savename)


savename_2 = 'secir_noagegroup_90days_I_based_MAPEs_two_axex'
lineplots_compartments_twoaxes(mape_per_day, mape_reversed_per_day, savename_2)

lineplots_pred_labels(pred_reversed, labels_reversed, num_plots=100)

savename = 'histogram_mape_per_run_simpleNN.png'
histogram_mape_per_run(mape_per_run, savename)
