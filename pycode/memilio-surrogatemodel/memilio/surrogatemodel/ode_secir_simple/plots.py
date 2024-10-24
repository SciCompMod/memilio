import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pickle

from memilio.simulation.osecir import InfectionState
import seaborn as sns

filename = 'data_secir_simple_90days_10k.pickle'

# import data
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
    os.path.dirname(os.path.realpath(path)))), 'data_paper')
# path_data = '/home/schm_a45/Documents/Code/memilio_test/memilio/pycode/memilio-surrogatemodel/memilio/data_paper'

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

def boxplot_inputs():
    file = open('/home/schm_a45/Documents/Code/memilio_test/memilio/pycode/memilio-surrogatemodel/memilio/data_paper/data_secir_simple_10k.pickle', 'rb')
    file_I = open('/home/schm_a45/Documents/Code/memilio_test/memilio/pycode/memilio-surrogatemodel/memilio/data_paper/data_secir_simple_30days_I_based_10k.pickle', 'rb')

    data = pickle.load(file)
    data = np.expm1(data['inputs'])
    data_w = pickle.load(file_I)
    data_w = np.expm1(data_w['inputs'])
    # data_w2 = pickle.load(file_w2)
    # data_w2 = np.expm1(data_w2['inputs'])

    # compartment_array = []
    # for compartment in InfectionState.values():
    #        compartment_array.append(compartment)
    # infectionstates = [str(compartment).split('.')[1] for compartment in compartment_array]

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
        nrows=4, ncols=2, sharey=False, figsize=(8, 10))
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

    # for ax, compartment, d, d_i, dw2 in zip(axes, index, np.asarray(data).transpose(),
    #                                        np.asarray(data_w).transpose(), np.asarray(data_w2).transpose()):

    for ax, compartment, d, d_i in zip(axes, infectionstates, np.asarray(data).transpose(),
                                       np.asarray(data_w).transpose()):

        d_df = pd.DataFrame(data=d)
        d_df = pd.melt(d_df.transpose(), var_name="Day")
        d_df['type'] = 'b'

        d_i_df = pd.DataFrame(data=d_i)
        d_i_df = pd.melt(d_i_df.transpose(), var_name="Day")
        d_i_df['type'] = 'b_i'

        # dw2_df = pd.DataFrame(data = dw2)
        # dw2_df = pd.melt(dw2_df.transpose(), var_name="Day")
        # dw2_df['type'] = 'w2'
        # df_all = pd.concat([d_df, dw1_df, dw2_df], ignore_index=True)
        df_all = pd.concat([d_df, d_i_df], ignore_index=True)

        sns.boxplot(ax=ax, x='Day', y='value', data=df_all,
                    hue='type', palette='Set1', width=0.8, legend='auto')
        ax.set_title(compartment, fontsize=10)
        ax.legend().set_visible(False)

        handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', ncol=3,
               bbox_to_anchor=(0.7, 1), frameon=False)
    plt.tight_layout()

    plt.savefig("boxplot_input_compartments_noagegroups_comparison.png")


def boxplot_inputs_single(file):

    data = pickle.load(file)
    data = np.expm1(data['inputs'])

    # compartment_array = []
    # for compartment in InfectionState.values():
    #        compartment_array.append(compartment)
    # infectionstates = [str(compartment).split('.')[1] for compartment in compartment_array]

    infectionstates = ['Susceptible', 'Exposed', 'InfectedNoSymptoms',
                       'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
        nrows=4, ncols=2, sharey=False, figsize=(8, 10))
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

    for ax, compartment, d in zip(axes, infectionstates, np.asarray(data).transpose()):

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

    plt.savefig("boxplot_input_compartments_noagegroups_I_based.png")


# lineplots_compartments(data['inputs'], data['labels'], 100)
# SINGLE_lineplot_compartments(data['inputs'], data['labels'], 100)

# boxplot_inputs()

# file = open('/home/schm_a45/Documents/Code/memilio_test/memilio/pycode/memilio-surrogatemodel/memilio/data_paper/data_secir_simple_30days_I_based_10k.pickle', 'rb')
# boxplot_inputs_single(file)


SINGLE_lineplot_compartments_log_and_nolog(
    data['inputs'], data['labels'], 100, 95)

# plt.savefig("/localdata1/gnn_paper_2024/images/no_agegroups/compartment_lines_noagegroups_30days_10k_paper_no" +
# str(i) + ".png")
