import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import os

colors = ['#1f77b4', '#2ca02c', '#ff7f0e']
linestyles = ['-', '--', '-.', ':']
fontsize_labels = 16
fontsize_legends = 12

models = ['Model C (ODE)', 'Model D (Graph-ODE)']


def plot_flops(name='number_flops'):
    fig, ax = plt.subplots()

    def flops_equation_based(x, eta):
        return (4*x**2+22*x)/eta

    def flops_graph_based(x, eta):
        return (43*x**2+23*x/eta)

    x = np.linspace(0, 400, 80)

    for idx, eta in enumerate([0.05, 0.1, 0.2, 0.5]):
        ax.plot(
            x, flops_equation_based(x, eta),
            linewidth=1.5, linestyle=linestyles[idx],
            color=colors[0],
            label=models[0] + ', $h=$' + str(eta))
        ax.plot(
            x, flops_graph_based(x, eta),
            linewidth=1.5, linestyle=linestyles[idx],
            color=colors[1],
            label=models[1] + ', $h=$' + str(eta))
    ax.set_ylim(bottom=0.)
    ax.set_xlim(left=0., right=400.)
    ax.set_xlabel('Number of regions', fontsize=fontsize_labels)
    ax.set_ylabel('Number of FLOP', fontsize=fontsize_labels)

    handles, labels = ax.get_legend_handles_labels()
    sorted_handles_labels = sorted(zip(handles, labels), key=lambda x: x[1])

    sorted_handles, sorted_labels = zip(*sorted_handles_labels)

    plt.tight_layout()
    ax.legend(sorted_handles, sorted_labels, fontsize=fontsize_legends)
    plt.grid(linestyle='--')

    plot_dir = os.path.join(os.path.dirname(__file__), '../Plots')
    plt.savefig(os.path.join(plot_dir, name), bbox_inches='tight', dpi=500)
    plt.close()


def compare_runtime_and_flops(files, name=''):
    fig, ax1 = plt.subplots()
    plt.grid(True, linestyle='--')

    for idx, file in enumerate(files):
        df = pd.read_json(file)

        ax1.plot(df["Regions"], df["Time"],
                 linestyle='--', marker='o', linewidth=1.2, label=models[idx])

    ax1.set_ylim(bottom=0.)
    ax1.set_xlim(left=0., right=400.)
    ax1.set_xlabel('Number of regions', fontsize=fontsize_labels)
    ax1.set_ylabel('Run time [seconds]', fontsize=fontsize_labels)

    ax2 = ax1.twinx()

    def flops_equation_based(x):
        return (4*x**2+22*x)*200

    def flops_graph_based(x):
        return (43*x**2+230*x)*20

    x = np.linspace(0, 400, 400)

    ax2.plot(x, flops_equation_based(x), linewidth=2,
             color='darkblue', label='FLOP Model C')
    ax2.plot(x, flops_graph_based(x), linewidth=2,
             color='#9C180D', label='FLOP Model D')
    ax2.set_ylabel('Number of FLOP', fontsize=fontsize_labels)
    ax2.set_ylim(bottom=0.)

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()

    handles = handles1 + handles2
    labels = labels1 + labels2

    plt.tight_layout()
    fig.legend(handles, labels, loc='upper left',
               bbox_to_anchor=(0.125, 0.925), ncols=2)

    plot_dir = os.path.join(os.path.dirname(__file__), '../Plots')
    plt.savefig(os.path.join(plot_dir, name), bbox_inches='tight', dpi=500)
    plt.close()


def compare_runtimes(files, name='', title='', models=[]):
    merged_df = pd.DataFrame()
    i = 0
    for file in files:
        df = pd.read_json(file)
        df = df.filter(items=['Regions', 'Time'])
        df.rename(columns={'Time': models[i]}, inplace=True)

        if merged_df.empty:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on='Regions', how='outer')
        i = i+1

    merged_df = merged_df.set_index('Regions')
    for column in merged_df.columns:
        plt.plot(merged_df.index, merged_df[column], label=column,
                 linestyle='--', marker='o', linewidth=1.2)
    plt.ylim(bottom=0.)
    plt.xlim(left=merged_df.index.min()-1, right=merged_df.index.max()+1)
    plt.xlabel('Number of regions', fontsize=fontsize_labels)
    plt.ylabel('Run time [seconds]', fontsize=fontsize_labels)
    plt.yticks(fontsize=fontsize_legends)
    plt.xticks(fontsize=fontsize_legends)
    plt.grid(True, linestyle='--')
    plt.legend(fontsize=fontsize_legends)
    plt.title(title)
    plt.tight_layout()

    plot_dir = os.path.join(os.path.dirname(__file__), '../Plots')
    plt.savefig(os.path.join(plot_dir, name), bbox_inches='tight', dpi=500)
    plt.close()


def plot_steps(files, name='', models=[]):
    for idx, file in enumerate(files):
        df = pd.read_json(file)
        df.set_index('Absolute tolerance', inplace=True)
        model_type = os.path.basename(file).split('_')[0]
        if model_type == 'ode':
            plt.plot(df, label=models[idx], color=colors[idx])
        else:
            plt.plot(df['Steps Hotspot'], color=colors[idx], label=models[idx])
            plt.plot(df['Steps other Regions'], color=colors[idx])
            plt.fill_between(
                df.index, df['Steps Hotspot'],
                df['Steps other Regions'],
                color=colors[idx],
                alpha=0.15)

        plt.ylim(bottom=10.)
        plt.xlim(left=df.index.min()/1.2, right=df.index.max()*1.2)
        plt.yticks(fontsize=fontsize_legends)
        plt.xticks(df.index, fontsize=fontsize_legends)

    plt.xscale('log')
    plt.gca().invert_xaxis()
    plt.ylabel('Number of steps', fontsize=fontsize_labels)
    plt.xlabel('Absolute tolerance', fontsize=fontsize_labels)
    plt.grid(True, linestyle='--')
    plt.legend(fontsize=fontsize_legends)

    plot_dir = os.path.join(os.path.dirname(__file__), '../Plots')
    plt.savefig(
        os.path.join(plot_dir, 'compare_steps.png'),
        bbox_inches='tight', dpi=500)
    plt.close()


def plot_quotient(files, name='', title='', models=[]):
    merged_df = pd.DataFrame()
    i = 0
    for file in files:
        df = pd.read_json(file)
        df = df.filter(items=['Regions', 'Time'])
        df.rename(columns={'Time': models[i]}, inplace=True)

        if merged_df.empty:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on='Regions', how='outer')
        i = i+1

    merged_df = merged_df.set_index('Regions')
    plt.plot(
        merged_df.index, merged_df[models[1]] / merged_df[models[0]],
        label='Quotient', linestyle='--', marker='o', linewidth=1.2)
    plt.ylim(bottom=0.)
    plt.xlim(left=merged_df.index.min()-1, right=merged_df.index.max()+1)
    plt.xlabel('Number of regions', fontsize=fontsize_labels)
    plt.ylabel('Run time [seconds]', fontsize=fontsize_labels)
    plt.yticks(fontsize=fontsize_legends)
    plt.xticks(fontsize=fontsize_legends)
    plt.grid(True, linestyle='--')
    plt.legend(fontsize=fontsize_legends)
    plt.title(title)
    plt.tight_layout()

    plot_dir = os.path.join(os.path.dirname(__file__), '../Plots')
    plt.savefig(os.path.join(plot_dir, name), bbox_inches='tight', dpi=500)
    plt.close()


if __name__ == "__main__":
    result_dir = os.path.join(os.path.dirname(__file__), '../results')

    ode_timing_euler = os.path.join(result_dir, 'ode_timing_euler.json')
    ode_timing_noage_euler = os.path.join(
        result_dir, 'ode_timing_noage_euler.json')
    graphbased_timing_euler = os.path.join(
        result_dir, 'graphbased_timing_euler.json')
    graphbased_timing_noage_euler = os.path.join(
        result_dir, 'graphbased_timing_noage_euler.json')

    ode_steps = os.path.join(result_dir, 'ode_steps.json')
    graphbased_steps = os.path.join(result_dir, 'graphbased_steps.json')

    timings_euler = [ode_timing_euler, graphbased_timing_euler]
    timings_euler_noage = [ode_timing_noage_euler,
                           graphbased_timing_noage_euler]

    compare_runtimes(
        timings_euler, name='compare_runtimes_euler', models=models)
    compare_runtimes(
        timings_euler_noage, name='compare_runtimes_euler_noage',
        models=models)
    plot_quotient(timings_euler, name='quotient',
                  title='Quotient', models=models)
    plot_quotient(timings_euler_noage, name='quotient_noage',
                  title='Quotient', models=models)
    compare_runtime_and_flops(
        timings_euler_noage, 'compare_runtimes_and_flops')
    plot_flops('number_flops')
    plot_steps([ode_steps, graphbased_steps], models=models)
