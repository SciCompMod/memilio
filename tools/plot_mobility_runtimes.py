import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import os

colors = ['#1f77b4', '#2ca02c', '#ff7f0e']
linestyles=['-', '--', '-.', ':']
fontsize_labels = 16
fontsize_legends = 12

models = ['Equation-based model', 'Graph-based model']

def plot_runtime(file, name=''):
    fig = plt.figure()
    df = pd.read_json(file)

    plt.plot(df["Regions"], df["Time"],
             linestyle='--', marker='o', linewidth=1.2)
    plt.ylim(bottom=df['Time'].min())
    plt.xlim(left=df["Regions"].min()-1, right=df["Regions"].max()+1)
    plt.xlabel('Number of regions', fontsize=fontsize_labels)
    plt.ylabel('Run time [seconds]', fontsize=fontsize_labels)
    plt.yticks(fontsize=fontsize_legends)
    plt.xticks(fontsize=fontsize_legends)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    plot_dir = os.path.join(os.path.dirname(__file__), '../Plots')
    if name is None:
        name = os.path.splitext(os.path.basename(file))[0]
    plt.savefig(os.path.join(plot_dir, name), bbox_inches='tight', dpi=500)
    plt.close()

def plot_flops(name='number_flops'):
    fig, ax = plt.subplots()

    def flops_equation_based(x, eta):
            return (4*x**2+22*x+1)/eta
    
    def flops_graph_based(x, eta):
            return (43*x**2+24*x/eta+2)*1
    
    x = np.linspace(0, 400, 80)

    
    for idx, eta in enumerate([0.05, 0.1, 0.2, 0.5]):
        ax.plot(x, flops_equation_based(x, eta), linewidth=1.5, color=colors[0], linestyle=linestyles[idx], label='Model C, $\eta=$'+ str(eta))
        ax.plot(x, flops_graph_based(x, eta), linewidth=1.5, color=colors[1], linestyle=linestyles[idx], label='Model D, $\eta=$'+ str(eta))
    ax.set_ylim(bottom=0.)
    ax.set_xlim(left=0., right=400.)
    ax.set_xlabel('Number of regions', fontsize=fontsize_labels)
    ax.set_ylabel('Number of FLOPs', fontsize=fontsize_labels)

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
    
    for file in files:
        df = pd.read_json(file)

        ax1.plot(df["Regions"], df["Time"],
                linestyle='--', marker='o', linewidth=1.2, label=file)

    ax1.set_ylim(bottom=0.)
    ax1.set_xlim(left=0., right=400.)
    ax1.set_xlabel('Number of regions', fontsize=fontsize_labels)
    ax1.set_ylabel('Run time [seconds]', fontsize=fontsize_labels)

    ax2 = ax1.twinx()

    def flops_equation_based(x):
            return (4*x**2+22*x+1)*200
    
    def flops_graph_based(x):
            return (43*x**2+240*x+2)*20
    
    x = np.linspace(0, 400, 400)
    
    ax2.plot(x, flops_equation_based(x), linestyle='--', linewidth=1.2)
    ax2.plot(x, flops_graph_based(x), linestyle='--', linewidth=1.2)
    ax2.set_ylabel('Number of FLOPs', fontsize=fontsize_labels)
    ax2.set_ylim(bottom=0.)

    plt.tight_layout()

    plot_dir = os.path.join(os.path.dirname(__file__), '../Plots')
    plt.savefig(os.path.join(plot_dir, name), bbox_inches='tight', dpi=500)
    plt.close()
    
def compare_runtimes(files, name='', title='', models=[]):
    merged_df = pd.DataFrame()
    i = 0
    for file in files:
        df = pd.read_json(file)
        df = df.filter(items=['Regions', 'Time'])
        # df.drop(thisFilter, inplace=True, axis=1)
        df.rename(columns={'Time': models[i]}, inplace=True)

        if merged_df.empty:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on='Regions', how='outer')
        i = i+1

    merged_df = merged_df.set_index('Regions')
    for column in merged_df.columns:
    #     plt.plot(merged_df['Regions'], column,
    #             linestyle='--', marker='o', linewidth=1.2)
        plt.plot(merged_df.index, merged_df[column], label=column,
             linestyle='--', marker='o', linewidth=1.2)
    plt.ylim(bottom=0.)
    plt.xlim(left=merged_df.index.min()-1, right=merged_df.index.max()+1)
    plt.xlabel('Number of regions', fontsize=fontsize_labels)
    plt.ylabel('Run time [seconds]', fontsize=fontsize_labels)
    plt.yticks(fontsize=fontsize_legends)
    plt.xticks(fontsize=fontsize_legends)
    plt.grid(True, linestyle='--')
    plt.legend()
    plt.title(title)
    plt.tight_layout()

    plot_dir = os.path.join(os.path.dirname(__file__), '../Plots')
    plt.savefig(os.path.join(plot_dir, name), bbox_inches='tight', dpi=500)
    plt.close()

if __name__ == "__main__":
    result_dir = os.path.join(os.path.dirname(__file__), '../results')

    result_equationbased_euler = os.path.join(result_dir, 'timing_equationbased_euler.json')
    result_equationbased_noage_euler = os.path.join(result_dir, 'timing_equationbased_noage_euler.json')
    result_graphbased_euler = os.path.join(result_dir, 'timing_graphbased_euler.json')
    result_graphbased_noage_euler = os.path.join(result_dir, 'timing_graphbased_noage_euler.json')

    results_euler = [result_equationbased_euler, result_graphbased_euler]
    results_euler_noage = [result_equationbased_noage_euler, result_graphbased_noage_euler]

    # compare_runtimes(results_euler, name='compare_runtimes_euler', models=models)
    # compare_runtimes(results_euler_noage, name='compare_runtimes_euler_noage', models=models)
    # compare_runtime_and_flops(results_euler_noage, 'compare_runtimes_and_flops')
    plot_flops()
