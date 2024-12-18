import matplotlib.pyplot as plt
import pandas as pd

import os
import json
import re

colors = ["tab:blue", "tab:orange", "tab:green",
          "tab:red", "tab:purple", "tab:brown"]
fontsize_labels = 16
fontsize_legends = 12

models = ['Equation-based model', 'Graph-based model']

def plot_runtime(file, name=''):
    fig = plt.figure()
    df = pd.read_json(file)

    plt.plot(df["Regions"], df["Time"],
             linestyle='--', marker='o', linewidth=1.2)
    plt.ylim(bottom=0.)
    plt.xlim(left=0., right=df["Regions"].max()+1)
    plt.xlabel('Number of regions', fontsize=fontsize_labels)
    plt.ylabel('Run time [seconds]', fontsize=fontsize_labels)
    plt.yticks(fontsize=fontsize_legends)
    plt.xticks(fontsize=fontsize_legends)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    plot_dir = os.path.join(os.path.dirname(__file__), '../Plots')
    name  = os.path.splitext(os.path.basename(file))[0]
    plt.savefig(os.path.join(plot_dir, name), bbox_inches='tight', dpi=500)
    
def compare_runtimes(files, name='', title='', models=[]):
    merged_df = pd.DataFrame()
    i = 0
    for file in files:
        df = pd.read_json(file)

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

    result_equationbased_start = os.path.join(result_dir, 'timing_equationbased_start.json')
    result_equationbased = os.path.join(result_dir, 'timing_equationbased.json')
    result_equationbased_O3 = os.path.join(result_dir, 'timing_equationbased_O3.json')
    result_equationbased_O2 = os.path.join(result_dir, 'timing_equationbased_O2.json')
    result_equationbased_O1 = os.path.join(result_dir, 'timing_equationbased_O1.json')
    result_equationbased_O0 = os.path.join(result_dir, 'timing_equationbased_O0.json')
    result_graphbased_start = os.path.join(result_dir, 'timing_graphbased_start.json')
    result_graphbased = os.path.join(result_dir, 'timing_graphbased.json')
    result_graphbased_smallsteps = os.path.join(result_dir, 'timing_graphbased_01steps.json')
    result_graphbased_unoptimized = os.path.join(result_dir, 'timing_graphbased_unoptimized.json')

    result_equationbased_mod4_0 = os.path.join(result_dir, 'timing_equationbased_mod4_0.json')
    result_equationbased_mod4_1 = os.path.join(result_dir, 'timing_equationbased_mod4_1.json')
    result_equationbased_mod4_2 = os.path.join(result_dir, 'timing_equationbased_mod4_2.json')
    result_equationbased_mod4_3 = os.path.join(result_dir, 'timing_equationbased_mod4_3.json')

    results_start = [result_equationbased_start, result_graphbased_start]
    results = [result_equationbased, result_graphbased, result_graphbased_smallsteps]
    results_unoptimized = [result_equationbased_O3, result_equationbased_O2, result_equationbased_O1, result_equationbased_O0]
    results_mod4 = [result_equationbased_mod4_0, result_equationbased_mod4_1, result_equationbased_mod4_2, result_equationbased_mod4_3]

    # plot_runtime(result_equationbased)
    # plot_runtime(result_graphbased)

    # compare_runtimes(results_start,name='compare_runtimes_start', title='Runtimes for Euler Method', models=models)
    compare_runtimes(results, name='compare_runtimes', title='Runtimes for Euler Method', models=['Equation-based model', 'Graph-based model', 'Graph-based model with dt=0.1'])
    compare_runtimes(results_unoptimized, name='compare_runtimes_unoptimized', title='Runtimes for Euler Method', models=['-O3', '-O2', '-O1', '-O0'])
    compare_runtimes(results_mod4, name='compare_runtimes_mod4', title='Runtimes for Euler Method', models=['%4=0', '%4=1', '%4=2', '%4=3'])
