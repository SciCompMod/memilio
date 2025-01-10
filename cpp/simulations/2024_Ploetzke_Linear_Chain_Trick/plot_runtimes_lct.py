#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Lena Ploetzke
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
"""@plot_runtimes_lct.py
Functions to create the plots for section 4.5 of the paper.

The results of time measurements to be plotted should be stored in a folder 'lct_runtimes' in the current directory as 
'.txt' files and should be named as used in the main() function.
Have a look at the README for an explanation of how to measure the run times. 
The run times are contained in the '.out' file of the job. Just copy the whole output in a '.txt' file.
The functionality to clean up the output is in this file.
"""
import os
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import json

colors = ["tab:blue", "tab:orange", "tab:green",
          "tab:red", "tab:purple", "tab:brown"]
fontsize_labels = 16
fontsize_legends = 12
plotfolder = 'Plots/Plots_runtime'


def extract_json_segments(input_file, output_file):
    """ Cleans up the output of the run time measurement using the get_runtimes_lct.sh 
    script and writes the relevant results in a '.json' file.

    @param[in] input_file: Path to the '.txt' file where the information with the run time measurements is stored.
        The relevant data should be contained in brackets {} within the file and should have the same
         format as one json entry. Information in the file which is not contained in brackets will be ignored.
    @param[in] output_file: Path to the '.json' file with the pure run time measurement data. 
            This file will be created in this function.
    """
    # Open file and read data.
    with open(input_file, 'r') as file:
        text = file.read()

    # Extract data in {}.
    segments = re.findall(r'\{(.*?)\}', text, re.DOTALL)

    # Convert segments into JSON objects.
    json_data = []
    for segment in segments:
        json_object = json.loads('{' + segment + '}')
        json_data.append(json_object)

    # Save JSON-objects in JSON file.
    with open(output_file, 'w') as json_file:
        json.dump(json_data, json_file, indent=4)


def plot_runtime(json_file, file_name=''):
    """ Creates a plot visualizing the run time measurements for models with different assumptions
         regarding the number of subcompartments stored in a json file.

    @param[in] json_file: Path to the json file containing run time measurements for models with different assumptions
         regarding the number of subcompartments.
    @param[in] file_name: The name of the file where the plot will be saved.
            If an empty string is provided, the plot will not be saved.
    """
    fig = plt.figure()
    df = pd.read_json(json_file)

    plt.plot(df["Subcompartments"], df["Time"],
             linestyle='--', marker='o', linewidth=1.2)
    plt.ylim(bottom=0.)
    plt.xlim(left=0., right=df["Subcompartments"].max()+1)
    plt.xlabel('Number of subcompartments', fontsize=fontsize_labels)
    plt.ylabel('Run time [seconds]', fontsize=fontsize_labels)
    plt.yticks(fontsize=fontsize_legends)
    plt.xticks(fontsize=fontsize_legends)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if file_name:
        fig.savefig(plotfolder+'/'+file_name+'.png',
                    bbox_inches='tight', dpi=500)
    plt.close()


def plot_runtime_and_steps(json_file, file_name=''):
    """ Creates a plot with two y-axis: One for the run time measurements and one for the required number of time steps
    for models with different assumptions regarding the number of subcompartments.

    @param[in] json_file: Path to the json file containing run time measurements and the time steps for models with 
    different assumptions regarding the number of subcompartments.
    @param[in] file_name: The name of the file where the plot will be saved.
            If an empty string is provided, the plot will not be saved.
    """
    fig, ax1 = plt.subplots()
    df = pd.read_json(json_file)

    # Run time at the left y-axis.
    ax1.plot(df["Subcompartments"], df["Time"],
             linestyle='--', marker='o', linewidth=1.3, color=colors[0], label="Run time")
    # Curve to compare the run time with the shape of a quadratic function.
    ax1.plot(df["Subcompartments"], df["Time"].max()/(df["Subcompartments"].max()**2) *
             df["Subcompartments"]**2, linewidth=1.2, linestyle='--', color='black', label=r"$\mathcal{O}((n_Z)^{2})$")
    ax1.set_xlabel('Number of subcompartments $n_Z$', fontsize=fontsize_labels)
    ax1.set_ylabel('Run time [seconds]',
                   fontsize=fontsize_labels, color=colors[0])
    ax1.set_ylim(bottom=0.)
    ax1.set_xlim(left=0., right=df["Subcompartments"].max()+1)
    ax1.tick_params(axis='y', labelcolor=colors[0], labelsize=fontsize_legends)

    # Second y-axis for steps.
    ax2 = ax1.twinx()
    ax2.plot(df["Subcompartments"], df["Steps"]-1,
             linestyle='--', marker='x', linewidth=1.2, color=colors[1], label="Steps")
    ax2.set_ylabel('Steps', fontsize=fontsize_labels, color=colors[1])
    ax2.tick_params(axis='y', labelcolor=colors[1], labelsize=fontsize_legends)
    ax2.set_ylim(bottom=0., top=df["Steps"].max()+10)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines2 + lines1, labels2 + labels1, fontsize=fontsize_labels-2)

    ax1.tick_params(axis='x', labelsize=fontsize_legends)
    ax1.grid(True, axis='x', linestyle='--', alpha=0.9)
    ax2.grid(True, which='both', axis='both', linestyle='--', alpha=0.9)
    plt.tight_layout()

    if file_name:
        fig.savefig(plotfolder+'/'+file_name+'.png',
                    bbox_inches='tight', dpi=500)
    plt.close()


def main():
    if not os.path.isdir(plotfolder):
        os.makedirs(plotfolder)

    # Simulation results should be stored in this folder.
    data_dir = os.path.join(os.path.dirname(
        __file__), "lct_runtimes")

    # Define which figures of the paper should be created.
    figures = ['left', 'center', 'right']
    if 'left' in figures:
        file_name = 'lct_runtime_subcompartments100'
        paths_to_file = os.path.join(data_dir, file_name)
        extract_json_segments(paths_to_file+'.txt', paths_to_file+'.json')
        plot_runtime(paths_to_file+'.json', file_name)
    if 'center' in figures:
        file_name = 'lct_runtime_subcompartments100_opt0'
        paths_to_file = os.path.join(data_dir, file_name)
        extract_json_segments(paths_to_file+'.txt', paths_to_file+'.json')
        plot_runtime(paths_to_file+'.json', file_name)
    if 'right' in figures:
        file_name = 'lct_runtime_subcompartments200_adaptive'
        paths_to_file = os.path.join(data_dir, file_name)
        extract_json_segments(paths_to_file+'.txt', paths_to_file+'.json')
        plot_runtime_and_steps(paths_to_file+'.json', file_name)


if __name__ == "__main__":
    main()
