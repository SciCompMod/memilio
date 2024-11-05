#############################################################################
# Copyright (C) 2020-2024 MEmilio
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

import os
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import json

colors = ["tab:blue", "tab:orange", "tab:green",
          "tab:red", "tab:purple", "tab:brown"]


def parse_valgrind_output(input_file, output_file):
    # Define relevant pattern.
    refs_pattern = re.compile(r"refs:\s+([0-9,]+)")
    d1_misses_pattern = re.compile(r"D1\s+misses:\s+([0-9,]+)")
    lld_misses_pattern = re.compile(r"LLd\s+misses:\s+([0-9,]+)")
    # List to store data.
    data = []
    # Open file.
    with open(input_file, 'r') as file:
        lines = file.readlines()

    refs = None
    d1_misses = None
    lld_misses = None
    index = 1

    # Go through all lines of the file.
    for line in lines:
        # Look if matching to refs pattern.
        refs_match = refs_pattern.search(line)
        if refs_match:
            # Remove , and store value.
            refs = refs_match.group(1).replace(',', '')

        # Look if matching to D1-Misses pattern.
        d1_misses_match = d1_misses_pattern.search(line)
        if d1_misses_match:
            d1_misses = d1_misses_match.group(1).replace(
                ',', '')

        # Look if matching to LLd-Misses pattern.
        lld_misses_match = lld_misses_pattern.search(line)
        if lld_misses_match:
            lld_misses = lld_misses_match.group(1).replace(
                ',', '')

        # If all relevant data for one index is found, store them and reset all values.
        if refs and d1_misses and lld_misses:
            data.append({
                'Subcompartments': index,
                'refs': refs,
                'D1_misses': d1_misses,
                'LLd_misses': lld_misses
            })
            refs, d1_misses, lld_misses = None, None, None
            index += 1

    # Save as JSON.
    with open(output_file, 'w') as json_file:
        json.dump(data, json_file, indent=4)


def extract_json_segments(input_file, output_file):
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


def plot_runtime_2d(name=''):
    fig = plt.figure()
    df = pd.read_json("lct_runtimes_1day.json")

    for agegroup in [1]:
        df_age = df[(df["Agegroups"] == agegroup)]
        plt.plot(df_age["Subcompartments"], df_age["Time"], linestyle='--', marker='o', linewidth=1.2,
                 label=str(agegroup), color=colors[agegroup-1])
    plt.ylim(bottom=0.)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13,
               title="Number of \nage groups:", title_fontsize=13)
    plt.xlabel('Number of subcompartments', fontsize=13)
    plt.ylabel('Run time [seconds]', fontsize=13)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/run_time_lct'+name+'.png', bbox_inches='tight', dpi=500)


def plot_runtime_2d_noage(jsonfilename, name=''):
    fig = plt.figure()
    df = pd.read_json(jsonfilename)

    df_age = df[(df["Agegroups"] == 1)]
    plt.plot(df_age["Subcompartments"], df_age["Time"],
             linestyle='--', marker='o', linewidth=1.2)
    plt.title(
        'Run time and cache miss rate with one age group for 20 simulation days')
    plt.ylim(bottom=0.)
    plt.xlim(left=0., right=51)
    plt.xlabel('Number of subcompartments', fontsize=13)
    plt.ylabel('Run time [seconds]', fontsize=13)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/run_time_lct'+name+'.png', bbox_inches='tight', dpi=500)


def plot_cachemisses(jsonfilename, savename=''):
    fig = plt.figure()
    df = pd.read_json(jsonfilename)

    plt.plot(df["Subcompartments"], df["LLd_misses"]/df["refs"]*100,
             linestyle='--', marker='o', linewidth=1.2)
    plt.title(
        'LL cache miss rate with one age group for ten simulation days')
    plt.ylim(bottom=0.)
    plt.xlim(left=0., right=51)
    plt.xlabel('Number of subcompartments', fontsize=13)
    plt.ylabel('Cache misses [%]', fontsize=13)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+savename+'.png', bbox_inches='tight', dpi=500)


def main():
    # run times
    jsonfilenameruntimes = 'runtimes/lct_runtimes_20day_onegroup.json'
    # extract_json_segments('runtimes/times_20days.txt', jsonfilenameruntimes)
    plot_runtime_2d_noage(jsonfilenameruntimes, '20days_onegroup')

    # Cache
    jsonfilenamecache = 'runtimes/valgrind_20day_onegroup.json'
    # parse_valgrind_output('runtimes/valgrind_20days.txt',jsonfilenamecache)
    # plot_cachemisses(jsonfilenamecache, 'LLcachemissrate')


if __name__ == "__main__":
    main()
