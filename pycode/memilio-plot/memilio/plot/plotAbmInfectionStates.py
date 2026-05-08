#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Sascha Korf
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

import sys
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py
from datetime import datetime
from scipy.ndimage import gaussian_filter1d


# Module for plotting number of agents per infection state and number of infected agents per location type from ABM results.
# This module provides functions to load and visualize infection states and
# location types from simulation results of the agent-based model (ABM) stored in HDF5 format.

# The used Loggers are:
# LogInfectionStatePerAgeGroup 
# LogInfectionPerLocationTypePerAgeGroup 
# The output of the loggers of several runs is stored in HDF5 files using mio::save_results in mio/io/result_io.h, see abm_history_object.cpp.

# Adjust these as needed.
state_labels = {
    1: 'Exposed',
    2: 'I_Asymp',
    3: 'I_Symp',
    4: 'I_Severe',
    5: 'I_Critical',
    7: 'Dead'
}

age_groups = ['Group1', 'Group2', 'Group3', 'Group4',
              'Group5', 'Group6', 'Total']

age_groups_dict = {
    'Group1': 'Ages 0-4',
    'Group2': 'Ages 5-14',
    'Group3': 'Ages 15-34',
    'Group4': 'Ages 35-59',
    'Group5': 'Ages 60-79',
    'Group6': 'Ages 80+',
    'Total': 'All Ages'
}

location_type_labels = {
    0: 'Home',
    1: 'School',
    2: 'Work',
    3: 'SocialEvent',
    4: 'BasicsShop',
    5: 'Hospital',
    6: 'ICU'
}


def load_h5_results(base_path, percentile):
    """ Reads HDF5 results for a given group and percentile.

    @param[in] base_path Path to results directory.
    @param[in] percentile Subdirectory for percentile (e.g. 'p50').
    @return Dictionary with data arrays. Keys are dataset names from the HDF5 file 
            (e.g., 'Time', 'Total', age group names like 'Group1', 'Group2', etc.).
            Values are numpy arrays containing the corresponding time series data.
    """
    # Try flat structure first (e.g., Results_p50.h5)
    file_path_flat = os.path.join(base_path, f"Results_{percentile}.h5")
    # Try subdirectory structure (e.g., p50/Results.h5)
    file_path_nested = os.path.join(base_path, percentile, "Results.h5")
    
    # Determine which file exists
    if os.path.exists(file_path_flat):
        file_path = file_path_flat
    elif os.path.exists(file_path_nested):
        file_path = file_path_nested
    else:
        raise FileNotFoundError(
            f"Could not find percentile results file. Tried:\n"
            f"  - {file_path_flat}\n"
            f"  - {file_path_nested}\n\n"
            f"Expected directory structure (option 1):\n"
            f"  {base_path}/\n"
            f"    Results_p05.h5\n"
            f"    Results_p25.h5\n"
            f"    Results_p50.h5\n"
            f"    Results_p75.h5\n"
            f"    Results_p95.h5\n\n"
            f"Or (option 2):\n"
            f"  {base_path}/\n"
            f"    p05/Results.h5\n"
            f"    p25/Results.h5\n"
            f"    p50/Results.h5\n"
            f"    p75/Results.h5\n"
            f"    p95/Results.h5\n"
        )
    
    with h5py.File(file_path, 'r') as f:
        data = {k: v[()] for k, v in f['0'].items()}
    return data


def plot_infections_loc_types_average(
        path_to_loc_types,
        start_date='2021-03-01',
        colormap='Set1',
        smooth_sigma=1,
        rolling_window=24,
        xtick_step=150):
    """ Plots rolling sum of new infections per 24 hours location type for the median run.

    @param[in] base_path Path to results directory.
    @param[in] start_date Start date as string.
    @param[in] colormap Matplotlib colormap.
    @param[in] smooth_sigma Sigma for Gaussian smoothing.
    @param[in] rolling_window Window size for rolling sum.
    @param[in] xtick_step Step size for x-axis ticks.
    """
    # Load data
    p50 = load_h5_results(path_to_loc_types, "p50")
    time = p50['Time']
    total_50 = p50['Total']

    plt.figure('Infection_location_types')
    plt.title(
        'Number of new infections per location type for the median run, rolling sum over 24 hours')
    color_plot = matplotlib.colormaps.get_cmap(colormap).colors

    # Check if data dimensions match expected location types
    num_cols = total_50.shape[1] if len(total_50.shape) > 1 else 1
    num_location_types = len(location_type_labels)
    
    if num_cols < num_location_types:
        print(f"Warning: Data has {num_cols} columns but {num_location_types} location types defined.")
        print(f"Only plotting first {num_cols} location types.")
    
    for idx, i in enumerate(location_type_labels.keys()):
        if i >= num_cols:
            break  # Skip if we don't have data for this location type
        color = color_plot[i % len(color_plot)] if i < len(
            color_plot) else "black"
        # Sum up every 24 hours, then smooth
        indexer = pd.api.indexers.FixedForwardWindowIndexer(
            window_size=rolling_window)
        y = pd.DataFrame(total_50[:, i]).rolling(
            window=indexer, min_periods=1).sum().to_numpy()
        y = y[0::rolling_window].flatten()
        y = gaussian_filter1d(y, sigma=smooth_sigma, mode='nearest')
        plt.plot(time[0::rolling_window], y, color=color, linewidth=2.5, label=location_type_labels[i])

    plt.legend()  # Auto-generate legend from plot labels
    _format_x_axis(time, start_date, xtick_step)
    plt.xlabel('Date')
    plt.ylabel('Number of individuals')


def plot_infection_states_results(
        path_to_infection_states,
        start_date='2021-03-01',
        colormap='Set1',
        xtick_step=150,
        show90=False
):
    """ Loads and plots infection state results. 

    @param[in] path_to_infection_states Path to results directory containing infection state data.
    @param[in] start_date Start date as string (YYYY-MM-DD format).
    @param[in] colormap Matplotlib colormap name.
    @param[in] xtick_step Step size for x-axis ticks.
    @param[in] show90 If True, plot 90% percentile (5% and 95%) in addition to 50% percentile.
    """

    # Load data
    p50 = load_h5_results(path_to_infection_states, "p50")
    p25 = load_h5_results(path_to_infection_states, "p25")
    p75 = load_h5_results(path_to_infection_states, "p75")
    time = p50['Time']
    total_50 = p50['Total']
    total_25 = p25['Total']
    total_75 = p75['Total']
    p05 = p95 = None
    total_05 = total_95 = None
    if show90:
        total_95 = load_h5_results(path_to_infection_states, "p95")
        total_05 = load_h5_results(path_to_infection_states, "p05")
        p95 = total_95['Total']
        p05 = total_05['Total']

    plot_infection_states_by_age_group(
        time, p50, p25, p75, colormap,
        p05_bs=total_05 if show90 else None,
        p95_bs=total_95 if show90 else None,
        show90=show90
    )
    plot_infection_states(time, total_50, total_25,
                          total_75, start_date, colormap, xtick_step,
                          y05=p05, y95=p95, show_90=show90)


def plot_infection_states(
        x, y50, y25, y75,
        start_date='2021-03-01',
        colormap='Set1',
        xtick_step=150,
        y05=None, y95=None, show_90=False):
    """ Plots infection states with percentile bands.

    @param[in] x Time array for x-axis.
    @param[in] y50 50th percentile data array.
    @param[in] y25 25th percentile data array.
    @param[in] y75 75th percentile data array.
    @param[in] start_date Start date as string (YYYY-MM-DD format).
    @param[in] colormap Matplotlib colormap name.
    @param[in] xtick_step Step size for x-axis ticks.
    @param[in] y05 5th percentile data array (optional).
    @param[in] y95 95th percentile data array (optional).
    @param[in] show_90 If True, plot 90% percentile bands in addition to 50% percentile.
    """

    plt.figure('Infection_states')

    plt.title('Infection states with 50% percentile')
    if show_90:
        plt.title('Infection states with 50% and 90% percentiles')

    color_plot = matplotlib.colormaps.get_cmap(colormap).colors

    states_plot = list(state_labels.keys())

    for i in states_plot:
        plt.plot(x, y50[:, i], color=color_plot[i],
                 linewidth=2.5, label=state_labels[i])
    # needs to be after the plot calls
    plt.legend([state_labels[i] for i in states_plot])
    for i in states_plot:
        plt.plot(x, y25[:, i], color=color_plot[i],
                 linestyle='dashdot', linewidth=1.2, alpha=0.7)
        plt.plot(x, y75[:, i], color=color_plot[i],
                 linestyle='dashdot', linewidth=1.2, alpha=0.7)
        plt.fill_between(x, y25[:, i], y75[:, i],
                         alpha=0.2, color=color_plot[i])
        # Optional: 90% percentile
        if show_90 and y05 is not None and y95 is not None:
            plt.plot(x, y05[:, i], color=color_plot[i],
                     linestyle='dashdot', linewidth=1.0, alpha=0.4)
            plt.plot(x, y95[:, i], color=color_plot[i],
                     linestyle='dashdot', linewidth=1.0, alpha=0.4)
            plt.fill_between(x, y05[:, i], y95[:, i],
                             # More transparent
                             alpha=0.25, color=color_plot[i])

    _format_x_axis(x, start_date, xtick_step)
    plt.xlabel('Date')
    plt.ylabel('Number of individuals')


def plot_infection_states_by_age_group(
    x, p50_bs, p25_bs, p75_bs, colormap='Set1',
    p05_bs=None, p95_bs=None, show90=False
):
    """ Plots infection states for each age group, with optional 90% percentile. 

    @param[in] x Time array for x-axis.
    @param[in] p50_bs Dictionary containing 50th percentile data for all age groups.
    @param[in] p25_bs Dictionary containing 25th percentile data for all age groups.
    @param[in] p75_bs Dictionary containing 75th percentile data for all age groups.
    @param[in] colormap Matplotlib colormap name.
    @param[in] p05_bs Dictionary containing 5th percentile data for all age groups (optional).
    @param[in] p95_bs Dictionary containing 95th percentile data for all age groups (optional).
    @param[in] show90 If True, plot 90% percentile bands in addition to 50% percentile.
    """

    # Dynamically detect available age groups from the data
    available_groups = [key for key in p50_bs.keys() if key.startswith('Group') or key == 'Total']
    # Sort to ensure Group1, Group2, ..., Total order
    available_groups = sorted(available_groups, key=lambda x: (x != 'Total', x))
    
    color_plot = matplotlib.colormaps.get_cmap(colormap).colors
    n_states = len(state_labels)
    fig, ax = plt.subplots(
        n_states, len(available_groups), constrained_layout=True, figsize=(20, 3 * n_states))

    for col_idx, group in enumerate(available_groups):
        y50 = p50_bs[group]
        y25 = p25_bs[group]
        y75 = p75_bs[group]
        y05 = p05_bs[group] if (show90 and p05_bs is not None) else None
        y95 = p95_bs[group] if (show90 and p95_bs is not None) else None
        
        # Get group label, using default if not in predefined dict
        group_label = age_groups_dict.get(group, group)
        
        for row_idx, (state_idx, label) in enumerate(state_labels.items()):
            _plot_state(
                ax[row_idx, col_idx], x, y50[:, state_idx], y25[:,
                                                                state_idx], y75[:, state_idx],
                color_plot[col_idx], f'#{label}, {group_label}',
                y05=y05[:, state_idx] if y05 is not None else None,
                y95=y95[:, state_idx] if y95 is not None else None,
                show90=show90
            )
            # The legend should say: solid line = median, dashed line = 25% and 75% perc. and if show90 is True, dotted line = 5%, 25%, 75%, 95% perc.
            perc_string = '25/75%' if not show90 else '5/25/75/95%'
            ax[row_idx, col_idx].legend(
                ['Median', f'{perc_string} perc.'],
                loc='upper left', fontsize=8)

            # Add y label for leftmost column
            if col_idx == 0:
                ax[row_idx, col_idx].set_ylabel('Number of individuals')

            # Add x label for bottom row
            if row_idx == n_states - 1:
                ax[row_idx, col_idx].set_xlabel('Time (days)')

    string_short = ' and 90%' if show90 else ''
    fig.suptitle(
        'Infection states per age group with 50%' + string_short + ' percentile',
        fontsize=16)


def _plot_state(ax, x, y50, y25, y75, color, title, y05=None, y95=None, show90=False):
    """ Helper to plot a single state with fill_between and optional 90% percentile. """
    ax.plot(x, y50, color=color, label='Median')
    ax.fill_between(x, y25, y75, alpha=0.5, color=color)
    if show90 and y05 is not None and y95 is not None:
        ax.plot(x, y05, color=color, linestyle='dotted',
                linewidth=1.0, alpha=0.4)
        ax.plot(x, y95, color=color, linestyle='dotted',
                linewidth=1.0, alpha=0.4)
        ax.fill_between(x, y05, y95, alpha=0.15, color=color)
    ax.tick_params(axis='y')
    ax.set_title(title)


def _format_x_axis(x, start_date, xtick_step):
    """ Helper to format x-axis as dates. """
    start = datetime.strptime(start_date, '%Y-%m-%d')
    xx = [start + pd.Timedelta(days=int(i)) for i in x]
    xx_str = [dt.strftime('%Y-%m-%d') for dt in xx]
    plt.gca().set_xticks(x[::xtick_step])
    plt.gca().set_xticklabels(xx_str[::xtick_step])
    plt.gcf().autofmt_xdate()


def main():
    """ Main function for CLI usage. """
    parser = argparse.ArgumentParser(
        description="Plot infection state and location type results.")
    parser.add_argument("--path-to-infection-states",
                        help="Path to infection states results")
    parser.add_argument("--path-to-loc-types",
                        help="Path to location types results")
    parser.add_argument("--start-date", type=str, default='2021-03-01',
                        help="Simulation start date (YYYY-MM-DD)")
    parser.add_argument("--colormap", type=str,
                        default='Set1', help="Matplotlib colormap")
    parser.add_argument("--xtick-step", type=int,
                        default=150, help="Step for x-axis ticks (usually hours)")
    parser.add_argument("--90percentile", action="store_true",
                        help="If set, plot 90% percentile as well")
    args = parser.parse_args()

    if not args.path_to_infection_states and not args.path_to_loc_types:
        print("Please provide a path to infection states or location types results.")
        return

    if args.path_to_infection_states:
        plot_infection_states_results(
            path_to_infection_states=args.path_to_infection_states,
            start_date=args.start_date,
            colormap=args.colormap,
            xtick_step=args.xtick_step,
            show90=True
        )
    
    if args.path_to_loc_types:
        plot_infections_loc_types_average(
            path_to_loc_types=args.path_to_loc_types,
            start_date=args.start_date,
            colormap=args.colormap,
            xtick_step=args.xtick_step)

    plt.show()


if __name__ == "__main__":
    main()
