#############################################################################
# Copyright (C) 2020-2024 MEmilio
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


# Module for plotting infection states and location types from ABM results.
# This module provides functions to load and visualize infection states and
# location types from simulation results stored in HDF5 format and are output
# by the MEmilio agent-based model (ABM).
# The used  Loggers are:
# struct LogInfectionStatePerAgeGroup : mio::LogAlways {
#     using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
#     /** 
#      * @brief Log the TimeSeries of the number of Person%s in an #InfectionState.
#      * @param[in] sim The simulation of the abm.
#      * @return A pair of the TimePoint and the TimeSeries of the number of Person%s in an #InfectionState.
#      */
#     static Type log(const mio::abm::Simulation& sim)
#     {
# 
#         Eigen::VectorXd sum = Eigen::VectorXd::Zero(
#             Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups()));
#         const auto curr_time = sim.get_time();
#         const auto persons   = sim.get_world().get_persons();
# 
#         // PRAGMA_OMP(parallel for)
#         for (auto i = size_t(0); i < persons.size(); ++i) {
#             auto& p = persons[i];
#             if (p.get_should_be_logged()) {
#                 auto index = (((size_t)(mio::abm::InfectionState::Count)) * ((uint32_t)p.get_age().get())) +
#                              ((uint32_t)p.get_infection_state(curr_time));
#                 // PRAGMA_OMP(atomic)
#                 sum[index] += 1;
#             }
#         }
#         return std::make_pair(curr_time, sum);
#     }
# };
# 
# struct LogInfectionPerLocationTypePerAgeGroup : mio::LogAlways {
#     using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
#     /** 
#      * @brief Log the TimeSeries of the number of Person%s in an #InfectionState.
#      * @param[in] sim The simulation of the abm.
#      * @return A pair of the TimePoint and the TimeSeries of the number of Person%s in an #InfectionState.
#      */
#     static Type log(const mio::abm::Simulation& sim)
#     {
# 
#         Eigen::VectorXd sum = Eigen::VectorXd::Zero(
#             Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups()));
#         auto curr_time     = sim.get_time();
#         auto prev_time     = sim.get_prev_time();
#         const auto persons = sim.get_world().get_persons();
# 
#         // PRAGMA_OMP(parallel for)
#         for (auto i = size_t(0); i < persons.size(); ++i) {
#             auto& p = persons[i];
#             if (p.get_should_be_logged()) {
#                 // PRAGMA_OMP(atomic)
#                 if ((p.get_infection_state(prev_time) != mio::abm::InfectionState::Exposed) &&
#                     (p.get_infection_state(curr_time) == mio::abm::InfectionState::Exposed)) {
#                     auto index = (((size_t)(mio::abm::LocationType::Count)) * ((uint32_t)p.get_age().get())) +
#                                  ((uint32_t)p.get_location().get_type());
#                     sum[index] += 1;
#                 }
#             }
#         }
#         return std::make_pair(curr_time, sum);
#     }
# };
# 
# The output of the loggers of several runs is stored in HDF5 files, with the memilio funciton mio::save_results in mio/io/result_io.h.


def load_h5_results(base_path, percentile):
    """ Reads HDF5 results for a given group and percentile.

    @param[in] base_path Path to results directory.
    @param[in] percentile Subdirectory for percentile (e.g. 'p50').
    @return Dictionary with data arrays.
    """
    file_path = os.path.join(base_path, percentile, "Results.h5")
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
    """ Plots rolling average infections per location type for the median run.

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
    plt.title('Infection per location type for the median run, rolling sum over 24 hours')
    color_plot = matplotlib.colormaps.get_cmap(colormap).colors
    # If you define further location types, you need to adjust this list
    states_plot = [0, 1, 2, 3, 4, 5 , 6] 
    legend_plot = ['Home', 'School', 'Work', 'SocialEvent', 'BasicsShop', 'Hospital', 'ICU']

    for idx, i in enumerate(states_plot):
        color = color_plot[i % len(color_plot)] if i < len(color_plot) else "black"
        # Sum up every 24 hours, then smooth
        indexer = pd.api.indexers.FixedForwardWindowIndexer(window_size=rolling_window)
        y = pd.DataFrame(total_50[:, i]).rolling(window=indexer, min_periods=1).sum().to_numpy()
        y = y[0::rolling_window].flatten()
        y = gaussian_filter1d(y, sigma=smooth_sigma, mode='nearest')
        plt.plot(time[0::rolling_window], y, color=color)

    plt.legend(legend_plot)
    _format_x_axis(time, start_date, xtick_step)
    plt.xlabel('Date')
    plt.ylabel('Number of individuals')
    plt.show()

def plot_infection_states_results(
        path_to_infection_states,
        start_date='2021-03-01',
        colormap='Set1',
        xtick_step=150):
    """ Loads and plots infection state results.

    @param[in] base_path Path to results directory.
    @param[in] start_date Start date as string.
    @param[in] colormap Matplotlib colormap.
    @param[in] xtick_step Step size for x-axis ticks.
    """
    # Load data
    p50 = load_h5_results(path_to_infection_states, "p50")
    p25 = load_h5_results(path_to_infection_states, "p25")
    p75 = load_h5_results(path_to_infection_states, "p75")
    time = p50['Time']
    total_50 = p50['Total']
    total_25 = p25['Total']
    total_75 = p75['Total']

    plot_infection_states_individual(time, p50, p25, p75, colormap)
    plot_infection_states(time, total_50, total_25, total_75, start_date, colormap, xtick_step)

def plot_infection_states(
        x, y50, y25, y75,
        start_date='2021-03-01',
        colormap='Set1',
        xtick_step=150):
    """ Plots infection states with percentiles.

    @param[in] x Time array.
    @param[in] y50 Median values.
    @param[in] y25 25th percentile values.
    @param[in] y75 75th percentile values.
    @param[in] start_date Start date as string.
    @param[in] colormap Matplotlib colormap.
    @param[in] xtick_step Step size for x-axis ticks.
    """
    plt.figure('Infection_states with 50% percentile')
    plt.title('Infection states with 50% percentile')
    color_plot = matplotlib.colormaps.get_cmap(colormap).colors
    states_plot = [1, 2, 3, 4, 5, 7]
    legend_plot = ['E', 'I_NSymp', 'I_Symp', 'I_Sev', 'I_Crit', 'Dead']

    for i in states_plot:
        plt.plot(x, y50[:, i], color=color_plot[i])
    plt.legend(legend_plot) # Needs to be done here, otherwise the percentage fill_between will not work correctly

    for i in states_plot:
        plt.fill_between(x, y50[:, i], y25[:, i], alpha=0.5, color=color_plot[i])
        plt.fill_between(x, y50[:, i], y75[:, i], alpha=0.5, color=color_plot[i])

    plt.legend(legend_plot)
    _format_x_axis(x, start_date, xtick_step)
    plt.xlabel('Time')
    plt.ylabel('Number of individuals')
    plt.show()

def plot_infection_states_individual(x, p50_bs, p25_bs, p75_bs, colormap='Set1'):
    """ Plots infection states for each age group.

    @param[in] x Time array.
    @param[in] p50_bs Median values by group.
    @param[in] p25_bs 25th percentile values by group.
    @param[in] p75_bs 75th percentile values by group.
    @param[in] colormap Matplotlib colormap.
    """
    age_groups = ['Group1', 'Group2', 'Group3', 'Group4', 'Group5', 'Group6', 'Total'] # Adjust as needed
    color_plot = matplotlib.colormaps.get_cmap(colormap).colors
    fig, ax = plt.subplots(6, len(age_groups), constrained_layout=True, figsize=(20, 9))

    for col_idx, group in enumerate(age_groups):
        y50 = p50_bs[group]
        y25 = p25_bs[group]
        y75 = p75_bs[group]
        # Infected no symptoms
        _plot_state(ax[0, col_idx], x, y50[:, 1], y25[:, 1], y75[:, 1], color_plot[col_idx], '#Infected_no_symptoms, Age' + str(group))
        # Infected symptoms
        _plot_state(ax[1, col_idx], x, y50[:, 2], y25[:, 2], y75[:, 2], color_plot[col_idx], '#Infected_symptoms, Age' + str(group))
        # Severe
        _plot_state(ax[2, col_idx], x, y50[:, 4], y25[:, 4], y75[:, 4], color_plot[col_idx], '#Severe, Age' + str(group))
        # Critical
        _plot_state(ax[3, col_idx], x, y50[:, 5], y25[:, 5], y75[:, 5], color_plot[col_idx], '#Critical, Age' + str(group))
        # Dead
        _plot_state(ax[4, col_idx], x, y50[:, 7], y25[:, 7], y75[:, 7], color_plot[col_idx], '#Dead, Age' + str(group))
        # Recovered
        _plot_state(ax[5, col_idx], x, y50[:, 6], y25[:, 6], y75[:, 6], color_plot[col_idx], '#Recovered, Age' + str(group))
    
    fig.suptitle('Infection states per age group with 50% percentile', fontsize=16)

    # We hide the Legend for the individual plots as it is too cluttered
    for ax_row in ax:
        for ax_col in ax_row:
            ax_col.legend().set_visible(False)

    plt.show()

def _plot_state(ax, x, y50, y25, y75, color, title):
    """ Helper to plot a single state with fill_between. """
    ax.set_xlabel('time (days)')
    ax.plot(x, y50, color=color, label='Median')
    ax.fill_between(x, y50, y25, alpha=0.5, color=color)
    ax.fill_between(x, y50, y75, alpha=0.5, color=color)
    ax.tick_params(axis='y')
    ax.set_title(title)
    ax.legend(['Simulation'])

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
    parser = argparse.ArgumentParser(description="Plot infection state and location type results.")
    parser.add_argument("--path-to-infection-states", help="Path to infection states results")
    parser.add_argument("--path-to-loc-types", help="Path to location types results")
    parser.add_argument("--start-date", type=str, default='2021-03-01', help="Simulation start date (YYYY-MM-DD)")
    parser.add_argument("--colormap", type=str, default='Set1', help="Matplotlib colormap")
    parser.add_argument("--xtick-step", type=int, default=150, help="Step for x-axis ticks")
    args = parser.parse_args()

    if  args.path_to_infection_states:
        plot_infection_states_results(
            args.path_to_infection_states,
            start_date=args.start_date,
            colormap=args.colormap,
            xtick_step=args.xtick_step)
    if args.path_to_loc_types:
        plot_infections_loc_types_average(
            args.path_to_loc_types,
            start_date=args.start_date,
            colormap=args.colormap,
            xtick_step=args.xtick_step)
        
    if not args.path_to_infection_states and not args.path_to_loc_types:
        print("Please provide a path to infection states or location types results.")
        sys.exit(1)
    plt.show()

    
if __name__ == "__main__":
    main()
