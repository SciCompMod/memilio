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

# The used  Loggers are:
# struct LogInfectionStatePerAgeGroup : mio::LogAlways {
#     using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
#     /**
#      * @brief Log the TimeSeries of the number of Person%s in an #InfectionState for every age group.
#      * @param[in] sim The simulation of the abm.
#      * @return A pair of the TimePoint and the TimeSeries of the number of Person%s in an #InfectionState for every age group.
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
#             auto index = (((size_t)(mio::abm::InfectionState::Count)) * ((uint32_t)p.get_age().get())) +
#                 ((uint32_t)p.get_infection_state(curr_time));
# // PRAGMA_OMP(atomic)
#              sum[index] += 1;
#         }
#         return std::make_pair(curr_time, sum);
#     }
# };
#
# struct LogInfectionPerLocationTypePerAgeGroup : mio::LogAlways {
#     using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
#     /**
#      * @brief Log the TimeSeries of the number of newly infected Person%s for each Location Type and each age.
#      * @param[in] sim The simulation of the abm.
#      * @return A pair of the TimePoint and the TimeSeries of newly infected Person%s for each Location Type and each age.
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
#                 // PRAGMA_OMP(atomic)
#                 if ((p.get_infection_state(prev_time) != mio::abm::InfectionState::Exposed) &&
#                     (p.get_infection_state(curr_time) == mio::abm::InfectionState::Exposed)) {
#                     auto index = (((size_t)(mio::abm::LocationType::Count)) * ((uint32_t)p.get_age().get())) +
#                                  ((uint32_t)p.get_location().get_type());
#                     sum[index] += 1;
#                 }
#         }
#         return std::make_pair(curr_time, sum);
#     }
# };
#
# The output of the loggers of several runs is stored in HDF5 files using mio::save_results in mio/io/result_io.h.

# Adjust these as needed.
state_labels = {
    2: 'Exposed',
    3: 'I_Asymp',
    4: 'I_Symp',
    5: 'I_Severe',
    6: 'I_Critical',
    7: 'Recovered',
    8: 'Dead'
}

state_groups = ['Group1', 'Group2', 'Group3', 'Group4',
                'Group5', 'Group6', 'Group7', 'Group8', 'Total']

color_plot = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']


def load_h5_results(base_path, percentile):
    """ Reads HDF5 results for a given group and percentile.

    @param[in] base_path Path to results directory.
    @param[in] percentile Subdirectory for percentile (e.g. 'p50').
    @return Dictionary with data arrays. Keys are dataset names from the HDF5 file
            (e.g., 'Time', 'Total', age group names like 'Group1', 'Group2', etc.).
            Values are numpy arrays containing the corresponding time series data.
    """
    file_path = os.path.join(base_path, "Results_" + percentile + ".h5")
    with h5py.File(file_path, 'r') as f:
        data = {k: v[()] for k, v in f['0'].items()}
    return data


def plot_infection_states(
        x, y50, y25, y75,
        path_to_infection_states,
        colormap='Set1',
        xtick_step=150,
        y05=None, y95=None, show_90=False):
    """ Plots infection states with percentile bands.

    @param[in] x Time array for x-axis.
    @param[in] y50 50th percentile data array.
    @param[in] y25 25th percentile data array.
    @param[in] y75 75th percentile data array.
    @param[in] path_to_infection_states Path to infection states directory.
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

    for i, group in enumerate(state_groups):
        if group != 'Total' and group != 'Group1':
            plt.plot(x, y50[group].flatten(), color=color_plot[i+1],
                     linewidth=2.5, label=state_labels[i+1])
    # needs to be after the plot calls
    plt.legend([state_labels[i] for i in states_plot])
    for i, group in enumerate(state_groups):
        if group != 'Total' and group != 'Group1':
            plt.plot(x, y25[group].flatten(), color=color_plot[i+1],
                     linestyle='dashdot', linewidth=1.2, alpha=0.7)
            plt.plot(x, y75[group].flatten(), color=color_plot[i+1],
                     linestyle='dashdot', linewidth=1.2, alpha=0.7)
            plt.fill_between(x, y25[group].flatten(), y75[group].flatten(),
                             alpha=0.2, color=color_plot[i+1])
            # Optional: 90% percentile
            if show_90 and y05 is not None and y95 is not None:
                plt.plot(x, y05[group].flatten(), color=color_plot[i+1],
                         linestyle='dashdot', linewidth=1.0, alpha=0.4)
                plt.plot(x, y95[group].flatten(), color=color_plot[i+1],
                         linestyle='dashdot', linewidth=1.0, alpha=0.4)
                plt.fill_between(x, y05[group].flatten(), y95[group].flatten(),
                                 # More transparent
                                 alpha=0.25, color=color_plot[i+1])

    _format_x_axis_days(x, xtick_step)
    plt.xlabel('Days')
    plt.ylabel('Number of individuals')

    # Save the plot to parent folder
    output_dir = os.path.dirname(path_to_infection_states)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "infection_states.png")
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")
    plt.close()


def plot_infection_states_results(
        path_to_infection_states,
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
    p05 = p95 = None
    total_05 = total_95 = None
    if show90:
        p95 = load_h5_results(path_to_infection_states, "p95")
        p05 = load_h5_results(path_to_infection_states, "p05")

    plot_infection_states(time, p50, p25,
                          p75, path_to_infection_states, colormap, xtick_step,
                          y05=p05, y95=p95, show_90=show90)


def _format_x_axis_days(x, xtick_step):
    """ Helper to format x-axis as days. """
    plt.gca().set_xticks(x[::xtick_step])
    plt.gca().set_xticklabels([int(day) for day in x[::xtick_step]])


def main():
    """ Main function for CLI usage. """
    parser = argparse.ArgumentParser(
        description="Plot infection state and location type results.")
    parser.add_argument("--path-to-infection-states",
                        help="Path to infection states results")
    parser.add_argument("--path-to-loc-types",
                        help="Path to location types results")
    parser.add_argument("--path-to-infection-states-2",
                        help="Path to second infection states results for comparison")
    parser.add_argument("--label1", type=str, default="Simulation 1",
                        help="Label for first simulation")
    parser.add_argument("--label2", type=str, default="Simulation 2",
                        help="Label for second simulation")
    parser.add_argument("--output-path", type=str, default=".",
                        help="Output path for comparison plots")
    parser.add_argument("--start-date", type=str, default='2021-03-01',
                        help="Simulation start date (YYYY-MM-DD)")
    parser.add_argument("--colormap", type=str,
                        default='Set1', help="Matplotlib colormap")
    parser.add_argument("--xtick-step", type=int,
                        default=1, help="Step for x-axis ticks (usually hours)")
    parser.add_argument("--s90percentile", action="store_false",
                        help="If set, plot 90% percentile as well")
    args = parser.parse_args()

    args.path_to_infection_states = "/Users/saschakorf/Nosynch/Arbeit/memilio/example_results"

    print("Arguments received:")
    print(f"Path to infection states: {args.path_to_infection_states}")
    print(f"Path to infection states 2: {args.path_to_infection_states_2}")
    print(f"Path to location types: {args.path_to_loc_types}")

    # Check if comparison mode is requested
    plot_infection_states_results(
        path_to_infection_states=args.path_to_infection_states,
        colormap=args.colormap,
        xtick_step=args.xtick_step,
        show90=args.s90percentile
    )

    if not args.path_to_infection_states and not args.path_to_loc_types:
        print("Please provide a path to infection states or location types results.")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("No arguments provided. Running in interactive mode.")
        main()
    else:
        print("Running in CLI mode with provided arguments.")
        main()
        sys.exit(0)
    sys.exit(1)
