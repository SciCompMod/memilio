#############################################################################
# Copyright (C) 2020-2024 German Aerospace Center (DLR-SC)
#
# Authors: Anna Wendler, Lena Ploetzke
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
import h5py
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from memilio.epidata import getDataIntoPandasDataFrame as gd


def plot_changepoint(files, fileending="", save_dir=""):
    """
    Plots the result of the changepoint simulation.

    @param[in] files Expects list of two files with ODE and IDE simulation results for flows, respectively, in this order. 
    @param[in] fileending Determines file ending of saved plot.
    @param[in] save_dir Directory where plot will be stored.
    """

    fig, ax = plt.subplots()
    legendplot = list(["ODE", "IDE"])

    # Define colors, we use helmholtzdarkblue, helmholtzclaim.
    colors = [(0, 40/255, 100/255), (20/255, 200/255, 255/255)]
    linestyles = ['-', '--']
    # Add results to plot-
    for file in range(len(files)):
        # Load data.
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]

        
        # As there should be only one Group, total is the simulation result.
        total = data['Total'][:, :]

        dates = data['Time'][:]

        timestep = np.diff(dates)[0]
        tmax = dates[-1]

        # Get indices where dates are >=0.
        indices = np.where(dates >= 0)

        # Plot data.        
        # ODE
        if file == 0:
            # Transform cumulative flows to absolute flows,
            # then transform from flows over time interval to flows at time points.
            ax.plot(dates[indices[0][1:]], np.diff(total[indices[0], 0])/np.diff(dates[indices[0]]), label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])

        # IDE
        elif file == 1:
            # Transform from flows over time interval to flows at time points.
            ax.plot(dates[1:], total[1:, 0]/np.diff(dates), label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])
                
        h5file.close()

        ax.set_xlim(left=0, right=tmax)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(fontsize=12)

    fig.supxlabel('Simulation time [days]')
    fig.supylabel('Daily new transmissions')
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=None, hspace=0.6)

    plt.tight_layout()

    # Save result.
    if save_dir != "":
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f"changepoint_{fileending}.png",
                    bbox_inches='tight', dpi=500)


if __name__ == '__main__':
    # Paths are valid if script is executed e.g. in memilio/cpp/simulations/IDE_paper
    # Path where simulation results (generated with ide_changepoints.cpp) are stored. 
    result_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/simulation_results/changepoints/")
    # Path where plots will be stored. 
    save_dir =  os.path.join(os.path.dirname(
        __file__), "../../..", "data/plots/changepoints/")

    plot_changepoint([os.path.join(result_dir, f"changepoint_ode_0.5_12_0.0100_flows"),
                     os.path.join(result_dir, f"changepoint_ide_0.5_12_0.0100_flows")],
                     fileending="0.5_12_0.0100", save_dir=save_dir)

    plot_changepoint([os.path.join(result_dir, f"changepoint_ode_2.0_12_0.0100_flows"),
                     os.path.join(result_dir, f"changepoint_ide_2.0_12_0.0100_flows")],
                     fileending="2.0_12_0.0100", save_dir=save_dir)
