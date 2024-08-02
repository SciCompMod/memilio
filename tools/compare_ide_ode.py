#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
"""@plot_results_lct_secir.py
functions to plot results of a simulation with a LCT SECIR model with subcompartments.
There is also a method to compare different results of different models.
The data to be plotted should be stored in a '../data/simulation_lct' folder as .h5 files.
Data could be generated eg by executing the file ./cpp/examples/lct_secir_compare.cpp.
"""

import h5py
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from memilio.epidata import getDataIntoPandasDataFrame as gd

# Define compartments
secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
              5: 'ICU', 6: 'Recovered', 7: 'Dead'}


def compare_results(files, legendplot, flows=True, fileending="", save=True, save_dir=None):
    """ Creates a 4x2 Plot with one subplot per compartment and one line per result one wants to compare.
    @param[in] files: paths of the files (without file extension .h5) with the simulation results that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] save: if save is True, the plot is saved in a folder named Plots.
    """
    if flows:

        secir_dict = {0: r"$\sigma_S^E$", 1: r"$\sigma_E^C$", 2: r"$\sigma_C^I$", 3: r"$\sigma_C^R$", 4: r"$\sigma_I^H$",
                      5: r"$\sigma_I^R$", 6: r"$\sigma_H^U$", 7: r"$\sigma_H^R$", 8: r"$\sigma_U^D$", 9: r"$\sigma_U^R$"}

        fig, axs = plt.subplots(5, 2, sharex='all', num='Compare files')
        num_plots = 10

    else:

        # Define compartments
        secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
                      5: 'ICU', 6: 'Recovered', 7: 'Dead'}

        fig, axs = plt.subplots(4, 2, sharex='all', num='Compare files')
        num_plots = 8

    # helmholtzdarkblue, helmholtzclaim
    colors = [(0, 40/255, 100/255), (20/255, 200/255, 255/255)]
    linestyles = ['-', '--']
    # add results to plot
    for file in range(len(files)):
        # load data
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]

        if flows:
            # As there should be only one Group, total is the simulation result
            if len(data['Total'][0]) == 10:
                total = data['Total'][:, :]
            elif len(data['Total'][0]) == 15:
                # in ODE: two compartments more are used which results in more flows;
                # throw out additional flows
                total = data['Total'][:, [0, 1, 2, 3, 6, 7, 10, 11, 13, 14]]
        else:
            if len(data['Total'][0]) == 8:
                # As there should be only one Group, total is the simulation result
                total = data['Total'][:, :]
            elif len(data['Total'][0]) == 10:
                # in ODE there are two compartments we don't use, throw these out
                total = data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]]

        dates = data['Time'][:]

        # plot data
        if flows:
            # ODE
            if file == 0:
                for i in range(num_plots):
                    axs[int(i/2), i % 2].plot(dates[1:],
                                              np.diff(total[:, i])/np.diff(dates), label=legendplot[file], color=colors[file], linestyle=linestyles[file])
            # IDE
            elif file == 1:
                for i in range(num_plots):
                    axs[int(i/2), i % 2].plot(dates[1:],
                                              total[1:, i]/np.diff(dates), label=legendplot[file], color=colors[file], linestyle=linestyles[file])
                min_date = dates[1]
        else:
            if file == 0:
                for i in range(num_plots):
                    axs[int(i/2), i % 2].plot(dates,
                                              total[:, i], label=legendplot[file], color=colors[file], linestyle=linestyles[file])
            elif file == 1:
                for i in range(num_plots):
                    axs[int(i/2), i % 2].plot(dates,
                                              total[:, i], label=legendplot[file], color=colors[file], linestyle=linestyles[file])
                min_date = dates[0]

        h5file.close()

    # define some characteristics of the plot
    for i in range(num_plots):
        axs[int(i/2), i % 2].set_title(secir_dict[i], fontsize=8)
        # axs[int(i/2), i % 2].set_ylim(bottom=10**3, top=5*10**3)
        axs[int(i/2), i % 2].set_xlim(left=-20, right=dates[-1])
        axs[int(i/2), i % 2].grid(True, linestyle='--')
        # axs[int(i/2), i % 2].legend(fontsize=8)
        axs[int(i/2), i % 2].ticklabel_format(axis='y',
                                              style='sci', scilimits=(0, 0))

    labels = ['ODE', 'IDE']
    fig.legend(labels, bbox_to_anchor=(0.1, -0.73, 0.8, 0.8),
               fancybox=False, shadow=False, ncol=1)  # bbox_to_anchor=(0.1, -0.73, 0.8, 0.8),

    fig.supxlabel(' Time')
    fig.supylabel('Number of persons')
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=None, hspace=0.6)
    # plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.3)

    # save result
    if save:
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f'ide_ode_{fileending}.png',
                    bbox_inches='tight', dpi=500)

    plt.close()

    # # plot legend separately
    # figsize = (2, 0.5)
    # fig_leg = plt.figure(figsize=figsize)
    # ax_leg = fig_leg.add_subplot(111)
    # # add the legend from the previous axes
    # ax_leg.legend(*axs[0, 0].get_legend_handles_labels(), loc='center', ncol=2)
    # # hide the axes frame and the x/y labels
    # ax_leg.axis('off')
    # fig_leg.savefig('plots/legend.png', dpi=500)


if __name__ == '__main__':
    # Path to simulation results

    dt_ode = '1e-1'
    dt_ide = '1e-1'

    setting = 2

    flows = True

    legendplot = list(["ODE", "IDE"])

    # Plot comparison of ODE and IDE models

    # -------------- Simulations based on LEOSS data ---------------

    # data_dir = os.path.join(os.path.dirname(
    #     __file__), "../results/fictional/leoss/")

    # compare_results([os.path.join(data_dir, f"fictional_ode_leoss_2.0_12_0.1000_flows"),
    #                  os.path.join(data_dir, f"fictional_ide_leoss_2.0_12_0.1000_flows")],
    #                 legendplot, flows=True, fileending="2.0_12_0.1000_flows", save=True, save_dir='plots/leoss/flows/')

    # # # # # # # # # Simulations based on Covasim data
    # data_dir = os.path.join(os.path.dirname(
    #     __file__), "../results/fictional/covasim/")

    # compare_results([os.path.join(data_dir, f"fictional_ode_covasim_0.5_12_0.1000_flows"),
    #                  os.path.join(data_dir, f"fictional_ide_covasim_0.5_12_0.1000_flows")],
    #                 legendplot, flows=True, fileending="0.5_12_0.1000_flows", save=True, save_dir='plots/covasim/flows/')

    # compare_results([os.path.join(data_dir, f"fictional_ode_covasim_0.0_12_0.1000_flows"),
    #                  os.path.join(data_dir, f"fictional_ide_covasim_0.0_12_0.1000_flows")],
    #                 legendplot, flows=True, fileending="0.0_12_0.1000_flows", save=True, save_dir='plots/covasim/flows/')

    # compare_results([os.path.join(data_dir, f"fictional_ode_covasim_0.0_12_0.1000_compartments"),
    #                  os.path.join(data_dir, f"fictional_ide_covasim_0.0_12_0.1000_compartments")],
    #                 legendplot, flows=False, fileending="0.0_12_0.1000_compartments", save=True, save_dir='plots/covasim/compartments/')

    # # # #  Real scenario
    data_dir = os.path.join(os.path.dirname(
        __file__), "../results/real/")

    # compare_results([os.path.join(data_dir, f"ode_2020-11-01_30_0.1000_flows"),
    #                  os.path.join(data_dir, f"ide_2020-11-01_30_0.1000_flows")],
    #                 legendplot, flows=True, fileending="2020-11-01_30_0.1000_flows", save=True, save_dir='plots/real/')

    # compare_results([os.path.join(data_dir, f"ode_2020-11-01_30_0.1000_compartments"),
    #                  os.path.join(data_dir, f"ide_2020-11-01_30_0.1000_compartments")],
    #                 legendplot, flows=False, fileending="2020-11-01_30_0.1000_compartments", save=True, save_dir='plots/real/')

    start_date = '2020-10-1'
    simulation_time = 45
    timestep = '0.1000'
    compare_results([os.path.join(data_dir, f"ode_{start_date}_{simulation_time}_{timestep}_flows"),
                     os.path.join(data_dir, f"ide_{start_date}_{simulation_time}_{timestep}_flows")],
                    legendplot, flows=True, fileending=f"{start_date}_{simulation_time}_{timestep}_flows", save=True, save_dir=f"plots/real/{start_date}/{simulation_time}/")

    compare_results([os.path.join(data_dir, f"ode_{start_date}_{simulation_time}_{timestep}_compartments"),
                     os.path.join(data_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")],
                    legendplot, flows=False, fileending=f"{start_date}_{simulation_time}_{timestep}_compartments", save=True, save_dir=f'plots/real/{start_date}/{simulation_time}/')

    # # # # #  Constant init
    # data_dir = os.path.join(os.path.dirname(
    #     __file__), "../results/fictional/constant_init/")

    # compare_results([os.path.join(data_dir, f"ode_constant_init_30_0.1000_flows"),
    #                  os.path.join(data_dir, f"ide_constant_init_30_0.1000_flows")],
    #                 legendplot, flows=True, fileending="30_0.1000_flows", save=True, save_dir='plots/constant_init/')

    # compare_results([os.path.join(data_dir, f"ode_constant_init_30_0.1000_compartments"),
    #                  os.path.join(data_dir, f"ide_constant_init_30_0.1000_compartments")],
    #                 legendplot, flows=False, fileending="30_0.1000_compartments", save=True, save_dir='plots/constant_init/')
