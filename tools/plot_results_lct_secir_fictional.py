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

"""@plot_results_lct_secir_fictional.py
Functions to plot and compare results of simulations with different kind of models,
eg LCT, IDE or ODE SECIR models without division in agegroups.
To compare with real data, use the script plot_results_lct_secir_real.py.

The data to be plotted should be stored in a '../data/simulation_lct' folder as .h5 files.
Data could be generated eg by executing the file ./cpp/examples/lct_secir_fictional_scenario.cpp.
"""

import h5py
import os
import matplotlib.pyplot as plt

import memilio.epidata.getDataIntoPandasDataFrame as gd

# Define compartments.
secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Carrier', 3: 'Infected', 4: 'Hospitalized',
              5: 'ICU', 6: 'Recovered', 7: 'Dead'}

# Define color and style to be used while plotting for different models to make plots consistent.
color_dict = {'ODE': '#1f77b4',
              'LCT3':'#2ca02c' ,
              'LCT10':'#ff7f0e' ,
              'LCT20': '#9467bd',
              'IDE3': '#e377c2',
              'IDE10': '#17becf',
              'GLCT3': '#bbf90f',
              'GLCT10':'#ef4026'
              }
linestyle_dict = {'ODE': 'solid',
                  'LCT3': 'solid',
                  'LCT10': 'solid',
                  'LCT20': 'solid',
                  'IDE3': 'dashdot',
                  'IDE10': 'dashdot',
                  'GLCT3': 'solid',
                  'GLCT10':'solid'
                  }

def compare_all_compartments(files, legendplot, filename_plot="compare_compartments"):
    """ Creates a 4x2 Plot with one subplot per compartment and one line per result one wants to compare.

    @param[in] files: paths of the files (without file extension .h5) with the simulation results 
        that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models).
        Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """
    fig, axs = plt.subplots(
        4, 2, sharex='all', num=filename_plot, tight_layout=True)

    # Add simulation results to plot.
    for file in range(len(files)):
        # Load data.
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result.
        total = data['Total'][:, :]
        if (total.shape[1] != 8):
            raise gd.DataError("Expected a different number of compartments.")
        # Plot result.
        if legendplot[file] in linestyle_dict:
            for i in range(8):
                axs[int(i/2), i % 2].plot(dates,
                                          total[:, i], label=legendplot[file], linewidth=1.2, 
                                          linestyle=linestyle_dict[legendplot[file]], 
                                          color=color_dict[legendplot[file]])
        else:
            for i in range(8):
                axs[int(i/2), i % 2].plot(dates,
                                          total[:, i], label=legendplot[file], linewidth=1.2)
        h5file.close()

    # Define some characteristics of the plot.
    for i in range(8):
        axs[int(i/2), i % 2].set_title(secir_dict[i], fontsize=8)
        axs[int(i/2), i % 2].set_xlim(left=0, right=dates[-1])
        axs[int(i/2), i % 2].grid(True, linestyle='--')
        axs[int(i/2), i % 2].tick_params(axis='y', labelsize=7)
        axs[int(i/2), i % 2].tick_params(axis='x', labelsize=7)
        # axs[int(i/2), i % 2].xaxis.set_ticks(np.arange(0, dates[-1]+1, 5))

    fig.supxlabel('Time (in days)', fontsize=9)

    lines, labels = axs[0, 0].get_legend_handles_labels()
    lgd = fig.legend(lines, labels, ncol=len(legendplot),  loc='outside lower center',
                     fontsize=10, bbox_to_anchor=(0.5, - 0.06), bbox_transform=fig.transFigure)

    plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.1)
    plt.subplots_adjust(bottom=0.09)

    # Save result.
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    fig.savefig('Plots/'+filename_plot+'.png',
                bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=500)


def plot_new_infections(files, ylim, legendplot, filename_plot="compare_new_infections"):
    """ Single plot to compare the incidence of different results. Incidence means the number of people leaving the 
        susceptible class per day.

    @param[in] files: paths of the files (without file extension .h5) with the simulation results
         that should be compared.
        Results should contain exactly 8 compartments (so use accumulated numbers for LCT models). 
        Names can be given in form of a list.
        One could compare results with eg different parameters or different models.
    @param[in] ylim: upper limit for the y-axis of the plot.
    @param[in] legendplot: list with names for the results that should be used for the legend of the plot.
    @param[in] filename_plot: name to use as the file name for the saved plot.
    """
    plt.figure(filename_plot)

    # Add simulation results to plot.
    for file in range(len(files)):
        # Load data.
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result.
        total = data['Total'][:, :]
        if (total.shape[1] != 8):
            raise gd.DataError(
                "Expected a different number of compartments.")
        incidence = (total[:-1, 0]-total[1:, 0])/(dates[1:]-dates[:-1])
        if (file == 0):
            incidence_ref = incidence
            print("Relative deviation of the new infections at time " +
                  str(dates[-1])+" from the results of "+legendplot[0])
        deviation = (incidence[-1]-incidence_ref[-1])/incidence_ref[-1]
        print(legendplot[file]+": "+str(deviation))
        # Plot result.
        if legendplot[file] in linestyle_dict:
            plt.plot(dates[1:], incidence, linewidth=1.2,
                     linestyle=linestyle_dict[legendplot[file]], color=color_dict[legendplot[file]])
        else:
            plt.plot(dates[1:], incidence, linewidth=1.2)

        h5file.close()

    plt.xlabel('Time (in days)', fontsize=16)
    # plt.xticks(np.arange(0, dates[-1]+1, 5))
    plt.yticks(fontsize=9)
    plt.ylabel('New infections per day', fontsize=14)
    plt.ylim(bottom=0, top=ylim)
    plt.xlim(left=0, right=dates[-1])
    plt.legend(legendplot, fontsize=14, framealpha=0.5)
    plt.grid(True, linestyle='--')
    plt.tight_layout()

    # Save result.
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    plt.savefig('Plots/'+filename_plot+'.png', bbox_inches='tight', dpi=500)


if __name__ == '__main__':
    # simulation results should be stored in folder "../data/simulation_lct".
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "data", "simulation_lct")

    case=3
    
    if case==1:
        ### Compare simulation results of different LCT models with a rise of R0 and a long simulation time.
        ## R0=4.0.
        plot_new_infections([os.path.join(data_dir, "riseR04long", "fictional_lct_4.0_1_long"),  
                            os.path.join(data_dir, "riseR04long", "fictional_lct_4.0_3_long"),  
                            os.path.join(data_dir, "riseR04long", "fictional_lct_4.0_10_long"), 
                            os.path.join(data_dir, "riseR04long", "fictional_lct_4.0_20_long")],
                            6e6, legendplot=list(["ODE", "LCT3", "LCT10", "LCT20"]),  
                            filename_plot="new_infections_lct_rise4long")
        compare_all_compartments([os.path.join(data_dir, "riseR04long", "fictional_lct_4.0_1_long"),  
                                os.path.join(data_dir, "riseR04long", "fictional_lct_4.0_3_long"), 
                                os.path.join(data_dir, "riseR04long", "fictional_lct_4.0_10_long"),
                                os.path.join(data_dir, "riseR04long", "fictional_lct_4.0_20_long")],
                                legendplot=list(["ODE", "LCT3", "LCT10", "LCT20"]), 
                                filename_plot="compartments_lct_rise4long")
        ## R0=2.0.
        plot_new_infections([os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_1_long"),  
                            os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_3_long"),  
                            os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_10_long"), 
                            os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_20_long")],
                            2e6, legendplot=list(["ODE", "LCT3", "LCT10", "LCT20"]),  
                            filename_plot="new_infections_lct_rise2long")
        compare_all_compartments([os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_1_long"),  
                                os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_3_long"),  
                                os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_10_long"), 
                                os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_20_long")],
                                legendplot=list(["ODE", "LCT3", "LCT10", "LCT20"]),   
                                filename_plot="compartments_lct_rise2long")
    elif case==2:
        ### Compare simulation results of different LCT and IDE models with a rise of R0 and a long simulation time.
        ## R0=2.0.
        plot_new_infections([os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_1_long"), 
                            os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_3_long"),  
                            os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_10_long"), 
                            os.path.join(data_dir, "riseR02long", "fictional_ide_2.0_3_long"),  
                            os.path.join(data_dir, "riseR02long", "fictional_ide_2.0_10_long")],
                            2e6, legendplot=list(["ODE", "LCT3", "LCT10", "IDE3", "IDE10"]),  
                            filename_plot="new_infections_ide_rise2long")
        compare_all_compartments([os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_1_long"), 
                                os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_3_long"),  
                                os.path.join(data_dir, "riseR02long", "fictional_lct_2.0_10_long"), 
                                os.path.join(data_dir, "riseR02long", "fictional_ide_2.0_3_long"),  
                                os.path.join(data_dir, "riseR02long", "fictional_ide_2.0_10_long")],
                                legendplot=list(["ODE", "LCT3", "LCT10", "IDE3", "IDE10"]),    
                                filename_plot="compartments_ide_rise2long")
    elif case==3:
        ### Compare simulation results of different LCT models with a rise of R0 to 2.
        plot_new_infections([os.path.join(data_dir, "riseR02short", "fictional_lct_2.0_1"), 
                            os.path.join(data_dir, "riseR02short", "fictional_lct_2.0_3"),  
                            os.path.join(data_dir, "riseR02short", "fictional_lct_2.0_10"), 
                            os.path.join(data_dir, "riseR02short", "fictional_lct_2.0_20")],
                            15500, legendplot=list(["ODE", "LCT3", "LCT10", "LCT20"]), 
                            filename_plot="new_infections_lct_rise2") 
        compare_all_compartments([os.path.join(data_dir, "riseR02short", "fictional_lct_2.0_1"), 
                            os.path.join(data_dir, "riseR02short", "fictional_lct_2.0_3"),  
                            os.path.join(data_dir, "riseR02short", "fictional_lct_2.0_10"), 
                            os.path.join(data_dir, "riseR02short", "fictional_lct_2.0_20")],
                             legendplot=list(["ODE", "LCT3", "LCT10", "LCT20"]), 
                            filename_plot="compartments_lct_rise2"   
                            )
    elif case==4:
        ### Compare simulation results of different LCT models with a drop of R0 to 0.5.
        plot_new_infections([os.path.join(data_dir, "dropR0short", "fictional_lct_0.5_1"), 
                            os.path.join(data_dir, "dropR0short", "fictional_lct_0.5_3"),  
                            os.path.join(data_dir, "dropR0short", "fictional_lct_0.5_10"), 
                            os.path.join(data_dir, "dropR0short", "fictional_lct_0.5_20")],
                            4100, legendplot=list(["ODE", "LCT3", "LCT10", "LCT20"]), 
                            filename_plot="new_infections_lct_drop0.5")
    elif case==5:
        ### Compare simulation results of different LCT, GLCT and IDE models with a rise of R0 and a long simulation time.
        ## R0=2.0.
        plot_new_infections([os.path.join(data_dir, "riseR02long", "fictional_glct_2.0_3_long"),  
                            os.path.join(data_dir, "riseR02long", "fictional_glct_2.0_10_long"),
                            os.path.join(data_dir, "riseR02long", "fictional_ide_2.0_3_long"),  
                            os.path.join(data_dir, "riseR02long", "fictional_ide_2.0_10_long")],
                            2e6, legendplot=list(["GLCT3", "GLCT10", "IDE3", "IDE10"]),  
                            filename_plot="new_infections_glct_rise2long")
        compare_all_compartments([os.path.join(data_dir, "riseR02long", "fictional_glct_2.0_3_long"),  
                            os.path.join(data_dir, "riseR02long", "fictional_glct_2.0_10_long"),
                            os.path.join(data_dir, "riseR02long", "fictional_ide_2.0_3_long"),  
                            os.path.join(data_dir, "riseR02long", "fictional_ide_2.0_10_long")],
                            legendplot=list(["GLCT3", "GLCT10", "IDE3", "IDE10"]),    
                            filename_plot="compartments_glct_rise2long")
    elif case==6:
        ### Compare simulation results of GLCT nd IDE models with a drop of R0 to 0.5.
        plot_new_infections([os.path.join(data_dir, "dropR0short", "fictional_glct_0.5_3"), 
                            os.path.join(data_dir, "dropR0short", "fictional_glct_0.5_10"),  
                            os.path.join(data_dir, "dropR0short", "fictional_ide_0.5_3"), 
                            os.path.join(data_dir, "dropR0short", "fictional_ide_0.5_10")],
                            4100, legendplot=list(["GLCT3", "GLCT10", "IDE3","IDE10"]), 
                            filename_plot="new_infections_glct_drop0.5")
        compare_all_compartments([os.path.join(data_dir, "dropR0short", "fictional_glct_0.5_3"), 
                            os.path.join(data_dir, "dropR0short", "fictional_glct_0.5_10"),  
                            os.path.join(data_dir, "dropR0short", "fictional_ide_0.5_3"), 
                            os.path.join(data_dir, "dropR0short", "fictional_ide_0.5_10")],
                            legendplot=list(["GLCT3", "GLCT10", "IDE3","IDE10"]), 
                            filename_plot="compartments_glct_drop0.5")
    elif case==7:
        ### Compare simulation results of GLCT and IDE models with a rise of R0 to 2.
        plot_new_infections([os.path.join(data_dir, "riseR02short", "fictional_glct_2.0_3"), 
                            os.path.join(data_dir, "riseR02short", "fictional_glct_2.0_10"),  
                            os.path.join(data_dir, "riseR02short", "fictional_ide_2.0_3"), 
                            os.path.join(data_dir, "riseR02short", "fictional_ide_2.0_10")],
                            15500, legendplot=list(["GLCT3", "GLCT10", "IDE3","IDE10"]), 
                            filename_plot="new_infections_glct_rise2.0")
        compare_all_compartments([os.path.join(data_dir, "riseR02short", "fictional_glct_2.0_3"), 
                            os.path.join(data_dir, "riseR02short", "fictional_glct_2.0_10"),  
                            os.path.join(data_dir, "riseR02short", "fictional_ide_2.0_3"), 
                            os.path.join(data_dir, "riseR02short", "fictional_ide_2.0_10")],
                            legendplot=list(["GLCT3", "GLCT10", "IDE3","IDE10"]), 
                            filename_plot="compartments_glct_rise2.0")
          
          
          
          