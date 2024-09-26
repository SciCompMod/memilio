import h5py
import os
import matplotlib.pyplot as plt

import memilio.epidata.getDataIntoPandasDataFrame as gd

# Define compartments.
secir_dict = {0: 'Susceptible', 1: 'Exposed', 2: 'Infected', 3: 'Recovered'}

# Define color and style to be used while plotting for different models to
# make plots consistent.
color_dict = {0: '#1f77b4',
              1: '#2ca02c'
              }
linestyle_dict = {"ODE SI": 'dashed',
                  "ODE Improved": 'dotted',
                  "Graph": 'dashdot'
                  }


def compare_all_compartments(
        files,
        legendplot,
        filename_plot="compare_compartments"):

    fig, axs = plt.subplots(
        2, 2, sharex='all', num=filename_plot, tight_layout=False)

    # Add simulation results to plot.
    for file in range(len(files)):
        # Load data.
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        number_regions = len(list(h5file.keys()))
        for region in range(number_regions):
            if (len(list(h5file[list(h5file.keys())[region]].keys())) > 3):
                data = h5file[list(h5file.keys())[region]]
                dates = data['Time'][:]

                number_regions = len(
                    list(h5file[list(h5file.keys())[region]].keys())) - 2
                for region in range(number_regions):
                    total = data['Group' + str(region + 1)][:, :]
                    if (total.shape[1] != 4):
                        raise gd.DataError(
                            "Expected a different number of compartments.")
                    # Plot result.
                    if legendplot[file] in linestyle_dict:
                        for i in range(4):
                            axs[int(i / 2),
                                i % 2].plot(dates,
                                            total[:,
                                                  i],
                                            label=legendplot[file] +
                                            " Region " + str(region),
                                            linewidth=1.2,
                                            linestyle=linestyle_dict[legendplot[file]],
                                            color=color_dict[region])
                    else:
                        for i in range(4):
                            axs[int(i / 2), i % 2].plot(dates, total[:, i],
                                                        label=legendplot[file], linewidth=1.2)
            else:
                data = h5file[list(h5file.keys())[region]]
                dates = data['Time'][:]
                # As there should be only one Group, total is the simulation
                # result.
                total = data['Total'][:, :]
                if (total.shape[1] != 4):
                    raise gd.DataError(
                        "Expected a different number of compartments.")
                # Plot result.
                if legendplot[file] in linestyle_dict:
                    for i in range(4):
                        axs[int(i / 2),
                            i % 2].plot(dates,
                                        total[:,
                                              i],
                                        label=legendplot[file] +
                                        " Region " + str(region),
                                        linewidth=1.2,
                                        linestyle=linestyle_dict[legendplot[file]],
                                        color=color_dict[region])
                else:
                    for i in range(4):
                        axs[int(i / 2), i % 2].plot(dates, total[:, i],
                                                    label=legendplot[file], linewidth=1.2)
        h5file.close()

    # Define some characteristics of the plot.
    for i in range(4):
        axs[int(i / 2), i % 2].set_title(secir_dict[i], fontsize=8)
        axs[int(i / 2), i % 2].set_xlim(left=0, right=dates[-1])
        axs[int(i / 2), i % 2].grid(True, linestyle='--')
        axs[int(i / 2), i % 2].tick_params(axis='y', labelsize=7)
        axs[int(i / 2), i % 2].tick_params(axis='x', labelsize=7)
        # axs[int(i/2), i % 2].xaxis.set_ticks(np.arange(0, dates[-1]+1, 5))

    fig.supxlabel('Time (in days)', fontsize=9)

    lines, labels = axs[0, 0].get_legend_handles_labels()
    lgd = fig.legend(lines, labels, ncol=len(legendplot), loc='outside lower center',
                     fontsize=10, bbox_to_anchor=(0.5, - 0.06), bbox_transform=fig.transFigure)

    plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.1)
    plt.subplots_adjust(bottom=0.09)

    # Save result.
    if not os.path.isdir('Plots'):
        os.makedirs('Plots')
    fig.savefig('Plots/' + filename_plot + '.png',
                bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=500)


def plot_new_infections(
        files,
        ylim,
        legendplot,
        filename_plot="compare_new_infections"):

    plt.figure(filename_plot)

    # Add simulation results to plot.
    for file in range(len(files)):
        # Load data.
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        number_regions = len(list(h5file.keys()))
        for region in range(number_regions):
            if (len(list(h5file[list(h5file.keys())[region]].keys())) > 3):
                data = h5file[list(h5file.keys())[region]]
                dates = data['Time'][:]

                number_regions = len(
                    list(h5file[list(h5file.keys())[region]].keys())) - 2
                for region_ in range(number_regions):
                    total = data['Group' + str(region_ + 1)][:, :]
                    if (total.shape[1] != 4):
                        raise gd.DataError(
                            "Expected a different number of compartments.")
                    incidence = (total[:-1, 0] - total[1:, 0]
                                 ) / (dates[1:] - dates[:-1])
                    # Plot result.
                    if legendplot[file] in linestyle_dict:
                        plt.plot(dates[1:],
                                 incidence,
                                 linewidth=1.2,
                                 linestyle=linestyle_dict[legendplot[file]],
                                 color=color_dict[region_])
                    else:
                        plt.plot(dates[1:], incidence, linewidth=1.2)
            else:
                data = h5file[list(h5file.keys())[region]]
                dates = data['Time'][:]
                # As there should be only one Group, total is the simulation
                # result.
                total = data['Total'][:, :]
                if (total.shape[1] != 4):
                    raise gd.DataError(
                        "Expected a different number of compartments.")
                incidence = (total[:-1, 0] - total[1:, 0]) / \
                    (dates[1:] - dates[:-1])
                # Plot result.
                if legendplot[file] in linestyle_dict:
                    plt.plot(dates[1:],
                             incidence,
                             linewidth=1.2,
                             linestyle=linestyle_dict[legendplot[file]],
                             color=color_dict[region])
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
    plt.savefig(
        'Plots/' +
        filename_plot +
        '.png',
        bbox_inches='tight',
        dpi=500)


if __name__ == '__main__':
    data_dir = os.path.join(os.path.dirname(__file__), "..", "cpp", "build")
    plot_new_infections([os.path.join(data_dir, "ode_result_standard"),
                         os.path.join(data_dir, "ode_result_improved"),
                         os.path.join(data_dir, "graph_result")],
                        2e3, legendplot=list(["ODE SI", "ODE Improved", "Graph"]))
    compare_all_compartments([os.path.join(data_dir, "ode_result_standard"),
                              os.path.join(data_dir, "ode_result_improved"),
                              os.path.join(data_dir, "graph_result")],
                             legendplot=list(["ODE SI", "ODE Improved", "Graph"]))
