import h5py
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from memilio.epidata import getDataIntoPandasDataFrame as gd


def plot_changepoint(files, legendplot, flows=True, fileending="", save=True, save_dir='plots/'):

    fig, ax = plt.subplots()

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
            total = data['Total'][:, :]
        else:
            if len(data['Total'][0]) == 8:
                # As there should be only one Group, total is the simulation result
                total = data['Total'][:, :]
            elif len(data['Total'][0]) == 10:
                # in ODE there are two compartments we don't use, throw these out
                total = data['Total'][:, [0, 1, 2, 4, 6, 7, 8, 9]]

        dates = data['Time'][:]

        timestep = np.diff(dates)[0]
        tmax = dates[-1]

        # get indices where dates are >=0
        indices = np.where(dates >= 0)
        # plot data
        if flows:
            # ODE
            if file == 0:
                # transform cumulative flows to flows absolute flows
                # then transform from flows over time interval to flows at time points
                ax.plot(dates[indices[0][1:]], np.diff(total[indices[0], 0])/np.diff(dates[indices[0]]), label=legendplot[file],
                        color=colors[file], linestyle=linestyles[file])

            # IDE
            elif file == 1:
                # transform from flows over time interval to flows at time points
                ax.plot(dates[1:], total[1:, 0]/np.diff(dates), label=legendplot[file],
                        color=colors[file], linestyle=linestyles[file])

                date_idx = -int(tmax/timestep)-2
                print(
                    f"New infections at {dates[date_idx]}: {total[date_idx,0]}")
                date_idx = -int(tmax/timestep)-1
                print(
                    f"New infections at {dates[date_idx]}: {total[date_idx,0]}")
                date_idx = -int(tmax/timestep)
                print(
                    f"New infections at {dates[date_idx]}: {total[date_idx,0]}")
        else:
            incidence = (total[:-1, 0]-total[1:, 0])/(dates[1:]-dates[:-1])
            ax.plot(dates[indices[0][1:]], incidence, label=legendplot[file],
                    color=colors[file], linestyle=linestyles[file])

        h5file.close()

        # ax.set_title(secir_dict[i], fontsize=8)
        # axs[int(i/2), i % 2].set_ylim(bottom=0)
        ax.set_xlim(left=0, right=tmax)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(fontsize=12)

    fig.supxlabel('Simulation time [days]')
    fig.supylabel('Daily new infections')
    plt.subplots_adjust(left=None, bottom=None, right=None,
                        top=None, wspace=None, hspace=0.6)

    plt.tight_layout()

    # save result
    if save:
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + f"changepoint_{fileending}.png",
                    bbox_inches='tight', dpi=500)


if __name__ == '__main__':
    legendplot = list(["ODE", "IDE"])
    # Path to simulation results
    data_dir = os.path.join(os.path.dirname(
        __file__), "..", "results/fictional/covasim/")

    plot_changepoint([os.path.join(data_dir, f"fictional_ode_covasim_0.5_12_0.0100_flows"),
                     os.path.join(data_dir, f"fictional_ide_covasim_0.5_12_0.0100_flows")],
                     legendplot, flows=True, fileending="0.5_12_0.0100", save=True, save_dir='plots/covasim/changepoints/')

    plot_changepoint([os.path.join(data_dir, f"fictional_ode_covasim_2.0_12_0.0100_flows"),
                     os.path.join(data_dir, f"fictional_ide_covasim_2.0_12_0.0100_flows")],
                     legendplot, flows=True, fileending="2.0_12_0.0100", save=True, save_dir='plots/covasim/changepoints/')
