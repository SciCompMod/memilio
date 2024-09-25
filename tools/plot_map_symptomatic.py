import datetime as dt
import os.path
import imageio

import numpy as np
import pandas as pd

import memilio.epidata.getPopulationData as gpd
import memilio.plot.plotMap as pm
from memilio.epidata import geoModificationGermany as geoger
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import SymLogNorm, Normalize
from matplotlib.colors import LinearSegmentedColormap


res_dir = "/localdata1/code/memilio/results_paper/transmission1.5/"
quartil = "50"
infections = False  # set true when we want to plot the transsmissions.
relative = True  # set true when we want to plot relative to the total population

if __name__ == '__main__':

    masks = ["0", "1_ffp2"]
    inf_mobility_node_mask_0 = []
    inf_mobility_node_mask_1_ffp2 = []
    for mode in masks:
        dir = res_dir + "mask_" + mode + "/flows_mb"
        files = os.listdir(dir)
        files = [x for x in files if x.startswith(
            "transmission_mobility_run_")]
        file_mobility = "/transmission_mobility_run_"
        file_local = "/total_transmission_local_run_"
        if not infections:
            file_mobility = "/symptomps_mobility_run_"
            file_local = "/total_symptomps_local_run_"

        num_files = len(files)
        progress = 0

        for indx in range(num_files):
            data = np.loadtxt(
                dir + file_local + str(indx) + ".txt", skiprows=1)
            infections_local = data[:, 1:]
            infections_local = np.diff(infections_local)

            data = np.loadtxt(
                dir + file_mobility + str(indx) + ".txt", skiprows=1)
            indices = data[:, 0].astype(int)
            values = data[:, 2:]
            np.add.at(infections_local, indices, values)

            if mode == "0":
                inf_mobility_node_mask_0.append(infections_local)

            else:
                inf_mobility_node_mask_1_ffp2.append(infections_local)

            progress += 1
            print("Progress: " + str(progress) + "/" + str(num_files))

    inf_mobility_node_mask_0 = np.array(inf_mobility_node_mask_0)
    inf_mobility_node_mask_0_quartils = np.percentile(
        inf_mobility_node_mask_0, [int(quartil)], axis=0)

    inf_mobility_node_mask_1_ffp2 = np.array(inf_mobility_node_mask_1_ffp2)
    inf_mobility_node_mask_1_ffp2_quartils = np.percentile(
        inf_mobility_node_mask_1_ffp2, [int(quartil)], axis=0)

    if relative:
        population = pd.read_json(
            'data/pydata/Germany/county_current_population.json')

        for i in range(0, population.shape[0]):
            inf_mobility_node_mask_0_quartils[0][i] = inf_mobility_node_mask_0_quartils[0][i] / \
                population['Population'][i] * 100000
            inf_mobility_node_mask_1_ffp2_quartils[0][i] = inf_mobility_node_mask_1_ffp2_quartils[0][i] / \
                population['Population'][i] * 100000

    min_val = np.min(np.concatenate([inf_mobility_node_mask_0_quartils,
                     inf_mobility_node_mask_1_ffp2_quartils]))
    max_val = np.max(np.concatenate([inf_mobility_node_mask_0_quartils,
                     inf_mobility_node_mask_1_ffp2_quartils]))

    norm = SymLogNorm(linthresh=1, linscale=0.7,
                      vmin=min_val, vmax=max_val)

    dir_plots = res_dir + "map"
    if infections:
        dir_plots = dir_plots + "/transmissions"
    else:
        dir_plots = dir_plots + "/infections"

    if relative:
        dir_plots = dir_plots + "/relative"
    if not os.path.exists(dir_plots):
        os.makedirs(dir_plots)

    colors = ["white", "green", "yellow", "red", "purple"]
    cmap = LinearSegmentedColormap.from_list("my_colormap", colors)

    # norm = Normalize(vmin=min_val, vmax=max_val)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Create a new figure for the colorbar
    cbar_fig, ax = plt.subplots(figsize=(8, 1))
    cbar = plt.colorbar(sm, orientation='horizontal', cax=ax)

    # cbar.ax.tick_params(labelsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(
        dir_plots, 'colorbar.png'), dpi=300)
    plt.clf()

    counties_considered = geoger.get_county_ids()
    for day in range(0, inf_mobility_node_mask_0_quartils.shape[-1]):

        day_mask_0 = inf_mobility_node_mask_0_quartils[0, :, day]
        day_mask_1_ffp2 = inf_mobility_node_mask_1_ffp2_quartils[0, :, day]

        df_mask_0 = pd.DataFrame(
            {'Region': counties_considered, 'Count': day_mask_0})
        df_mask_1_ffp2 = pd.DataFrame(
            {'Region': counties_considered, 'Count': day_mask_1_ffp2})

        print("day: " + str(day))
        for mode in masks:
            filename = "mask_" + mode + "_day_" + str(day)
            df = df_mask_0
            if mode == "1_ffp2":
                df = df_mask_1_ffp2
                filename = "mask_" + mode + "_day_" + str(day)

            pm.plot_map(norm,
                        df, scale_colors=[min_val, max_val],
                        legend=['', ''],
                        title='Day' + str(day),
                        plot_colorbar=False,
                        output_path=dir_plots,
                        fig_name=filename, dpi=300,
                        outercolor='white')
