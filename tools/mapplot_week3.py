import datetime as dt
import os.path
import imageio
import cv2

import numpy as np
import pandas as pd

import memilio.epidata.getPopulationData as gpd
import memilio.plot.plotMap as pm
from memilio.epidata import geoModificationGermany as geoger
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import SymLogNorm, LinearSegmentedColormap
from tqdm.auto import tqdm

from PIL import Image


if __name__ == '__main__':
    path_results = "/localdata1/test/memilio/test"

    plot_flows = False
    percentile = "p50"

    # create dir plots in path_results
    plot_path = os.path.join(path_results, "plots")
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)

    interventions = ["No_intervention", "20p_reduc", "40p_reduc", "60p_reduc"]
    # alle zahlen von 0 bis 30
    days = [i for i in range(14, 14 + 31)]

    #  choose compartments, we want to plot
    # TODO: Missing for flows
    timms = [27, 28]
    deads = [24, 25, 26]
    susceptible = [0, 1, 23]
    infected = [5, 6, 7, 8, 9, 10, 11, 12, 13,
                14, 15, 16, 17, 18, 19, 20, 21, 22]
    infected_symp = [11, 12, 13, 14, 15, 16]
    symptomatic = [11, 12, 13, 14, 15, 16]  # 2, 3, 4,
    icu = [20, 21, 22]
    exposed = [2, 3, 4]
    hospitalized = [17, 18, 19]

    get_max_val = 0

    # progress bar
    total_iterations = len(interventions) * len(days)

    progress_bar = tqdm(total=total_iterations)

    for intervention in interventions:
        path = os.path.join(
            path_results, intervention, "run_0",  percentile, "Results")

        files_input = {
            'Data set 1': path}
        file_format = 'h5'
        age_groups = {0: '0-4', 1: '5-14', 2: '15-34',
                      3: '35-59', 4: '60-79', 5: '80+'}
        if len(age_groups) == 6:
            filter_age = None
        else:
            if file_format == 'json':
                filter_age = [val for val in age_groups.values()]
            else:
                filter_age = ['Group' + str(key)
                              for key in age_groups.keys()]

        relative = False
        num_days = 31

        min_val = 0
        max_val = 350

        # max_val_2 = 200

        norm = SymLogNorm(linthresh=1, linscale=0.7,
                          vmin=min_val, vmax=max_val)

        # Own colorbar
        colors = ["green", "yellow", "red", "purple"]
        cmap = LinearSegmentedColormap.from_list(
            "my_colormap", colors)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        # Create a new figure for the colorbar
        cbar_fig, ax = plt.subplots(figsize=(8, 1))
        cbar = plt.colorbar(sm, orientation='horizontal', cax=ax)
        cbar.ax.tick_params(labelsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(
            plot_path, 'colorbar.png'), dpi=300)
        plt.clf()

        for day in days:

            i = 0
            for file in files_input.values():
                # MEmilio backend hdf5 example

                df = pm.extract_data(
                    file, region_spec=None, column=None, date=day,
                    filters={'Group': filter_age,
                             'InfectionState': infected_symp},
                    file_format=file_format)

                # all entries to numeric
                df = df.apply(pd.to_numeric, errors='coerce')
                dfs_all = df

                if relative:

                    try:
                        population = pd.read_json(
                            'data/pydata/Germany/county_current_population.json')
                    # pandas>1.5 raise FileNotFoundError instead of ValueError
                    except (ValueError, FileNotFoundError):
                        print(
                            "Population data was not found. Download it from the internet.")
                        population = gpd.get_population_data(
                            read_data=False, file_format=file_format,
                            out_folder='data/pydata/Germany/', no_raw=True,
                            split_gender=False, merge_eisenach=True)

                    # For fitting of different age groups we need format ">X".
                    age_group_values = list(age_groups.values())
                    age_group_values[-1] = age_group_values[-1].replace(
                        '80+', '>79')
                    # scale data
                    df = pm.scale_dataframe_relative(
                        df, age_group_values, population)

            filename = intervention + "_day_" + str(day-14)

            dfs_all = dfs_all.apply(pd.to_numeric, errors='coerce')

            dfs_all_sorted = dfs_all.sort_values(by='Region')
            dfs_all_sorted = dfs_all_sorted.reset_index(drop=True)

            # max_val_2 = dfs_all_sorted['Count'].max()
            # if max_val_2 > get_max_val:
            #     get_max_val = max_val_2

            pm.plot_map(norm,
                        dfs_all_sorted, scale_colors=[
                            min_val, max_val],
                        legend=['', ''],
                        title='Total Individuals hosptialized day ' +
                        str(day - 14),
                        plot_colorbar=False,
                        output_path=plot_path,
                        fig_name=filename, dpi=300,
                        outercolor='white')

            progress_bar.update()
    progress_bar.close()
    # print("max value: ", get_max_val)

    # create single gif for each intervention
    path_colorbar = "/localdata1/test/memilio/test/plots/colorbar.png"
    colorbar = imageio.v2.imread(path_colorbar)

    num_days = 31
    interventions_gif = ["No_intervention",
                         "20p_reduc", "40p_reduc", "60p_reduc"]
    for intervention in interventions_gif:
        frames = []
        for day in range(0, num_days):
            filename = intervention + '_day_' + str(day)
            filename = filename + '.png'
            output_path = "/localdata1/test/memilio/test/plots/"
            image = imageio.v2.imread(os.path.join(output_path, filename))

            # Resize colorbar to match the width of the image
            # Ensure colorbar is not wider than the main image
            if colorbar.shape[1] > image.shape[1]:
                # Resize colorbar to match the width of the image
                colorbar = cv2.resize(
                    colorbar, (image.shape[1], colorbar.shape[0]))

            # Calculate padding sizes
            padding_left = (image.shape[1] - colorbar.shape[1]) // 2
            padding_right = image.shape[1] - colorbar.shape[1] - padding_left

            # Resize colorbar to match the width of the image
            colorbar_resized = np.pad(
                colorbar, ((0, 0), (padding_left, padding_right), (0, 0)), mode='constant')

            # Append colorbar under the image
            image = np.concatenate((image, colorbar_resized), axis=0)

            frames.append(image)

        imageio.mimsave(output_path + "/gifs/" + intervention + ".gif",  # output gif
                        frames,          # array of input frames
                        duration=1)

    for day in range(0, num_days):
        images = []
        for intervention in interventions_gif:
            filename = intervention + '_day_' + str(day) + '.png'
            output_path = "/localdata1/test/memilio/test/plots/"
            image = imageio.v2.imread(os.path.join(output_path, filename))
            images.append(image)
        # Stack images horizontally
        frame = np.concatenate(images, axis=1)

        # Ensure colorbar is not wider than the main image
        if colorbar.shape[1] > frame.shape[1]:
            # Resize colorbar to match the width of the image
            colorbar = cv2.resize(
                colorbar, (frame.shape[1], colorbar.shape[0]))

        # Calculate padding sizes
        padding_left = (frame.shape[1] - colorbar.shape[1]) // 2
        padding_right = frame.shape[1] - colorbar.shape[1] - padding_left

        # Resize colorbar to match the width of the image
        colorbar_padded = np.pad(
            colorbar, ((0, 0), (padding_left, padding_right), (0, 0)), mode='constant')

        # Append colorbar under the frame, centered
        frame = np.concatenate((frame, colorbar_padded), axis=0)

        # Ensure all frames have the same shape
        if frames:
            if frame.shape != frames[0].shape:
                # Resize frame to match the shape of the first frame
                frame = cv2.resize(
                    frame, (frames[0].shape[1], frames[0].shape[0]))

        frames.append(frame)

    imageio.mimsave(output_path + "/gifs/" + intervention +
                    "combined.gif", frames, duration=1)
