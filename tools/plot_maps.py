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
from matplotlib.colors import SymLogNorm, LinearSegmentedColormap, Normalize
from tqdm.auto import tqdm

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def plot_risk_map(path_results, path_plots, days, percentile):
    if not os.path.exists(path_plots):
        os.makedirs(path_plots)
    path_risk_results = os.path.join(path_results, "risk")
    plot_maps(path_risk_results, path_plots, compartments=[0], percentile=percentile,
              days=days, min_val=0, max_val=1, filename="risk_map", relative=False, age_groups={0: '0-4'})


def plot_icu_map(path_results, path_plots, days, percentile, max_val, icu_cap):
    path = os.path.join(path_plots)
    if not os.path.exists(path):
        os.makedirs(path)
    plot_maps(path_results, path, compartments=[7], percentile=percentile,
              days=days, min_val=0, max_val=max_val, filename="icu_map", relative=True, icu_cap=icu_cap)


def plot_flows(path_results, path_plots, days, percentile, mode, fn, symptomatic=True):
    flows_indx = [0]
    if symptomatic:
        flows_indx = [2, 4]
        fn += "_symptomatic"
    for mode in modes:
        path_plots_mode = os.path.join(path_plots, mode)
        if not os.path.exists(path_plots_mode):
            os.makedirs(path_plots_mode)
        plot_maps(os.path.join(path_results, mode, "flows"), path_plots_mode, compartments=flows_indx, percentile=percentile,
                  days=days, min_val=1, max_val=3000, filename=fn, relative=True, flows=True)


def plot_r0_map(path_results, path_plots, days, percentile, mode):
    for mode in modes:
        path_r = os.path.join(path_results, mode, 'r0')
        path_p = os.path.join(path_plots, mode)
        if not os.path.exists(path_r):
            os.makedirs(path_r)
        if not os.path.exists(path_p):
            os.makedirs(path_p)
        plot_maps(path_r, path_p, compartments=[0], percentile=percentile,
                  days=days, min_val=0, max_val=50, filename="r0_map", relative=False)
    x = 1


def create_colorbar(path_plots, norm, title):
    colors = ["green", "yellow", "red", "purple"]
    cmap = LinearSegmentedColormap.from_list("my_colormap", colors)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Create a new figure for the colorbar
    cbar_fig, ax = plt.subplots(figsize=(8, 1))
    cbar = plt.colorbar(sm, orientation='horizontal', cax=ax)
    cbar.ax.tick_params(labelsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(path_plots, title + '_colorbar.png'), dpi=300)
    plt.clf()


def plot_maps(path_results, path_plots, compartments, percentile, days, min_val, max_val, filename="data", relative=True, age_groups={0: '0-4', 1: '5-14', 2: '15-34',
                                                                                                                                      3: '35-59', 4: '60-79', 5: '80+'}, flows=False, path_results2="", icu_cap=1):
    progress_bar = tqdm(total=len(days))

    path = os.path.join(path_results, percentile, "Results")

    get_max_val = 0

    files_input = {
        'Data set 1': path}

    if path_results2 != "":
        path2 = os.path.join(path_results2, percentile, "Results")
        files_input['Data set 2'] = path2

    file_format = 'h5'

    if len(age_groups) == 6:
        filter_age = None
    else:
        if file_format == 'json':
            filter_age = [val for val in age_groups.values()]
        else:
            filter_age = ['Group' + str(key + 1)
                          for key in age_groups.keys()]

    norm = SymLogNorm(linthresh=1, linscale=0.7,
                      vmin=min_val, vmax=max_val)
    # linear norm
    # norm = Normalize(vmin=min_val, vmax=max_val)
    create_colorbar(path_plots, norm, filename)

    for day in days:
        i = 0
        for file in files_input.values():
            # MEmilio backend hdf5 example

            df = pm.extract_data(
                file, region_spec=None, column=None, date=day,
                filters={'Group': filter_age,
                         'InfectionState': compartments},
                file_format=file_format)

            # all entries to numeric
            df = df.apply(pd.to_numeric, errors='coerce')

            if flows:
                # if flows are plotted, the data is commulative
                # format to daily data by getting the data from the day before
                if day > 0:
                    df_previous = pm.extract_data(
                        file, region_spec=None, column=None, date=day-1,
                        filters={'Group': filter_age,
                                 'InfectionState': compartments},
                        file_format=file_format)
                    df['Count'] = df['Count'] - df_previous['Count']

                    # check if all values are positive
                    if df['Count'].min() < 0:
                        print("Negative values in data.")
                        df['Count'] = df['Count'].clip(lower=0)

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

                # overwrite data in col 'Count' with data from 'Count (rel)'
                df['Count'] = df['Count (rel)'] * 100_000 / icu_cap
                # del old column 'Count (rel)'
                df = df.drop(columns=['Count (rel)'])

            if i == 0:
                dfs_all = pd.DataFrame(df.iloc[:, 0])

            dfs_all[df.columns[-1] + ' ' + str(i)] = df[df.columns[-1]]
            i += 1

        fn = filename + "_day_" + str(day)

        dfs_all = dfs_all.apply(pd.to_numeric, errors='coerce')

        dfs_all_sorted = dfs_all.sort_values(by='Region')
        dfs_all_sorted = dfs_all_sorted.reset_index(drop=True)

        # or dfs_all_sorted['Count 1'].max() > get_max_val:
        if dfs_all_sorted['Count 0'].max() > get_max_val:
            if path_results2 == "":
                get_max_val = dfs_all_sorted['Count 0'].max()
            else:
                get_max_val = max(
                    dfs_all_sorted['Count 0'].max(), dfs_all_sorted['Count 1'].max())

        pm.plot_map(norm=norm, data=dfs_all_sorted,
                    scale_colors=[min_val, max_val],
                    legend=['', ''],
                    title='',
                    plot_colorbar=False,
                    output_path=path_plots,
                    fig_name=fn,
                    dpi=300,
                    outercolor='white')

        progress_bar.update()
    progress_bar.close()
    print("max value: ", get_max_val)


# /localdata1/code_2024/memilio/results/october_2020_original_params/kmin_0.200000_kmax_0.800000/


if __name__ == '__main__':
    modes = ["FeedbackDamping"]
    path_cwd = os.getcwd()
    icu_cap = [6, 9, 12, 15]
    cap_indx = 1
    # results/fixed_damping_kmin_0.300000/ClassicDamping/mse_428223962312.262817
    # path_results = os.path.join(
    #     path_cwd, "results", "ICUCap_" + str(icu_cap[cap_indx]) + ".000000", 'kmin_0.000000_kmax_0.700000')

    # /localdata1/code_2024/memilio/results/bremen/kmin_1.000000_kmax_1.000000/FeedbackDamping/contacts

    # /localdata1/code_2024/memilio/results/ICUCap_9.000000/rho_1.080000/BlendingFactorRegional_0.000000/kmin_0.000000_kmax_0.360000

    # /localdata1/code_2024/memilio/results/ICUCap_9.000000/rho_1.000000/BlendingFactorRegional_0.700000/kmin_0.000000_kmax_0.360000/FeedbackDamping/risk

    path_results = os.path.join(
        path_cwd, "results/new_regional_def", "ICUCap_9.000000", "rho_1.000000", "BlendingFactorRegional_0.700000", "kmin_0.000000_kmax_0.360000", "FeedbackDamping")
    path_plots = os.path.join(
        path_cwd, "plots", "ICUCap_" + str(icu_cap[cap_indx]) + ".000000", "rp_kmin_0.000000_kmax_0.360000", "new")
    path_data = os.path.join(path_cwd, "data")
    num_days = 198

    percentile = "p50"
    days = list(range(0, num_days + 1, 20))

    days.append(num_days)

    # plot_risk_map(path_results, path_plots, days, percentile)
    plot_icu_map(path_results, path_plots, days, percentile,
                 10, icu_cap[cap_indx])
    regional_proximity = True
    kmax = [x / 100 for x in range(0, 101, 10)]
    loc = 'ICU_Proximity' if regional_proximity else "ICU_Federal_State"
    # path_plots = os.path.join(
    #     path_plots, loc)
    # max_val = [430, 350, 250, 150, 90, 60, 45, 35, 25, 20, 18]
    # if regional_proximity:
    #     max_val = [430, 320, 200, 130, 80, 60, 40, 30, 25, 20, 16]
    # max_val_indx = 0
    # for km in kmax:
    #     km_formatted = f"{km:.2f}"
    #     res_dir = "results"
    #     if regional_proximity:
    #         res_dir += "/new_regional_def"
    #     plot_icu_map(os.path.join(
    #         path_cwd, "results/new_regional_def", "ICUCap_9.000000/rho_1.000000/BlendingFactorRegional_0.000000", f'kmin_0.000000_kmax_{km_formatted}0000'),
    #         os.path.join(path_plots, f"MAP_kmax_{km_formatted}_" + loc), days, percentile, ["FeedbackDamping"],  os.path.join(
    #         path_cwd, "results/new_regional_def", "ICUCap_9.000000/rho_1.000000/BlendingFactorRegional_0.700000", f'kmin_0.000000_kmax_{km_formatted}0000'), max_val[max_val_indx])
    #     max_val_indx += 1
    # plot_r0_map(path_results, path_plots, np.arange(
    #     0, 101, 5), percentile, modes)
    # plot_icu_map(path_results, path_plots, days, percentile, modes)

    # plot_flows(path_results, path_plots, days,
    #            percentile, modes, "Daily_Infections")
