import matplotlib.pyplot as plt
import datetime as dt
import os.path
import imageio

import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

from memilio.epidata import modifyDataframeSeries as mdfs
import memilio.epidata.getPopulationData as gpd
import memilio.plot.plotMap as pm
from memilio.epidata import geoModificationGermany as geoger
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import SymLogNorm, LinearSegmentedColormap, Normalize
from tqdm.auto import tqdm
from datetime import datetime, timedelta
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks
import h5py

import seaborn as sns

start_date = "2020-10-01"
total_pop = 83278910.0
opacity = 0.15
lineWidth = 2.
fontsize = 14
legendsize = 14
ticks = 14
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors_state_risk = sns.color_palette("tab20")
color_map = 'viridis'


def read_total_results_h5(path, comp, group_key='Total'):
    f = h5py.File(path, 'r')
    group = f['0']
    total = group[group_key][()]
    comp_simulated = np.sum(total[:, comp], axis=1)
    f.close()

    return comp_simulated


def read_results_h5(path, comp, group_key='Total'):
    with h5py.File(path, 'r') as f:
        keys = list(f.keys())
        res = np.zeros((len(keys), f[keys[0]][group_key].shape[0]))
        for i, key in enumerate(keys):
            group = f[key]
            total = group[group_key][()]
            if total.shape[1] > 1:
                comp_simulated = np.sum(total[:, comp], axis=1)
            else:
                comp_simulated = total[:, comp]
            res[i] = comp_simulated

    return res


def read_county_results_h5(path, comp, group_key='Total'):
    data = {}
    with h5py.File(path, 'r') as f:
        for key in f.keys():
            group = f[key]
            total = group[group_key][()]
            comp_simulated = np.sum(total[:, comp], axis=1)
            data[key] = comp_simulated
    df = pd.DataFrame(data)
    return df


def get_pop(file_format='h5'):
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

    return population


def create_folder_if_not_exists(path):
    if not os.path.exists(path):
        os.makedirs(path)


def get_dirs_starting_with_x(x, path_results):
    def get_kmax_value(dir_name):
        return float(dir_name.split('_')[3])
    dirs = [d for d in os.listdir(path_results) if os.path.isdir(
        os.path.join(path_results, d)) and d.startswith(x)]
    return sorted(dirs, key=get_kmax_value)


def read_data(path, target_indx, percentile, flows):
    df = read_county_results_h5(os.path.join(
        path, percentile, "Results.h5"), target_indx[0], group_key='Total')
    if flows:
        df = df.diff().iloc[1:].reset_index(drop=True)
    return df


def plot_and_save(fig, path, filename):
    plt.tight_layout()
    fig.savefig(os.path.join(path, filename), dpi=600)
    plt.close(fig)


def get_risk_total(path):
    data_risk = read_results_h5(path, 0, 'Group1')
    pop = get_pop()
    total_pop = pop['Population'].sum()
    # iterate over all rows from data_risk and scale the risk with the population
    risk_total = np.zeros(data_risk.shape[1])
    for i in range(data_risk.shape[0]):
        risk_total += data_risk[i, :] * \
            pop['Population'].values[i] / total_pop
    return risk_total


def plot_peak_values_splitted(data, path_plots, modes, target_indx, percentile="p50", flows=True, plot_type='day', title='Peak Value', vertical=False, log_scale=False):
    if len(modes) > 1:
        print("Only one mode is allowed for the grid peak plot.")
        return

    plot_dir = os.path.join(path_plots, "peaks", title)
    create_folder_if_not_exists(plot_dir)

    num_subplots = 10
    allowed_kmax = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]  # , 0.8, 0.9]

    # delete all entries with kmax not in allowed_kmax. data is a df
    data = data[data['kmax'].isin(allowed_kmax)]

    # Adjusting the width
    subplot_width = 1.4  # ~88mm
    subplot_height = 1.0  # Keeping height fixed

    num_rows = int(len(allowed_kmax) / 2)
    num_cols = 2

    fig_width = num_cols * subplot_width
    fig_height = num_rows * subplot_height
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(num_rows, num_cols, figure=fig)
    tiks_size = 5
    all_peaks = []
    for fig_index in range(len(allowed_kmax)):

        if plot_type == 'day':
            row = fig_index % num_rows
            col = fig_index // num_rows
            ax = fig.add_subplot(gs[row, col])
            peaks_each_day = data[data['kmax'] == allowed_kmax[fig_index]]
            peaks_each_day = peaks_each_day.explode('peak_value')['day'].reset_index(
                drop=True)
            peaks_each_day = peaks_each_day.value_counts().sort_index()
            full_index = pd.RangeIndex(start=0, stop=201, step=1)
            peaks_each_day = peaks_each_day.reindex(full_index, fill_value=0)
            ax.plot(peaks_each_day, color=colors[0], linewidth=lineWidth)
            all_peaks.append(peaks_each_day)
            # if row == 2 and col == 0:
            # ax.set_ylabel('Number of Peaks', fontsize=fontsize)

        else:
            peaks_each_day = data[data['kmax'] == allowed_kmax[fig_index]].explode(
                'peak_value').reset_index(drop=True)
            # delete all cols except day and peak_value
            peaks_each_day = peaks_each_day[['day', 'peak_value']]

            # Creating the boxplot
            row = fig_index % num_rows
            col = fig_index // num_rows
            ax = fig.add_subplot(gs[row, col])
            flierprops = dict(marker='')
            peaks_each_day.boxplot(
                column="peak_value",
                by="day",
                grid=False,
                positions=peaks_each_day['day'].unique(),
                flierprops=flierprops,
                ax=ax
            )
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title('')
            # ax.set_xlim(80, 180)
            # ax.set_xticks(range(80, 180, 20))
            # ax.set_xticklabels(range(80, 180, 20), fontsize=tiks_size)
            plt.suptitle("")
            all_peaks.append(peaks_each_day)

            # if row == 2 and col == 0:
            #     ax.set_ylabel(title + ' per 100,000k', fontsize=fontsize)

            if row == 4:
                ax.set_xlabel('Days', fontsize=8)

        # ax.set_title(
        #     rf'$\psi_{{max}}$: {allowed_kmax[fig_index]:.1f}', fontsize=fontsize)

    global_max = float('-inf')
    global_min = float('inf')
    if plot_type == 'day':
        global_min = min(min(peaks) for peaks in all_peaks if peaks.size)
        global_max = max(max(peaks)
                         for peaks in all_peaks if peaks.size) + 10
        global_min = 0
        global_max = 100
    else:
        # Iterate over each DataFrame in all_peaks to find the global max and min
        for df in all_peaks:
            local_max = df['peak_value'].max()
            local_min = df['peak_value'].min()
            global_max = max(global_max, local_max)
            global_min = min(global_min, local_min)
        global_min = 1
        global_max = 1000  # global_max + 100

    count_ax = 0
    xlim_min = 0
    xlim_max = 140
    for ax in fig.get_axes():
        if count_ax > num_subplots:
            break
        ax.grid()
        if log_scale:
            ax.set_yscale('symlog')
        ax.set_ylim(global_min, global_max)
        ax.tick_params(axis='x', labelsize=tiks_size)
        ax.tick_params(axis='y', labelsize=tiks_size)
        ax.set_xlim(xlim_min, xlim_max)
        ax.set_xticks(range(xlim_min, xlim_max, 20))
        ax.set_xticklabels(range(xlim_min, xlim_max, 20), fontsize=tiks_size)
        if vertical:
            if ax == fig.get_axes()[int(len(allowed_kmax) / 2) - 1] or ax == fig.get_axes()[int(len(allowed_kmax) - 1)]:
                ax.set_xlabel('Days', fontsize=8)

        count_ax += 1

    fn = f'splitted_peaks_grid_{plot_type}.png'
    if log_scale:
        fn = f'splitted_peaks_grid_{plot_type}_log_scale.png'
    plot_and_save(
        fig, plot_dir, fn)


def plot_peak_values(path_results, path_plots, modes, target_indx, percentile="p50", flows=True, plot_type='kmin', title='Peak Value', vertical=False, dir_value='kmin', peak_prominence=1):
    if len(modes) > 1:
        print("Only one mode is allowed for the grid peak plot.")
        return

    plot_dir = os.path.join(path_plots, "peaks", title)
    create_folder_if_not_exists(plot_dir)

    allowed_kmax = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    kmin_dirs = get_dirs_starting_with_x(dir_value, path_results)
    all_kmins = 0
    if dir_value == 'kmin':
        all_kmins = sorted(set(float(d.split("_")[1]) for d in kmin_dirs))
    else:
        all_kmins = sorted(set(float(d.split("_")[-1]) for d in kmin_dirs))

    num_days = 0

    for index, kmin in enumerate(all_kmins):
        if dir_value == 'fixed':
            if kmin > all_kmins[0]:
                continue

        num_subplots = round(11 - index)
        if dir_value == 'fixed':
            num_subplots = 11
        subplot_width = 8  # Breite eines einzelnen Subplots
        subplot_height = 3  # Höhe eines einzelnen Subplots

        if vertical:
            fig_height = num_subplots * subplot_height
            fig_width = subplot_width
            fig = plt.figure(figsize=(fig_width, fig_height))
            gs = GridSpec(num_subplots, 1, figure=fig)
        else:
            fig_width = num_subplots * subplot_width
            fig_height = subplot_height
            fig = plt.figure(figsize=(fig_width, fig_height))
            gs = GridSpec(1, num_subplots, figure=fig)

        all_peaks = []
        all_risks = []

        for run in kmin_dirs:
            if dir_value == 'kmin':
                if kmin != float(run.split("_")[1]):
                    continue
            kmax = float(run.split("_")[3])
            if kmax not in allowed_kmax:
                continue

            path = os.path.join(path_results, run,
                                modes[0], "flows" if flows else "")
            df = read_data(path, target_indx, percentile, flows)

            if num_days == 0:
                num_days = df.shape[0]

            # smooth the data from the df and find the peaks
            def get_peak(df):
                smoothed = df.rolling(window=7, min_periods=1).mean()
                # iterate over each column and find the peaks
                num_counties_no_peak = 0
                peaks = []
                for col in smoothed.columns:
                    peak_indices = find_peaks(
                        smoothed[col].values, prominence=peak_prominence, distance=21)[0]
                    if len(peak_indices) == 0:
                        peak_indices = find_peaks(
                            smoothed[col].values, prominence=0, distance=21)[0]
                    peak_values = smoothed[col].values[peak_indices]
                    if len(peak_indices) == 0:
                        num_counties_no_peak += 1
                        continue

                    peaks.append((peak_indices, peak_values))
                if num_counties_no_peak > 0:
                    print("Number of counties without peak: ",
                          num_counties_no_peak)
                return peaks

            # peaks contains the peak indices and the peak values for each county
            peaks = get_peak(df)

            if plot_type == 'kmin':
                # peaks = np.bincount(df.idxmax().values, minlength=df.shape[0])
                # calculate the peaks per day out of peaks
                peaks_each_day = np.zeros(num_days)
                for peak_indices, _ in peaks:
                    for peak_index in peak_indices:
                        peaks_each_day[peak_index] += 1
                all_peaks.append(peaks_each_day)
                if vertical:
                    ax = fig.add_subplot(gs[len(all_peaks) - 1, 0])
                else:
                    ax = fig.add_subplot(gs[0, len(all_peaks) - 1])
                ax.plot(peaks_each_day, color=colors[0], linewidth=lineWidth)
                if len(all_peaks) == 1:
                    ax.set_ylabel('Number of Peaks')
            else:
                peak_indices = df.idxmax().values
                # iterate over peaks and scale the peak_value
                for _, p_values in peaks:
                    for i in range(len(p_values)):
                        p_values[i] = p_values[i] / \
                            get_pop()['Population'].values[i] * 100_000

                # create a  empty dict with the peak values for each day
                peak_values = {i: [] for i in range(df.shape[0])}
                for idx, val in peaks:
                    for i, v in zip(idx, val):
                        peak_values[i].append(v)

                # convert the dict to a dict to use sns.boxplot
                df_long = pd.DataFrame(
                    [{'Day': day, 'Peak Value': val} for day, vals in peak_values.items() for val in vals])
                all_peaks.extend(peak_values.values())
                all_days = pd.DataFrame({'Day': range(num_days)})
                df_long = pd.merge(all_days, df_long, on='Day', how='left')
                if vertical:
                    ax = fig.add_subplot(gs[len(all_peaks) // len(df) - 1, 0])
                else:
                    ax = fig.add_subplot(gs[0, len(all_peaks) // len(df) - 1])
                sns.boxplot(x='Day', y='Peak Value', data=df_long, ax=ax)
                ax.set_ylabel(title + ' per 100,000k')

            if modes[0] == 'FeedbackDamping':
                risk = get_risk_total(os.path.join(
                    path_results, run, "FeedbackDamping", "risk", percentile, "Results.h5"))
                all_risks.append(risk)

            kmax = float(run.split("_")[3])
            ax.set_title(
                rf'$\psi_{{min}}$: {kmin}, $\psi_{{max}}$: {kmax:.1f}', fontsize=fontsize)

            if dir_value == 'fixed':
                ax.set_title(
                    f'kmin: { float(run.split("_")[-1])}', fontsize=fontsize)
        global_min = 0
        global_max = 0
        if plot_type == 'kmin':
            global_min = min(min(peaks) for peaks in all_peaks if peaks.size)
            global_max = max(max(peaks)
                             for peaks in all_peaks if peaks.size) + 10
        else:
            global_min = min(min(peaks) for peaks in all_peaks if peaks)
            global_max = max(max(peaks)
                             for peaks in all_peaks if peaks) + 10

        count_ax = 0
        for ax in fig.get_axes():
            if count_ax > num_subplots:
                break
            # grid
            ax.grid()
            # log scale
            ax.set_yscale('symlog')
            ax.set_ylim(global_min, global_max)
            ax.set_xticks(range(0, num_days, 50))
            ax.tick_params(axis='x', labelsize=ticks)
            ax.tick_params(axis='y', labelsize=ticks)
            if not vertical:
                ax.set_xlabel('')
                if ax != fig.get_axes()[0]:
                    ax.set_yticklabels([])
                    ax.set_ylabel('')
                if ax == fig.get_axes()[len(fig.get_axes()) // 2]:
                    ax.set_xlabel('Days', fontsize=fontsize)
            else:
                ax.set_xlabel('')
                ax.set_ylabel('')
                if ax == fig.get_axes()[len(fig.get_axes()) // 2]:
                    ax.set_ylabel(title + ' per 100,000k', fontsize=fontsize)
                if ax == fig.get_axes()[len(fig.get_axes()) - 1]:
                    ax.set_xlabel('Days', fontsize=fontsize)

            count_ax += 1

        plot_and_save(
            fig, plot_dir, f'peaks_grid_{plot_type}_{dir_value}_{kmin}.pdf')


def get_peak_data(path_results, blending_fact, modes, target_indx, percentile="p50", flows=True, title='Peak Value', dir_type='peak', filter=None, peak_prominence=1):
    if len(modes) > 1:
        print("Only one mode is allowed for the grid peak plot.")
        return

    population = get_pop()

    data = {'b_fact': [], 'kmax': [], 'day': [], 'peak_value': []}
    for b_fact in blending_fact:
        path_results_blended = os.path.join(
            path_results, f"BlendingFactorRegional_{b_fact:.2f}" + "0000")

        kmin_dirs = get_dirs_starting_with_x('kmin', path_results_blended)
        all_kmax = sorted(set(float(d.split("_")[-1]) for d in kmin_dirs))

        num_days = 0

        for index, kmax in enumerate(all_kmax):
            path = os.path.join(path_results_blended, kmin_dirs[index],
                                modes[0], "flows" if flows else "")
            df = read_data(path, target_indx, percentile, flows)

            # columns as integers
            df.columns = df.columns.astype(int)
            # sort the columns
            df = df.reindex(sorted(df.columns), axis=1)

            # apply filter if not None
            if filter is not None:
                # Flatten the filter list to match the column names
                filter_flat = [str(item[0]) for item in filter]
                df = df[filter_flat]

            if num_days == 0:
                num_days = df.shape[0]

            # smooth the data from the df and find the peaks
            def get_peak(df, population=population):
                smoothed = df.rolling(window=7, min_periods=1).mean()
                # iterate over each column and find the peaks
                num_counties_no_peak = 0
                peaks = []
                for col in smoothed.columns:
                    peak_indices = find_peaks(
                        smoothed[col].values, prominence=peak_prominence, distance=21)[0]
                    # if len(peak_indices) == 0:
                    #     peak_indices = find_peaks(
                    #         smoothed[col].values, prominence=0, distance=21)[0]
                    peak_values = smoothed[col].values[peak_indices]
                    if len(peak_indices) == 0:
                        num_counties_no_peak += 1
                        continue

                    # scale the peak values per 100_000 inhabitants
                    # search which index from population['ID_County] is equal to col
                    idx = population[population['ID_County'] == col].index[0]
                    pop = population['Population'].values[idx]
                    peak_values = peak_values / pop * 100_000

                    peaks.append((peak_indices, peak_values))
                if num_counties_no_peak > 0:
                    print("Number of counties without peak: ",
                          num_counties_no_peak)
                return peaks

            # peaks contains the peak indices and the peak values for each county
            if dir_type == 'peak':
                peaks = get_peak(df)

                peak_values = {i: [] for i in range(df.shape[0])}
                for idx, val in peaks:
                    for i, v in zip(idx, val):
                        peak_values[i].append(v)

                # delete all empty entries
                peak_values = {k: v for k,
                               v in peak_values.items() if len(v) > 0}

            elif dir_type == 'max_val':
                # get max value and index for each column
                max_vals_df = df.max()
                max_indices_df = df.idxmax()

                # scale max_vals_df per 100_000
                max_vals_df = max_vals_df / \
                    population['Population'].values * 100_000

                peak_values = {i: [] for i in range(df.shape[0])}
                for i in range(len(max_vals_df)):
                    peak_values[max_indices_df.values[i]].append(
                        max_vals_df.values[i])

                # delete all empty entries
                peak_values = {k: v for k,
                               v in peak_values.items() if len(v) > 0}
            else:
                print("dir_type not known")
                return

            for day, values in peak_values.items():
                data['b_fact'].append(b_fact)
                data['kmax'].append(kmax)
                data['day'].append(day)
                data['peak_value'].append(values)

    # Convert dictionary to DataFrame
    df_peaks = pd.DataFrame(data)

    return df_peaks


def plot_peak_time_kmax_variance(path_results, path_plots, blending_fact, modes, target_indx, rhos, percentile="p50", flows=True, title='Peak Value'):

    # filter_counties = pd.read_csv(
    #     '/localdata1/code_2024/memilio/tools/border_counties.txt')
    # filter_counties = filter_counties.values.tolist()
    filter_counties = None

    for rho in rhos:
        path_rho = os.path.join(
            path_results, f"rho_{rho:.6f}")

        data = []
        if title == 'Peak Value':
            data = get_peak_data(path_rho, blending_fact,
                                 modes, target_indx, percentile, flows, title, 'peak', filter=filter_counties)
        elif title == 'Max Value':
            data = get_peak_data(path_rho, blending_fact,
                                 modes, target_indx, percentile, flows, title, 'max_val', filter=filter_counties)
        else:
            print('Title not supported')
            return

        # Process data for 'kmin'
        all_b_fact_kmin = data['b_fact'].unique()
        all_kmax_kmin = data['kmax'].unique()
        # delete 0.36 from all_kmax_kmin
        all_kmax_kmin = all_kmax_kmin[all_kmax_kmin != 0.36]
        data_variance_kmin = {'b_fact': [], 'kmax': [],
                              'mean': [], 'std': [], 'variance': []}

        for b_fact in all_b_fact_kmin:
            for kmax in all_kmax_kmin:
                peak_data = data[(data['b_fact'] == b_fact) & (
                    data['kmax'] == kmax)]
                # i want a list of all peak_days. So  just the col 'day'. However, if there are n peak_values, the day should be n times in the list
                peak_day_list = []
                for i in range(len(peak_data)):
                    entry = peak_data.iloc[i]
                    for _ in range(len(entry['peak_value'])):
                        peak_day_list.append(entry['day'])
                # to np array
                peak_day_list = np.array(peak_day_list)
                mean_peak = peak_day_list.mean()
                std_peak = peak_day_list.std()
                var_peak = peak_day_list.var()

                data_variance_kmin['b_fact'].append(b_fact)
                data_variance_kmin['kmax'].append(kmax)
                data_variance_kmin['mean'].append(mean_peak)
                data_variance_kmin['std'].append(std_peak)
                data_variance_kmin['variance'].append(var_peak)

        # Process data for 'val'
        # TODO: Noch prüfen, obs jetzt korrekt ist!!
        data_variance_val = {'b_fact': [], 'kmax': [],
                             'mean': [], 'std': [], 'variance': []}

        for b_fact in all_b_fact_kmin:
            for kmax in all_kmax_kmin:
                peak_data = data[(data['b_fact'] == b_fact)
                                 & (data['kmax'] == kmax)]
                peak_data = peak_data[peak_data['peak_value'].apply(
                    lambda x: len(x) > 0)].reset_index(drop=True)
                peak_data = peak_data.explode('peak_value')
                peak_data['peak_value'] = peak_data['peak_value'].astype(float)

                mean_peak = peak_data['peak_value'].values.mean()
                std_peak = peak_data['peak_value'].values.std()
                var_peak = peak_data['peak_value'].values.var()

                data_variance_val['b_fact'].append(b_fact)
                data_variance_val['kmax'].append(kmax)
                data_variance_val['mean'].append(mean_peak)
                data_variance_val['std'].append(var_peak)
                data_variance_val['variance'].append(var_peak)

        # Plot and save each heatmap separately
        size_ticks = 18

       # transform b_fact in interval [0, 100]
        max_b_fact = 0.7
        data_variance_kmin['b_fact'] = [
            round((b_fact * 100 / max_b_fact)) for b_fact in data_variance_kmin['b_fact']
        ]
        data_variance_val['b_fact'] = [
            round((b_fact * 100 / max_b_fact)) for b_fact in data_variance_val['b_fact']
        ]

        # same for kmax
        data_variance_kmin['kmax'] = [
            round((kmax * 100 / 1.0)) for kmax in data_variance_kmin['kmax']
        ]
        data_variance_val['kmax'] = [
            round((kmax * 100 / 1.0)) for kmax in data_variance_val['kmax']
        ]

        results_df_val = pd.DataFrame(data_variance_val)
        results_df_kmin = pd.DataFrame(data_variance_kmin)

        # delete all entries with kmax == 0, 2, or 4
        results_df_val = results_df_val[results_df_val['kmax'] != 0]
        results_df_val = results_df_val[results_df_val['kmax'] != 2]
        # results_df_val = results_df_val[results_df_val['kmax'] != 4]

        results_df_kmin = results_df_kmin[results_df_kmin['kmax'] != 0]
        results_df_kmin = results_df_kmin[results_df_kmin['kmax'] != 2]
        # results_df_kmin = results_df_kmin[results_df_kmin['kmax'] != 4]

        xlabel_title = 'Impact Federal State'  # 'Impact Regional Proximity'  #

        statistical_measures = ['mean', 'std']
        for measure in statistical_measures:

            cmap = 'viridis' if measure == 'mean' else 'plasma'
            # Plot for 'kmin'
            log_scale = False
            fig_kmin, ax_kmin = plt.subplots(
                figsize=(10, 8))  # Increase figure size
            # results_df_kmin = pd.DataFrame(data_variance_kmin)
            heatmap_data_kmin = results_df_kmin.pivot_table(
                values=measure, index='kmax', columns='b_fact', aggfunc=lambda x: x)
            # invert the y-axis
            heatmap_data_kmin = heatmap_data_kmin.iloc[::-1]
            # save heatmap data to csv
            # heatmap_data_kmin.to_csv(os.path.join(
            #     path_plots, f'{title}_{measure}_heatmap_kmin_{rho}_{xlabel_title}.csv'))
            sns.heatmap(heatmap_data_kmin, cmap=cmap, ax=ax_kmin, norm=SymLogNorm(
                linthresh=1) if log_scale else None, cbar=True)
            # cbar = ax_kmin.collections[0].colorbar
            # cbar.locator = MaxNLocator(nbins=10)  # Erlaubt bis zu 10 Ticks
            # cbar.update_ticks()
            # ax_kmin.set_xlabel(xlabel_title, fontsize=fontsize)
            # ax_kmin.set_ylabel(rf'$\psi_{{max}}$', fontsize=fontsize)
            # deactive x and y labels
            ax_kmin.set_xlabel('')
            ax_kmin.set_ylabel('')
            ax_kmin.tick_params(axis='x', labelsize=size_ticks)
            ax_kmin.tick_params(axis='y', labelsize=size_ticks)

            # increase ticks size colorbar
            cbar = ax_kmin.collections[0].colorbar
            cbar.ax.tick_params(labelsize=size_ticks)

            # aktuell werden die labels in allowed_labels_x  angezeigt. DIe labels sollen aber durch new_labels_x ersetzt werden
            # Set the ticks explicitly to match the allowed positions
            # ax_kmin.set_xticks(transformed_current_labels)

            # Set the new labels for these ticks
            # ax_kmin.set_xticklabels([str(label)
            #                         for label in transformed_current_labels])

            allowed_labels_x = [0.0, 14, 28, 42, 56, 70]
            # new_labels_x = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
            allowed_labels = [0,  10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
            for n, label in enumerate(ax_kmin.xaxis.get_ticklabels()):
                if float(label.get_text()) not in allowed_labels:
                    label.set_visible(False)

            for n, label in enumerate(ax_kmin.yaxis.get_ticklabels()):
                if float(label.get_text()) not in allowed_labels:
                    label.set_visible(False)

            plt.yticks(rotation=0)  # Keep the labels horizontal
            plt.xticks(rotation=0)

            plt.savefig(os.path.join(
                path_plots, f'{title}_{measure}_heatmap_kmin_{rho}.pdf'), dpi=600)
            plt.clf()
            plt.close(fig_kmin)

            # Plot for 'val'
            log_scale = False
            fig_val, ax_val = plt.subplots(
                figsize=(10, 8))  # Increase figure size
            # results_df_val = pd.DataFrame(data_variance_val)
            heatmap_data_val = results_df_val.pivot_table(
                values=measure, index='kmax', columns='b_fact', aggfunc=lambda x: x)
            # save heatmap data to csv
            heatmap_data_val.to_csv(os.path.join(
                path_plots, f'{title}_{measure}_heatmap_val_{rho}_{xlabel_title}.csv'))
            # invert the y-axis
            heatmap_data_val = heatmap_data_val.iloc[::-1]
            sns.heatmap(heatmap_data_val, cmap='plasma', ax=ax_val, norm=SymLogNorm(
                linthresh=1) if log_scale else None, cbar=True)
            # cbar = ax_val.collections[0].colorbar
            # num_ticks = 10
            # ticks = np.linspace(cbar.vmin, cbar.vmax, num_ticks)
            # cbar.set_ticks(ticks)
            # ax_val.set_xlabel(xlabel_title, fontsize=fontsize)
            # ax_val.set_ylabel(rf'$\psi_{{max}}$', fontsize=fontsize)
            # deactive x and y labels
            ax_val.set_xlabel('')
            ax_val.set_ylabel('')
            ax_val.tick_params(axis='x', labelsize=size_ticks)
            ax_val.tick_params(axis='y', labelsize=size_ticks)

            # deactive x and y labels
            ax_kmin.set_xlabel('')
            ax_kmin.set_ylabel('')

            # increase ticks size colorbar
            cbar = ax_val.collections[0].colorbar
            cbar.ax.tick_params(labelsize=size_ticks)

            for n, label in enumerate(ax_val.xaxis.get_ticklabels()):
                if float(label.get_text()) not in allowed_labels:
                    label.set_visible(False)

            for n, label in enumerate(ax_val.yaxis.get_ticklabels()):
                if float(label.get_text()) not in allowed_labels:
                    label.set_visible(False)

            plt.yticks(rotation=0)  # Keep the labels horizontal
            plt.xticks(rotation=0)

            plt.savefig(os.path.join(
                path_plots, f'{title}_{measure}_heatmap_val_{rho}.pdf'), dpi=600)
            plt.clf()
            plt.close(fig_val)

    return 0


def plot_peaks_all_blending_factors(path_results, path_plots, blending_fact, log_scale=True, modes=["FeedbackDamping"]):
    icu_compartment = [[7]]
    data = get_peak_data(path_results, blending_fact,
                         modes, icu_compartment, "p50", False, 'Peak Value', 'peak')
    for b_fact in blending_fact:
        data_bfact = data[data['b_fact'] == b_fact]
        b_fact_add = str(b_fact) + "00000"
        if len(b_fact_add) == 9:
            b_fact_add = b_fact_add[:8]
        path_plot_b_fact = os.path.join(
            path_plots, "BlendingFactorRegional_" + b_fact_add)
        create_folder_if_not_exists(path_plot_b_fact)
        plot_peak_values_splitted(data_bfact, path_plot_b_fact, [
            "FeedbackDamping"], icu_compartment, plot_type='day', title='ICU Occupancy', flows=False, vertical=True, log_scale=log_scale)
        plot_peak_values_splitted(data_bfact, path_plot_b_fact, [
            "FeedbackDamping"], icu_compartment, plot_type='val', title='ICU Occupancy', flows=False, vertical=True, log_scale=log_scale)


def violin_plot_all_kmax(path_results, path_plots, modes, blending_fact, kmaxs, target_indx, percentile="p50", flows=False, title='Peak Value', dir_type='kmin', log_scale=True, rho=1.0):
    # Retrieve peak data using the get_peak_data function
    data = get_peak_data(path_results, blending_fact,
                         modes, target_indx, percentile, flows, title, dir_type)

    # Filter data for specific blending factors and kmax values
    data = data[data['b_fact'].isin(blending_fact)]
    data = data[data['kmax'].isin(kmaxs)]

    ticks_size = 8

    for b_fact in blending_fact:
        peak_data = data[data['b_fact'] == b_fact]

        # explode peak_value
        peak_data = peak_data.explode('peak_value')

        # Reduce the dimension of peak_value
        peak_data = peak_data[peak_data['peak_value']
                              > 0.0].reset_index(drop=True)
        peak_data = peak_data.loc[peak_data.index.repeat(
            peak_data['peak_value'])].copy().reset_index(drop=True)
        peak_data.drop(columns=['peak_value'], inplace=True)

        # gehe durch alle Einträge und gucke ob ein eintrag mit day kleiner 10 existiert, derren peak_value > 0 ist. Setze diesen auf 0
        # for i in range(len(peak_data)):
        #     if peak_data.loc[i, 'day'] < 10:
        #         peak_data.loc[i, 'day'] = 0

        plt.figure(figsize=(4., 4.))
        sns.violinplot(x='kmax', y='day', data=peak_data,
                       inner="box", linewidth=1, palette='Blues', width=0.8)

        ax = plt.gca()
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10))

        # i want more ticks on the y-axis
        # ax.yaxis.set_major_locator(MaxNLocator(nbins=10))

        # set ticks to size ticks_size
        ax.tick_params(axis='x', labelsize=ticks_size)
        ax.tick_params(axis='y', labelsize=ticks_size)

        # Set labels and title
        # Using raw string for latex
        plt.xlabel(r'$\psi_{\max}$', fontsize=fontsize)
        plt.ylabel('Day', fontsize=fontsize)
        # title = 'Impact Regional Proximity'  # 'Impact Federal State'  #
        # plt.title(title + f': {b_fact}', fontsize=18)

        # grid
        plt.grid(axis='y', linestyle='--', alpha=0.6)

        # Adjust layout to prevent cutting off
        plt.tight_layout()

        # set ylim between 50 und 190
        plt.ylim(30, 210)
        # y axis log scale
        # plt.yscale('log')

        # Save the plot as a PDF with tight bounding box
        pdf_filename = os.path.join(
            path_plots, f'violin_plot_{b_fact}_rho_{rho}.pdf')
        plt.savefig(pdf_filename, format='pdf', bbox_inches='tight')

        # Clear the plot after saving to avoid overlap when looping
        plt.close()


def plot_overall_contacts_and_icu_overload(list_scenario, path_plots, icu_comp, labels, log_scale=False, kmax="", percentile="p50"):
    population = get_pop()
    data_risk = {}
    data_icu = {}
    data_icu_total = {}
    for scenario in list_scenario:
        ############### contact reduction ###############
        path_results_risk = os.path.join(
            scenario, "risk", percentile, "Results.h5")
        data = {}
        with h5py.File(path_results_risk, 'r') as f:
            for key in f.keys():
                group = f[key]
                total = group['Group1'][()]
                data[key] = np.sum(total[:], axis=1)
        df = pd.DataFrame(data)

        # transform risk to contact reduction
        ContactReductionMax = float(kmax)
        ContactReductionMin = 0.0
        EpsilonContacts = 0.1

        def calculate_reduc_fac(perceived_risk_contacts):
            reduc_fac = (
                (ContactReductionMax - ContactReductionMin) *
                EpsilonContacts *
                np.log(np.exp(perceived_risk_contacts / EpsilonContacts) + 1.0) +
                ContactReductionMin
            )
            return reduc_fac

        # transform every value in df to contact reduction and get mean, min, max per day (row)
        df = df.applymap(calculate_reduc_fac)
        df['mean'] = df.mean(axis=1)
        df['min'] = df.min(axis=1)
        df['max'] = df.max(axis=1)
        # df = df[['mean', 'min', 'max']]
        data_risk[scenario] = df

        ######################### ICU overload #########################
        path_results_icu = os.path.join(
            scenario, percentile, "Results.h5")
        data = {}
        with h5py.File(path_results_icu, 'r') as f:
            for key in f.keys():
                group = f[key]
                total = group['Total'][()]
                icu_occ = total[1:, icu_comp[0]]
                # search pop in key in population
                population_key = population[population['ID_County'] == int(
                    key)]['Population'].reset_index(drop=True)[0]
                # calc local icu capacity
                icu_cap_per_100k = 9
                icu_cap_key = icu_cap_per_100k * population_key / 100_000
                # calculate overload
                overload = icu_occ - icu_cap_key
                # set all negative values to 0
                overload[overload < 0] = 0
                data[key] = overload.flatten()
        df_icu = pd.DataFrame(data)
        df_icu['sum'] = df_icu.sum(axis=1)
        # df_icu = df_icu[['sum']]
        data_icu[scenario] = df_icu

    # get total icu occupancy
    total_pop = population['Population'].sum()
    for scenario in list_scenario:
        path_results_icu = os.path.join(
            scenario, percentile, "Results_sum.h5")
        data = read_total_results_h5(path_results_icu, icu_comp).flatten()
        data = data / total_pop * 100_000
        data_icu_total[scenario] = data

    # Create a single figure with two subplots (side by side)
    # fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 6), sharex=False)
    fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(16, 6), sharex=False)

    # Plot Contact Reduction (left plot)
    for indx, scenario in enumerate(list_scenario):
        label_scenario = labels[indx]
        ax1.plot(
            data_risk[scenario]['mean'],
            label=label_scenario,
            color=colors[indx],
            linestyle='-', linewidth=lineWidth
        )
    # ax1.set_title('Contact Reduction', fontsize=fontsize)
    ax1.set_ylabel('Contact Reduction', fontsize=fontsize)
    # ax1.set_xlabel('Day', fontsize=fontsize)
    # ax1.legend()
    ax1.grid(linestyle='--', alpha=0.6)
    ax1.tick_params(axis='x', labelsize=ticks)
    ax1.tick_params(axis='y', labelsize=ticks)
    if log_scale:
        ax1.set_yscale('log')

    # Plot ICU Overload (middle plot)
    # for indx, scenario in enumerate(list_scenario):
    #     label_scenario = labels[indx]
    #     ax2.plot(
    #         data_icu[scenario]['sum'],
    #         label=label_scenario,
    #         color=colors[indx],
    #         linestyle='-', linewidth=lineWidth
    #     )
    # ax2.set_title('ICU Overload', fontsize=fontsize)
    # ax2.set_ylabel('ICU Overflow', fontsize=fontsize)
    # ax2.set_xlabel('Day', fontsize=fontsize)
    # # ax2.legend(fontsize=legendsize, loc='upper left')
    # # save legend from ax2 in separate file
    # fig_legend = plt.figure(figsize=(2, 1))
    # ax_legend = fig_legend.add_subplot(111)
    # ax_legend.legend(*ax2.get_legend_handles_labels(), fontsize=legendsize)
    # ax_legend.axis('off')
    # fig_legend.savefig(os.path.join(
    #     path_plots, f'diff_scales_kmax_{kmax}_legend.pdf'), bbox_inches='tight')
    # ax2.legend().set_visible(False)
    # ax2.grid(linestyle='--', alpha=0.6)
    # if log_scale:
    #     ax2.set_yscale('log')
    # ax2.tick_params(axis='x', labelsize=ticks)
    # ax2.tick_params(axis='y', labelsize=ticks)

    # Plot Total ICU Occupancy (right plot)
    for indx, scenario in enumerate(list_scenario):
        label_scenario = labels[indx]
        ax3.plot(
            data_icu_total[scenario],
            label=label_scenario,
            color=colors[indx],
            linestyle='-', linewidth=lineWidth
        )
    # ax3.set_title('Total ICU Occupancy', fontsize=14)
    ax3.set_ylabel('ICU Occupancy per 100,000', fontsize=fontsize)
    # ax3.set_xlabel('Day', fontsize=12)
    # ax3.legend()
    ax3.grid(linestyle='--', alpha=0.6)
    if log_scale:
        ax3.set_yscale('log')

    # ylim von 0.1 bis 10
    # ax3.set_ylim(0.1, 100)

    ax3.tick_params(axis='x', labelsize=ticks)
    ax3.tick_params(axis='y', labelsize=ticks)

    # Adjust layout and save to a single file
    plt.tight_layout()
    plt.savefig(os.path.join(path_plots, f'diff_scales_kmax_{kmax}.pdf'))
    plt.close(fig)


if __name__ == '__main__':
    path_cwd = os.getcwd()
    dirs_rp = ["results/Sep_new_regional_def",
               "results/Nov_new_regional_def", "results/Dez_new_regional_def"]

    dirs_fs = ["results/sim_september",
               "results/sim_november", "results/sim_december"]
    plots_add = ["/september", "/november", "/december"]

    for dir_rp, dir_fs, plot_add in zip(dirs_rp, dirs_fs, plots_add):
        icu_compartment = [[7]]

        # skip first iteration and second
        if dir_rp == "results/Sep_new_regional_def" or dir_rp == "results/Nov_new_regional_def":
            continue

        path_plots = os.path.join(
            path_cwd, "plots" + plot_add, "ICUCap_9.000000")

        # ---------- FIGURE 1 ---------------
        plot_peaks_all_blending_factors(
            os.path.join(path_cwd, dir_rp, "ICUCap_9.000000", "rho_1.000000"), path_plots, [0.3], log_scale=True)

        # # ---------- FIGURE 2 ---------------
        # all_kmax = [x / 100 for x in range(0, 105, 4)]
        # rhos = [0.92, 1.0, 1.08]
        # for rho in rhos:
        #     path_results_rho = os.path.join(
        #         path_cwd, dir_fs, "ICUCap_9.000000", f"rho_{rho:.6f}")
        #     violin_plot_all_kmax(path_results_rho, path_plots, ["FeedbackDamping"], [0.3], all_kmax, icu_compartment,
        #                          flows=False, title='Peak Value', dir_type='peak', log_scale=True, rho=rho)

        # # ---------- FIGURE 3 ---------------
        # plot_peak_time_kmax_variance(os.path.join(path_cwd, dir_rp, "ICUCap_9.000000"), path_plots, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], ["FeedbackDamping"],
        #                              icu_compartment, [1.0], flows=False)

        # ---------- FIGURE 4 ---------------
        # kmaxs_scenarios = ['0.36', '0.80']
        # labels = [
        #     'Regional prox.', 'Federal state', 'National']
        # for kmax_scenario in kmaxs_scenarios:
        #     list_scenario = [
        #         # os.path.join(path_cwd, "results/new_regional_def_735", "ICUCap_9.000000", "rho_1.000000", "BlendingFactorRegional_0.700000", f"kmin_0.000000_kmax_{kmax}0000", "FeedbackDamping"),
        #         os.path.join(
        #             path_cwd, dir_rp, "ICUCap_9.000000", "rho_1.000000", "BlendingFactorRegional_0.700000", f"kmin_0.000000_kmax_{kmax_scenario}0000", "FeedbackDamping"),
        #         os.path.join(
        #             path_cwd, dir_fs, "ICUCap_9.000000", "rho_1.000000", "BlendingFactorRegional_0.700000", f"kmin_0.000000_kmax_{kmax_scenario}0000", "FeedbackDamping"),
        #         os.path.join(
        #             path_cwd, dir_fs, "ICUCap_9.000000", "rho_1.000000", "BlendingFactorRegional_0.000000", f"kmin_0.000000_kmax_{kmax_scenario}0000", "FeedbackDamping")]
        #     plot_overall_contacts_and_icu_overload(list_scenario, path_plots, [
        #                                            [7]], labels, log_scale=False, kmax=kmax_scenario)
