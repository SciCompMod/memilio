import os
import matplotlib.pyplot as plt
import datetime as dt
import os.path
import imageio

import numpy as np
import pandas as pd
import copy

import matplotlib.dates as mdates
from memilio.epidata import modifyDataframeSeries as mdfs
import memilio.epidata.getPopulationData as gpd
import memilio.plot.plotMap as pm
from memilio.epidata import geoModificationGermany as geoger
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import SymLogNorm, LinearSegmentedColormap
from tqdm.auto import tqdm
from datetime import datetime, timedelta
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks
import h5py

import seaborn as sns


# sns.set_style("darkgrid")

start_date = "2021-10-01"
num_days = 701
total_pop = 83278910.0
opacity = 0.15
lineWidth = 4
fontsize = 20
legendsize = 8
ticks = 12
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors_state_risk = sns.color_palette("tab20")

scenario = "amelag_single_1_april"
plot_dir = "/localdata1/pure/memilio/plots/" + scenario


def read_total_results_h5(path, group_key='Total'):
    f = h5py.File(path, 'r')
    group = f['0']
    total = group[group_key][()]
    f.close()
    return total


def compare_results(results_dir, compare_data, compartments, percentile='p50'):
    # Read data
    shift = False
    days_shift = 250

    percentiles = ["p50", "p75", "p95"]

    # Read simulated data for p25, p50, p75
    data_simulated_p25 = read_total_results_h5(
        os.path.join(results_dir, percentiles[0], "Results_sum.h5"))
    data_simulated_p25 = pd.DataFrame(data_simulated_p25)
    data_simulated_p25 = data_simulated_p25[compartments]
    data_simulated_p25["sum"] = data_simulated_p25.sum(axis=1)
    if shift:
        data_simulated_p25 = data_simulated_p25.iloc[days_shift:].reset_index(
            drop=True)
        # data_simulated_p25 = data_simulated_p25.drop(
        #     data_simulated_p25.index[100:160]).reset_index(drop=True)
        # data_simulated_p25 = data_simulated_p25.rolling(
        #     window=15, min_periods=1).mean()

    data_simulated_p50 = read_total_results_h5(
        os.path.join(results_dir, percentiles[1], "Results_sum.h5"))
    data_simulated_p50 = pd.DataFrame(data_simulated_p50)
    data_simulated_p50 = data_simulated_p50[compartments]
    data_simulated_p50["sum"] = data_simulated_p50.sum(axis=1)
    if shift:
        data_simulated_p50 = data_simulated_p50.iloc[days_shift:].reset_index(
            drop=True)
        # data_simulated_p50 = data_simulated_p50.drop(
        #     data_simulated_p50.index[120:160]).reset_index(drop=True)
        # data_simulated_p50 = data_simulated_p50.rolling(
        #     window=15, min_periods=1).mean()

    data_simulated_p75 = read_total_results_h5(
        os.path.join(results_dir, percentiles[2], "Results_sum.h5"))
    data_simulated_p75 = pd.DataFrame(data_simulated_p75)
    data_simulated_p75 = data_simulated_p75[compartments]
    data_simulated_p75["sum"] = data_simulated_p75.sum(axis=1)
    if shift:
        data_simulated_p75 = data_simulated_p75.iloc[days_shift:].reset_index(
            drop=True)
        # data_simulated_p75 = data_simulated_p75.drop(
        #     data_simulated_p75.index[120:160]).reset_index(drop=True)
        # data_simulated_p75 = data_simulated_p75.rolling(
        #     window=15, min_periods=1).mean()

    # read vvs data
    path_vvs = "/localdata1/pure/memilio/results/amelag_single_vvs/fac_variant_3.500000_red_fact_exposed_0.400000_t_immunity_1.000000_t_waning_immunity_365.000000/"
    data_simulated_p25_vvs = read_total_results_h5(
        os.path.join(path_vvs, percentiles[0], "Results_sum.h5"))
    data_simulated_p25_vvs = pd.DataFrame(data_simulated_p25_vvs)
    data_simulated_p25_vvs = data_simulated_p25_vvs[compartments]
    data_simulated_p25_vvs["sum"] = data_simulated_p25_vvs.sum(axis=1)
    if shift:
        data_simulated_p25_vvs = data_simulated_p25_vvs.iloc[days_shift:].reset_index(
            drop=True)

    data_simulated_p50_vvs = read_total_results_h5(
        os.path.join(path_vvs, percentiles[1], "Results_sum.h5"))
    data_simulated_p50_vvs = pd.DataFrame(data_simulated_p50_vvs)
    data_simulated_p50_vvs = data_simulated_p50_vvs[compartments]
    data_simulated_p50_vvs["sum"] = data_simulated_p50_vvs.sum(axis=1)
    if shift:
        data_simulated_p50_vvs = data_simulated_p50_vvs.iloc[days_shift:].reset_index(
            drop=True)

    data_simulated_p75_vvs = read_total_results_h5(
        os.path.join(path_vvs, percentiles[2], "Results_sum.h5"))
    data_simulated_p75_vvs = pd.DataFrame(data_simulated_p75_vvs)
    data_simulated_p75_vvs = data_simulated_p75_vvs[compartments]
    data_simulated_p75_vvs["sum"] = data_simulated_p75_vvs.sum(axis=1)
    if shift:
        data_simulated_p75_vvs = data_simulated_p75_vvs.iloc[days_shift:].reset_index(
            drop=True)

    # Read comparison data
    data_compare = pd.read_csv(compare_data, sep="\t")
    # delete all rows with datum smaller start_date
    data_compare = data_compare[data_compare['datum'] >= start_date].reset_index(
        drop=True)
    dates = data_compare[['datum']]
    data_compare = data_compare[['loess_vorhersage']]

    end_date = '2022-06-01'  # Enddatum im Format YYYY-MM-DD
    days_to_st = pd.date_range(start=start_date, end=end_date)

    days_to_st_series = pd.Series(days_to_st, name='datum')
    days_to_st_series = pd.to_datetime(days_to_st_series)
    days_to_st_series = days_to_st_series.dt.to_pydatetime()

    extendes_dates = copy.deepcopy(dates['datum'].values)
    extendes_dates = np.append(days_to_st_series, extendes_dates)

    extendes_dates = pd.DataFrame(extendes_dates, columns=['datum'])

    # format dates['datum'] to format YYYY-MM-DD
    extendes_dates['datum'] = pd.to_datetime(extendes_dates['datum'])
    dates['datum'] = pd.to_datetime(dates['datum'])

    # read RKI cases
    path_rki = "/localdata1/pure/memilio/data/pydata/Germany/cases_all_germany_ma7.json"
    df_rki = pd.read_json(path_rki)
    df_rki = df_rki.loc[(df_rki['Date'] >= start_date)]
    df_rki = df_rki.reset_index(drop=True)
    cases_rki = df_rki['Confirmed'].values
    # data is accumulated, so we need to calculate the daily cases
    cases_rki = np.diff(cases_rki, prepend=0)
    cases_rki = cases_rki[1:]

    # Truncate data to match length
    days_between_start_and_end = (pd.to_datetime(
        end_date) - pd.to_datetime(start_date)).days
    num_days = min(len(data_compare), len(
        data_simulated_p25) - days_between_start_and_end)
    data_compare = data_compare.iloc[:num_days]
    data_simulated_p25 = data_simulated_p25[:len(extendes_dates)]
    data_simulated_p50 = data_simulated_p50[:len(extendes_dates)]
    data_simulated_p75 = data_simulated_p75[:len(extendes_dates)]
    data_simulated_p25_vvs = data_simulated_p25_vvs[:len(extendes_dates)]
    data_simulated_p50_vvs = data_simulated_p50_vvs[:len(extendes_dates)]
    data_simulated_p75_vvs = data_simulated_p75_vvs[:len(extendes_dates)]
    dates = dates[:num_days]
    extendes_dates = extendes_dates[:len(data_simulated_p25)]
    cases_rki = cases_rki[:len(extendes_dates)]

    # Adjust comparison data scale
    data_compare['loess_vorhersage'] = data_compare['loess_vorhersage'] / \
        data_compare['loess_vorhersage'].iloc[0] * \
        data_simulated_p50['sum'].iloc[0]

    # Compute and print error
    error = np.linalg.norm(
        data_compare['loess_vorhersage'] - data_simulated_p50['sum'])
    print(f"Error: {error:.2e}")

    # Convert date column to datetime format
    dates['datum'] = pd.to_datetime(dates['datum'])

    # Plotting
    fig, ax1 = plt.subplots(figsize=(20, 10))

    # Plot simulation results
    ax1.set_xlabel('Date', fontsize=fontsize)
    ax1.set_ylabel('Cases', color='black', fontsize=fontsize)
    line2, = ax1.plot(extendes_dates['datum'], data_simulated_p25["sum"],
                      color=colors[0], linestyle='--', linewidth=lineWidth)
    line3, = ax1.plot(extendes_dates['datum'], data_simulated_p50["sum"],
                      color=colors[0], label='Simulated symptomatic cases (with waning)', linewidth=lineWidth)
    line4, = ax1.plot(extendes_dates['datum'], data_simulated_p75["sum"],
                      color=colors[0], linestyle='--', linewidth=lineWidth)
    ax1.fill_between(extendes_dates['datum'], data_simulated_p25["sum"],
                     data_simulated_p75["sum"], color=colors[0], alpha=opacity)

    line2_vvs, = ax1.plot(extendes_dates['datum'], data_simulated_p25_vvs["sum"],
                          color=colors[1], linestyle='--', linewidth=lineWidth)
    line3_vvs, = ax1.plot(extendes_dates['datum'], data_simulated_p50_vvs["sum"],
                          color=colors[1], label='Simulated symptomatic cases (without waning)', linewidth=lineWidth)
    line4_vvs, = ax1.plot(extendes_dates['datum'], data_simulated_p75_vvs["sum"],
                          color=colors[1], linestyle='--', linewidth=lineWidth)
    ax1.fill_between(extendes_dates['datum'], data_simulated_p25_vvs["sum"],
                     data_simulated_p75_vvs["sum"], color=colors[1], alpha=opacity)

    ax1.tick_params(axis='y', labelcolor=colors[0], labelsize=ticks)
    # ax1.set_ylim([1, 6000000])  # Set limits to avoid log scale errors
    # ax1.set_yscale('log')  # Set y-axis to log scale
    ax1.set_yscale('symlog', linthresh=1)
    ax1.set_ylim([0, 1e8])

    # Plot RKI data on ax1
    dates['datum'] = pd.to_datetime(dates['datum'])
    line_rki, = ax1.plot(extendes_dates['datum'], cases_rki, label='Reported Cases',
                         color=colors[2], linewidth=lineWidth)

    # Plot Amelag data on the second y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('Predicted Viral Load (LOESS Regression)',
                   color='black', fontsize=fontsize)
    line1, = ax2.plot(dates['datum'], data_compare['loess_vorhersage'],
                      label="AMELAG Data", color=colors[3], linewidth=lineWidth)
    ax2.tick_params(axis='y', labelcolor=colors[3], labelsize=ticks)
    ax2.set_yscale('log')  # Set second y-axis to log scale
    ax2.set_ylim([10000, 5e6])

    # Set date format for the x-axis
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%B %Y'))
    ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=6))

    # x ticks size
    # plt.xticks(fontsize=25)

    # Add grid
    ax1.grid(True)

    # Set font size of x-axis labels
    plt.setp(ax1.get_xticklabels(), fontsize=20)

    # Create a global legend
    lines = [line1, line3, line3_vvs, line_rki]
    labels = [line.get_label() for line in lines]
    fig.legend(lines, labels, loc='lower right',  bbox_to_anchor=(0.95, 0.1),
               fontsize=fontsize, frameon=True)

    fig.tight_layout()
    plt.savefig(os.path.join(plot_dir, "compare_amelag.png"))
    plt.clf()

    return 0


def plot_icu(results_dir):
    quartils = ["05", "25", "50", "75", "95"]

    icu_sim_m_false = []
    end_date = '2024-07-24'

    for quartil in quartils:
        # masks false
        data = read_total_results_h5(
            results_dir + "/p" + quartil + "/Results_sum.h5")
        data = data[:, [20, 21, 22]]
        data_summed = np.sum(data, axis=1)
        icu_sim_m_false.append(data_summed)
        # [20, 21, 22]))

    path_json_icu = "data//pydata//Germany//germany_divi_ma7.json"
    df = pd.read_json(path_json_icu)
    dates = pd.date_range(start=start_date, end=end_date)

    filtered_df = df.loc[(df['Date'] >= start_date) &
                         (df['Date'] <= end_date)]

    max_length = min(len(icu_sim_m_false[0]), len(filtered_df['Date']))
    icu_sim_m_false = [icu_sim_m_false[i][:max_length]
                       for i in range(len(icu_sim_m_false))]

    plt.plot(dates, icu_sim_m_false[1],
             color=colors[0], linewidth=lineWidth, linestyle='--')
    plt.plot(dates, icu_sim_m_false[2],
             color=colors[0], label='No Masks', linewidth=lineWidth, linestyle='-')
    plt.plot(dates, icu_sim_m_false[3],
             color=colors[0], linewidth=lineWidth, linestyle='--')
    plt.fill_between(
        dates, icu_sim_m_false[1],
        icu_sim_m_false[3],
        alpha=0.2, color=colors[0])

    plt.plot(filtered_df['Date'], filtered_df['ICU'],
             label='ICU people reported', color=colors[1], linewidth=lineWidth)

    dates = pd.date_range(start=start_date, end=end_date)[
        : len(icu_sim_m_false)]
    weeks = pd.date_range(start=start_date, end=end_date, freq='90D')
    plt.xticks(weeks, fontsize=15, rotation=45)

    plt.xlabel('Date', fontsize=18)
    plt.ylabel('ICU Occupancy', fontsize=18)
    # plt.yscale('log')
    plt.xticks(rotation=45)
    plt.grid(True)
    # plt.yscale('log')
    plt.tight_layout()
    plt.legend(fontsize=11)
    plt.savefig(plot_dir + "/icu_germany.png")
    plt.clf()

    return 0


def check_einzelstandorte(path_file):
    data = pd.read_csv(path_file, sep="\t")

    # Get unique locations and convert dates
    locations = data['standort'].unique()
    data['datum'] = pd.to_datetime(data['datum'])

    # Split locations into chunks of 20
    location_chunks = [locations[i:i + 20]
                       for i in range(0, len(locations), 20)]

    # Determine the number of subplots needed
    num_chunks = len(location_chunks)
    rows = (num_chunks + 1) // 2  # Arrange plots in two columns

    # Create the figure with subplots
    fig, axes = plt.subplots(rows, 2, figsize=(20, 5 * rows))
    axes = axes.flatten()  # Flatten axes array for easy indexing

    for i, chunk in enumerate(location_chunks):
        ax = axes[i]

        # Plot data for each location in the current chunk
        for location in chunk:
            data_location = data.loc[data['standort'] == location]
            ax.plot(data_location['datum'],
                    data_location['loess_vorhersage'], label=location)

        # Label axes and format dates
        ax.set_xlabel('Date')
        ax.set_ylabel('Amelag Data')
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %Y'))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
        ax.grid(True)
        ax.set_yscale('log')  # Set log scale for y-axis
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8)

    plt.tight_layout()
    fig.suptitle('Amelag Data by Location', fontsize=16, y=1.02)

    # Save the figure as one image
    plot_filename = os.path.join(
        plot_dir, "einzelstandorte", "compare_amelag_einzelstandorte.png")
    plt.savefig(plot_filename, bbox_inches='tight')

    return data


def check_einzelstandorte_all_in_one(path_file):
    data = pd.read_csv(path_file, sep="\t")
    data['datum'] = pd.to_datetime(data['datum'])
    locations = data['standort'].unique()
    dates = data['datum'].unique()

    for location in locations:
        data_location = data.loc[data['standort'] == location]
        plt.plot(data_location['datum'],
                 data_location['loess_vorhersage'], label=location)

    # label x and y axis
    plt.xlabel('Date', fontsize=fontsize)
    plt.ylabel('Amelag Data', fontsize=fontsize)

    # Set date format for the x-axis
    ax = plt.gca()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %Y'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))

    # size ticks
    plt.xticks(fontsize=ticks)

    # Add grid
    plt.grid(True)

    # Set font size of x-axis labels
    plt.setp(ax.get_xticklabels(), fontsize=10)

    # log scale
    plt.yscale('log')

    plt.savefig(os.path.join(
        plot_dir, "einzelstandorte/compare_amelag_einzelstandorte.png"))

    plt.clf()

    return data


def einzelstandorte_per_state(path_file):
    # Define the mapping of State IDs to State Names
    state_dict = {
        1: 'Schleswig-Holstein',
        2: 'Hamburg',
        3: 'Lower Saxony',
        4: 'Bremen',
        5: 'North Rhine-Westphalia',
        6: 'Hesse',
        7: 'Rhineland-Palatinate',
        8: 'Baden-Württemberg',
        9: 'Bavaria',
        10: 'Saarland',
        11: 'Berlin',
        12: 'Brandenburg',
        13: 'Mecklenburg-Western Pomerania',
        14: 'Saxony',
        15: 'Saxony-Anhalt',
        16: 'Thuringia',
    }

    abbreviation_to_id = {
        'SH': 1,  # Schleswig-Holstein
        'HH': 2,  # Hamburg
        'NI': 3,  # Lower Saxony
        'HB': 4,  # Bremen
        'NW': 5,  # North Rhine-Westphalia
        'HE': 6,  # Hesse
        'RP': 7,  # Rhineland-Palatinate
        'BW': 8,  # Baden-Württemberg
        'BY': 9,  # Bavaria
        'SL': 10,  # Saarland
        'BE': 11,  # Berlin
        'BB': 12,  # Brandenburg
        'MV': 13,  # Mecklenburg-Western Pomerania
        'SN': 14,  # Saxony
        'ST': 15,  # Saxony-Anhalt
        'TH': 16,  # Thuringia
    }

    # Load the main data file
    data = pd.read_csv(path_file, sep="\t")
    data['datum'] = pd.to_datetime(data['datum'])

    # lösche alle spalten aus datum , bundesland, und loess_vorhersage
    data = data[['datum', 'standort', 'bundesland', 'loess_vorhersage']]

    # Create a new column for the State ID
    data['State ID'] = data['bundesland'].map(abbreviation_to_id)

    # Split states into two groups
    state_groups = [list(state_dict.keys())[:8], list(state_dict.keys())[8:]]

    for idx, states in enumerate(state_groups):
        # Create the figure with subplots for the current group of states
        fig, axes = plt.subplots(4, 2, figsize=(20, 20))
        axes = axes.flatten()  # Flatten axes array for easy indexing

        for i, state in enumerate(states):
            ax = axes[i]  # Select the appropriate subplot
            state_data = data[data['State ID'] == state]
            locations = state_data['standort'].unique()

            # Plot data for each location within the current state
            for location in locations:
                data_location = state_data[state_data['standort'] == location]
                ax.plot(data_location['datum'],
                        data_location['loess_vorhersage'], label=location, linewidth=lineWidth)

            # Label axes and format dates
            ax.set_title(state_dict[state], fontsize=fontsize)
            ax.set_xlabel('Date', fontsize=fontsize)
            ax.set_ylabel('Predicted Viral Load (LOESS Regression)',
                          fontsize=10)
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %Y'))
            ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
            ax.grid(True)
            ax.set_yscale('log')  # Set log scale for y-axis
            ax.legend(loc='upper left', bbox_to_anchor=(
                1, 1), fontsize=legendsize)
            ax.tick_params(axis='both', which='major', labelsize=ticks)
            ax.set_xlim([dt.datetime(2022, 10, 1), dt.datetime(2024, 6, 1)])

        plt.tight_layout()
        fig.suptitle(
            f'AMELAG Data by Federal State (Part {idx + 1})', fontsize=fontsize, y=1.02)

        plot_filename = os.path.join(
            plot_dir, "einzelstandorte", f"compare_amelag_einzelstandorte_by_state_part_{idx + 1}.png")
        plt.savefig(plot_filename, bbox_inches='tight')
        plt.clf()

    # # Aggregate loess_vorhersage by state and date
    # agg_data = data.groupby(['State ID', 'datum']).agg(
    #     {'loess_vorhersage': 'sum'}).reset_index()

    # # Plot aggregated data the same way as above
    # fig, axes = plt.subplots(rows, 2, figsize=(20, 5 * rows))
    # axes = axes.flatten()  # Flatten axes array for easy indexing
    # for i, state in enumerate(state_dict):
    #     ax = axes[i]
    #     state_data = agg_data[agg_data['State ID'] == state]
    #     ax.plot(state_data['datum'],
    #             state_data['loess_vorhersage'], linewidth=lineWidth)
    #     ax.set_title(state_dict[state], fontsize=fontsize)
    #     ax.set_xlabel('Date', fontsize=fontsize)
    #     ax.set_ylabel('AMELAG Data', fontsize=fontsize)
    #     ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %Y'))
    #     ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
    #     ax.grid(True)
    #     ax.set_yscale('log')
    #     ax.tick_params(axis='both', which='major', labelsize=ticks)

    # plt.tight_layout()
    # fig.suptitle('Aggregated AMELAG Data by Federal State',
    #              fontsize=fontsize, y=1.02)
    # plot_filename = os.path.join(
    #     plot_dir, "einzelstandorte", "compare_amelag_einzelstandorte_aggregated_by_state.png")
    # plt.savefig(plot_filename, bbox_inches='tight')
    # plt.clf()


if __name__ == '__main__':
    # check if plot_dir exists, if not create it
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    path_cwd = os.getcwd()
    model_type = "amelag_single" + "_vvs"
    path_results = path_cwd + \
        "/results/amelag_single_1_april/fac_variant_3.500000_red_fact_exposed_0.400000_t_immunity_60.000000_t_waning_immunity_365.000000"
    path_compare = path_cwd + "/tools/amelag_aggregierte_kurve.tsv"
    path_compare_einzel = path_cwd + "/tools/amelag_einzelstandorte.tsv"
    icu_comp = [20, 21, 22]
    inf_symptomatic = [11, 12, 13, 14, 15, 16]
    # check_einzelstandorte(path_compare_einzel)
    # check_einzelstandorte_all_in_one(path_compare_einzel)
    einzelstandorte_per_state(path_compare_einzel)
    # compare_results(path_results, path_compare, inf_symptomatic)
    # plot_icu(path_results)
