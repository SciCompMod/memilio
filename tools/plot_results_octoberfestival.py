import json
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import matplotx

import memilio.epidata.getPopulationData as gpd
import memilio.plot.plotMap as pm
from memilio.epidata import geoModificationGermany as geoger
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import SymLogNorm, LinearSegmentedColormap

from scipy.ndimage import gaussian_filter

import seaborn as sns

# sns.set_style("darkgrid")


path = "/localdata1/code/memilio/DigiHero_SARSCoV2_Data.csv"

path_pop_json = "/localdata1/code/memilio/data/pydata/Germany/county_population.json"

scenario = "octoberfest1_1"
res_dir = "/localdata1/pure/memilio/results/" + scenario + "/secirts"


path_results_rki = "results_paper/Results_rki_sum.h5"
plot_dir = "/localdata1/pure/memilio/plots/" + scenario

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

start_date = '2022-10-10'
num_days = 83

opacity = 0.15
lineWidth = 3
fontsize = 28
tickssize = 12
figsize = (13, 10)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

start_date_obj = datetime.strptime(start_date, '%Y-%m-%d')
end_date_obj = start_date_obj + timedelta(days=num_days)
end_date = end_date_obj.strftime('%Y-%m-%d')


def read_results_h5(path, comp, filter_id,  group_key='Total'):
    with h5py.File(path, 'r') as f:
        keys = list(f.keys())
        res = []
        for i, key in enumerate(keys):
            # filter is a id (f.e. 9). Discard all keys not starting with the filter_id
            if not key.startswith(filter_id):
                continue
            group = f[key]
            total = group[group_key][()]
            if total.shape[1] > 1:
                comp_simulated = np.sum(total[:, comp], axis=1)
            else:
                comp_simulated = total[:, comp]
            res.append(comp_simulated)

    return res


# def read_results_h5(path, comp, group_key='Total'):
#     with h5py.File(path, 'r') as f:
#         keys = list(f.keys())
#         res = []
#         for i, key in enumerate(keys):
#             group = f[key]
#             total = group[group_key][()]
#             if total.shape[1] > 1:
#                 comp_simulated = np.sum(total[:, comp], axis=1)
#             else:
#                 comp_simulated = total[:, comp]
#             res.append(comp_simulated)

#     return res


def read_results_summed(res_dir, filename, masks, comp):

    f = h5py.File(res_dir + masks + filename, 'r')

    # Get the HDF5 group; key needs to be a group name from above
    group = f['0']

    # This assumes group[some_key_inside_the_group] is a dataset,
    # and returns a np.array:
    # time = group['Time'][()]
    total = group['Total'][()]
    # group1 = group['Group1'][()]

    comp_simulated = np.sum(total[:, comp], axis=1)

    # After you are done
    f.close()

    return comp_simulated[:num_days]


def read_results_quartil(path, comp, quartil):

    path_quartil = path + "//p" + str(quartil) + "//Results_sum.h5"
    f = h5py.File(path_quartil, 'r')

    # Get the HDF5 group; key needs to be a group name from above
    group = f['0']

    # This assumes group[some_key_inside_the_group] is a dataset,
    # and returns a np.array:
    # time = group['Time'][()]
    total = group['Total'][()]
    # group1 = group['Group1'][()]

    comp_simulated = np.sum(total[:, comp], axis=1)

    # After you are done
    f.close()

    return comp_simulated[:num_days]


def plot_immunity_levels(quartil="50"):

    naive = [0, 2, 5, 8, 11, 14, 17, 20, 24, 27]
    pi = [1, 3, 6, 9, 12, 15, 18, 21, 25]
    ii = [4, 7, 10, 13, 16, 19, 22, 23, 26, 28]

    modes = ["0", "1_ffp2"]
    for mode in modes:
        N_sim = read_results_summed(
            res_dir, "/mask_" + str(mode) + "/p" +
            quartil + "/Results_sum.h5", "",
            naive)
        PI_sim = read_results_summed(
            res_dir, "/mask_" + str(mode) + "/p" +
            quartil + "/Results_sum.h5", "",
            pi)
        II_sim = read_results_summed(
            res_dir, "/mask_" + str(mode) + "/p" +
            quartil + "/Results_sum.h5", "",
            ii)

        total = N_sim + PI_sim + II_sim
        N_sim = N_sim / total
        PI_sim = PI_sim / total
        II_sim = II_sim / total

        N_sim = N_sim[:num_days-1]
        PI_sim = PI_sim[:num_days-1]
        II_sim = II_sim[:num_days-1]

        # all days between start and end date
        dates = pd.date_range(start=start_date, end=end_date)[:num_days-1]

        width = 1.0
        plt.bar(dates, N_sim, width=width, label='Naive')
        plt.bar(dates, PI_sim, bottom=N_sim,
                width=width, label='Partial Immunity')
        plt.bar(dates, II_sim, bottom=N_sim + PI_sim,
                width=width, label='Improved Immunity')
        plt.xlabel('Date', fontsize=15)
        plt.ylabel('Percentage of population', fontsize=15)
        # plt.yscale('log')
        plt.title('Immunity levels distribution during simulation', fontsize=15)
        plt.xticks(rotation=45)
        plt.grid(True)
        plt.xlim(dates[0], dates[-1])
        plt.tight_layout()
        plt.legend(fontsize=10)
        plt.savefig(plot_dir + "/immunity_levels_masks_" + str(mode) + ".png")
        plt.clf()


def plot_susceptible(quartil="50"):
    # masks false
    SN_sim_m_false = read_results_summed(
        res_dir, "/masks_0/p" + quartil + "/Results_sum.h5", "", [0])
    SPI_sim_m_false = read_results_summed(
        res_dir, "/masks_0/p" + quartil + "/Results_sum.h5", "", [1])
    SII_sim_m_false = read_results_summed(
        res_dir, "/masks_0/p" + quartil + "/Results_sum.h5", "", [23])

    SN_sim_m_true = read_results_summed(
        res_dir, "/mask_1/p" + quartil + "/Results_sum.h5", "", [0])
    SPI_sim_m_true = read_results_summed(
        res_dir, "/mask_1/p" + quartil + "/Results_sum.h5", "", [1])
    SII_sim_m_true = read_results_summed(
        res_dir, "/mask_1/p" + quartil + "/Results_sum.h5", "", [23])

    # all days between start and end date
    dates = pd.date_range(start=start_date, end=end_date)[:num_days-1]

    # kürze alle auf num_days - 1
    SN_sim_m_false = SN_sim_m_false[:num_days-1]
    SPI_sim_m_false = SPI_sim_m_false[:num_days-1]
    SII_sim_m_false = SII_sim_m_false[:num_days-1]
    SN_sim_m_true = SN_sim_m_true[:num_days-1]
    SPI_sim_m_true = SPI_sim_m_true[:num_days-1]
    SII_sim_m_true = SII_sim_m_true[:num_days-1]

    # Plotten der Todesfälle über die Tage
    plt.plot(
        dates,
        SN_sim_m_false,
        label='Suscebtible Naive, no masks transport')
    plt.plot(
        dates,
        SPI_sim_m_false,
        label='Suscebtible PI, no masks transport')
    plt.plot(
        dates,
        SII_sim_m_false,
        label='Suscebtible II, no masks transport')
    plt.plot(
        dates,
        SN_sim_m_true,
        label='Suscebtible Naive, with masks transport')
    plt.plot(
        dates,
        SPI_sim_m_true,
        label='Suscebtible PI, with masks transport')
    plt.plot(
        dates,
        SII_sim_m_true,
        label='Suscebtible II, with masks transport')
    plt.xlabel('Date', fontsize=15)
    plt.ylabel('Number of individuals', fontsize=15)
    # plt.yscale('log')
    plt.title('Susceptibles', fontsize=15)
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(fontsize=10)
    plt.savefig(plot_dir + "/susceptibles.png")
    plt.clf()

    return 0


def plot_deaths(comulate=True):

    quartils = ["05", "25", "50", "75", "95"]

    death_sim_m_false = []
    death_sim_m_true = []
    death_sim_m_true_ffp2 = []

    for quartil in quartils:
        # masks false
        death_sim_m_false.append(read_results_summed(
            res_dir, "/masks_0/p" + quartil + "/Results_sum.h5", "",
            list(range(24, 27))))

        # masks true
        death_sim_m_true.append(read_results_summed(
            res_dir, "/mask_1/p" + quartil + "/Results_sum.h5", "",
            list(range(24, 27))))

        # ff2 masks
        death_sim_m_true_ffp2.append(read_results_summed(
            res_dir, "/mask_1_ffp2/p" + quartil + "/Results_sum.h5", "",
            list(range(24, 27))))

    # read rki excel data
    results_rki_extrapolated = read_results_summed(
        path_results_rki, "", "", [24])
    results_rki_extrapolated = np.diff(results_rki_extrapolated)

    # death_0108 = 143979
    # death_3009 = 149948
    # diff_death = (death_3009 - death_0108) / 60
    # results_rki_reported_weekly = np.full(60, diff_death)

    results_rki_extrapolated = results_rki_extrapolated[:num_days]

    # all days between start and end date
    dates = pd.date_range(start=start_date, end=end_date)[:num_days]

    dates_rki = pd.date_range(start=start_date, end=end_date)[
        :results_rki_extrapolated.shape[0]]

    for i in range(1, 4):
        if i == 2:
            linestyle = '-'
        else:
            linestyle = '--'

        plt.plot(dates, death_sim_m_false[i],
                 label='No Masks P' + quartils[i], linestyle=linestyle, color='red')
        plt.plot(dates, death_sim_m_true[i],
                 label='Mask Surgical P' + quartils[i], linestyle=linestyle, color='green')
        plt.plot(dates, death_sim_m_true_ffp2[i],
                 label='Mask  FFP2 P' + quartils[i], linestyle=linestyle, color='blue')

    # plt.fill_between(
    #     dates, death_sim_m_false[1],
    #     death_sim_m_false[3],
    #     alpha=0.2, color='lightblue')
    # plt.fill_between(
    #     dates, death_sim_m_true[1],
    #     death_sim_m_true[3],
    #     alpha=0.2, color='lightblue')
    # plt.fill_between(
    #     dates, death_sim_m_true_ffp2[1],
    #     death_sim_m_true_ffp2[3],
    #     alpha=0.2, color='lightblue')
    plt.plot(
        dates_rki,
        results_rki_extrapolated,
        label='reported Deaths extrapolated')

    plt.xlabel('Date', fontsize=15)
    plt.ylabel('Deaths', fontsize=15)
    # plt.yscale('log')
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.xlim(dates[0], dates[-1])
    plt.tight_layout()
    plt.legend(fontsize=10)
    if comulate:
        plt.savefig(plot_dir + "/cumulated_deaths.png")
    else:
        plt.savefig(plot_dir + "/daily_deaths.png")
    plt.clf()

    return 0


def plot_icu_global():
    quartils = ["05", "25", "50", "75", "95"]

    icu_sim_m_false = []

    for quartil in quartils:
        # masks false
        icu_sim_m_false.append(read_results_summed(
            res_dir, "/masks_0/p" + quartil + "/Results_sum.h5", "",
            [20, 21, 22]))

    path_json_icu = "data//pydata//Germany//germany_divi_ma7.json"
    df = pd.read_json(path_json_icu)
    dates = pd.date_range(start=start_date, end=end_date)[:num_days]

    filtered_df = df.loc[(df['Date'] >= start_date) &
                         (df['Date'] <= end_date)][:num_days]

    for i in range(1, 4):
        if i == 2:
            linestyle = '-'
        else:
            linestyle = '--'

        plt.plot(dates, icu_sim_m_false[i],
                 label='No Masks P' + quartils[i], linestyle=linestyle, color=colors[0])

    plt.fill_between(
        dates, icu_sim_m_false[1],
        icu_sim_m_false[3],
        alpha=0.2, color=colors[0])

    plt.plot(filtered_df['Date'], filtered_df['ICU'],
             label='ICU people reported', color=colors[1])
    plt.xlabel('Date', fontsize=18)
    plt.ylabel('People ICU', fontsize=18)
    # plt.yscale('log')
    plt.title('Number of people located in ICU', fontsize=20)
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.yscale('log')
    plt.tight_layout()
    plt.legend(fontsize=11)
    plt.savefig(plot_dir + "/icu.png")
    plt.clf()

    return 0


def plot_icu_bavaria():
    quartils = ["05", "25", "50", "75", "95"]

    icu_sim_m_false = []

    for quartil in quartils:
        # results for bayern
        icu_sim_m_false.append(read_results_h5(
            res_dir + "/masks_0/p" + quartil + "/Results.h5", [20, 21, 22], "9"))

    icu_sim_summed = np.sum(icu_sim_m_false, axis=1)

    # for i in range(len(icu_sim_summed)):
    #     # multipliziere alle außer den ersten wert mit 0.8
    #     icu_sim_summed[i][:] = icu_sim_summed[i][:] * 3.0

    path_json_icu = "data//pydata//Germany//county_divi.json"
    df = pd.read_json(path_json_icu)

    dates = pd.date_range(start=start_date, end=end_date)[:num_days]

    filtered_df = df[
        (df['ID_County'].astype(str).str.startswith('9')) &
        (df['Date'] >= start_date) &
        (df['Date'] <= end_date)
    ]
    daily_icu_sum = filtered_df.groupby('Date')['ICU'].sum().reset_index()

    # for i in range(1, 4):
    #     if i == 2:
    #         linestyle = '-'
    #     else:
    #         linestyle = '--'

    #     plt.plot(dates, icu_sim_summed[i][:num_days],
    #              label='No Masks P' + quartils[i], linestyle=linestyle, color=colors[0], linewidth=lineWidth)

    plt.plot(dates, icu_sim_summed[1][:num_days],
             linestyle='--', color=colors[0], linewidth=lineWidth)
    plt.plot(dates, icu_sim_summed[2][:num_days],
             color=colors[0], label='Simulated', linewidth=lineWidth, linestyle='-')
    plt.plot(dates, icu_sim_summed[3][:num_days],
             linestyle='--', color=colors[0], linewidth=lineWidth)

    plt.fill_between(
        dates, icu_sim_summed[1][:num_days],
        icu_sim_summed[3][:num_days],
        alpha=0.2, color=colors[0])

    plt.plot(daily_icu_sum['Date'], daily_icu_sum['ICU'],
             label='Reported', color=colors[1],  linewidth=lineWidth)

    plt.xlabel('Date', fontsize=18)
    plt.ylabel('ICU Occupancy', fontsize=18)
    # plt.title('Number of people located in ICU', fontsize=20)
    plt.xticks(rotation=45)

    # Set date format on x-axis
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%B %Y'))
    plt.gca().xaxis.set_major_locator(
        mdates.MonthLocator(interval=1))  # Show every month

    plt.xticks(fontsize=tickssize)
    # plt.ylim(0, 1000)

    plt.grid(True)
    plt.tight_layout()
    plt.legend(fontsize=15)
    plt.savefig(plot_dir + "/icu_bavaria.png")
    plt.clf()

    return 0


def plot_sheet():

    modes = ["0"]

    quartils = ["05", "25", "50", "75", "95"]
    for mode in modes:
        sim_exposed, sim_INS, sim_Isy, sim_hosp, sim_icu, sim_sus, sim_dead, sim_ti = [
        ], [], [], [], [], [], [], []
        path_quartils = res_dir + "/masks_" + mode
        for quartil in quartils:
            sim_exposed.append(read_results_quartil(
                path_quartils, [2, 3, 4], quartil))
            sim_INS.append(read_results_quartil(
                path_quartils, list(range(5, 8)), quartil))
            sim_Isy.append(read_results_quartil(
                path_quartils, list(range(11, 14)), quartil))
            sim_hosp.append(read_results_quartil(
                path_quartils, list(range(17, 20)), quartil))
            sim_icu.append(read_results_quartil(
                path_quartils, list(range(20, 23)), quartil))
            sim_sus.append(read_results_quartil(
                path_quartils, [0, 1, 23], quartil))
            sim_dead.append(read_results_quartil(
                path_quartils, list(range(24, 27)), quartil))
            sim_ti.append(read_results_quartil(
                path_quartils, [27, 28], quartil))

        fig, axs = plt.subplots(4, 2, figsize=(8.27, 11.69))
        for i in range(1, 4):
            if i == 2:
                linestyle = '-'
                color = "blue"
            else:
                linestyle = '--'
                color = "blue"
            axs[0, 0].plot(range(num_days), sim_exposed[i],
                           label="P" + quartils[i], linestyle=linestyle, color=color)
            axs[0, 1].plot(range(num_days), sim_INS[i],
                           label="P" + quartils[i], linestyle=linestyle, color=color)
            axs[1, 0].plot(range(num_days), sim_Isy[i],
                           label="P" + quartils[i], linestyle=linestyle, color=color)
            axs[1, 1].plot(range(num_days), sim_hosp[i],
                           label="P" + quartils[i], linestyle=linestyle, color=color)
            axs[2, 0].plot(range(num_days), sim_icu[i],
                           label="P" + quartils[i], linestyle=linestyle, color=color)
            axs[2, 1].plot(range(num_days), sim_sus[i],
                           label="P" + quartils[i], linestyle=linestyle, color=color)
            axs[3, 0].plot(range(num_days), sim_dead[i],
                           label="P" + quartils[i], linestyle=linestyle, color=color)
            axs[3, 1].plot(range(num_days), sim_ti[i],
                           label="P" + quartils[i], linestyle=linestyle, color=color)

        axs[0, 0].set_title("Exposed")
        axs[0, 1].set_title("Infected no symptoms")
        axs[1, 0].set_title("Infected symptomatic")
        axs[1, 1].set_title("Severe")
        axs[2, 0].set_title("ICU")
        axs[2, 1].set_title("Susceptible")
        axs[3, 0].set_title("Deaths")
        axs[3, 1].set_title("Temporary immunity")

        axs[0, 0].set_yscale('log')
        axs[0, 1].set_yscale('log')
        axs[1, 0].set_yscale('log')
        axs[1, 1].set_yscale('log')
        axs[3, 1].set_yscale('log')

        fillcolor = "blue"
        axs[0, 0].fill_between(range(num_days), sim_exposed[1],
                               sim_exposed[3], alpha=0.2, color=fillcolor)
        axs[0, 1].fill_between(range(num_days), sim_INS[1],
                               sim_INS[3], alpha=0.2, color=fillcolor)
        axs[1, 0].fill_between(range(num_days), sim_Isy[1],
                               sim_Isy[3], alpha=0.2, color=fillcolor)
        axs[1, 1].fill_between(range(num_days), sim_hosp[1],
                               sim_hosp[3], alpha=0.2, color=fillcolor)
        axs[2, 0].fill_between(range(num_days), sim_icu[1],
                               sim_icu[3], alpha=0.2, color=fillcolor)
        axs[2, 1].fill_between(range(num_days), sim_sus[1],
                               sim_sus[3], alpha=0.2, color=fillcolor)
        axs[3, 0].fill_between(range(num_days), sim_dead[1],
                               sim_dead[3], alpha=0.2, color=fillcolor)
        axs[3, 1].fill_between(range(num_days), sim_ti[1],
                               sim_ti[3], alpha=0.2, color=fillcolor)

        for ax in axs.flat:
            ax.set(xlabel="Tage", ylabel="Werte")
            ax.grid(True)
        plt.tight_layout()

        plt.savefig(plot_dir + "/combined_plots_mask" + mode + ".png")
        plt.clf()


def plot_total_infections(daily_transmissions=True):
    # load data digihero
    pop_germany = 84607000
    data_digihero = pd.read_csv(path)
    start_indx = data_digihero.index[data_digihero['Date']
                                     == start_date].tolist()

    end_indx = data_digihero.index[data_digihero['Date'] == end_date].tolist()
    data_digihero = data_digihero.iloc[start_indx[0]:end_indx[0]+1]
    cases_digihero = data_digihero['N_Inf'] * \
        pop_germany / data_digihero['N']

    #     flows S->E: 0,17,33
    # flows Ins->Isy: 3,5,20,22, 36,38
    flow_indx = [0, 17, 33]
    if not daily_transmissions:
        flow_indx = [3, 5, 20, 22, 36, 38]

    # dates = pd.date_range(
    #     start=data_digihero['Date'][0], end=data_digihero['Date'][1155])
    # plt.plot(dates, cases_digihero,
    #          label='DigiHero')
    # weeks = pd.date_range(
    #     start=data_digihero['Date'][0], end=data_digihero['Date'][1155], freq='30D')

    # plt.xticks(weeks, fontsize=15, rotation=90)
    # plt.show()

    # read case data rki
    path_rki_data = "/localdata1/code/memilio/data/pydata/Germany/cases_all_germany.json"
    df = pd.read_json(path_rki_data)
    df = df.loc[(df['Date'] >= start_date) & (df['Date'] <= end_date)]
    cases_rki = np.diff(df['Confirmed'])

    masks = ["0"]  # "0", "1", "1_ffp2"]
    inf_mobility_node_masks_0 = []
    inf_mobility_node_mask_1_ffp2 = []
    for mode in masks:
        dir_mobility = res_dir + "/masks_" + mode + "/mobility/flows/"
        data_p25 = read_results_quartil(dir_mobility, flow_indx, "25")
        data_p50 = read_results_quartil(dir_mobility, flow_indx, "50")
        data_p75 = read_results_quartil(dir_mobility, flow_indx, "75")

        dir_local = res_dir + "/masks_" + mode + "/flows/"
        data_p25 += read_results_quartil(dir_local, flow_indx, "25")
        data_p50 += read_results_quartil(dir_local, flow_indx, "50")
        data_p75 += read_results_quartil(dir_local, flow_indx, "75")

        if mode == "0":
            inf_mobility_node_masks_0.append(data_p25)
            inf_mobility_node_masks_0.append(data_p50)
            inf_mobility_node_masks_0.append(data_p75)

        else:
            inf_mobility_node_mask_1_ffp2.append(data_p25)
            inf_mobility_node_mask_1_ffp2.append(data_p50)
            inf_mobility_node_mask_1_ffp2.append(data_p75)

    # values are cumulated. Calculate the difference
    inf_mobility_node_masks_0 = np.array(inf_mobility_node_masks_0)
    inf_mobility_node_masks_0 = np.diff(inf_mobility_node_masks_0, axis=1)

    inf_mobility_node_mask_1_ffp2 = np.array(inf_mobility_node_mask_1_ffp2)
    inf_mobility_node_mask_1_ffp2 = np.diff(
        inf_mobility_node_mask_1_ffp2, axis=1)

    inf_mobility_node_masks_0 = np.array(inf_mobility_node_masks_0)
    inf_mobility_node_masks_0_quartils = np.percentile(
        inf_mobility_node_masks_0, [5, 25, 50, 75, 95], axis=0)

    # inf_mobility_node_mask_1 = np.array(inf_mobility_node_mask_1)
    # inf_mobility_node_mask_1_quartils = np.percentile(
    #     inf_mobility_node_mask_1, [5, 25, 50, 75, 95], axis=0)

    inf_mobility_node_mask_1_ffp2 = np.array(inf_mobility_node_mask_1_ffp2)
    inf_mobility_node_mask_1_ffp2_quartils = np.percentile(
        inf_mobility_node_mask_1_ffp2, [5, 25, 50, 75, 95], axis=0)

    # all days between start and end date
    dates = pd.date_range(start=start_date, end=end_date)[
        :inf_mobility_node_mask_1_ffp2_quartils.shape[1]]
    weeks = pd.date_range(start=start_date, end=end_date, freq='7D')

    # no masks
    plt.plot(
        dates, inf_mobility_node_masks_0_quartils[1], color=colors[0], linestyle='--', linewidth=lineWidth)
    plt.plot(
        dates, inf_mobility_node_masks_0_quartils[2], color=colors[0], label='No face masks', linewidth=lineWidth)
    plt.plot(
        dates, inf_mobility_node_masks_0_quartils[3], color=colors[0], linestyle='--', linewidth=lineWidth)
    plt.fill_between(
        dates, inf_mobility_node_masks_0_quartils[1], inf_mobility_node_masks_0_quartils[3], color=colors[0], alpha=opacity)

    # ffp2 masks
    plt.plot(
        dates, inf_mobility_node_mask_1_ffp2_quartils[1], color=colors[2], linestyle='--', linewidth=lineWidth)
    plt.plot(
        dates, inf_mobility_node_mask_1_ffp2_quartils[2], color=colors[2], label='FFP2 masks', linewidth=lineWidth)
    plt.plot(
        dates, inf_mobility_node_mask_1_ffp2_quartils[3], color=colors[2], linestyle='--', linewidth=lineWidth)
    plt.fill_between(
        dates, inf_mobility_node_mask_1_ffp2_quartils[1], inf_mobility_node_mask_1_ffp2_quartils[3], color=colors[2], alpha=opacity, linewidth=lineWidth)

    # official RKI Reporting
    plt.plot(dates, cases_rki[:inf_mobility_node_mask_1_ffp2_quartils.shape[1]],
             label='Reported Infections RKI', color=colors[1], linewidth=lineWidth)

    # digihero data
    # plt.plot(dates, cases_digihero[:inf_mobility_node_mask_1_ffp2_quartils.shape[1]],
    #          label='Extrapolated DigiHero Reports ', color=colors[3], linewidth=lineWidth)
    plt.xlabel('Date', fontsize=18)
    if daily_transmissions:
        plt.ylabel('Total Infections', fontsize=18)
    else:
        plt.ylabel('Daily Symptomatic', fontsize=18)
    # plt.yscale('log')
    plt.xticks(weeks, fontsize=tickssize, rotation=45)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(fontsize=11)
    if daily_transmissions:
        plt.savefig(plot_dir + "/total_infections.png")
    else:
        plt.savefig(plot_dir + "/total_symptoms.png")

    plt.clf()


def plot_total_infections_bavaria(daily_transmissions=True):
    # load data digihero
    pop_germany = 13080000
    data_digihero = pd.read_csv(path)
    start_indx = data_digihero.index[data_digihero['Date']
                                     == start_date].tolist()
    end_indx = data_digihero.index[data_digihero['Date'] == end_date].tolist()
    # data_digihero = data_digihero.iloc[start_indx[0]:end_indx[0]+1]
    # cases_digihero = data_digihero['N_Inf'] * pop_germany / data_digihero['N']
    # cases_digihero = cases_digihero.rolling(window=7, min_periods=1).mean()

    # Define flow indices based on daily_transmissions flag
    flow_indx = [0, 17, 33] if daily_transmissions else [3, 5, 20, 22, 36, 38]

    # read case data rki
    path_rki_data = "/localdata1/pure/memilio/data/pydata/Germany/cases_all_county_ma7.json"
    df = pd.read_json(path_rki_data)
    # delete all entries in 'ID_County' that are not 9xxxxx (Bavaria)
    df = df[df['ID_County'].astype(str).str.startswith('9')]
    # sum up the cases for all counties in Bavaria
    df = df.groupby('Date').sum().reset_index()
    df = df.loc[(df['Date'] >= start_date) & (df['Date'] <= end_date)]
    cases_rki = np.diff(df['Confirmed'])

    # path_rki_data_ma7 = '/localdata1/code/memilio/data/pydata/Germany/cases_all_germany_ma7.json'
    # df_ma7 = pd.read_json(path_rki_data_ma7)
    # df_ma7 = df_ma7.loc[(df_ma7['Date'] >= start_date)
    #                     & (df_ma7['Date'] <= end_date)]
    # cases_rki_ma7 = np.diff(df_ma7['Confirmed'])

    masks = ["0"]
    inf_mobility_node_masks_0 = []
    for mode in masks:
        dir_mobility = res_dir + "/masks_" + mode + "/mobility/flows/"
        data_p25 = read_results_h5(
            dir_mobility + "/p25/Results.h5", flow_indx, "9")
        data_p50 = read_results_h5(
            dir_mobility + "/p50/Results.h5", flow_indx, "9")
        data_p75 = read_results_h5(
            dir_mobility + "/p75/Results.h5", flow_indx, "9")

        dir_local = res_dir + "/masks_" + mode + "/flows/"
        data_p25 += read_results_h5(dir_local +
                                    "/p25/Results.h5", flow_indx, "9")
        data_p50 += read_results_h5(dir_local +
                                    "/p50/Results.h5", flow_indx, "9")
        data_p75 += read_results_h5(dir_local +
                                    "/p75/Results.h5", flow_indx, "9")

        if mode == "0":
            inf_mobility_node_masks_0.extend([data_p25, data_p50, data_p75])

    # Sum and differentiate the simulation data
    icu_sim_summed = np.sum(inf_mobility_node_masks_0, axis=1)
    inf_mobility_node_masks_0 = np.array(icu_sim_summed)
    inf_mobility_node_masks_0 = np.diff(inf_mobility_node_masks_0, axis=1)

    # cases_rki[0] should be the same as the first value in the simulation data
    cases_rki[0] = inf_mobility_node_masks_0[1][0]

    # Generate dates and weeks for plotting
    dates = pd.date_range(start=start_date, end=end_date)[
        :inf_mobility_node_masks_0.shape[1]]

    # Plotting data
    plt.plot(dates, inf_mobility_node_masks_0[0],
             color=colors[0], linestyle='--', linewidth=lineWidth)
    plt.plot(dates, inf_mobility_node_masks_0[1], color=colors[0],
             label='Simulated symptomatic Cases', linewidth=lineWidth)
    plt.plot(dates, inf_mobility_node_masks_0[2],
             color=colors[0], linestyle='--', linewidth=lineWidth)
    plt.fill_between(
        dates, inf_mobility_node_masks_0[0], inf_mobility_node_masks_0[2], color=colors[0], alpha=opacity)

    # Official RKI Reporting
    plt.plot(dates, cases_rki[:inf_mobility_node_masks_0.shape[1]],
             label='Reported Cases', color=colors[1], linewidth=lineWidth)
    # plotte ma7 daten mit gestrichelter linie
    # plt.plot(dates, cases_rki_ma7[:inf_mobility_node_masks_0.shape[1]],
    #          label='Reported Infections RKI MA7', color=colors[1], linestyle='--', linewidth=lineWidth)

    # DigiHero data
    # plt.plot(dates, cases_digihero[:inf_mobility_node_masks_0.shape[1]],
    #          label='Extrapolated DigiHero Reports ', color=colors[3], linewidth=lineWidth)

    plt.xlabel('Date', fontsize=18)
    if daily_transmissions:
        plt.ylabel('Total Infections', fontsize=18)
    else:
        plt.ylabel('Cases', fontsize=18)

    # Format x-axis to display 'Month Year' and increase tick font size
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%B %Y'))
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    plt.xticks(fontsize=tickssize, rotation=45)

    plt.grid(True)
    plt.tight_layout()
    plt.legend(fontsize=15)
    # log scale for y-axis
    plt.yscale('log')
    if daily_transmissions:
        plt.savefig(plot_dir + "/total_infections_bavaria.png")
    else:
        plt.savefig(plot_dir + "/total_symptoms_bavaria.png")

    plt.clf()


def map_plot_bavaria(res_dir, relative=False, rki_data=False):

    # timms = [27, 28]
    # deads = [24, 25, 26]
    # susceptible = [0, 1, 23]
    # infected = [5, 6, 7, 8, 9, 10, 11,
    #             12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    symptomatic = [11, 12, 13, 14, 15, 16]  # 2, 3, 4,
    # icu = [20, 21, 22]
    # exposed = [2, 3, 4]
    files_input = {
        'Data set 1': res_dir}
    file_format = 'h5'
    # Define age groups which will be considered through filtering
    # Keep keys and values as well as its assignment constant, remove entries
    # if only part of the population should be plotted or considered, e.g., by
    # setting:
    # age_groups = {1: '5-14', 2: '15-34'}
    age_groups = {0: '0-4', 1: '5-14', 2: '15-34',
                  3: '35-59', 4: '60-79', 5: '80+'}
    if len(age_groups) == 6:
        filter_age = None
    else:
        if file_format == 'json':
            filter_age = [val for val in age_groups.values()]
        else:
            filter_age = ['Group' + str(key) for key in age_groups.keys()]

    # results_bavaria = read_results_h5(
    #     res_dir + ".h5", symptomatic, "9")

    # results_total =  read_results_h5(
    #     res_dir + ".h5", symptomatic)
    # get the max and min value global in results_bavaria
    min_val = 0.001  # np.min(np.concatenate(results_bavaria))
    max_val = 0.005  # np.max(np.concatenate(results_bavaria))

    calc_min_val = 0.5
    calc_max_val = 0

    norm = SymLogNorm(linthresh=1, linscale=0.7,
                      vmin=min_val, vmax=max_val)

    # Create a ScalarMappable with the same Norm and Colormap as your data
    colors = ["green", "yellow", "red", "purple"]
    cmap = LinearSegmentedColormap.from_list("my_colormap", colors)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Create a new figure for the colorbar
    cbar_fig, ax = plt.subplots(figsize=(8, 1))
    cbar = plt.colorbar(sm, orientation='horizontal', cax=ax)

    # Set ticks every 0.01
    # cbar.set_ticks([0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012])
    # cbar.set_ticks([0, 0.0005, 0.001, 0.0015, 0.002])

    # Format ticks as standard decimal notation
    cbar.ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    # Optional: Adjust tick label size
    cbar.ax.tick_params(labelsize=10)
    plt.tight_layout()

    # Save the colorbar
    fn_colorbar = os.path.join(plot_dir, 'colorbar.png')
    if rki_data:
        fn_colorbar = os.path.join(plot_dir, 'colorbar_rki.png')
    plt.savefig(fn_colorbar, dpi=900)
    plt.clf()

    for day in range(20, 61, 20):

        i = 0
        for file in files_input.values():
            # MEmilio backend hdf5 example
            df = pm.extract_data(
                file, region_spec=None, column=None, date=day,
                filters={'Group': filter_age, 'InfectionState': symptomatic},
                file_format=file_format)
            if rki_data:
                path_rki_data = "/localdata1/pure/memilio/data/pydata/Germany/cases_all_county_ma7.json"
                df = pd.read_json(path_rki_data)
                # delete all entries in 'ID_County' that are not 9xxxxx (Bavaria)
                # df = df[df['ID_County'].astype(str).str.startswith('9')]
                # calculate start_date + day
                start_date_dt = datetime.strptime(start_date, '%Y-%m-%d')
                date_cur = start_date_dt + timedelta(days=day)
                date_cur = date_cur.strftime('%Y-%m-%d')

                date_before_curr = start_date_dt + timedelta(days=day - 1)
                date_before_curr = date_before_curr.strftime('%Y-%m-%d')
                df_day_before = df.loc[(df['Date'] == date_before_curr)].reset_index(
                    drop=True)
                df = df.loc[(df['Date'] == date_cur)].reset_index(drop=True)
                # df data is accumulated. Calculate the difference using df_day_before
                df['Confirmed'] = df['Confirmed'] - df_day_before['Confirmed']
                # delete all columns except 'ID_County' and 'Confirmed'
                df = df[['ID_County', 'Confirmed']]
                # rename ID_County to Region
                df = df.rename(columns={'ID_County': 'Region'})
                # rename Confirmed to Count
                df = df.rename(columns={'Confirmed': 'Count'})

            # all entries to numeric
            df = df.apply(pd.to_numeric, errors='coerce')

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

                # overwrite column Count with Count (rel)
                df['Count'] = df['Count (rel)']
                df = df.drop(columns=['Count (rel)'])

        # all entries to numeric
        df = df.apply(pd.to_numeric, errors='coerce')
        dfs_all = df

        filename = 'day_' + str(day)
        if rki_data:
            filename += '_rki'

        dfs_all = dfs_all.apply(pd.to_numeric, errors='coerce')

        dfs_all_sorted = dfs_all.sort_values(by='Region')
        dfs_all_sorted = dfs_all_sorted.reset_index(drop=True)

        # delete entires in there Region does not start with 9 (Bavaria)
        dfs_all_sorted = dfs_all_sorted[dfs_all_sorted['Region'].astype(
            str).str.startswith('9')]

        # suche nach dem index von münchen und multipliziere Count mit 1.5
        id_munich = 9162
        dfs_all_sorted.loc[dfs_all_sorted['Region']
                           == id_munich, 'Count'] *= 1.2

        if (calc_max_val < dfs_all_sorted['Count'].max()):
            calc_max_val = dfs_all_sorted['Count'].max()

        if (calc_min_val > dfs_all_sorted['Count'].min()):
            calc_min_val = dfs_all_sorted['Count'].min()

        print("day: " + str(day))

        pm.plot_map(norm,
                    dfs_all_sorted, scale_colors=[min_val, max_val],
                    legend=['', ''],
                    title='',
                    plot_colorbar=False,
                    output_path=plot_dir,
                    fig_name=filename, dpi=900,
                    outercolor='white')
        plt.clf()
    print("min_val: " + str(calc_min_val))
    print("max_val: " + str(calc_max_val))


def map_plot_bavaria_diff_t0(res_dir, relative=False, rki_data=False):

    # timms = [27, 28]
    # deads = [24, 25, 26]
    # susceptible = [0, 1, 23]
    # infected = [5, 6, 7, 8, 9, 10, 11,
    #             12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    symptomatic = [11, 12, 13, 14, 15, 16]  # 2, 3, 4,
    # icu = [20, 21, 22]
    # exposed = [2, 3, 4]
    files_input = {
        'Data set 1': res_dir}
    file_format = 'h5'
    # Define age groups which will be considered through filtering
    # Keep keys and values as well as its assignment constant, remove entries
    # if only part of the population should be plotted or considered, e.g., by
    # setting:
    # age_groups = {1: '5-14', 2: '15-34'}
    age_groups = {0: '0-4', 1: '5-14', 2: '15-34',
                  3: '35-59', 4: '60-79', 5: '80+'}
    if len(age_groups) == 6:
        filter_age = None
    else:
        if file_format == 'json':
            filter_age = [val for val in age_groups.values()]
        else:
            filter_age = ['Group' + str(key) for key in age_groups.keys()]

    # results_bavaria = read_results_h5(
    #     res_dir + ".h5", symptomatic, "9")

    # results_total =  read_results_h5(
    #     res_dir + ".h5", symptomatic)
    # get the max and min value global in results_bavaria
    min_val = 0  # np.min(np.concatenate(results_bavaria))
    max_val = 1  # np.max(np.concatenate(results_bavaria))

    calc_max_val = 0
    calc_min_val = 100000

    norm = SymLogNorm(linthresh=1, linscale=0.7,
                      vmin=min_val, vmax=max_val)

    # Create a ScalarMappable with the same Norm and Colormap as your data
    colors = ["green", "yellow", "red", "purple"]
    cmap = LinearSegmentedColormap.from_list("my_colormap", colors)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Create a new figure for the colorbar
    cbar_fig, ax = plt.subplots(figsize=(8, 1))
    cbar = plt.colorbar(sm, orientation='horizontal', cax=ax)

    # Set ticks every 0.01
    # cbar.set_ticks([0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012])
    # cbar.set_ticks([0, 0.0005, 0.001, 0.0015, 0.002])

    # Format ticks as standard decimal notation
    cbar.ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    # Optional: Adjust tick label size
    cbar.ax.tick_params(labelsize=10)
    plt.tight_layout()

    # Save the colorbar
    fn_colorbar = os.path.join(plot_dir, 'diff_colorbar.png')
    if rki_data:
        fn_colorbar = os.path.join(plot_dir, 'diff_colorbar_rki.png')
    plt.savefig(fn_colorbar, dpi=900)
    plt.clf()

    data_t0 = []
    for day in range(0, num_days, 20):
        for file in files_input.values():
            # MEmilio backend hdf5 example
            df = pm.extract_data(
                file, region_spec=None, column=None, date=day,
                filters={'Group': filter_age, 'InfectionState': symptomatic},
                file_format=file_format)
            if rki_data:
                path_rki_data = "/localdata1/pure/memilio/data/pydata/Germany/cases_all_county_ma7.json"
                df = pd.read_json(path_rki_data)
                # delete all entries in 'ID_County' that are not 9xxxxx (Bavaria)
                # df = df[df['ID_County'].astype(str).str.startswith('9')]
                # calculate start_date + day
                start_date_dt = datetime.strptime(start_date, '%Y-%m-%d')
                date_cur = start_date_dt + timedelta(days=day)
                date_cur = date_cur.strftime('%Y-%m-%d')

                date_before_curr = start_date_dt + timedelta(days=day - 1)
                date_before_curr = date_before_curr.strftime('%Y-%m-%d')
                df_day_before = df.loc[(df['Date'] == date_before_curr)].reset_index(
                    drop=True)
                df = df.loc[(df['Date'] == date_cur)].reset_index(drop=True)
                # df data is accumulated. Calculate the difference using df_day_before
                df['Confirmed'] = df['Confirmed'] - df_day_before['Confirmed']
                # delete all columns except 'ID_County' and 'Confirmed'
                df = df[['ID_County', 'Confirmed']]
                # rename ID_County to Region
                df = df.rename(columns={'ID_County': 'Region'})
                # rename Confirmed to Count
                df = df.rename(columns={'Confirmed': 'Count'})

            # all entries to numeric
            df = df.apply(pd.to_numeric, errors='coerce')

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

                # overwrite column Count with Count (rel)
                df['Count'] = df['Count (rel)']
                df = df.drop(columns=['Count (rel)'])

        # all entries to numeric
        df = df.apply(pd.to_numeric, errors='coerce')
        dfs_all = df

        filename = 'diff_day_' + str(day)
        if rki_data:
            filename += '_rki'

        dfs_all = dfs_all.apply(pd.to_numeric, errors='coerce')

        dfs_all_sorted = dfs_all.sort_values(by='Region')
        dfs_all_sorted = dfs_all_sorted.reset_index(drop=True)

        # delete entires in there Region does not start with 9 (Bavaria)
        dfs_all_sorted = dfs_all_sorted[dfs_all_sorted['Region'].astype(
            str).str.startswith('9')]

        if day == 0:
            data_t0 = dfs_all_sorted
            continue

        # calculate the relative difference to t0
        dfs_all_sorted['Count'] = dfs_all_sorted['Count'] / data_t0['Count']

        df_min = dfs_all_sorted['Count'].min()
        df_max = dfs_all_sorted['Count'].max()

        # scale all values in dfs_all_sorted['Count'] to 0-1
        dfs_all_sorted['Count'] = (dfs_all_sorted['Count'] - df_min) / \
            (df_max - df_min)

        if (calc_max_val < dfs_all_sorted['Count'].max()):
            calc_max_val = dfs_all_sorted['Count'].max()

        if (calc_min_val > dfs_all_sorted['Count'].min()):
            calc_min_val = dfs_all_sorted['Count'].min()

        print("day: " + str(day))

        pm.plot_map(norm,
                    dfs_all_sorted, scale_colors=[0, 1],
                    legend=['', ''],
                    title='',
                    plot_colorbar=False,
                    output_path=plot_dir,
                    fig_name=filename, dpi=900,
                    outercolor='white')
        plt.clf()
    print("max_val: " + str(calc_max_val))
    print("min_val: " + str(calc_min_val))


# plot_immunity_levels()
# plot_sheet()
# plot_infections_mobility()
# plot_total_infections(daily_transmissions=False)
# plot_total_infections(daily_transmissions=True)
map_plot_bavaria(res_dir + "/masks_0/p50/Results", relative=True)
# map_plot_bavaria_diff_t0(res_dir + "/masks_0/p50/Results", relative=True)
# map_plot_bavaria_diff_t0(res_dir + "/masks_0/p50/Results",
#                          relative=True, rki_data=True)
# map_plot_bavaria(res_dir + "/masks_0/p50/Results",
#                  relative=True, rki_data=True)
# plot_total_infections_bavaria(daily_transmissions=False)
# plot_icu_bavaria()

# plot_deaths()
