import json
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import matplotx

from scipy.ndimage import gaussian_filter

import seaborn as sns

sns.set_style("darkgrid")


path = "/localdata1/code/memilio/DigiHero_SARSCoV2_Data.csv"

path_pop_json = "/localdata1/code/memilio/data/pydata/Germany/county_population.json"

# transmission_delta = True
scenario = "new_parameters_adjicu_01"
# scenario = "octoberfest1_2"
# scenario = "amelag"
res_dir = "/localdata1/pure/memilio/results/" + scenario + "/secirts"
# if transmission_delta:
#     res_dir += "transmission1.5/"


path_results_rki = "results_paper/Results_rki_sum.h5"
plot_dir = "/localdata1/pure/memilio/plots/" + scenario

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

start_date = '2022-08-01'
num_days = 83

opacity = 0.15
lineWidth = 3
fontsize = 28
figsize = (13, 10)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

start_date_obj = datetime.strptime(start_date, '%Y-%m-%d')
end_date_obj = start_date_obj + timedelta(days=num_days)
end_date = end_date_obj.strftime('%Y-%m-%d')


def read_results(res_dir, filename, masks, comp):

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

    f.close()

    return comp_simulated[:num_days]


def plot_immunity_levels(quartil="50"):

    naive = [0, 2, 5, 8, 11, 14, 17, 20, 24, 27]
    pi = [1, 3, 6, 9, 12, 15, 18, 21, 25]
    ii = [4, 7, 10, 13, 16, 19, 22, 23, 26, 28]

    modes = ["0", "1_ffp2"]
    for mode in modes:
        N_sim = read_results(
            res_dir, "/mask_" + str(mode) + "/p" +
            quartil + "/Results_sum.h5", "",
            naive)
        PI_sim = read_results(
            res_dir, "/mask_" + str(mode) + "/p" +
            quartil + "/Results_sum.h5", "",
            pi)
        II_sim = read_results(
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
    SN_sim_m_false = read_results(
        res_dir, "/mask_0/p" + quartil + "/Results_sum.h5", "", [0])
    SPI_sim_m_false = read_results(
        res_dir, "/mask_0/p" + quartil + "/Results_sum.h5", "", [1])
    SII_sim_m_false = read_results(
        res_dir, "/mask_0/p" + quartil + "/Results_sum.h5", "", [23])

    SN_sim_m_true = read_results(
        res_dir, "/mask_1/p" + quartil + "/Results_sum.h5", "", [0])
    SPI_sim_m_true = read_results(
        res_dir, "/mask_1/p" + quartil + "/Results_sum.h5", "", [1])
    SII_sim_m_true = read_results(
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
        death_sim_m_false.append(read_results(
            res_dir, "/mask_0/p" + quartil + "/Results_sum.h5", "",
            list(range(24, 27))))

        # masks true
        death_sim_m_true.append(read_results(
            res_dir, "/mask_1/p" + quartil + "/Results_sum.h5", "",
            list(range(24, 27))))

        # ff2 masks
        death_sim_m_true_ffp2.append(read_results(
            res_dir, "/mask_1_ffp2/p" + quartil + "/Results_sum.h5", "",
            list(range(24, 27))))

    # read rki excel data
    results_rki_extrapolated = read_results(
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


def plot_icu():
    quartils = ["05", "25", "50", "75", "95"]

    icu_sim_m_false = []
    # icu_sim_m_true = []
    icu_sim_m_true_ffp2 = []

    for quartil in quartils:
        # masks false
        icu_sim_m_false.append(read_results(
            res_dir, "/masks_0/p" + quartil + "/Results_sum.h5", "",
            [20, 21, 22]))

        # masks true
        # icu_sim_m_true.append(read_results(
        #     res_dir, "/masks_1/p" + quartil + "/Results_sum.h5", "",
        # [20, 21, 22]))

        # ff2 masks
        icu_sim_m_true_ffp2.append(read_results(
            res_dir, "/masks_1_ffp2/p" + quartil + "/Results_sum.h5", "",
            [20, 21, 22]))

    path_json_icu = "data//pydata//Germany//germany_divi_ma7.json"
    df = pd.read_json(path_json_icu)
    dates = pd.date_range(start=start_date, end=end_date)[:num_days]

    filtered_df = df.loc[(df['Date'] >= start_date) &
                         (df['Date'] <= end_date)][:num_days]

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

    plt.plot(dates, icu_sim_m_true_ffp2[1],
             color=colors[2], linewidth=lineWidth, linestyle='--')
    plt.plot(dates, icu_sim_m_true_ffp2[2],
             color=colors[2], label='FFP2 Masks', linewidth=lineWidth, linestyle='-')
    plt.plot(dates, icu_sim_m_true_ffp2[3],
             color=colors[2], linewidth=lineWidth, linestyle='--')
    plt.fill_between(
        dates, icu_sim_m_true_ffp2[1],
        icu_sim_m_true_ffp2[3],
        alpha=0.2, color=colors[2])

    plt.plot(filtered_df['Date'], filtered_df['ICU'],
             label='ICU Occupancy reported', color=colors[1], linewidth=lineWidth)

    dates = pd.date_range(start=start_date, end=end_date)[
        : len(icu_sim_m_true_ffp2)]
    weeks = pd.date_range(start=start_date, end=end_date, freq='7D')
    plt.xticks(weeks, fontsize=15, rotation=45)

    plt.xlabel('Date', fontsize=18)
    plt.ylabel('ICU Occupancy', fontsize=18)
    # plt.yscale('log')
    plt.xticks(rotation=45)
    plt.grid(True)
    # plt.yscale('log')
    plt.tight_layout()
    plt.legend(fontsize=11)
    plt.savefig(plot_dir + "/icu_germany_masks.png")
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
    cases_digihero = data_digihero['N_Inf'] * pop_germany / data_digihero['N']

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

    masks = ["0", "1_ffp2"]  # "0", "1", "1_ffp2"]
    inf_mobility_node_mask_0 = []
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
            inf_mobility_node_mask_0.append(data_p25)
            inf_mobility_node_mask_0.append(data_p50)
            inf_mobility_node_mask_0.append(data_p75)

        else:
            inf_mobility_node_mask_1_ffp2.append(data_p25)
            inf_mobility_node_mask_1_ffp2.append(data_p50)
            inf_mobility_node_mask_1_ffp2.append(data_p75)

    # values are cumulated. Calculate the difference
    inf_mobility_node_mask_0 = np.array(inf_mobility_node_mask_0)
    inf_mobility_node_mask_0 = np.diff(inf_mobility_node_mask_0, axis=1)

    inf_mobility_node_mask_1_ffp2 = np.array(inf_mobility_node_mask_1_ffp2)
    inf_mobility_node_mask_1_ffp2 = np.diff(
        inf_mobility_node_mask_1_ffp2, axis=1)

    inf_mobility_node_mask_0 = np.array(inf_mobility_node_mask_0)
    inf_mobility_node_mask_0_quartils = np.percentile(
        inf_mobility_node_mask_0, [5, 25, 50, 75, 95], axis=0)

    # inf_mobility_node_mask_1 = np.array(inf_mobility_node_mask_1)
    # inf_mobility_node_mask_1_quartils = np.percentile(
    #     inf_mobility_node_mask_1, [5, 25, 50, 75, 95], axis=0)

    inf_mobility_node_mask_1_ffp2 = np.array(inf_mobility_node_mask_1_ffp2)
    inf_mobility_node_mask_1_ffp2_quartils = np.percentile(
        inf_mobility_node_mask_1_ffp2, [5, 25, 50, 75, 95], axis=0)

    # all days between start and end date
    dates = pd.date_range(start=start_date, end=end_date)[
        : inf_mobility_node_mask_1_ffp2_quartils.shape[1]]
    weeks = pd.date_range(start=start_date, end=end_date, freq='7D')

    # no masks
    plt.plot(
        dates, inf_mobility_node_mask_0_quartils[1], color=colors[0], linestyle='--', linewidth=lineWidth)
    plt.plot(
        dates, inf_mobility_node_mask_0_quartils[2], color=colors[0], label='No face masks', linewidth=lineWidth)
    plt.plot(
        dates, inf_mobility_node_mask_0_quartils[3], color=colors[0], linestyle='--', linewidth=lineWidth)
    plt.fill_between(
        dates, inf_mobility_node_mask_0_quartils[1], inf_mobility_node_mask_0_quartils[3], color=colors[0], alpha=opacity)

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
    plt.plot(dates, cases_digihero[:inf_mobility_node_mask_1_ffp2_quartils.shape[1]],
             label='Extrapolated DigiHero Reports ', color=colors[3], linewidth=lineWidth)
    plt.xlabel('Date', fontsize=18)
    if daily_transmissions:
        plt.ylabel('Total Infections', fontsize=18)
    else:
        plt.ylabel('Daily Symptomatic', fontsize=18)
    # plt.yscale('log')
    plt.xticks(weeks, fontsize=15, rotation=45)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(fontsize=11)
    if daily_transmissions:
        plt.savefig(plot_dir + "/total_infections.png")
    else:
        plt.savefig(plot_dir + "/total_symptoms.png")

    plt.clf()


def plot_waning(quartil="50"):
    masks = ["0", "1"]
    inf_mobility_node_mask_0_pi = []
    inf_mobility_node_mask_0_ii = []
    inf_mobility_node_mask_1_pi = []
    inf_mobility_node_mask_1_ii = []
    quartils = ["05", "25", "50", "75", "95"]
    try:
        indx_quartil = quartils.index(quartil)
    except ValueError:
        print("Quartil not found")
        return 0
    for mode in masks:
        dir = res_dir + "/mask_" + mode + "/flows_mb"
        files = os.listdir(dir)
        num_files = len(files)

        if mode == "0":
            for indx in range(num_files):
                df = pd.read_csv(
                    dir + "/results_run_" + str(indx) + ".txt", sep=" ")
                inf_mobility_node_mask_0_pi.append(
                    df['Waning_Partial_Immunity'].values[:-1])
                inf_mobility_node_mask_0_ii.append(
                    df['Waning_Improved_Immunity'].values[:-1])

        else:
            for indx in range(num_files):
                df = pd.read_csv(
                    dir + "/results_run_" + str(indx) + ".txt", sep=" ")
                inf_mobility_node_mask_1_pi.append(
                    df['Waning_Partial_Immunity'].values[:-1])
                inf_mobility_node_mask_1_ii.append(
                    df['Waning_Improved_Immunity'].values[:-1])

    inf_mobility_node_mask_0_pi = np.array(inf_mobility_node_mask_0_pi)
    inf_mobility_node_mask_0_pi_quartils = np.percentile(
        inf_mobility_node_mask_0_pi, [5, 25, 50, 75, 95], axis=0)
    inf_mobility_node_mask_0_ii = np.array(inf_mobility_node_mask_0_ii)
    inf_mobility_node_mask_0_ii_quartils = np.percentile(
        inf_mobility_node_mask_0_ii, [5, 25, 50, 75, 95], axis=0)

    inf_mobility_node_mask_1_pi = np.array(inf_mobility_node_mask_1_pi)
    inf_mobility_node_mask_1_pi_quartils = np.percentile(
        inf_mobility_node_mask_1_pi, [5, 25, 50, 75, 95], axis=0)
    inf_mobility_node_mask_1_ii = np.array(inf_mobility_node_mask_1_ii)
    inf_mobility_node_mask_1_ii_quartils = np.percentile(
        inf_mobility_node_mask_1_ii, [5, 25, 50, 75, 95], axis=0)

    inf_mobility_node_mask_0_pi_quartils = np.diff(
        inf_mobility_node_mask_0_pi_quartils, axis=1)
    inf_mobility_node_mask_0_ii_quartils = np.diff(
        inf_mobility_node_mask_0_ii_quartils, axis=1)
    inf_mobility_node_mask_1_pi_quartils = np.diff(
        inf_mobility_node_mask_1_pi_quartils, axis=1)
    inf_mobility_node_mask_1_ii_quartils = np.diff(
        inf_mobility_node_mask_1_ii_quartils, axis=1)

    dates = pd.date_range(start=start_date, end=end_date)[
        : inf_mobility_node_mask_1_ii_quartils.shape[1]]

    plt.plot(dates, inf_mobility_node_mask_0_pi_quartils[indx_quartil],
             label='Waning partial immunity, no masks in transport')
    plt.plot(dates, inf_mobility_node_mask_1_pi_quartils[indx_quartil],
             label='Waning partial immunity, masks in transport')
    plt.plot(dates, inf_mobility_node_mask_0_ii_quartils[indx_quartil],
             label='Waning improved immunity, no masks in transport')
    plt.plot(dates, inf_mobility_node_mask_1_ii_quartils[indx_quartil],
             label='Waning improved immunity, masks in transport')
    plt.xlabel('Date', fontsize=18)
    plt.ylabel('Immunity Wanings', fontsize=18)
    # plt.yscale('log')
    plt.title(
        'Number of persons whose immunity wanes during the simulation',
        fontsize=10)
    plt.xticks(rotation=45)
    plt.grid(True)
    # plt.yscale('log')
    plt.tight_layout()
    plt.legend(fontsize=11)
    plt.savefig(plot_dir + "/wanings.png")
    plt.clf()

    return 0


def plot_infections_mobility(quartil="50"):
    commuter = np.loadtxt(
        "/localdata1/code/memilio/test/commuter_migration_with_locals.txt")

    path_tt = "/localdata1/code/memilio/test/travel_times_pathes.txt"
    file = "transmission_mobility_run_"
    travel_times = np.loadtxt(path_tt)
    masks = ["0", "1_ffp2"]
    flow_indx = [0, 17, 33]
    #     flow_indx = [3, 5, 20, 22, 36, 38]  # flows to symptomatic infected
    inf_mobility_node_mask_0 = []
    # inf_mobility_node_mask_1 = []
    inf_mobility_node_mask_1_ffp2 = []
    for mode in masks:
        dir = res_dir + "/masks_" + mode + "/mobility/flows/"
        data_p25 = read_results_quartil(dir, flow_indx, "25")
        data_p50 = read_results_quartil(dir, flow_indx, "50")
        data_p75 = read_results_quartil(dir, flow_indx, "75")

        if mode == "0":
            inf_mobility_node_mask_0.append(data_p25)
            inf_mobility_node_mask_0.append(data_p50)
            inf_mobility_node_mask_0.append(data_p75)

        else:
            inf_mobility_node_mask_1_ffp2.append(data_p25)
            inf_mobility_node_mask_1_ffp2.append(data_p50)
            inf_mobility_node_mask_1_ffp2.append(data_p75)

    # values are cumulated. Calculate the difference
    inf_mobility_node_mask_0 = np.array(inf_mobility_node_mask_0)
    inf_mobility_node_mask_0 = np.diff(inf_mobility_node_mask_0, axis=1)

    inf_mobility_node_mask_1_ffp2 = np.array(inf_mobility_node_mask_1_ffp2)
    inf_mobility_node_mask_1_ffp2 = np.diff(
        inf_mobility_node_mask_1_ffp2, axis=1)

    dates = pd.date_range(start=start_date, end=end_date)[
        : inf_mobility_node_mask_0.shape[1]]
    weeks = pd.date_range(start=start_date, end=end_date, freq='7D')

    # plot summed transmissions
    plt.plot(
        dates, inf_mobility_node_mask_0[0], color=colors[0], linestyle='--', linewidth=lineWidth)
    plt.plot(
        dates, inf_mobility_node_mask_0[1], color=colors[0], label='No face masks', linewidth=lineWidth)
    plt.plot(
        dates, inf_mobility_node_mask_0[2], color=colors[0], linestyle='--', linewidth=lineWidth)
    plt.fill_between(
        dates, inf_mobility_node_mask_0[0], inf_mobility_node_mask_0[2], color=colors[0], alpha=opacity)

    plt.plot(
        dates, inf_mobility_node_mask_1_ffp2[0], color=colors[2], linestyle='--', linewidth=lineWidth)
    plt.plot(
        dates, inf_mobility_node_mask_1_ffp2[1], color=colors[2], label='FFP2 masks', linewidth=lineWidth)
    plt.plot(
        dates, inf_mobility_node_mask_1_ffp2[2], color=colors[2], linestyle='--', linewidth=lineWidth)
    plt.fill_between(
        dates, inf_mobility_node_mask_1_ffp2[0], inf_mobility_node_mask_1_ffp2[2], color=colors[2], alpha=opacity)

    plt.gca().xaxis.set_major_locator(
        mdates.DayLocator(interval=7))
    plt.gca().xaxis.set_major_formatter(
        mdates.DateFormatter('%Y-%m-%d'))
    plt.xlabel('Date', fontsize=20)
    plt.ylabel('Daily Transmissions', fontsize=20)
    plt.xticks(weeks, fontsize=15, rotation=45)
    plt.xlim(dates[0], dates[-1])
    plt.grid(True)
    plt.yscale('log')
    # # adjust logarithmic ticks
    # plt.yticks([10, 100, 1000, 10000, 100000, 1000000, 10000000])
    plt.tight_layout()
    plt.legend(fontsize=11)
    plt.savefig(plot_dir + "/transmissions_mobility_daily.png")
    plt.clf()

    # for no masks plot the 3 different commuting types
    # plt.plot(dates, inf_mobility_node_mask_0_quartils[1][:, 0],
    #          label='Inner county')
    # plt.plot(dates, inf_mobility_node_mask_0_quartils[1][:, 1],
    #          label='Under 1 hour travel time')
    # plt.plot(dates, inf_mobility_node_mask_0_quartils[1][:, 2],
    #          label='Over 1 hour travel time')
    # plt.gca().xaxis.set_major_locator(
    #     mdates.DayLocator(interval=7))  # Jeden 7.Tag ein Tick
    # plt.gca().xaxis.set_major_formatter(
    #     mdates.DateFormatter('%Y-%m-%d'))  # Format des Datums
    # plt.xlabel('Date', fontsize=20)
    # plt.ylabel('Daily Transmissions', fontsize=20)
    # plt.xticks(weeks, rotation=45, fontsize=15)
    # plt.xlim(dates[0], dates[-1])
    # plt.grid(True)
    # plt.yscale('log')
    # plt.tight_layout()
    # plt.legend(fontsize=11)
    # plt.savefig(plot_dir + "/transmissions_daily_no_masks.png")
    # plt.clf()

    return 0


# plot_immunity_levels()
# plot_sheet()
# plot_infections_mobility()
# plot_total_infections(daily_transmissions=False)
# plot_total_infections(daily_transmissions=True)
plot_icu()
# plot_susceptible()

# plot_deaths_quartils()
# plot_waning()
# plot_deaths(False)

# plot_deaths()
