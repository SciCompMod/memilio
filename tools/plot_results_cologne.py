import json
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from datetime import datetime, timedelta
import memilio.plot.plotMap as pm
from memilio.epidata import geoModificationGermany as geoger

# plt.style.use('ggplot')

old_model = True

# res_dir_old = "/localdata1/code/memilio/results_paper/mask_0_cologne_old"
res_dir_old = "/localdata1/code/memilio/results_paper/mask_0_cologne_old/"

# res_dir_new = "/localdata1/code/memilio/results_paper/mask_0_cologne"
res_dir_new = "/localdata1/code/memilio/results_paper/mask_0_cologne/"


path_results_rki = "/localdata1/code/memilio/results_paper/Results_rki_sum.h5"
plot_dir = "/localdata1/code/memilio/results_paper/Plots"
start_date = '2022-08-01'
num_days = 90

fontsize_title = 25

start_date_obj = datetime.strptime(start_date, '%Y-%m-%d')
end_date_obj = start_date_obj + timedelta(days=num_days)
end_date = end_date_obj.strftime('%Y-%m-%d')

opacity = 0.15
lineWidth = 2
fontsize = 28
figsize = (13, 10)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


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

    # After you are done
    f.close()

    return comp_simulated[:num_days]


# sum results except cologne
def read_results_quartil_cologne(path, comp, quartil):

    path_quartil = path + "//p" + str(quartil) + "//Results.h5"
    f = h5py.File(path_quartil, 'r')

    id_cologne = "5315"

    # Get the HDF5 group; key needs to be a group name from above
    group = f[id_cologne]

    # This assumes group[some_key_inside_the_group] is a dataset,
    # and returns a np.array:
    # time = group['Time'][()]
    total = group['Total'][()]
    # group1 = group['Group1'][()]

    comp_simulated = np.sum(total[:, comp], axis=1)

    # After you are done
    f.close()

    return comp_simulated[:num_days]


def read_results_quartil_except_cologne(path, comp, quartil):

    path_quartil = path + "//p" + str(quartil) + "//Results.h5"
    f = h5py.File(path_quartil, 'r')

    infected = np.zeros(num_days + 1)

    # iterate over all keys
    for key in f.keys():
        if key != "5315":
            group = f[key]
            total = group['Total'][()]
            comp_simulated = np.sum(total[:, comp], axis=1)
            infected += comp_simulated

    # After you are done
    f.close()

    return infected


def plot_total_pop():
    total_pop = []
    res_dir = "/localdata1/code/memilio/results_paper/mask_0_cologne/"

    quartils = ["05", "25", "50", "75", "95"]
    if old_model:
        all_comps = [x for x in range(0, 27)]
    else:
        all_comps = [x for x in range(0, 29)]
    for quartil in quartils:
        total_pop.append(read_results_quartil(
            res_dir, all_comps, quartil))

    num_days = len(total_pop[0])

    print("Anzahl Tage = " + str(num_days))

    # # prüfe ob die Summe der Populationen gleich bleibt
    # for i in range(1, num_days):
    #     if total_pop[0][i] != total_pop[0][i-1]:
    #         print("Summe der Populationen ist nicht gleich. Diff ) = " +
    #               str(total_pop[0][i] - total_pop[0][i-1]))

    # print("Summe der Populationen ist gleich")


def plot_sheet():

    modes = ["0"]

    quartils = ["05", "25", "50", "75", "95"]
    for mode in modes:
        sim_exposed, sim_INS, sim_Isy, sim_hosp, sim_icu, sim_sus, sim_dead, sim_ti = [
        ], [], [], [], [], [], [], []
        path_quartils = res_dir
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
            if not old_model:
                sim_ti.append(read_results_quartil(
                    path_quartils, [27, 28], quartil))

        num_days = len(sim_exposed[0])

        fig, axs = plt.subplots(4, 2, figsize=(8.27, 11.69))

        # Plotten der Ergebnisse für jede Simulation in einem Subplot

        # passe die plots so an, dass in jedem plot jeweils alle quartile geplottet werden
        for i in range(5):
            axs[0, 0].plot(range(num_days), sim_exposed[i],
                           label="P" + quartils[i])
            axs[0, 1].plot(range(num_days), sim_INS[i],
                           label="P" + quartils[i])
            axs[1, 0].plot(range(num_days), sim_Isy[i],
                           label="P" + quartils[i])
            axs[1, 1].plot(range(num_days), sim_hosp[i],
                           label="P" + quartils[i])
            axs[2, 0].plot(range(num_days), sim_icu[i],
                           label="P" + quartils[i])
            axs[2, 1].plot(range(num_days), sim_sus[i],
                           label="P" + quartils[i])
            axs[3, 0].plot(range(num_days), sim_dead[i],
                           label="P" + quartils[i])
            if not old_model:
                axs[3, 1].plot(range(num_days), sim_ti[i],
                               label="P" + quartils[i])

        axs[0, 0].set_title("Exposed")
        axs[0, 1].set_title("Infected no symptoms")
        axs[1, 0].set_title("Infected symptomatic")
        axs[1, 1].set_title("Severe")
        axs[2, 0].set_title("ICU")
        axs[2, 1].set_title("Susceptible")
        axs[3, 0].set_title("Deaths")
        if not old_model:
            axs[3, 1].set_title("Temporary immunity")

        axs[0, 0].set_yscale('log')
        axs[0, 1].set_yscale('log')
        axs[1, 0].set_yscale('log')
        axs[1, 1].set_yscale('log')
        axs[3, 1].set_yscale('log')

        axs[0, 0].fill_between(range(num_days), sim_exposed[0],
                               sim_exposed[4], alpha=0.2)
        axs[0, 1].fill_between(range(num_days), sim_INS[0],
                               sim_INS[4], alpha=0.2)
        axs[1, 0].fill_between(range(num_days), sim_Isy[0],
                               sim_Isy[4], alpha=0.2)
        axs[1, 1].fill_between(range(num_days), sim_hosp[0],
                               sim_hosp[4], alpha=0.2)
        axs[2, 0].fill_between(range(num_days), sim_icu[0],
                               sim_icu[4], alpha=0.2)
        axs[2, 1].fill_between(range(num_days), sim_sus[0],
                               sim_sus[4], alpha=0.2)
        axs[3, 0].fill_between(range(num_days), sim_dead[0],
                               sim_dead[4], alpha=0.2)
        if not old_model:
            axs[3, 1].fill_between(range(num_days), sim_ti[0],
                                   sim_ti[4], alpha=0.2)

        for ax in axs.flat:
            ax.set(xlabel="Tage", ylabel="Werte")
            ax.grid(True)

        plt.tight_layout()

        if old_model:
            plt.savefig(plot_dir + "/cologne_combined_plots_mask" +
                        mode + "_old_model.png")
        else:
            plt.savefig(
                plot_dir + "/cologne_combined_plots_mask" + mode + ".png")


def plot_infections():
    path_old = "/localdata1/code/memilio/results_paper/mask_0_cologne_old"
    path_new = "/localdata1/code/memilio/results_paper/mask_0_cologne"
    quartils = "50"
    infected = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    res_old = read_results_quartil(path_old, infected, quartils)
    res_new = read_results_quartil(path_new, infected, quartils)

    plt.plot(res_old,
             label='Exisiting method')
    plt.plot(res_new,
             label='Introduced method')
    plt.xlabel('Date', fontsize=18)
    plt.ylabel('Infected Individuals', fontsize=18)
    plt.xticks(np.arange(0, 21, 5))
    plt.xlim([0, 20])
    plt.grid(True)
    plt.tight_layout()
    plt.legend(fontsize=11)
    plt.savefig(plot_dir + "/cologne_sim_infections.png")
    plt.clf()


def read_county_data(path, day):
    quartils = "50"
    infected = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                12, 13, 14, 15, 16, 17, 18, 19, 20, 21]

    files_input = {
        'Data set 1': path}
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

    for file in files_input.values():
        # MEmilio backend hdf5 example

        df = pm.extract_data(
            file, region_spec=None, column=None, date=day,
            filters={'Group': filter_age, 'InfectionState': infected},
            file_format=file_format)
    return df


def plot_infections_neighboors(quartil="50"):
    counties_considered = geoger.get_county_ids()
    infected = [5, 6, 7, 8, 9, 10, 11,
                12, 13, 14, 15, 16, 17, 18, 19, 20, 21]  # 2, 3, 4,

    infected_cologne_old = read_results_quartil_except_cologne(
        res_dir_old, infected, quartil)
    infected_cologne_new = read_results_quartil_except_cologne(
        res_dir_new, infected, quartil)

    # old_p25 = read_results_quartil_except_cologne(
    #     res_dir_old, infected, "25")
    # old_p50 = read_results_quartil_except_cologne(
    #     res_dir_old, infected, "50")
    # old_p75 = read_results_quartil_except_cologne(
    #     res_dir_old, infected, "75")

    # new_p25 = read_results_quartil_except_cologne(
    #     res_dir_new, infected, "25")
    # new_p50 = read_results_quartil_except_cologne(
    #     res_dir_new, infected, "50")
    # new_p75 = read_results_quartil_except_cologne(
    #     res_dir_new, infected, "75")

    # Plot for Connected Counties
    linewidth = 3.0
    fontsize_axes = 16

    plt.figure(figsize=(10, 6))

    dates = pd.date_range(start=start_date, end=end_date)[
        :num_days+1]
    # erst den ersten und danach jeden 7. tag
    weeks = pd.date_range(start=start_date, end=end_date, freq='7D')

    # plt.plot(
    #     dates, old_p25, color=colors[0], linestyle='--', linewidth=lineWidth)
    # plt.plot(
    #     dates, old_p50, color=colors[0], label='Existing method', linewidth=lineWidth)
    # plt.plot(
    #     dates, old_p75, color=colors[0], linestyle='--', linewidth=lineWidth)
    # plt.fill_between(
    #     dates, old_p25, old_p75, color=colors[0], alpha=opacity)

    # plt.plot(
    #     dates, new_p25, color=colors[2], linestyle='--', linewidth=lineWidth)
    # plt.plot(
    #     dates, new_p50, color=colors[2], label='Novel method', linewidth=lineWidth)
    # plt.plot(
    #     dates, new_p75, color=colors[2], linestyle='--', linewidth=lineWidth)
    # plt.fill_between(
    #     dates, new_p25, new_p75, color=colors[2], alpha=opacity)

    # plt.plot(dates, infected_cologne_old,
    #          label='Existing method', linewidth=linewidth)
    # plt.plot(dates, infected_cologne_new,
    #          label='Novel method', linewidth=linewidth)
    plt.xlabel('Days', fontsize=18)
    plt.ylabel('Infected Individuals', fontsize=18)
    # plt.xticks(np.arange(0, num_days+1, 5), fontsize=fontsize_axes)
    plt.yticks(fontsize=fontsize_axes)
    # plt.xlim([1, 10])
    plt.xticks(weeks, fontsize=15, rotation=45)
    plt.title('Infected people in all Counties except Cologne',
              fontsize=fontsize_title)
    plt.grid(True)
    # plt.yscale('log')
    plt.tight_layout()
    plt.legend(fontsize=18)
    plt.savefig(plot_dir + "/cologne_sim_infections_not_cologne.png")
    plt.clf()


def plot_infections_cologne(quartil="50"):
    infected = [5, 6, 7, 8, 9, 10, 11,
                12, 13, 14, 15, 16, 17, 18, 19, 20, 21]  # 2, 3, 4,
    infected_cologne_old = read_results_quartil_cologne(
        res_dir_old, infected, quartil)
    infected_cologne_new = read_results_quartil_cologne(
        res_dir_new, infected, "25")

    # id_cologne = 161
    # infected_cologne_old = []
    # infected_cologne_new = []

    # for day in range(0, 21):
    #     df_old = read_county_data(path_old, day)
    #     infected_cologne_old.append(df_old.loc[id_cologne]["Count"])

    #     df_new = read_county_data(path_new, day)
    #     infected_cologne_new.append(df_new.loc[id_cologne]["Count"])

    # Plot for Connected Counties
    linewidth = 3.0
    fontsize_axes = 16

    dates = pd.date_range(start=start_date, end=end_date)[
        :num_days]
    # erst den ersten und danach jeden 7. tag
    weeks = pd.date_range(start=start_date, end=end_date, freq='7D')

    plt.figure(figsize=(10, 6))

    old_p25 = read_results_quartil_cologne(
        res_dir_old, infected, "25")
    old_p50 = read_results_quartil_cologne(
        res_dir_old, infected, "50")
    old_p75 = read_results_quartil_cologne(
        res_dir_old, infected, "75")

    new_p25 = read_results_quartil_cologne(
        res_dir_new, infected, "25")
    new_p50 = read_results_quartil_cologne(
        res_dir_new, infected, "50")
    new_p75 = read_results_quartil_cologne(
        res_dir_new, infected, "75")

    plt.plot(
        dates, old_p25, color=colors[0], linestyle='--', linewidth=lineWidth)
    plt.plot(
        dates, old_p50, color=colors[0], label='Existing method', linewidth=lineWidth)
    plt.plot(
        dates, old_p75, color=colors[0], linestyle='--', linewidth=lineWidth)
    plt.fill_between(
        dates, old_p25, old_p75, color=colors[0], alpha=opacity)

    plt.plot(
        dates, new_p25, color=colors[2], linestyle='--', linewidth=lineWidth)
    plt.plot(
        dates, new_p50, color=colors[2], label='Novel method', linewidth=lineWidth)
    plt.plot(
        dates, new_p75, color=colors[2], linestyle='--', linewidth=lineWidth)
    plt.fill_between(
        dates, new_p25, new_p75, color=colors[2], alpha=opacity)

    # plt.plot(dates, infected_cologne_old[:num_days],
    #          label='Existing method', linewidth=linewidth)
    # plt.plot(dates, infected_cologne_new[:num_days],
    #          label='Novel method', linewidth=linewidth)
    plt.xlabel('Days', fontsize=18)
    plt.xticks(weeks, fontsize=15, rotation=45)
    plt.ylabel('Infected Individuals', fontsize=18)
    # plt.xticks(np.arange(0,  num_days+1, 5), fontsize=fontsize_axes)
    plt.yticks(fontsize=fontsize_axes)
    plt.grid(True)
    plt.title('Infected people in Cologne', fontsize=fontsize_title)
    plt.xlim([dates[0], dates[-1]])
    plt.tight_layout()
    # log scale
    # plt.yscale('log')
    plt.legend(fontsize=18)
    # plt.show()
    plt.savefig(plot_dir + "/cologne_sim_infections_cologne_log.png")


# plot_total_pop()
# plot_sheet()
plot_infections_neighboors()
# plot_infections_cologne()
# plot_infections()
