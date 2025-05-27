# Python script to analyze bs runs
# input is a bs run folder with the following structure:
# bs_run_folder has a txt file for each bs run
# each txt file has a line for each time step
# each line has a column for each compartment as well as the timestep
# each column has the number of individuals in that compartment
# the first line of each txt file is the header

import sys
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import h5py
from datetime import datetime
from matplotlib.dates import DateFormatter

fontsize = 20


def plot_infections_loc_types_avarage(path):
    # 50-percentile
    f_p50 = h5py.File(
        path+"/infection_per_location_type_per_age_group/0/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]

    # 25-percentile
    f_p25 = h5py.File(
        path+"/infection_per_location_type_per_age_group/0/p05/Results.h5", 'r')
    p25_bs = f_p25['0']
    total_25 = p25_bs['Total'][()]

    # 75-percentile
    f_p75 = h5py.File(
        path + "/infection_per_location_type_per_age_group/0/p95/Results.h5", 'r')
    p75_bs = f_p75['0']
    total_75 = p75_bs['Total'][()]

    time = p50_bs['Time'][()]

    plot_infection_per_location_type_mean(
        time, total_50, total_25, total_75)

def plot_infection_per_location_type_mean(x, y50, y25, y75):

    plt.figure('Infection_location_types')
    plt.title('At which location type an infection happened, avaraged over all runs')

    color_plot = matplotlib.colormaps.get_cmap('Set1').colors

    states_plot = [0, 1, 2, 3, 4, 10]
    legend_plot = ['Home', 'School', 'Work',
                   'SocialEvent', 'BasicsShop','Event']

    for i in states_plot:
        # rolling average#
        ac_color = color_plot[i%len(color_plot)]
        if(i > len(color_plot)):
            ac_color = "black"
        
        # we need to sum up every 24 hours
        indexer = pd.api.indexers.FixedForwardWindowIndexer(window_size=24)
        np_y50 = pd.DataFrame(y50[:, i]).rolling(window=indexer, min_periods=1).sum().to_numpy()
        np_y50=np_y50[0::24].flatten()
        # now smoothen this with a gaussian filter
        np_y50 = gaussian_filter1d(np_y50, sigma=1, mode='nearest')
        
        plt.plot(x[0::24], np_y50, color=ac_color)

    plt.legend(legend_plot)

    # currently the x axis has the values of the time steps, we need to convert them to dates and set the x axis to dates
    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in x]
    xx = [xx[i].strftime('%Y-%m-%d') for i in range(len(xx))]
    # but just take every 10th date to make it more readable
    plt.gca().set_xticks(x[::150])
    plt.gca().set_xticklabels(xx[::150])
    plt.gcf().autofmt_xdate()

    plt.xlabel('Date')
    plt.ylabel('Number of individuals')
    plt.show()

def plot_infection_states_results(path):
    # 50-percentile
    f_p50 = h5py.File(
        path+"/infection_state_per_age_group/0/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]
    # 25-percentile
    f_p25 = h5py.File(
        path+"/infection_state_per_age_group/0/p05/Results.h5", 'r')
    p25_bs = f_p25['0']
    total_25 = p25_bs['Total'][()]
    # 75-percentile
    f_p75 = h5py.File(
        path + "/infection_state_per_age_group/0/p95/Results.h5", 'r')
    p75_bs = f_p75['0']
    total_75 = p75_bs['Total'][()]

    time = p50_bs['Time'][()]

    plot_infection_states_individual(
        time, p50_bs, p25_bs, p75_bs)
    plot_infection_states(time, total_50, total_25, total_75)

def plot_infection_states(x, y50, y25, y75, y_real=None):
    plt.figure('Infection_states')
    plt.title('Infection states')

    color_plot = matplotlib.colormaps.get_cmap('Set1').colors

    states_plot = [1, 2, 3, 4, 5, 7]
    legend_plot = ['E', 'I_NSymp', 'I_Symp',
                   'I_Sev', 'I_Crit', 'Dead', 'Sm. re. pos.']

    for i in states_plot:
        plt.plot(x, y50[:, i], color=color_plot[i])


    plt.legend(legend_plot)
    for i in states_plot:
        plt.fill_between(x, y50[:, i], y25[:, i],
                         alpha=0.5, color=color_plot[i])
        plt.fill_between(x, y50[:, i], y75[:, i],
                         alpha=0.5, color=color_plot[i])

    # currently the x axis has the values of the time steps, we need to convert them to dates and set the x axis to dates
    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in x]
    xx = [xx[i].strftime('%Y-%m-%d') for i in range(len(xx))]
    # but just take every 10th date to make it more readable
    plt.gca().set_xticks(x[::150])
    plt.gca().set_xticklabels(xx[::150])
    plt.gcf().autofmt_xdate()

    plt.xlabel('Time')
    plt.ylabel('Number of individuals')
    plt.show()

def plot_infection_states_individual(x, p50_bs, p25_bs, p75_bs):


    age_group_access = ['Group1', 'Group2', 'Group3',
                        'Group4', 'Group5', 'Group6', 'Total']

    color_plot = matplotlib.colormaps.get_cmap('Set1').colors

    fig, ax = plt.subplots(6, len(age_group_access), constrained_layout=True)
    fig.set_figwidth(20)
    fig.set_figheight(9)
    for j, count in zip(age_group_access, range(len(age_group_access))):
        y50 = p50_bs[j][()]
        y25 = p25_bs[j][()]
        y75 = p75_bs[j][()]


        # infeced no symptoms
        ax_infected_no_symptoms = ax[0, count]
        ax_infected_no_symptoms.set_xlabel('time (days)')
        ax_infected_no_symptoms.plot(
            x, y50[:, 1], color=color_plot[count], label='y50')
        ax_infected_no_symptoms.fill_between(
            x, y50[:, 1], y25[:, 1], alpha=0.5, color=color_plot[count])
        ax_infected_no_symptoms.fill_between(
            x, y50[:, 1], y75[:, 1], alpha=0.5, color=color_plot[count])
        ax_infected_no_symptoms.tick_params(axis='y')
        ax_infected_no_symptoms.title.set_text(
            '#Infected_no_symptoms, Age{}'.format(j))
        ax_infected_no_symptoms.legend(['Simulation'])

        # Infected_symptoms
        ax_infected_symptoms = ax[1, count]
        ax_infected_symptoms.set_xlabel('time (days)')
        ax_infected_symptoms.plot(
            x, y50[:, 2], color=color_plot[count], label='y50')
        ax_infected_symptoms.fill_between(
            x, y50[:, 2], y25[:, 2], alpha=0.5, color=color_plot[count])
        ax_infected_symptoms.fill_between(
            x, y50[:, 2], y75[:, 2], alpha=0.5, color=color_plot[count])
        ax_infected_symptoms.tick_params(axis='y')
        ax_infected_symptoms.title.set_text(
            '#Infected_symptoms, Age{}'.format(j))
        ax_infected_symptoms.legend(['Simulation'])

        # Severe
        ax_severe = ax[2, count]
        ax_severe.set_xlabel('time (days)')
        ax_severe.plot(x, y50[:, 4], color=color_plot[count], label='y50')
        ax_severe.fill_between(
            x, y50[:, 4], y25[:, 4], alpha=0.5, color=color_plot[count])
        ax_severe.fill_between(
            x, y50[:, 4], y75[:, 4], alpha=0.5, color=color_plot[count])
        ax_severe.tick_params(axis='y')
        ax_severe.title.set_text('#Severe, Age{}'.format(j))
        ax_severe.legend(['Simulation'])

        # Critical
        ax_critical = ax[3, count]
        ax_critical.set_xlabel('time (days)')
        ax_critical.plot(x, y50[:, [5]], color=color_plot[count], label='y50')
        ax_critical.fill_between(
            x, y50[:, 5], y25[:, 5], alpha=0.5, color=color_plot[count])
        ax_critical.fill_between(
            x, y50[:, 5], y75[:, 5], alpha=0.5, color=color_plot[count])
        ax_critical.tick_params(axis='y')
        ax_critical.title.set_text('#Critical, Age{}'.format(j))
        ax_critical.legend(['Simulation'])

        # Dead
        ax_dead = ax[4, count]
        ax_dead.set_xlabel('time (days)')
        ax_dead.plot(x, y50[:, [7]], color=color_plot[count], label='y50')
        ax_dead.fill_between(x, y50[:, 7], y25[:, 7],
                             alpha=0.5, color=color_plot[count])
        ax_dead.fill_between(x, y50[:, 7], y75[:, 7],
                             alpha=0.5, color=color_plot[count])
        ax_dead.tick_params(axis='y')
        ax_dead.title.set_text('#Dead, Age{}'.format(j))
        ax_dead.legend(['Simulation'])

        # Recovered
        ax_dead = ax[5, count]
        ax_dead.set_xlabel('time (days)')
        ax_dead.plot(x, y50[:, [6]], color=color_plot[count], label='y50')
        ax_dead.fill_between(x, y50[:, 6], y25[:, 6],
                             alpha=0.5, color=color_plot[count])
        ax_dead.fill_between(x, y50[:, 6], y75[:, 6],
                             alpha=0.5, color=color_plot[count])
        ax_dead.tick_params(axis='y')
        ax_dead.title.set_text('#Recovered, Age{}'.format(j))
        ax_dead.legend(['Simulation'])

    # fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()


if __name__ == "__main__":
    path = ""

    if (len(sys.argv) > 1):
        n_runs = sys.argv[1]
    else:
        n_runs = len([entry for entry in os.listdir(path)
                     if os.path.isfile(os.path.join(path, entry))])
    
    plot_infection_states_results(path)
    plot_infections_loc_types_avarage(path)