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
from scipy.ndimage import gaussian_filter1d


def plot_infections_loc_types_avarage(path):
    # 50-percentile
    f_p50 = h5py.File(
        path+"/infection_per_location_type/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]

    # 25-percentile
    f_p25 = h5py.File(
        path+"/infection_per_location_type/p25/Results.h5", 'r')
    p25_bs = f_p25['0']
    total_25 = p25_bs['Total'][()]

    # 75-percentile
    f_p75 = h5py.File(
        path + "/infection_per_location_type/p75/Results.h5", 'r')
    p75_bs = f_p75['0']
    total_75 = p75_bs['Total'][()]

    time = p50_bs['Time'][()]

    plot_infection_per_location_type_mean(
        time, total_50, total_25, total_75)
    
def plot_infection_per_location_type_mean(x, y50, y25, y75):

    plt.figure('Infection_states_location_types')
    plt.title('Infection states per location types avaraged over all runs')

    color_plot = matplotlib.colormaps.get_cmap('Set1').colors

    states_plot = [0,1,2,3,4]
    legend_plot = ['Home', 'School', 'Work', 'SocialEvent', 'BasicsShop', 'Hospital', 'ICU']

    for i in states_plot:
        #rolling average
        plt.plot(x, gaussian_filter1d(pd.DataFrame(y50[:, i]).rolling(24*3, min_periods=1).sum(), sigma=15), color=color_plot[i])

    plt.legend(legend_plot)

    for i in states_plot:
        y50_smoothed= gaussian_filter1d(pd.DataFrame(y50[:, i]).rolling(24*3, min_periods=1).sum(), sigma=15).flatten()
        y25_smoothed= gaussian_filter1d(pd.DataFrame(y25[:, i]).rolling(24*3, min_periods=1).sum(), sigma=15).flatten()
        y75_smoothed= gaussian_filter1d(pd.DataFrame(y75[:, i]).rolling(24*3, min_periods=1).sum(), sigma=15).flatten()
        plt.fill_between(x, y50_smoothed,  y25_smoothed ,
                         alpha=0.5, color=color_plot[i])
        plt.fill_between(x, y50_smoothed,  y75_smoothed,
                         alpha=0.5, color=color_plot[i])
        
    #currently the x axis has the values of the time steps, we need to convert them to dates and set the x axis to dates
    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in x]
    xx = [xx[i].strftime('%Y-%m-%d') for i in range(len(xx))]
    #but just take every 10th date to make it more readable
    plt.gca().set_xticks(x[::150])
    plt.gca().set_xticklabels(xx[::150])
    plt.gcf().autofmt_xdate()


    plt.xlabel('Time')
    plt.ylabel('Number of individuals')
    plt.show()

def plot_infectoin_states_results(path):
    # 50-percentile
    f_p50 = h5py.File(
        path+"/infection_state_per_age_group/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]
    # 25-percentile
    f_p25 = h5py.File(
        path+"/infection_state_per_age_group/p25/Results.h5", 'r')
    p25_bs = f_p25['0']
    total_25 = p25_bs['Total'][()]
    # 75-percentile
    f_p75 = h5py.File(
        path + "/infection_state_per_age_group/p75/Results.h5", 'r')
    p75_bs = f_p75['0']
    total_75 = p75_bs['Total'][()]

    time = p50_bs['Time'][()]

    # real world
    f_real = h5py.File(
        path + "/Results_rki.h5", 'r')
    real_bs = f_real['3101']
    total_real = real_bs['Total'][()]
    

    plot_infection_states_individual(
        time, p50_bs, p25_bs, p75_bs, real_bs)
    plot_infection_states(time, total_50, total_25, total_75, total_real)

def plot_infection_states(x, y50, y25, y75, y_real=None):
    plt.figure('Infection_states')
    plt.title('Infection states')

    color_plot = matplotlib.colormaps.get_cmap('Set1').colors

    states_plot = [1, 2, 3, 4, 5, 7]
    legend_plot = ['E', 'I_NSymp', 'I_Symp', 'I_Sev', 'I_Crit', 'Dead', 'Sm. rep. Sympt.']

    for i in states_plot:
        plt.plot(x, y50[:, i], color=color_plot[i])

        #plot real data
    if y_real is not None:
        x_real = np.linspace(0, y_real.shape[0]-1, y_real.shape[0])
        plt.plot(x_real, y_real[:,4], '.', color='tab:red')

    plt.legend(legend_plot)

   
    # plt.show()

    

    for i in states_plot:
        plt.fill_between(x, y50[:, i], y25[:, i],
                         alpha=0.5, color=color_plot[i])
        plt.fill_between(x, y50[:, i], y75[:, i],
                         alpha=0.5, color=color_plot[i])
        
       
    #currently the x axis has the values of the time steps, we need to convert them to dates and set the x axis to dates
    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in x]
    xx = [xx[i].strftime('%Y-%m-%d') for i in range(len(xx))]
    #but just take every 10th date to make it more readable
    plt.gca().set_xticks(x[::150])
    plt.gca().set_xticklabels(xx[::150])
    plt.gcf().autofmt_xdate()


    plt.xlabel('Time')
    plt.ylabel('Number of individuals')
    plt.show()

def plot_infection_states_individual(x, p50_bs, p25_bs, p75_bs, real_bs):

    y_real_total = real_bs['Total'][()]
    x_real = np.linspace(0, y_real_total.shape[0]-1, y_real_total.shape[0])

    age_group_access = ['Group1', 'Group2', 'Group3', 'Group4', 'Group5', 'Group6', 'Total']

    color_plot = matplotlib.colormaps.get_cmap('Set1').colors

    fig, ax = plt.subplots(3, len(age_group_access), constrained_layout = True)
    fig.set_figwidth(20)
    fig.set_figheight(9)
    for j, count in zip(age_group_access, range(len(age_group_access))):
        y50 = p50_bs[j][()]
        y25 = p25_bs[j][()]
        y75 = p75_bs[j][()]
        y_real = real_bs[j][()]
        
        # Severe
        ax_severe = ax[0, count]
        ax_severe.set_xlabel('time (days)')
        ax_severe.plot(x, y50[:, 4], color=color_plot[count], label='y50')
        ax_severe.plot(x_real, y_real[:, 6], '.', color=color_plot[count], label='y_real')
        ax_severe.fill_between(x, y50[:, 4], y25[:, 4], alpha=0.5, color=color_plot[count])
        ax_severe.fill_between(x, y50[:, 4], y75[:, 4], alpha=0.5, color=color_plot[count])
        ax_severe.tick_params(axis='y')
        ax_severe.title.set_text('Severe, Age{}'.format(j))
        ax_severe.legend(['Simulation','Real'])

        # Critical
        ax_critical = ax[1, count]
        ax_critical.set_xlabel('time (days)')
        ax_critical.plot(x, y50[:, [5]], color=color_plot[count], label='y50')
        ax_critical.plot(x_real, y_real[:, [7]], '.', color=color_plot[count], label='y_real')
        ax_critical.fill_between(x, y50[:, 5], y25[:, 5], alpha=0.5, color=color_plot[count])
        ax_critical.fill_between(x, y50[:, 5], y75[:, 5], alpha=0.5, color=color_plot[count])
        ax_critical.tick_params(axis='y')
        ax_critical.title.set_text('Critical, Age{}'.format(j))
        ax_critical.legend(['Simulation','Real'])

        # Dead
        ax_dead = ax[2, count]
        ax_dead.set_xlabel('time (days)')
        ax_dead.plot(x, y50[:, [7]], color=color_plot[count], label='y50')
        ax_dead.plot(x_real, y_real[:, [9]], '.', color=color_plot[count], label='y_real')
        ax_dead.fill_between(x, y50[:, 7], y25[:, 7], alpha=0.5, color=color_plot[count])
        ax_dead.fill_between(x, y50[:, 7], y75[:, 7], alpha=0.5, color=color_plot[count])
        ax_dead.tick_params(axis='y')
        ax_dead.title.set_text('Dead, Age{}'.format(j))
        ax_dead.legend(['Simulation','Real'])
    
    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()




if __name__ == "__main__":
    # path = "/Users/david/Documents/HZI/memilio/data/results_last_run"
    # path = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/results_last_run"
    path = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/cluster_results/1/results_last_run"
    # path = r"C:\Users\korf_sa\Documents\rep\data\results_last_run"

    if (len(sys.argv) > 1):
        n_runs = sys.argv[1]
    else:
        n_runs = len([entry for entry in os.listdir(path)
                     if os.path.isfile(os.path.join(path, entry))])
    plot_infectoin_states_results(path)
    plot_infections_loc_types_avarage(path)
