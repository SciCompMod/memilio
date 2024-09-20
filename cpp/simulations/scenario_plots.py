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
from scipy.signal import savgol_filter

fontsize = 18

def plot_scenario(path, folder, folder_high, folder_enough):

    # four plots with 3 lines for each scenario and another one for the real data
    # we do one plot with 3 subplots first two plots in upper left and right are 
    # 1. daily new infections and 2. cumulative infections
    # the third plot in the lower left is the number of tests per day
    # the fourth plot in the lower right is the number of cumulative positive tests per day
    # the x-axis is the time for a date 

    # load the data
    locations = [0, 1, 2, 3, 4, 10]
    legend_plot = ['Home', 'School', 'Work',
                    'SocialEvent', 'BasicsShop','Event']
    
    # folders to iterate over
    folders = [folder, folder_high, folder_enough]
    # for folder in [folder_normal, folder_high, folder_enough]:
    # we plot both in separate plots
    fig, axs = plt.subplots(3, 2, figsize=(19, 12))
    fig.subplots_adjust(hspace=0.5)

    inf_p50_normal = h5py.File(
            path+"/infection_per_location_type_per_age_group/"+folder+"/p50/Results.h5", 'r')
    time = inf_p50_normal['0']['Time'][()] 
    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in time]
    xx = [xx[i].strftime('%m-%d') for i in range(len(xx))]
    xx = xx[1:][::24][:90]

    #label string for each folder
    labels = ['Standart Testing', 'High Testing Scenario', 'Very High Testing Scenario']

    for folder in folders:
        inf_p50_normal = h5py.File(
            path+"/infection_per_location_type_per_age_group/"+folder+"/p50/Results.h5", 'r')
        inf_p75_normal = h5py.File(
            path+"/infection_per_location_type_per_age_group/"+folder+"/p75/Results.h5", 'r')
        inf_p25_normal = h5py.File(
            path+"/infection_per_location_type_per_age_group/"+folder+"/p25/Results.h5", 'r')
        p50_bs_normal = inf_p50_normal['0']['Total'][()] 
        p25_bs_normal = inf_p25_normal['0']['Total'][()] 
        p75_bs_normal = inf_p75_normal['0']['Total'][()] 

        tests_p50_normal =  h5py.File(
            path+"/test_per_location_type_per_age_group/"+folder+"/p50/Results.h5", 'r')
        tests_p25_normal =  h5py.File(
            path+"/test_per_location_type_per_age_group/"+folder+"/p25/Results.h5", 'r')
        tests_p75_normal =  h5py.File(
            path+"/test_per_location_type_per_age_group/"+folder+"/p75/Results.h5", 'r')
        p50_bs_tests_normal = tests_p50_normal['0']['Total'][()]
        p25_bs_tests_normal = tests_p25_normal['0']['Total'][()]
        p75_bs_tests_normal = tests_p75_normal['0']['Total'][()]


        test_p_pos_p50_normal = h5py.File(
            path+"/new_detected_infections/"+folder+"/p50/Results.h5", 'r')
        test_p_pos_p25_normal = h5py.File(
            path+"/new_detected_infections/"+folder+"/p25/Results.h5", 'r')
        test_p_pos_p75_normal = h5py.File(
            path+"/new_detected_infections/"+folder+"/p75/Results.h5", 'r')
        p50_bs_test_p_pos_normal = test_p_pos_p50_normal['0']['Total'][()].flatten()
        p25_bs_test_p_pos_normal = test_p_pos_p25_normal['0']['Total'][()].flatten()
        p75_bs_test_p_pos_normal = test_p_pos_p75_normal['0']['Total'][()].flatten()

        deaths_p50_normal = h5py.File(
            path+"/infection_state_per_age_group/"+folder+"/p50/Results.h5", 'r')
        deaths_p25_normal = h5py.File(
            path+"/infection_state_per_age_group/"+folder+"/p25/Results.h5", 'r')
        deaths_p75_normal = h5py.File(
            path+"/infection_state_per_age_group/"+folder+"/p75/Results.h5", 'r')
        p50_bs_deaths_normal = deaths_p50_normal['0']['Total'][()][:, 7][::24][0:90]-deaths_p50_normal['0']['Total'][()][:, 7][::24][0]
        p25_bs_deaths_normal = deaths_p25_normal['0']['Total'][()][:, 7][::24][0:90]-deaths_p25_normal['0']['Total'][()][:, 7][::24][0]
        p75_bs_deaths_normal = deaths_p75_normal['0']['Total'][()][:, 7][::24][0:90]-deaths_p75_normal['0']['Total'][()][:, 7][::24][0]

        reproduction_p50_normal = h5py.File(
            path+"/estimated_reproduction_number/"+folder+"/p50/Results.h5", 'r')
        reproduction_p25_normal = h5py.File(
            path+"/estimated_reproduction_number/"+folder+"/p25/Results.h5", 'r')
        reproduction_p75_normal = h5py.File(
            path+"/estimated_reproduction_number/"+folder+"/p75/Results.h5", 'r')
        p50_bs_reproduction_normal = reproduction_p50_normal['0']['Total'][()][::24][0:90].flatten()
        p25_bs_reproduction_normal = reproduction_p25_normal['0']['Total'][()][::24][0:90].flatten()
        p75_bs_reproduction_normal = reproduction_p75_normal['0']['Total'][()][::24][0:90].flatten()




        # total_50_positive = np.sum(p50_bs_test_p_pos_normal, axis=1)
        total_50_positive = np.cumsum(p50_bs_test_p_pos_normal, axis=0)
        # total_50_positive = total_50_positive[::24]
        total_50_positive = total_50_positive[0:91] # we still need to take the difference to get the daily amount
        total_50_positive = np.diff(total_50_positive, axis=0).flatten()

        # total_25_positive = np.sum(p25_bs_test_p_pos_normal, axis=1)
        total_25_positive = np.cumsum(p25_bs_test_p_pos_normal, axis=0)
        # total_25_positive = total_25_positive[::24]
        total_25_positive = total_25_positive[0:91] # we still need to take the difference to get the daily amount
        total_25_positive = np.diff(total_25_positive, axis=0).flatten()

        # total_75_positive = np.sum(p75_bs_test_p_pos_normal, axis=1)
        total_75_positive = np.cumsum(p75_bs_test_p_pos_normal, axis=0)
        # total_75_positive = total_75_positive[::24]
        total_75_positive = total_75_positive[0:91] # we still need to take the difference to get the daily amount
        total_75_positive = np.diff(total_75_positive, axis=0).flatten()

        total_50_done = np.sum(p50_bs_tests_normal, axis=1)
        total_50_done = np.cumsum(total_50_done, axis=0)
        total_50_done = total_50_done[::24]
        total_50_done = total_50_done[0:91] # we still need to take the difference to get the daily amount
        total_50_done = np.diff(total_50_done, axis=0).flatten()
        # we smooth this with a gaussian filter

        total_25_done = np.sum(p25_bs_tests_normal, axis=1)
        total_25_done = np.cumsum(total_25_done, axis=0)
        total_25_done = total_25_done[::24]
        total_25_done = total_25_done[0:91] # we still need to take the difference to get the daily amount
        total_25_done = np.diff(total_25_done, axis=0).flatten()
        # we smooth this with a gaussian filter

        total_75_done = np.sum(p75_bs_tests_normal, axis=1)
        total_75_done = np.cumsum(total_75_done, axis=0)
        total_75_done = total_75_done[::24]
        total_75_done = total_75_done[0:91] # we still need to take the difference to get the daily amount
        total_75_done = np.diff(total_75_done, axis=0).flatten()
        # we smooth this with a gaussian filter


        p50_bs_normal_all_locations = np.zeros(len(p50_bs_normal))
        p25_bs_normal_all_locations = np.zeros(len(p25_bs_normal))
        p75_bs_normal_all_locations = np.zeros(len(p75_bs_normal))
        for location in locations:
            p50_bs_normal_all_locations += p50_bs_normal[:, location]
            p25_bs_normal_all_locations += p25_bs_normal[:, location]
            p75_bs_normal_all_locations += p75_bs_normal[:, location]

        
        cum_inf_normal_50 = np.cumsum(p50_bs_normal_all_locations)
        cum_inf_normal_50 = cum_inf_normal_50[::24]
        cum_inf_normal_50 = cum_inf_normal_50[:91]
        new_inf_normal_50 = np.diff(cum_inf_normal_50)

        cum_inf_normal_25 = np.cumsum(p25_bs_normal_all_locations)
        cum_inf_normal_25 = cum_inf_normal_25[::24]
        cum_inf_normal_25 = cum_inf_normal_25[:91]
        new_inf_normal_25 = np.diff(cum_inf_normal_25)

        cum_inf_normal_75 = np.cumsum(p75_bs_normal_all_locations)
        cum_inf_normal_75 = cum_inf_normal_75[::24]
        cum_inf_normal_75 = cum_inf_normal_75[:91]
        new_inf_normal_75 = np.diff(cum_inf_normal_75)



        # first plot
        axs[0, 0].plot(xx, new_inf_normal_50, label='Normal')
        axs[0, 0].fill_between(xx, new_inf_normal_25, new_inf_normal_75, alpha=0.5)
        
        
        # second plot
        axs[0, 1].plot(xx, cum_inf_normal_50[0:90], label=labels[folders.index(folder)])
        axs[0, 1].fill_between(xx, cum_inf_normal_25[0:90], cum_inf_normal_75[0:90], alpha=0.5)
        
        # third plot
        axs[1, 0].plot(xx, total_50_done, label='Normal')
        axs[1, 0].fill_between(xx, total_25_done, total_75_done, alpha=0.5)
        # make x axis start at 0 and end at maximum of the data
        axs[1, 0].set_ylim([0, max(total_75_done)*1.5])
        
        # fourth plot
        axs[1, 1].plot(xx[:89], total_50_positive, label='Normal')
        axs[1, 1].fill_between(xx[:89], total_25_positive, total_75_positive, alpha=0.5)

        # plot dead in last row 
        axs[2, 0].plot(xx, p50_bs_deaths_normal, label='Normal')
        axs[2, 0].fill_between(xx, p25_bs_deaths_normal, p75_bs_deaths_normal, alpha=0.5)

        # plot reproduction number
        axs[2, 1].plot(xx, p50_bs_reproduction_normal, label='Normal')
        axs[2, 1].fill_between(xx, p25_bs_reproduction_normal, p75_bs_reproduction_normal, alpha=0.5)
        axs[2, 1].set_ylim([0, max(p75_bs_reproduction_normal)*1.5])


    
  

    axs[0, 0].set_title('Daily New Infections', fontsize=fontsize+5)
    axs[0, 0].set_ylabel('Number of infections', fontsize=fontsize-2)
    axs[0, 0].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[0, 0].tick_params(axis='both', which='minor', labelsize=fontsize-4)

    axs[0, 1].set_title('Cumulative Infections', fontsize=fontsize+5)
    axs[0, 1].set_ylabel('Number of infections', fontsize=fontsize-2)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=fontsize-4)

    axs[1, 0].set_title('Daily Tests', fontsize=fontsize+5)
    axs[1, 0].set_ylabel('Number of tests', fontsize=fontsize-2)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=fontsize-4)

    axs[1, 1].set_title('Daily New Detected Infections', fontsize=fontsize+5)
    axs[1, 1].set_ylabel('Number of infections', fontsize=fontsize-2)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=fontsize-4)

    axs[2, 0].set_title('Cumulative Deaths', fontsize=fontsize+5)
    axs[2, 0].set_ylabel('Number of deaths' , fontsize=fontsize-2)
    axs[2, 0].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[2, 0].tick_params(axis='both', which='minor', labelsize=fontsize-4)

    axs[2, 1].set_title('Reproduction Number', fontsize=fontsize+5)
    axs[2, 1].set_ylabel('Reproduction number', fontsize=fontsize-2)
    axs[2, 1].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[2, 1].tick_params(axis='both', which='minor', labelsize=fontsize-4)

    # We still need one legend for all plots, we plot these in the first plot
    axs[0, 1].legend(loc='upper left', fontsize=fontsize-4)


    for ax in axs.flat:
        ax.set_xticks(xx[::10])
        ax.set_xticklabels(xx[::10])
        ax.tick_params(axis='x')
    plt.savefig('scenario_plot.jpg', dpi=300)





if __name__ == "__main__":
    path_to_scenario_1 = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/cluster_results/vorlaufige_ergebnisse/results_2024-08-28211441" # more tesring 
    path_to_scenario_2 = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/cluster_results/vorlaufige_ergebnisse/results_2024-08-28113737" # no lockdown but much testing
    path_to_scenario_3 = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/"
    path_to_main_data = "0"
    path_to_high_testing_data = "1"
    path_to_enough_testing_data = "2"

    # plot_scenario(path_to_scenario_1, path_to_main_data, path_to_high_testing_data, path_to_enough_testing_data)
    plot_scenario(path_to_scenario_2, path_to_main_data, path_to_high_testing_data, path_to_enough_testing_data)
    # plot_scenario(path_to_scenario_3, path_to_main_data, path_to_high_testing_data, path_to_enough_testing_data)
