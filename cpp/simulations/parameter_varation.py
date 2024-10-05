import sys
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import colormaps as cmaps
import h5py
from datetime import datetime
from matplotlib.dates import DateFormatter
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
import seaborn as sns


fontsize = 20

def plot_of_cumuative_infections(path, index=0):
    # plot cumulative infections
    f_p50 = h5py.File(
        path+"/infection_per_location_type_per_age_group/"+str(index)+"/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]
    # we sum over all locations
    total_50 = np.sum(total_50, axis=1)
    total_50 = np.cumsum(total_50, axis=0)
    total_50 = total_50[::24]
    total_50 = total_50[0:91] # we still need to take the difference to get the daily amount
    # we smooth this with a gaussian filter
    return total_50

def plot_number_of_tests(path, index):
    f_p50_done = h5py.File(
        path+"/test_per_location_type_per_age_group/"+str(index)+"/p50/Results.h5", 'r')
    p50_bs_done = f_p50_done['0']
    total_50_done = p50_bs_done['Total'][()]
    #same for the done tests
    total_50_done = np.sum(total_50_done, axis=1)
    total_50_done = np.cumsum(total_50_done, axis=0)
    total_50_done = total_50_done[::24]
    total_50_done = total_50_done[0:91] # we still need to take the difference to get the daily amount
    total_50_done = np.diff(total_50_done, axis=0).flatten()
    # we smooth this with a gaussian filter
    total_50_done = gaussian_filter1d(total_50_done, sigma=1, mode='nearest')
    return total_50_done
     
def plot_positive_tests(path, index):
    f_p50_done = h5py.File(
        path+"/positive_test_per_location_type_per_age_group/"+str(index)+"/p50/Results.h5", 'r')
    p50_bs_done = f_p50_done['0']
    total_50_done = p50_bs_done['Total'][()]
    #same for the done tests
    total_50_done = np.sum(total_50_done, axis=1)
    total_50_done = np.cumsum(total_50_done, axis=0)
    total_50_done = total_50_done[::24]
    total_50_done = total_50_done[0:91] # we still need to take the difference to get the daily amount
    total_50_done = np.diff(total_50_done, axis=0).flatten()
    # we smooth this with a gaussian filter
    total_50_done = gaussian_filter1d(total_50_done, sigma=1, mode='nearest')
    return total_50_done

def get_maximum_dead(path, index):
    deaths_p50_normal = h5py.File(
            path+"/infection_state_per_age_group/"+str(index)+"/p50/Results.h5", 'r')
    p50_bs_deaths_normal = deaths_p50_normal['0']['Total'][()][:, 7][::24][0:90]-deaths_p50_normal['0']['Total'][()][:, 7][::24][0]
    return np.max(p50_bs_deaths_normal)

def get_maximum_hospitalized(path, index):
    hospitalized_p50_normal = h5py.File(
            path+"/infection_state_per_age_group/"+str(index)+"/p50/Results.h5", 'r')
    p50_bs_hospitalized_normal = hospitalized_p50_normal['0']['Total'][()][:, 5][::24][0:90]
    return np.max(p50_bs_hospitalized_normal)

def get_maximum_cum_infected(path, index):
    infected_p50_normal = h5py.File(
            path+"/infection_per_location_type_per_age_group/"+str(index)+"/p50/Results.h5", 'r')
    locations = [0, 1, 2, 3, 4, 10]
    p50_bs_infected_normal = infected_p50_normal['0']['Total'][()]
    p50_bs_normal_all_locations = np.zeros(len(p50_bs_infected_normal))
    for location in locations:
        p50_bs_normal_all_locations += p50_bs_infected_normal[:, location]
    cum_inf_normal_50 = np.cumsum(p50_bs_normal_all_locations)
    return np.max(cum_inf_normal_50)

def get_maximum_daily_infections(path, index):
    infected_p50_normal = h5py.File(
            path+"/infection_per_location_type_per_age_group/"+str(index)+"/p50/Results.h5", 'r')
    locations = [0, 1, 2, 3, 4, 10]
    p50_bs_infected_normal = infected_p50_normal['0']['Total'][()]
    p50_bs_normal_all_locations = np.zeros(len(p50_bs_infected_normal))
    for location in locations:
        p50_bs_normal_all_locations += p50_bs_infected_normal[:, location]
    cum_inf_normal_50 = np.cumsum(p50_bs_normal_all_locations)
    daily_inf = np.diff(cum_inf_normal_50)
    return np.max(daily_inf)


if __name__ == "__main__":
    
    # path_data = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/cluster_results/vorlaufige_ergebnisse/results_2024-08-28120457"
    # path_data = r"C:\Users\korf_sa\Documents\memilio\data\cluster_results\final_results\results_2024-09-20233035"
    # variable_1 = "Test likelihood with symptoms"
    # variable_2 = "Ratio for asymptomatic agents to test"
    # values_1 = [0.012, 0.024, 0.036, 0.048, 0.06]
    # values_2 = ["1/2","1/5","1/8","1/11","1/14"]

    # path_data = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/cluster_results/vorlaufige_ergebnisse/results_2024-08-28120954"
    # path_data = r"C:\Users\korf_sa\Documents\memilio\data\cluster_results\final_results\results_2024-09-20233103"
    # variable_1 = "Quarantine length (days)"
    # variable_2 = "Quarantine efficiency"
    # values_1 = [2,5,8,11,14]
    # values_2 = [0.0,0.25,0.5,0.75,1.0]

    # path_data = r"C:\Users\korf_sa\Documents\memilio\data\cluster_results\final_results\results_2024-09-29002004"
    # variable_1 = "Test likelihood with symptoms"
    # variable_2 = "Quarantine length (days)"
    # values_1 = [0.012, 0.024, 0.036, 0.048, 0.06]
    # values_2 = [2,5,8,11,14]

    path_data = r"C:\Users\korf_sa\Documents\memilio\data\cluster_results\final_results\results_2024-10-01225821"
    variable_1 = "Test likelihood with symptoms"
    variable_2 = "Mask Efficiency"
    values_1 = [0.012, 0.024, 0.036, 0.048, 0.06]
    values_2 = [0.2, 0.225, 0.25, 0.275, 0.3]

    # we want to have a Grid Plot with the size of the values
    # in each grid we plot a certain plot which we can define
    # on the general structure of the plot we we write the values of the parameters
    # we want to plot in the title of the plot 
    # they all should have the same y-axis, so we can compare them, the highest value should be the maximum of all plots
    
    # fig = plt.figure('Parameter Variation', figsize=(19, 10))
    # gs = fig.add_gridspec(len(values_1), len(values_2), hspace=0.4, wspace=0.4)
    # # fig.suptitle('Parameter Variation', fontsize=fontsize)
    # fig.show()
    # highest_value = 0
    # for i in range(len(values_1)):
    #     for j in range(len(values_2)):
    #         index = (i)*len(values_2)+(j)
    #         plot_data = plot_of_cumuative_infections(path_data, index)
    #         if np.max(plot_data) > highest_value:
    #             highest_value = np.max(plot_data)

    # for i, value_1 in enumerate(values_1):
    #     for j, value_2 in enumerate(values_2):
    #         index = (i)*len(values_2)+(j)
    #         ax = fig.add_subplot(gs[i, j])
    #         ax.set_title(f'{variable_1}={value_1}, {variable_2}={value_2}')
    #         plot_data = plot_of_cumuative_infections(path_data, index)
    #         ax.set_ylim(0, highest_value)
    #         ax.plot(plot_data, label='Cumulative Infections')
    #         ax.legend()
    # plt.show()


    # fig = plt.figure('Parameter Variation', figsize=(19, 10))
    # gs = fig.add_gridspec(len(values_1), len(values_2), hspace=0.4, wspace=0.4)
    # fig.suptitle('Parameter Variation', fontsize=16)
    # fig.show()
    # highest_value = 0
    # for i in range(len(values_1)):
    #     for j in range(len(values_2)):
    #         index = (i)*len(values_2)+(j)
    #         plot_data = plot_number_of_tests(path_data, index)
    #         if np.max(plot_data) > highest_value:
    #             highest_value = np.max(plot_data)
    # for i, value_1 in enumerate(values_1):
    #     for j, value_2 in enumerate(values_2):
    #         index = (i)*len(values_2)+(j)   
    #         ax = fig.add_subplot(gs[i, j])
    #         ax.set_title(f'{variable_1}={value_1}, {variable_2}={value_2}')
    #         ax.set_ylim(0, highest_value)
    #         plot_data = plot_number_of_tests(path_data, index)
    #         ax.plot(plot_data, label='Number of Tests')
    #         ax.legend()
    # plt.show()
    # fig = plt.figure('Parameter Variation', figsize=(19, 10))
    # gs = fig.add_gridspec(len(values_1), len(values_2), hspace=0.4, wspace=0.4)
    # fig.suptitle('Parameter Variation', fontsize=16)
    # fig.show()
    # highest_value = 0
    # for i in range(len(values_1)):
    #     for j in range(len(values_2)):
    #         index = (i)*len(values_2)+(j)
    #         plot_data = plot_positive_tests(path_data, index)
    #         if np.max(plot_data) > highest_value:
    #             highest_value = np.max(plot_data)
    # for i, value_1 in enumerate(values_1):
    #     for j, value_2 in enumerate(values_2):
    #         index = (i)*len(values_2)+(j)
            
    #         ax = fig.add_subplot(gs[i, j])
    #         ax.set_title(f'{variable_1}={value_1}, {variable_2}={value_2}')
    #         ax.set_ylim(0, highest_value)
    #         plot_data = plot_positive_tests(path_data, index)
    #         ax.plot(plot_data, label='Positive Tests')
    #         ax.legend()
    # plt.show()

    # we want to have four plots in one figure
    # in each plot for each value of the parameter we want to plot a big square in a color which is defined by the value of the interesting value in that plot
    # first plot that value is the maximum of the cumulative infections
    # second plot that value is the maximum of the daily infections
    # third plot that value is the maximum of persons dead
    # fourth plot that value is the maximum of persons in the hospital

    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    # space between the upper and lower plots
    plt.subplots_adjust(hspace=0.3, wspace=0.2)
    fig.show()
    sns.set(font_scale=1.2)
    # first plot
    value_for_values = np.zeros((len(values_1), len(values_2)))
    for i, value_1 in enumerate(values_1):
        for j, value_2 in enumerate(values_2):
            index = (i)*len(values_2)+(j)
            value_for_values[i, j] = get_maximum_cum_infected(path_data, index)
    # permute the values, so that the plot is correct
    # we need a seaborn heatmap for the colorbar

    sns.heatmap(value_for_values, ax=axs[0, 0], cmap='viridis', annot=True, fmt=".0f")
    axs[0, 0].set_title('Cumulative Infections', fontsize=fontsize+4)
    # x and y axis labels
    # axs[0, 0].invert_xaxis()
    axs[0, 0].set_xticklabels(values_2)
    axs[0, 0].set_yticklabels(values_1)
    axs[0, 0].set_ylabel(variable_1, fontsize=fontsize-4)
    axs[0, 0].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[0, 0].tick_params(axis='both', which='minor', labelsize=fontsize-4)
    



    
    # second plot
    value_for_values = np.zeros((len(values_1), len(values_2)))
    for i, value_1 in enumerate(values_1):
        for j, value_2 in enumerate(values_2):
            index = (i)*len(values_2)+(j)
            value_for_values[i, j] = get_maximum_dead(path_data, index)
    # permute the values, so that the plot is correct
    # we need a seaborn heatmap for the colorbar
    sns.heatmap(value_for_values, ax=axs[0, 1], cmap='viridis', annot=True, fmt=".0f")
    # invert x axis
    # axs[0, 1].invert_xaxis()
    axs[0, 1].set_title('Deaths', fontsize=fontsize+4)
    # x and y axis labels
    axs[0, 1].set_xticklabels(values_2)
    axs[0, 1].set_yticklabels(values_1)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=fontsize-4)

    
    

    # third plot
    value_for_values = np.zeros((len(values_1), len(values_2)))
    for i, value_1 in enumerate(values_1):
        for j, value_2 in enumerate(values_2):
            index = (i)*len(values_2)+(j)
            value_for_values[i, j] = get_maximum_hospitalized(path_data, index)
    # permute the values, so that the plot is correct
    # we need a seaborn heatmap for the colorbar
    sns.heatmap(value_for_values, ax=axs[1, 0], cmap='viridis', annot=True, fmt=".0f")
    axs[1, 0].set_title('Maximum Hospitalized Persons', fontsize=fontsize+4)
    # x and y axis labels
    axs[1, 0].set_xticklabels(values_2)
    axs[1, 0].set_yticklabels(values_1)
    axs[1, 0].set_ylabel(variable_1, fontsize=fontsize-4)
    axs[1, 0].set_xlabel(variable_2, fontsize=fontsize-4)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=fontsize-4)
    
    # axs[1, 0].invert_xaxis()


    # fourth plot
    value_for_values = np.zeros((len(values_1), len(values_2)))
    for i, value_1 in enumerate(values_1):
        for j, value_2 in enumerate(values_2):
            index = (i)*len(values_2)+(j)
            value_for_values[i, j] = get_maximum_daily_infections(path_data, index)
    # permute the values, so that the plot is correct
    # we need a seaborn heatmap for the colorbar
    sns.heatmap(value_for_values, ax=axs[1, 1], cmap='viridis', annot=True, fmt=".0f")
    axs[1, 1].set_title('Maximum Daily Infections', fontsize=fontsize+4)
    # x and y axis labels
    axs[1, 1].set_xticklabels(values_2)
    axs[1, 1].set_yticklabels(values_1)
    axs[1, 1].set_xlabel(variable_2, fontsize=fontsize-4)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=fontsize-4)

    
    # axs[1, 1].invert_xaxis()



    plt.savefig('parameter_variation.png', dpi=300)


        
