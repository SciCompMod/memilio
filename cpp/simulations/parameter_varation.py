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

def plot_of_cumuative_infections(path, index=0):
    # plot cumulative infections
    f_p50 = h5py.File(
        path+"/infection_per_location_type_per_age_group/"+str(index)+"/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]

    time = p50_bs['Time'][()]
    time = time[::24]
    time = time[0:90]


    states_plot = [0, 1, 2, 3, 4, 10]
    np_y50_total = np.zeros(90)
    for i in states_plot:
        # we need to sum up every 24 hours
        indexer = pd.api.indexers.FixedForwardWindowIndexer(window_size=24)
        np_y50 = pd.DataFrame(total_50[:, i]).rolling(window=indexer, min_periods=1).sum().to_numpy()
        np_y50=np_y50[0::24].flatten()[0:90]
        np_y50_total += np_y50
    return np_y50_total

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

if __name__ == "__main__":
    
    path_data = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/results_2024-08-15140128"
    variable_1 = "test_likelihood_symptom"
    variable_2 = "test_likelihood_asymptom"
    values_1 = [0.01, 0.03, 0.05]
    values_2 = [5, 10, 30]

    # we want to have a Grid Plot with the size of the values
    # in each grid we plot a certain plot which we can define
    # on the general structure of the plot we we write the values of the parameters
    # we want to plot in the title of the plot 
    # they all should have the same y-axis, so we can compare them, the highest value should be the maximum of all plots
    
    fig = plt.figure('Parameter Variation', figsize=(1, 1))
    gs = fig.add_gridspec(len(values_1), len(values_2), hspace=0.4, wspace=0.4)
    fig.suptitle('Parameter Variation', fontsize=16)
    fig.show()
    highest_value = 0
    for i in range(len(values_1)):
        for j in range(len(values_2)):
            index = (i)*len(values_2)+(j)
            plot_data = plot_of_cumuative_infections(path_data, index)
            if np.max(plot_data) > highest_value:
                highest_value = np.max(plot_data)

    for i, value_1 in enumerate(values_1):
        for j, value_2 in enumerate(values_2):
            index = (i)*len(values_2)+(j)
            ax = fig.add_subplot(gs[i, j])
            ax.set_title(f'{variable_1}={value_1}, {variable_2}={value_2}')
            plot_data = plot_of_cumuative_infections(path_data, index)
            ax.set_ylim(0, highest_value)
            ax.plot(plot_data, label='Cumulative Infections')
            ax.legend()
    plt.show()


    fig = plt.figure('Parameter Variation', figsize=(1, 1))
    gs = fig.add_gridspec(len(values_1), len(values_2), hspace=0.4, wspace=0.4)
    fig.suptitle('Parameter Variation', fontsize=16)
    fig.show()
    highest_value = 0
    for i in range(len(values_1)):
        for j in range(len(values_2)):
            index = (i)*len(values_2)+(j)
            plot_data = plot_number_of_tests(path_data, index)
            if np.max(plot_data) > highest_value:
                highest_value = np.max(plot_data)
    for i, value_1 in enumerate(values_1):
        for j, value_2 in enumerate(values_2):
            index = (i)*len(values_2)+(j)
            
            ax = fig.add_subplot(gs[i, j])
            ax.set_title(f'{variable_1}={value_1}, {variable_2}={value_2}')
            ax.set_ylim(0, highest_value)
            plot_data = plot_number_of_tests(path_data, index)
            ax.plot(plot_data, label='Number of Tests')
            ax.legend()
    plt.show()
    fig = plt.figure('Parameter Variation', figsize=(1, 1))
    gs = fig.add_gridspec(len(values_1), len(values_2), hspace=0.4, wspace=0.4)
    fig.suptitle('Parameter Variation', fontsize=16)
    fig.show()
    highest_value = 0
    for i in range(len(values_1)):
        for j in range(len(values_2)):
            index = (i)*len(values_2)+(j)
            plot_data = plot_positive_tests(path_data, index)
            if np.max(plot_data) > highest_value:
                highest_value = np.max(plot_data)
    for i, value_1 in enumerate(values_1):
        for j, value_2 in enumerate(values_2):
            index = (i)*len(values_2)+(j)
            
            ax = fig.add_subplot(gs[i, j])
            ax.set_title(f'{variable_1}={value_1}, {variable_2}={value_2}')
            ax.set_ylim(0, highest_value)
            plot_data = plot_positive_tests(path_data, index)
            ax.plot(plot_data, label='Positive Tests')
            ax.legend()

    plt.show()



    
