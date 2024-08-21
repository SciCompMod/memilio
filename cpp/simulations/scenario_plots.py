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

def plot_cumulative_infections(path, folder_normal, folder_high, folder_enough):

    # four plots with 3 lines for each scenario and another one for the real data
    # we do one plot with 3 subplots first two plots in upper left and right are 
    # 1. daily new infections and 2. cumulative infections
    # the third plot in the lower left is the number of tests per day
    # the fourth plot in the lower right is the number of positive tests per day

    





if __name__ == "__main__":
    path_to_data = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/cluster_results/6/"
    path_to_main_data = "results_2024-08-21084018"
    path_to_high_testing_data = "results_2024-08-21084018"
    path_to_enough_testing_data = "results_2024-08-21084018"

    if (len(sys.argv) > 1):
        n_runs = sys.argv[1]
    else:
        n_runs = len([entry for entry in os.listdir(path)
                     if os.path.isfile(os.path.join(path, entry))])
    plot(path)
