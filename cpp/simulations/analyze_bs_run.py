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


def main(path, n_runs):
    # get first_file in folder
    # first_file = os.listdir(path)[0]
    # file_path = os.path.join(path, first_file)
    # read in txt file
    # df = pd.read_csv(file_path, delim_whitespace=True)
    # convert to numpy array
    # df_np = df.to_numpy()
    # get the number of rows and columns
    # num_rows = df_np.shape[0]
    # num_cols = df_np.shape[1]
    # get the number of compartments
    # num_compartments = num_cols - 1
    # get the number of time steps
    # num_time_steps = num_rows-1
    # get the compartment names
    # compartment_names = df.columns[1:]
    # get the time steps
    # time_steps = df_np[:, 0]

    # get number of files in folder
    # num_files = len([entry for entry in os.listdir(path)])
    # read in each txt file and convert to numpy array
    # df_np_3d = np.empty((num_rows, num_cols, n_runs))
    # print(os.listdir(path))
    for file in os.listdir(path):
        file_path = os.path.join(path, file)
        # read in txt file

        if file.startswith("infection_per_location_type.txt"):
            df = pd.read_csv(file_path, delim_whitespace=True)
            plot_infection_per_location_type(df)
        if file.startswith("infection_per_age_group.txt"):
            df = pd.read_csv(file_path, delim_whitespace=True)
            # plot_infection_per_age_group(df)
        # if file.startswith("run_"):
            # convert to numpy array
        #    df_np = df.to_numpy()
            # attach to array
        #    df_np_3d[:, :, i] = df_np
        #    plot_mean_and_std(df_np_3d)


def plot_infection_per_location_type(df):
    # Calculate moving average for all location types
    df['Home'] = gaussian_filter1d(df.Home.rolling(24*3, min_periods=1).sum(), sigma=15)
    df['School'] = gaussian_filter1d(df.School.rolling(24*3, min_periods=1).sum(), sigma=15)
    df['Work'] =  gaussian_filter1d(df.Work.rolling(24*3, min_periods=1).sum(), sigma=15)
    df['SocialEvent'] =  gaussian_filter1d(df.SocialEvent.rolling(24*3, min_periods=1).sum(), sigma=15)
    df['BasicsShop'] =  gaussian_filter1d(df.BasicsShop.rolling(24*3, min_periods=1).sum(), sigma=15)

    df.plot(x='Time', y=['Home', 'School','Work',  'SocialEvent', 'BasicsShop'], figsize=(10, 6), title='Infections stratified according to location type')
   
   # x axis should be titled Time in days
    plt.xlabel('Time (days)')
    plt.ylabel('New infections')
    



    # Subplots of individual location types
    # fig, axs = plt.subplots(5, 2, constrained_layout=True)
    # df.plot(x='Time', y='MA_Home', color='tab:blue', ax=axs[0, 0])
    # df.plot(x='Time', y='MA_Work', color='tab:green', ax=axs[0, 1])
    # df.plot(x='Time', y='MA_School', color='tab:red', ax=axs[1, 0])
    # df.plot(x='Time', y='MA_SocialEvent', color='tab:orange', ax=axs[1, 1])
    # df.plot(x='Time', y='MA_BasicsShop', color='tab:purple', ax=axs[2, 0])
    # df.plot(x='Time', y='MA_Hospital', color='tab:brown', ax=axs[2, 1])
    # df.plot(x='Time', y='MA_ICU', color='tab:pink', ax=axs[3, 0])
    # df.plot(x='Time', y='MA_Car', color='tab:gray', ax=axs[3, 1])
    # df.plot(x='Time', y='MA_PublicTransport', color='tab:olive', ax=axs[4, 0])
    # df.plot(x='Time', y='MA_Cemetery', color='tab:cyan', ax=axs[4, 1])

    plt.show()


def plot_infection_per_age_group(df):
    df.plot(x='Time', y=['0_to_4', '5_to_14', '15_to_34',
            '35_to_59', '60_to_79', '80_plus'], figsize=(10, 6))

    plt.show()


def plot_results(path):
    # median / 50-percentile
    f_p50 = h5py.File(
        path+"/infection_state_per_age_group/p50/Results.h5", 'r')

    # Get the HDF5 group; key needs to be a group name from above
    # only one graph node saved
    p50_bs = f_p50['0']

    # This assumes group[some_key_inside_the_group] is a dataset,
    # and returns a np.array:
    time = p50_bs['Time'][()]
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

    # real world
    f_real = h5py.File(
        path + "/Results_rki.h5", 'r')
    real_bs = f_real['3101']
    total_real = real_bs['Total'][()]
    

    plot_infection_states_individual(
        time, p50_bs, p25_bs, p75_bs, real_bs)
    plot_infection_states(time, total_50, total_25, total_75, total_real)

    # After you are done
    f_p50.close()
    f_p25.close()
    f_p75.close()
    f_real.close()

def plot_infection_states(x, y50, y25, y75, y_real=None):
    plt.figure('Infection_states')
    plt.title('Infection states')

    color_plot = cmx.get_cmap('Set1').colors

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

    color_plot = cmx.get_cmap('Set1').colors

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


def plot_infections_per_age_group(path):
    f = h5py.File(
        path+"/infection_per_age_group/p50/Results.h5", 'r')

    # Get the HDF5 group; key needs to be a group name from above
    group = f['0']

    # This assumes group[some_key_inside_the_group] is a dataset,
    # and returns a np.array:
    time = group['Time'][()]

    plt.figure('Age_Group')
    for g in ['Group' + str(n) for n in range(1, 7)]:
        gr = group[g][()]
        plt.plot(time, gr)

    plt.legend(['0-4', '5-14', '15-34', '35-59', '60-79', '80+'])

    # After you are done
    f.close()


def plot_mean_and_std(Y):

    x_plot = Y[:, 0, 0]
    compartments = Y[:, 1:, 1:]
    # average value
    compartments_avg = np.mean(compartments, axis=2)
    # plot average
    for i in range(compartments_avg.shape[1]):
        plt.plot(x_plot, compartments_avg[:, i])

    # plt.plot(x_plot,compartments_avg)
    # legend
    plt.legend(['S', 'E', 'I_NS', 'I_Sy', 'I_Sev', 'I_Crit', 'R', 'D'])
    plt.show()
    # standard deviation
    # compartments_std = np.std(compartments,axis=2)
    # plt.plot(x_plot,compartments_avg + compartments_std)
    # plt.plot(x_plot,compartments_avg - compartments_std)
    # plt.show()


if __name__ == "__main__":
    # path to results
    # path = "/Users/david/Documents/HZI/memilio/data/results"
    # path = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/results/cluster/results"
    path = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/results_last_run"
    # path = r"C:\Users\korf_sa\Documents\rep\data\results"
    # path = r"C:\Users\korf_sa\Documents\rep\data\results_cluster\results"
    if (len(sys.argv) > 1):
        n_runs = sys.argv[1]
    else:
        n_runs = len([entry for entry in os.listdir(path)
                     if os.path.isfile(os.path.join(path, entry))])
    plot_results(path)
    main(path, n_runs)
