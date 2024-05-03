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
    df['MA_Home'] = df.Home.rolling(10, min_periods=1).mean()
    df['MA_Work'] = df.Work.rolling(10, min_periods=1).mean()
    df['MA_School'] = df.School.rolling(10, min_periods=1).mean()
    df['MA_SocialEvent'] = df.SocialEvent.rolling(10, min_periods=1).mean()
    df['MA_BasicsShop'] = df.BasicsShop.rolling(10, min_periods=1).mean()
    df['MA_Hospital'] = df.Hospital.rolling(10, min_periods=1).mean()
    df['MA_ICU'] = df.ICU.rolling(10, min_periods=1).mean()
    df['MA_Car'] = df.Car.rolling(10, min_periods=1).mean()
    df['MA_PublicTransport'] = df.PublicTransport.rolling(
        10, min_periods=1).mean()
    df['MA_Cemetery'] = df.Cemetery.rolling(10, min_periods=1).mean()

    df.plot(x='Time', y=['MA_Home', 'MA_Work', 'MA_School', 'MA_SocialEvent', 'MA_BasicsShop',
            'MA_Hospital', 'MA_ICU', 'MA_Car', 'MA_PublicTransport', 'MA_Cemetery'], figsize=(10, 6))

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
    f = h5py.File(
        path+"/infection_state_per_age_group/p50/Results.h5", 'r')

    # Get the HDF5 group; key needs to be a group name from above
    group = f['0']

    # This assumes group[some_key_inside_the_group] is a dataset,
    # and returns a np.array:
    time = group['Time'][()]
    total_50 = group['Total'][()]

    # After you are done
    f.close()

    # 05-percentile
    f = h5py.File(
        path+"/infection_state_per_age_group/p25/Results.h5", 'r')
    group = f['0']
    total_25 = group['Total'][()]
    f.close()

    # 95-percentile
    f = h5py.File(
        path + "/infection_state_per_age_group/p75/Results.h5", 'r')
    group = f['0']
    total_75 = group['Total'][()]
    f.close()

    # real world
    # TODO
    f = h5py.File(
        path + "/Results_rki.h5", 'r')
    group = f['3101']
    total_real = group['Total'][()]
    f.close()

    plot_infection_states_individual(
        time, total_50, total_25, total_75, total_real)
    plot_infection_states(time, total_50, total_25, total_75)


def plot_infection_states(x, y50, y25, y75):
    plt.figure('Infection_states')
    plt.title('Infection states')
    color_plot = cmx.get_cmap('Set1').colors

    states_plot = [1, 2, 3, 4, 5, 7]
    legend_plot = ['E', 'I_NS', 'I_S', 'I_Sev', 'I_Crit', 'Dead']

    for i in states_plot:
        plt.plot(x, y50[:, i], color=color_plot[i])

    plt.legend(legend_plot)

    for i in states_plot:
        plt.fill_between(x, y50[:, i], y25[:, i],
                         alpha=0.5, color=color_plot[i])
        plt.fill_between(x, y50[:, i], y75[:, i],
                         alpha=0.5, color=color_plot[i])


def plot_infection_states_individual(x, y50, y25, y75, y_real):
    # plt.figure('Infection_states_dead')
    # # plt.plot(x, y50[:,[5,7]])
    x_real = np.linspace(0, y_real.shape[0]-1, y_real.shape[0])
    # plt.plot(x_real, y_real[:,[2]])
    # plt.legend(['I_Crit', 'Dead'])

    fig, ax = plt.subplots(3, 1)

    # Severe
    ax[0].set_xlabel('time (s)')
    ax[0].plot(x, y50[:, [4]], color='tab:red', label='y50')
    ax[0].plot(x_real, y_real[:, [6]], color='tab:blue', label='y_real')
    ax[0].tick_params(axis='y')
    ax[0].title.set_text('Severe')
    ax[0].legend()
    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    # Critical
    ax[1].set_xlabel('time (s)')
    ax[1].plot(x, y50[:, [5]], color='tab:red', label='y50')
    ax[1].plot(x_real, y_real[:, [7]], color='tab:blue', label='y_real')
    ax[1].tick_params(axis='y')
    ax[1].title.set_text('Critical')
    ax[1].legend()
    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    # Dead
    ax[2].set_xlabel('time (s)')
    ax[2].plot(x, y50[:, [7]], color='tab:red', label='y50')
    ax[2].plot(x_real, y_real[:, [9]], color='tab:blue', label='y_real')
    ax[2].tick_params(axis='y')
    ax[2].title.set_text('Dead')
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    ax[2].legend()
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
    path = "./results"
    if (len(sys.argv) > 1):
        n_runs = sys.argv[1]
    else:
        n_runs = len([entry for entry in os.listdir(path)
                     if os.path.isfile(os.path.join(path, entry))])
    plot_results(path)
    main(path, n_runs)
