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
from scipy.signal import savgol_filter

fontsize = 20


def plot_dead(path):
    # we will have a seperate plot the cumulative infected individuals, cumulative symptomatic individuals and cumulative dead individual
    # we need to load the data
    f_p50 = h5py.File(
        path+"/infection_state_per_age_group/0/p50/Results.h5", 'r')
    p50_bs = f_p50['0']

    # do the same for 25 and 75 percentile
    f_p25 = h5py.File(
        path+"/infection_state_per_age_group/0/p25/Results.h5", 'r')
    p25_bs = f_p25['0']

    f_p75 = h5py.File(
        path+"/infection_state_per_age_group/0/p75/Results.h5", 'r')
    p75_bs = f_p75['0']

    # do the same for 05 and 95 percentile
    f_p05 = h5py.File(
        path+"/infection_state_per_age_group/0/p05/Results.h5", 'r')
    p05_bs = f_p05['0']

    f_p95 = h5py.File(
        path+"/infection_state_per_age_group/0/p95/Results.h5", 'r')
    p95_bs = f_p95['0']

    age_group_access = ['Group1', 'Group2', 'Group3',
                        'Group4', 'Group5', 'Group6', 'Total']

    # we need the real data json file cases_all_county_age
    df_abb = pd.read_json(
        path+"/../../../pydata/Germany/cases_all_county_age_ma1.json")

    # we just need the columns cases and date
    # we need to offset the dates by 19 day
    df_abb['Date'] = df_abb['Date'] + pd.DateOffset(days=18)
    # we need just the dates bewteen 2021-03-01 and 2021-06-01
    df_abb = df_abb[(df_abb['Date'] >= '2021-03-01') &
                    (df_abb['Date'] <= '2021-06-01')]
    # we just need the cases with id 3101
    df_abb = df_abb[df_abb['ID_County'] == 3101]
    # df_abb['Deaths'] = np.round(df_abb[['Deaths']].to_numpy())

    # we need the amount of dead persons for each age group: These are A00-A04, A05-A14, A15-A34, A35-A59, A60-A79, A80+
    age_groups = ['A00-A04', 'A05-A14', 'A15-A34', 'A35-A59', 'A60-A79', 'A80+']
    age_grous_string = ['Age 0-4', 'Age 5-14', 'Age 15-34', 'Age 35-59', 'Age 60-79', 'Age 80+']
    # we need to sum up the amount of dead persons for each age group

    # we want the deaths for the age groups
    df_abb = df_abb[['Date', 'Deaths', 'Age_RKI']]
    # we want a plot with 2 rows. Second row has a plot with each age group and the simulated and real dead persons
    # First row has the cumulative dead persons
    fig = plt.figure('Deaths')
    fig.set_figwidth(20)
    fig.set_figheight(9)
    gs = fig.add_gridspec(2,6)

    # we need the cumulative dead persons
    ax = fig.add_subplot(gs[0, :])
    df_total_dead = df_abb.groupby('Date').sum()[0:90]
    y_real = df_total_dead['Deaths'].to_numpy()
    # we need to substract the first value from the rest
    y_real = y_real - y_real[0]

    y_sim = p50_bs['Total'][()][:, 7][::24][0:90]
    y_sim = y_sim - y_sim[0]

    y_sim25 = p25_bs['Total'][()][:, 7][::24][0:90]
    y_sim25 = y_sim25 - y_sim25[0]

    y_sim75 = p75_bs['Total'][()][:, 7][::24][0:90]
    y_sim75 = y_sim75 - y_sim75[0]

    y_sim05 = p05_bs['Total'][()][:, 7][::24][0:90]
    y_sim05 = y_sim05 - y_sim05[0]
    
    y_sim95 = p95_bs['Total'][()][:, 7][::24][0:90]
    y_sim95 = y_sim95 - y_sim95[0]



    # we calculate the RMSE
    rmse_dead = np.sqrt(((y_real- y_sim)**2).mean())
    # we need to plot the cumulative dead persons from the real world and from the simulation
   
    ax.plot(df_total_dead.index, y_sim, color='tab:blue',label='Simulated deaths')
    ax.plot(df_total_dead.index, y_real, 'v',color='tab:red', linewidth=4, label='Extrapolated deaths from reported infection case data')
    ax.fill_between(df_total_dead.index, y_sim75, y_sim25,
                            alpha=0.5, color='tab:blue', label='50% Confidence interval')
    ax.fill_between(df_total_dead.index,y_sim95, y_sim05,
                            alpha=0.25, color='tab:blue', label='90% Confidence interval')
    # ax.text(0.25, 0.8, 'RMSE: '+str(float("{:.2f}".format(rmse_dead))), horizontalalignment='center',
    #         verticalalignment='center', transform=plt.gca().transAxes, color='pink', fontsize=15)
    ax.set_label('Number of individuals')
    ax.set_title('Cumulative Deaths', fontsize=fontsize)
    ax.set_ylabel('Number of individuals', fontsize=fontsize-8)
    ax.legend(fontsize=fontsize-8)

    # now for each age group
    for i, age_group in zip(range(6), age_group_access):
        ax = fig.add_subplot(gs[1, i])
        # we need the amount of dead persons for each age group 
        df_abb_age_group = df_abb[df_abb['Age_RKI'] == age_groups[i]][0:90]
        y_real =  np.round(df_abb_age_group['Deaths'].to_numpy())
        # we need to plot the dead persons from the real world and from the simulation
        ax.plot(df_abb_age_group['Date'], y_real-y_real[0], color='tab:red')
        ax.plot(df_abb_age_group['Date'], p50_bs[age_group_access[i]][()][:, 7][::24][0:90]-p50_bs[age_group_access[i]][()][:, 7][::24][0], color='tab:blue')
        ax.fill_between(df_abb_age_group['Date'], p75_bs[age_group_access[i]][()][:, 7][::24][0:90]-p75_bs[age_group_access[i]][()][:, 7][::24][0], p25_bs[age_group_access[i]][()][:, 7][::24][0:90]-p25_bs[age_group_access[i]][()][:, 7][::24][0],
                            alpha=0.5, color='tab:blue')
        ax.set_title('Deaths, '+age_grous_string[i])
        ax.set_ybound(lower=0)
        ax.set_xticks(df_abb_age_group['Date'][::50])
        ax.tick_params(axis='both', which='major', labelsize=fontsize-10)
        ax.tick_params(axis='both', which='minor', labelsize=fontsize-10)
        if i == 0:
            ax.set_ylabel('Number of individuals',fontsize=fontsize-8)
            ax.set_ybound(upper=1)
    
    plt.show()
   
def plot_icu(path):
    
    df_abb = pd.read_json(path+"/../../../pydata/Germany/county_divi.json")

    perc_of_critical_in_icu_age = [0.55,0.55,0.55,0.56,0.54,0.46]
    perc_of_critical_in_icu=0.55

    age_group_access = ['Group1', 'Group2', 'Group3',
                        'Group4', 'Group5', 'Group6', 'Total']


    # we just need the columns ICU_low and ICU_hig
    df_abb = df_abb[['ID_County', 'ICU', 'Date']]

    df_abb = df_abb[df_abb['ID_County'] == 3101]
    # we need just the dates bewteen 2021-03-01 and 2021-06-01
    df_abb = df_abb[(df_abb['Date'] >= '2021-03-01') &
                    (df_abb['Date'] <= '2021-06-01')]

    # we plot this against this the Amount of persons in the ICU from our model
    f_p50 = h5py.File(
        path+"/infection_state_per_age_group/0/p50/Results.h5", 'r')
    total_50 = f_p50['0']['Total'][()][::24][0:90]

    total_50_age = f_p50['0'][age_group_access[0]][()]
    for i in range(6):
              total_50_age += f_p50['0'][age_group_access[i]][()]*perc_of_critical_in_icu_age[i]
    total_50_age = total_50_age[::24][0:90]


     # we plot this against this the Amount of persons in the ICU from our model
    f_p75 = h5py.File(
        path+"/infection_state_per_age_group/0/p75/Results.h5", 'r')
    # total_75 = f_p75['0']['Total'][()][::24][0:90]
    total_75_age = f_p75['0'][age_group_access[0]][()]
    for i in range(6):
        total_75_age += f_p75['0'][age_group_access[i]][()]*perc_of_critical_in_icu_age[i]
    total_75_age = total_75_age[::24][0:90]

    # same with 25 percentile
    f_p25 = h5py.File(
        path+"/infection_state_per_age_group/0/p25/Results.h5", 'r')
    # total_25 = f_p25['0']['Total'][()][::24][0:90]
    total_25_age = f_p25['0'][age_group_access[0]][()]
    for i in range(6):
        total_25_age += f_p25['0'][age_group_access[i]][()]*perc_of_critical_in_icu_age[i]
    total_25_age = total_25_age[::24][0:90]

    # same with 05 and 95 percentile
    f_p05 = h5py.File(
        path+"/infection_state_per_age_group/0/p05/Results.h5", 'r')
    # total_05 = f_p05['0']['Total'][()][::24][0:90]
    total_05_age = f_p05['0'][age_group_access[0]][()]
    for i in range(6):
        total_05_age += f_p05['0'][age_group_access[i]][()]*perc_of_critical_in_icu_age[i]
    total_05_age = total_05_age[::24][0:90]

    f_p95 = h5py.File(
        path+"/infection_state_per_age_group/0/p95/Results.h5", 'r')
    # total_95 = f_p95['0']['Total'][()][::24][0:90]
    total_95_age = f_p95['0'][age_group_access[0]][()]
    for i in range(6):
        total_95_age += f_p95['0'][age_group_access[i]][()]*perc_of_critical_in_icu_age[i]
    total_95_age = total_95_age[::24][0:90]

    
    ICU_Simulation_one_percentile = np.floor(total_50[:, 5]*perc_of_critical_in_icu)
    ICU_Simulation = np.round(total_50_age[:, 5])
    ICU_Simulation75 = np.round(total_75_age[:, 5])
    ICU_Simulation25 = np.round(total_25_age[:, 5])
    ICU_Simulation05 = np.round(total_05_age[:, 5])
    ICU_Simulation95 = np.round(total_95_age[:, 5])
    ICU_Real = df_abb['ICU'][0:90]

    #smooth the data
    # ICU_Real = gaussian_filter1d(ICU_Real, sigma=1, mode='nearest')
    # ICU_Simulation = gaussian_filter1d(ICU_Simulation, sigma=1, mode='nearest')



    # we calculate the RMSE
    rmse_ICU = np.sqrt(((ICU_Real - ICU_Simulation_one_percentile)**2).mean())

    # plot the ICU beds and the ICU beds taken
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    fig.set_figwidth(12)
    fig.set_figheight(9)
    # we plot the ICU_low and the ICU_high
    ax.plot(df_abb['Date'][0:90], ICU_Real,'x', color='tab:red', linewidth=10, label='Data')
    ax.plot(df_abb['Date'][0:90], ICU_Simulation, color='tab:blue', label='Simulation')
    # ax.plot(df_abb['Date'][0:90], ICU_Simulation_one_percentile, color='tab:green', label='Simulated ICU beds')
    ax.fill_between(df_abb['Date'][0:90],ICU_Simulation75, ICU_Simulation25,
                         alpha=0.5, color='tab:blue', label='50% Confidence interval')
    ax.fill_between(df_abb['Date'][0:90],ICU_Simulation05, ICU_Simulation95,
                         alpha=0.25, color='tab:blue', label='90% Confidence interval')
    

    # we also write the rmse
    # ax.text(0.25, 0.8, 'RMSE: '+str(float("{:.2f}".format(rmse_ICU))), horizontalalignment='center',
    #         verticalalignment='center', transform=plt.gca().transAxes, color='pink', fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=fontsize-4)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize-4)
    ax.set_ylabel('Occupied ICU beds', fontsize=fontsize)
    ax.set_title('ICU beds', fontsize=fontsize+4)
    ax.legend(fontsize=fontsize-4)
    plt.show()




if __name__ == "__main__":
    path = ""

    if (len(sys.argv) > 1):
        n_runs = sys.argv[1]
    else:
        n_runs = len([entry for entry in os.listdir(path)
                     if os.path.isfile(os.path.join(path, entry))])
    
    # plot_icu(path)
    # plot_dead(path)
