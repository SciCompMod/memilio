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
    
    # We also print which percentage of the population was infected total, for the 50 percent interval
    # we need to sum up just the last timepoint
    # and divide it by the populat

    plt.legend(legend_plot)

    # currently the x axis has the values of the time steps, we need to convert them to dates and set the x axis to dates
    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in x]
    xx = [xx[i].strftime('%Y-%m-%d') for i in range(len(xx))]
    # but just take every 10th date to make it more readable
    plt.gca().set_xticks(x[::150])
    plt.gca().set_xticklabels(xx[::150])
    plt.gcf().autofmt_xdate()
    # these are matplotlib.patch.Patch properties
    # place a text box in upper left in axes coords

    plt.title('Infection location types')
    plt.gca().set_ylim(bottom=0)
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

    # plot_infection_states_individual(
    #     time, p50_bs, p25_bs, p75_bs)
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
        
    # We also print which percentage of the population was infected total, for the 50 percent interval
    # we need to sum up just the last timepoint

    population = 1000
    # we need to sum up just the last timepoint
    total_infected = y50[-1, 1] + y50[-1, 2] + y50[-1, 3] + \
        y50[-1, 4] + y50[-1, 5] + y50[-1, 7]
    # and divide it by the population
    total_infected = total_infected / population * 100
    # and print it
    
        

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


    time = p50_bs['Time'][()]
    time = time[::24]
    time = time[0:90]


    # we need the amount of dead persons for each age group: These are A00-A04, A05-A14, A15-A34, A35-A59, A60-A79, A80+
    age_groups = ['A00-A04', 'A05-A14', 'A15-A34', 'A35-A59', 'A60-A79', 'A80+']
    age_grous_string = ['Age 0-4', 'Age 5-14', 'Age 15-34', 'Age 35-59', 'Age 60-79', 'Age 80+']
    # we need to sum up the amount of dead persons for each age group


    # we want a plot with 2 rows. Second row has a plot with each age group and the simulated and real dead persons
    # First row has the cumulative dead persons
    fig = plt.figure('Deaths')
    fig.set_figwidth(20)
    fig.set_figheight(9)
    gs = fig.add_gridspec(2,6)

    # we need the cumulative dead persons
    ax = fig.add_subplot(gs[0, :])
    # we need to substract the first value from the rest

    y_sim = p50_bs['Total'][()][:, 7][::24][0:90]
    y_sim = y_sim - y_sim[0]

    y_sim25 = p25_bs['Total'][()][:,7][::24][0:90]
    y_sim25 = y_sim25 - y_sim25[0]

    y_sim75 = p75_bs['Total'][()][:,7][::24][0:90]
    y_sim75 = y_sim75 - y_sim75[0]

    y_sim05 = p05_bs['Total'][()][:,7][::24][0:90]
    y_sim05 = y_sim05 - y_sim05[0]
    
    y_sim95 = p95_bs['Total'][()][:, 7][::24][0:90]
    y_sim95 = y_sim95 - y_sim95[0]



   
    # we need to plot the cumulative dead persons from the real world and from the simulation
   
    ax.plot(time, y_sim, color='tab:blue',label='Simulated deaths')
    ax.fill_between(time, y_sim75, y_sim25, alpha=0.5, color='tab:blue', label='50% Confidence interval')
    ax.fill_between(time, y_sim95, y_sim05, alpha=0.25, color='tab:blue', label='90% Confidence interval')
    # ax.text(0.25, 0.8, 'RMSE: '+str(float("{:.2f}".format(rmse_dead))), horizontalalignment='center',
    #         verticalalignment='center', transform=plt.gca().transAxes, color='pink', fontsize=15)
    ax.set_label('Number of individuals')
    ax.set_title('Cumulative Deaths', fontsize=fontsize)
    ax.set_ylabel('Number of individuals', fontsize=fontsize-8)
    ax.legend(fontsize=fontsize-8)

    # # now for each age group
    # for i, age_group in zip(range(6), age_group_access):
    #     ax = fig.add_subplot(gs[1, i])
    #     # we need the amount of dead persons for each age group 
    #     df_abb_age_group = df_abb[df_abb['Age_RKI'] == age_groups[i]][0:90]
    #     y_real =  np.round(df_abb_age_group['Deaths'].to_numpy())
    #     # we need to plot the dead persons from the real world and from the simulation
    #     ax.plot(df_abb_age_group['Date'], y_real-y_real[0], color='tab:red')
    #     ax.plot(df_abb_age_group['Date'], p50_bs[age_group_access[i]][()][:, 7][::24][0:90]-p50_bs[age_group_access[i]][()][:, 7][::24][0], color='tab:blue')
    #     ax.fill_between(df_abb_age_group['Date'], p75_bs[age_group_access[i]][()][:, 7][::24][0:90]-p75_bs[age_group_access[i]][()][:, 7][::24][0], p25_bs[age_group_access[i]][()][:, 7][::24][0:90]-p25_bs[age_group_access[i]][()][:, 7][::24][0],
    #                         alpha=0.5, color='tab:blue')
    #     ax.set_title('Deaths, '+age_grous_string[i])
    #     ax.set_ybound(lower=0)
    #     ax.set_xticks(df_abb_age_group['Date'][::50])
    #     ax.tick_params(axis='both', which='major', labelsize=fontsize-10)
    #     ax.tick_params(axis='both', which='minor', labelsize=fontsize-10)
    #     if i == 0:
    #         ax.set_ylabel('Number of individuals',fontsize=fontsize-8)
    #         ax.set_ybound(upper=1)
    
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

def plot_tests(path):

    df_abb = pd.read_excel(
        path+"/pydata/Germany/SARS-CoV-2-PCR-Testungen_in_Deutschland.xlsx")
    # in the week row the format is YYYY-WX where X is the week number and Y is the year number
    # We need the week number and the year number and just take 2021-W9 to 2021-W21
    # we just Take the rows where YYYY is 2021 and X is between 9 and 21
    df_abb['Year'] = df_abb['date'].str.split('-').str[0]
    df_abb['Week'] = df_abb['date'].str.split('-').str[1]
    df_abb['Week'] = df_abb['Week'].str.split('W').str[1]
    df_abb['Week'] = df_abb['Week'].astype(int)
    df_abb['Year'] = df_abb['Year'].astype(int)
    # we just take the rows where the year is 2021 and the week is between 9 and 21
    df_abb = df_abb[(df_abb['Year'] == 2021) & (
        df_abb['Week'] >= 9) & (df_abb['Week'] <= 21)]
    # We just need the columns tests total, tests accumulated, tests positive
    df_abb = df_abb[['date', 'tests_total', 'tests_total_accumulated',
                     'tests_positive', 'tests_positive_accumulated', 'tests_positive_ratio']]
    # we assumethe tests get distributed over the week, so we divide the tests by 7 for each day, but we also need to infer the amount of tests in brunswick
    # as brunswick has around 250.000 persons, we take the 250.000/80.000.000 = 1/320 of the tests
    df_abb[['tests_total', 'tests_total_accumulated', 'tests_positive', 'tests_positive_accumulated']] = df_abb[[
        'tests_total', 'tests_total_accumulated', 'tests_positive', 'tests_positive_accumulated']]/320

    # we model this the following way
    # we know the amount of people who PCR test themselves positive and in general the amount of people who test themselves positive
    # as well as the amount of people who test themselves through other means
    # now we assume: the majority of people who are symptomatic (I_Isymp, I_Sev, I_Crit) will test themselves (80%)
    # the majority of people who are asymptomatic (I_NSymp) will not test themselves (8%)
    # and there is an amount of recently infected, which test themselves netherless (1%)

    # first thing we can do, is that we PCR test an amount of
    # we divide the persons into 2 groups: asymptomatic and symptomatic
    # we assume that the symptomatic persons will test themselves with an x time higher probability than the asymptomatic persons

    ratio_testing_symptomatic_vs_asymptomatic = 0.1

    f_p50 = h5py.File(
        path+"/results_last_run/infection_state_per_age_group/0/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]
    total_50 = total_50[::24]
    total_50 = total_50[0:90]

    # we make a np array with the amount of symptomatic and asymptomatic persons
    # first we calculate the amount of positive tests, and assume PCR is perfect
    # the amount of tests done on each day is;
    PCR_tests = df_abb['tests_total'].to_numpy()
    PCR_tests = np.repeat(PCR_tests, 7)
    PCR_tests = PCR_tests/7
    PCR_tests = PCR_tests[0:90]

    PCR_tests_positive = df_abb['tests_positive'].to_numpy()
    PCR_tests_positive = np.repeat(PCR_tests_positive, 7)
    PCR_tests_positive = PCR_tests_positive/7
    PCR_tests_positive = PCR_tests_positive[0:90]

    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in range(90)]
    xx = [xx[i].strftime('%Y-%m-%d') for i in range(len(xx))]

    # plot these
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    fig.set_figwidth(20)
    fig.set_figheight(9)
    # we plot the tests positive and the real cases
    ax.plot(xx, PCR_tests, color='tab:red')
    ax.plot(xx, PCR_tests_positive, color='tab:blue')
    ax.set_xlabel('Date')

    # The amount of persons, who do tests on  a day is:
    PCR_tests_symptomatic = np.zeros(90)
    PCR_tests_asymptomatic = np.zeros(90)
    PCR_tests_symptomatic = PCR_tests * \
        (1/(1+(1/ratio_testing_symptomatic_vs_asymptomatic)))
    PCR_tests_asymptomatic = PCR_tests - PCR_tests_symptomatic

    # the real amount of positive tested persons is:
    # likelihood of being poisiitve is (E+I_NSymp)/(S+E+I_NSymp)
    lik_being_positive_asymptomatic = (
        total_50[:, 1]+total_50[:, 2])/(total_50[:, 0]+total_50[:, 1]+total_50[:, 2])

    tests_positive = (lik_being_positive_asymptomatic *
                      PCR_tests_asymptomatic)+PCR_tests_symptomatic
    # we need to plot this
    ax.plot(xx, tests_positive, color='tab:green')
    ax.set_ylabel('Number of tests')
    ax.title.set_text(
        'Tests positive PCR and real confirmedcases from any source')
    ax.legend(['Tests', 'Tests positive', 'Tests positive inferred'])
    plt.show()

    # new plot which is showing the amount of symptomatic persons
    sympt_persons = total_50[:, 3]+total_50[:, 4]+total_50[:, 5]
    # we plot this
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    fig.set_figwidth(20)
    fig.set_figheight(9)
    # we plot the tests positive and the real cases
    ax.plot(xx, sympt_persons, color='tab:red')
    ax.set_xlabel('Date')
    ax.set_ylabel('Number of symptomatic persons')
    ax.title.set_text('Symptomatic persons')
    ax.legend(['Symptomatic persons'])
    plt.show()

    # also the amount of asymptomatic persons with the amount of persons infected
    asympt_persons = total_50[:, 2]+total_50[:, 1]+total_50[:, 0]
    asympt_positive = total_50[:, 2]+total_50[:, 1]
    # we plot this
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    fig.set_figwidth(20)
    fig.set_figheight(9)
    # we plot the tests positive and the real cases
    ax.plot(xx, asympt_persons, color='tab:red')
    ax.plot(xx, asympt_positive, color='tab:blue')
    ax.set_xlabel('time (days)')
    ax.set_ylabel('Number of asymptomatic persons')
    ax.title.set_text('Asymptomatic persons')
    ax.legend(['Asymptomatic persons', 'Asymptomatic persons positive'])
    plt.show()

def calc_positive_tests_overall(infection_states, sensitivity, specificity, r_sns, lt_sympt):

    lt_asympt = lt_sympt/r_sns
    inferred_positive_tests_sympt = (
        infection_states[:, 3]*lt_sympt+infection_states[:, 4]*lt_sympt+infection_states[:, 5]*lt_sympt)*sensitivity
    # asymptomatic persons
    inferred_positive_tests_asympt = (infection_states[:, 0]*lt_asympt)*(
        1-specificity)+((infection_states[:, 1]+infection_states[:, 2])*lt_asympt)*sensitivity
    return inferred_positive_tests_sympt+inferred_positive_tests_asympt, inferred_positive_tests_sympt, inferred_positive_tests_asympt

def plot_estimated_reproduction_number(path):
    f_p50 = h5py.File(
        path+"/estimated_reproduction_number/0/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]
    total_50 = total_50[::24]
    total_50 = total_50[0:90].flatten()
    # we smooth this with a gaussian filter
    total_50 = gaussian_filter1d(total_50, sigma=1, mode='nearest')
    time = p50_bs['Time'][()]
    time = time[::24]
    time = time[0:90]

    # we plot this
    # we plot the tests positive and the real cases
    plt.plot(time, total_50, color='tab:red')
    plt.xlabel('time (days)')
    plt.ylabel('Estimated reproduction number')
    plt.title('Estimated reproduction number')
    plt.show()

def plot_cumulative_detected_infections(path):

    df_abb = pd.read_json(
        path+"/../../../pydata/Germany/cases_all_county_repdate_ma1.json")
    # we need the 
    df_abb = df_abb[['Date', 'Confirmed', 'ID_County']]
    df_abb = df_abb[(df_abb['Date'] >= '2021-03-01') & (df_abb['Date'] <= '2021-06-01')]
    df_abb = df_abb[df_abb['ID_County'] == 3101]
    df_substract = np.round(df_abb['Confirmed'][0:1])
    df_abb =   np.round(df_abb['Confirmed'][0:90])
    df_abb = df_abb - df_substract.values[0]
    df_diff = df_abb.diff()
    # we also want to plot the amount of new detected infections on a 7 day rolling average
    df_abb_rolling = df_diff.rolling(window=7, min_periods=1, center=True).mean()
 
    


    f_p50 = h5py.File(
        path+"/cumulative_detected_infections/0/p50/Results.h5", 'r')
    total_50 = np.round(f_p50['0']['Total'][()][::24][0:90].flatten())
    # we smooth this with a gaussian filter

    f_p95 = h5py.File(
        path+"/cumulative_detected_infections/0/p95/Results.h5", 'r')
    total_95 = np.round(f_p95['0']['Total'][()][::24][0:90].flatten())

    f_p05 = h5py.File(
        path+"/cumulative_detected_infections/0/p05/Results.h5", 'r')
    total_05 = np.round(f_p05['0']['Total'][()][::24][0:90].flatten())

    f_p25 = h5py.File(
        path+"/cumulative_detected_infections/0/p25/Results.h5", 'r')
    total_25 = np.round(f_p25['0']['Total'][()][::24][0:90].flatten())

    f_p75 = h5py.File(
        path+"/cumulative_detected_infections/0/p75/Results.h5", 'r')
    total_75 = np.round(f_p75['0']['Total'][()][::24][0:90].flatten())


    # we do the same for the new infecitons (same as above but folder   new_detected_infections)
    f_p50_diff = h5py.File(
        path+"/new_detected_infections/0/p50/Results.h5", 'r')
    total_50_diff = np.round(f_p50_diff['0']['Total'][()].flatten())

    f_p95_diff = h5py.File(
        path+"/new_detected_infections/0/p95/Results.h5", 'r')
    total_95_diff = np.round(f_p95_diff['0']['Total'][()].flatten())
    
    f_p05_diff = h5py.File(
        path+"/new_detected_infections/0/p05/Results.h5", 'r')
    total_05_diff = np.round(f_p05_diff['0']['Total'][()].flatten())

    f_p25_diff = h5py.File(
        path+"/new_detected_infections/0/p25/Results.h5", 'r')
    total_25_diff = np.round(f_p25_diff['0']['Total'][()].flatten())

    f_p75_diff = h5py.File(
        path+"/new_detected_infections/0/p75/Results.h5", 'r')
    total_75_diff = np.round(f_p75_diff['0']['Total'][()].flatten())

    # we need to sum every 24 entries to get the daily amount
    total_50_diff = np.cumsum(total_50_diff, axis=0)
    total_50_diff = np.diff(total_50_diff)[0:90]

    total_95_diff = np.cumsum(total_95_diff, axis=0)
    total_95_diff = np.diff(total_95_diff)[0:90]

    total_05_diff = np.cumsum(total_05_diff, axis=0)
    total_05_diff = np.diff(total_05_diff)[0:90]

    total_75_diff = np.cumsum(total_75_diff, axis=0)
    total_75_diff = np.diff(total_75_diff)[0:90]

    total_25_diff = np.cumsum(total_25_diff, axis=0)
    total_25_diff = np.diff(total_25_diff)[0:90]





    # we smooth this with a gaussian filter
    time = f_p50['0']['Time'][()]
    time = time[::24]
    time = time[0:90]

    # we calculate the RMSE
    rmse_detected = np.sqrt(((df_abb - total_50)**2).mean())



    # we plot this
    # we plot the tests positive and the real cases
    fig = plt.figure('Cumulative detected infections', constrained_layout=True)
    fig.set_figwidth(20)
    fig.set_figheight(9)
    plt.plot(time, total_50, color='tab:blue', label='Simulation')
    plt.plot(time, df_abb, 'x', color='tab:red', label='Data', linewidth=4)
    plt.fill_between(time, total_75, total_25,
                            alpha=0.5, color='tab:blue', label='90% Confidence interval')
    plt.fill_between(time, total_95, total_05,
                            alpha=0.25, color='tab:blue', label='50% Confidence interval')

    
    plt.ylabel('Cumulative detected infections', fontsize=fontsize-8)
    plt.title('Cumulative detected infections', fontsize=fontsize)
    plt.legend(fontsize=fontsize-8)
    # currently the x axis has the values of the time steps, we need to convert them to dates and set the x axis to dates
    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in time]
    xx = [xx[i].strftime('%Y-%m-%d') for i in range(90)]
    # but just take every 10th date to make it more readable
    plt.gca().set_xticks(time[::14])
    plt.gca().set_xticklabels(xx[::14])
    plt.gcf().autofmt_xdate()
    #rmse
    # plt.text(0.25, 0.8, 'RMSE: '+str(float("{:.2f}".format(rmse_detected))), horizontalalignment='center',
    #         verticalalignment='center', transform=plt.gca().transAxes, color='pink', fontsize=15)
    
    plt.show()

    test_p_pos_p50_normal = h5py.File(
        path+"/positive_test_per_location_type_per_age_group/"+"0"+"/p50/Results.h5", 'r')
    p50_bs_test_p_pos_normal = test_p_pos_p50_normal['0']['Total'][()]

    # we need to sum every 24 entries to get the daily amount
    total_50_positive = np.sum(p50_bs_test_p_pos_normal, axis=1)
    total_50_positive = np.cumsum(total_50_positive, axis=0)
    total_50_positive = total_50_positive[::24]
    total_50_positive = total_50_positive[0:90] # we still need to take the difference to get the daily amount
    total_50_positive = np.diff(total_50_positive, axis=0).flatten()



    # also the amount of new detected infections
    fig = plt.figure('New detected infections', constrained_layout=True)
    fig.set_figwidth(20)
    fig.set_figheight(9)
    # we plot the tests positive and the real cases
    plt.plot(time[0:89], total_50_diff, color='tab:blue', label='Simulation')
    # we dont plot the real curve as a line but as x points and not every day but every 2nd day
    plt.plot(time[0:89], df_diff[0:89], 'x', color='tab:red', label='Data', linewidth=4)
    plt.plot(time[0:89], total_50_positive, color='tab:green', label='Data 7 day rolling average', linewidth=3, linestyle='dashed')
    # also the rolling average
    plt.plot(time[0:89], df_abb_rolling[0:89], color='tab:red', label='Data 7 day rolling average', linewidth=3, linestyle='dashed')
    plt.fill_between(time[0:89], total_75_diff, total_25_diff,
                            alpha=0.5, color='tab:blue', label='90% Confidence interval')   
    plt.fill_between(time[0:89], total_95_diff, total_05_diff,
                            alpha=0.25, color='tab:blue', label='50% Confidence interval')
  
    plt.ylabel('New detected infections', fontsize=fontsize-8)
    plt.title('Daily detected infections', fontsize=fontsize)
    plt.legend(fontsize=fontsize-5)
    # currently the x axis has the values of the time steps, we need to convert them to dates and set the x axis to dates
    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in time]
    xx = [xx[i].strftime('%Y-%m-%d') for i in range(90)]
    # but just take every 10th date to make it more readable
    plt.gca().set_xticks(time[0:89][::14])
    plt.gca().set_xticklabels(xx[::14])
    plt.gcf().autofmt_xdate()

    plt.show()

def plot_positive_and_done_test(path):
    f_p50_positive = h5py.File(
        path+"/positive_test_per_location_type_per_age_group/0/p50/Results.h5", 'r')
    p50_bs_positive = f_p50_positive['0']
    total_50_positive = p50_bs_positive['Total'][()]

    f_p50_done = h5py.File(
        path+"/test_per_location_type_per_age_group/0/p50/Results.h5", 'r')
    p50_bs_done = f_p50_done['0']
    total_50_done = p50_bs_done['Total'][()]

    time = p50_bs_positive['Time'][()][::24][0:90]

    # weas one entry is one hour we take the sum every 24 entries to get the daily amount, we do this with cumsum
    # first we need to sum up over all age groups
    total_50_positive = np.sum(total_50_positive, axis=1)
    total_50_positive = np.cumsum(total_50_positive, axis=0)
    total_50_positive = total_50_positive[::24]
    total_50_positive = total_50_positive[0:91] # we still need to take the difference to get the daily amount
    total_50_positive = np.diff(total_50_positive, axis=0).flatten()
    # we smooth this with a gaussian filter
    total_50_positive = gaussian_filter1d(total_50_positive, sigma=1, mode='nearest')

    #same for the done tests
    total_50_done = np.sum(total_50_done, axis=1)
    total_50_done = np.cumsum(total_50_done, axis=0)
    total_50_done = total_50_done[::24]
    total_50_done = total_50_done[0:91] # we still need to take the difference to get the daily amount
    total_50_done = np.diff(total_50_done, axis=0).flatten()
    # we smooth this with a gaussian filter
    total_50_done = gaussian_filter1d(total_50_done, sigma=2, mode='nearest')

    # we plot this
    # we plot the tests positive and the real cases
    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in range(90)]
    xx = [xx[i].strftime('%Y-%m-%d') for i in range(len(xx))]
    plt.gca().set_xticks(time[::5])
    plt.gca().set_xticklabels(xx[::5])
    plt.gcf().autofmt_xdate()

    plt.plot(xx, total_50_positive, color='tab:green')
    plt.plot(xx, total_50_done, color='tab:red')
    plt.xlabel('time (days)')
    plt.ylabel('Number of tests')
    plt.legend(['Positive tests', 'Done tests'])
    plt.title('Positive and done tests')
    plt.show()


def plot_fitting_plots(path):

    # We want to have the fitting for the four things we fittet against: Cumulative deaths, ICU, Cumulative detected infections and new detected infections
    # readin of the data
    db_abb_deaths = pd.read_json(
        path+"/../../../pydata/Germany/cases_all_county_age_ma1.json")

    db_abb_icu = pd.read_json(
        path+"/../../../pydata/Germany/county_divi.json")

    db_abb_detected = pd.read_json(
        path+"/../../../pydata/Germany/cases_all_county_repdate_ma1.json")

    # we want to plot the plots in a 4x1 grid and we want to have the same x axis for all of them
    # the first plot has also the legend
    # in each plot we want to plot the real data and the simulated data, as well as the confidence intervals 50% and 90%

    fig, axs = plt.subplots(2, 2, figsize=(25, 15), constrained_layout=True)

    # we want to plot the cumulative deaths
    # we need the cumulative dead persons
    df_total_dead = db_abb_deaths
    df_total_dead['Date'] = df_total_dead['Date']+pd.DateOffset(days=18)
    df_total_dead = df_total_dead[(df_total_dead['Date'] >= '2021-03-01') &
                    (df_total_dead['Date'] <= '2021-06-01')]
    df_total_dead = df_total_dead[df_total_dead['ID_County'] == 3101]
    df_total_dead = df_total_dead.groupby('Date').sum()[0:90]
    deaths_real = df_total_dead['Deaths'].to_numpy()
    mse_death = deaths_real
    deaths_real = deaths_real[0:90] - deaths_real[0]
    

    # simulation deaths and confidence intervals
    f_p50_deaths = h5py.File(
        path+"/infection_state_per_age_group/0/p50/Results.h5", 'r')
    total_50_deaths = f_p50_deaths['0']['Total'][()][:, 7][::24][0:90]  
    mse_sim_dead = total_50_deaths
    total_50_deaths = total_50_deaths - total_50_deaths[0]

    f_p75_deaths = h5py.File(
        path+"/infection_state_per_age_group/0/p75/Results.h5", 'r')
    total_75_deaths = f_p75_deaths['0']['Total'][()][:, 7][::24][0:90]
    total_75_deaths = total_75_deaths - total_75_deaths[0]

    f_p25_deaths = h5py.File(
        path+"/infection_state_per_age_group/0/p25/Results.h5", 'r')
    total_25_deaths = f_p25_deaths['0']['Total'][()][:, 7][::24][0:90]
    total_25_deaths = total_25_deaths - total_25_deaths[0]

    f_p05_deaths = h5py.File(
        path+"/infection_state_per_age_group/0/p05/Results.h5", 'r')
    total_05_deaths = f_p05_deaths['0']['Total'][()][:, 7][::24][0:90]
    total_05_deaths = total_05_deaths - total_05_deaths[0]

    f_p95_deaths = h5py.File(
        path+"/infection_state_per_age_group/0/p95/Results.h5", 'r')
    total_95_deaths = f_p95_deaths['0']['Total'][()][:, 7][::24][0:90]
    total_95_deaths = total_95_deaths - total_95_deaths[0]

    # we plot the deaths
    axs[0,0].plot(df_total_dead.index, deaths_real, 'v', color='tab:red',linewidth=4, label='Extrapolated deaths from reported infection case data')
    axs[0,0].plot(df_total_dead.index, total_50_deaths, color='tab:blue', label='Simulated deaths')
    axs[0,0].fill_between(df_total_dead.index, total_75_deaths, total_25_deaths,
                            alpha=0.5, color='tab:blue', label='50% Confidence interval')   
    axs[0,0].fill_between(df_total_dead.index, total_95_deaths, total_05_deaths,
                            alpha=0.25, color='tab:blue', label='90% Confidence interval')
    axs[0,0].set_ylabel('Cumulative deaths', fontsize=fontsize)
    axs[0,0].set_title('Cumulative deaths', fontsize=fontsize+4)
    axs[0,0].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[0,0].tick_params(axis='both', which='minor', labelsize=fontsize-4)
    axs[0,0].legend(fontsize=fontsize-4, loc='upper left')

    # we want to plot the ICU
    
    perc_of_critical_in_icu_age = [0.55,0.55,0.55,0.56,0.54,0.46]
    age_group_access = ['Group1', 'Group2', 'Group3',
                        'Group4', 'Group5', 'Group6', 'Total']
    df_abb_icu = db_abb_icu[['ID_County', 'ICU', 'Date']]
    df_abb_icu = df_abb_icu[df_abb_icu['ID_County'] == 3101]
    df_abb_icu = df_abb_icu[(df_abb_icu['Date'] >= '2021-03-01') &
                    (df_abb_icu['Date'] <= '2021-06-01')]
    ICU_real = df_abb_icu['ICU'].to_numpy()

    # we need to get the simulated ICU values
    f_p50_icu = h5py.File(
        path+"/infection_state_per_age_group/0/p50/Results.h5", 'r')
    total_50_icu = f_p50_icu['0']['Total'][()]
    total_50_icu = total_50_icu[::24][0:90]
    total_50_age = f_p50_icu['0'][age_group_access[0]][()]
    for i in range(6):
              total_50_age += f_p50_icu['0'][age_group_access[i]][()]*perc_of_critical_in_icu_age[i]
    total_50_age = total_50_age[::24][0:90]

    # we plot this against this the Amount of persons in the ICU from our model
    f_p75_icu = h5py.File(
        path+"/infection_state_per_age_group/0/p75/Results.h5", 'r')
    # total_75 = f_p75['0']['Total'][()][::24][0:90]
    total_75_age = f_p75_icu['0'][age_group_access[0]][()]
    for i in range(6):
        total_75_age += f_p75_icu['0'][age_group_access[i]][()]*perc_of_critical_in_icu_age[i]
    total_75_age = total_75_age[::24][0:90]

    # same with 25 percentile
    f_p25_icu = h5py.File(
        path+"/infection_state_per_age_group/0/p25/Results.h5", 'r')
    # total_25 = f_p25['0']['Total'][()][::24][0:90]
    total_25_age = f_p25_icu['0'][age_group_access[0]][()]
    for i in range(6):
        total_25_age += f_p25_icu['0'][age_group_access[i]][()]*perc_of_critical_in_icu_age[i]
    total_25_age = total_25_age[::24][0:90]

    # same with 05 and 95 percentile
    f_p05_icu = h5py.File(
        path+"/infection_state_per_age_group/0/p05/Results.h5", 'r')
    # total_05 = f_p05['0']['Total'][()][::24][0:90]
    total_05_age = f_p05_icu['0'][age_group_access[0]][()]
    for i in range(6):
        total_05_age += f_p05_icu['0'][age_group_access[i]][()]*perc_of_critical_in_icu_age[i]
    total_05_age = total_05_age[::24][0:90]

    f_p95_icu = h5py.File(
        path+"/infection_state_per_age_group/0/p95/Results.h5", 'r')
    # total_95 = f_p95['0']['Total'][()][::24][0:90]
    total_95_age = f_p95_icu['0'][age_group_access[0]][()]
    for i in range(6):
        total_95_age += f_p95_icu['0'][age_group_access[i]][()]*perc_of_critical_in_icu_age[i]
    total_95_age = total_95_age[::24][0:90]

    ICU_Simulation = np.round(total_50_age[:, 5])
    ICU_Simulation75 = np.round(total_75_age[:, 5])
    ICU_Simulation25 = np.round(total_25_age[:, 5])
    ICU_Simulation05 = np.round(total_05_age[:, 5])
    ICU_Simulation95 = np.round(total_95_age[:, 5])

    # we plot the ICU beds and the ICU beds taken
    axs[1,0].plot(df_abb_icu['Date'][0:90], ICU_real[0:90], 'v', color='tab:red', linewidth=10, label='Reported ICU beds taken')
    axs[1,0].plot(df_abb_icu['Date'][0:90], ICU_Simulation[0:90], color='tab:blue', label='Simulation')
    axs[1,0].fill_between(df_abb_icu['Date'][0:90], ICU_Simulation75, ICU_Simulation25,
                            alpha=0.5, color='tab:blue', label='50% Confidence interval')
    axs[1,0].fill_between(df_abb_icu['Date'][0:90], ICU_Simulation05, ICU_Simulation95, 
                            alpha=0.25, color='tab:blue', label='90% Confidence interval')
    axs[1,0].set_ylabel('Occupied ICU beds', fontsize=fontsize)
    axs[1,0].set_title('ICU beds', fontsize=fontsize+4)
    axs[1,0].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[1,0].tick_params(axis='both', which='minor', labelsize=fontsize-4)
    axs[1,0].legend(fontsize=fontsize-4)

    # we want to plot the cumulative detected infections
    df_total_detected = db_abb_detected[['Date', 'Confirmed', 'ID_County']]
    df_total_detected = df_total_detected[(df_total_detected['Date'] >= '2021-03-01') &
                    (df_total_detected['Date'] <= '2021-06-01')]
    df_total_detected = df_total_detected[df_total_detected['ID_County'] == 3101]
    df_substract = np.round(df_total_detected['Confirmed'][0:1])
    df_total_detected = np.round(df_total_detected['Confirmed'][0:90])
    df_total_detected = df_total_detected - df_substract.values[0]
    df_detected_diff = df_total_detected.diff()
    df_detected_rolling_diff = df_detected_diff.rolling(window=7, min_periods=1, center=True).mean()

    # simulation deaths and confidence intervals
    f_p50_detected = h5py.File(
        path+"/cumulative_detected_infections/0/p50/Results.h5", 'r')
    total_50_detected = np.round(f_p50_detected['0']['Total'][()][:, 0][::24][0:90].flatten())

    f_p75_detected = h5py.File(
        path+"/cumulative_detected_infections/0/p75/Results.h5", 'r')
    total_75_detected = np.round(f_p75_detected['0']['Total'][()][:, 0][::24][0:90].flatten())
    
    f_p25_detected = h5py.File(
        path+"/cumulative_detected_infections/0/p25/Results.h5", 'r')
    total_25_detected = np.round(f_p25_detected['0']['Total'][()][:, 0][::24][0:90].flatten())

    f_p05_detected = h5py.File(
        path+"/cumulative_detected_infections/0/p05/Results.h5", 'r')
    total_05_detected = np.round(f_p05_detected['0']['Total'][()][:, 0][::24][0:90].flatten())

    f_p95_detected = h5py.File(
        path+"/cumulative_detected_infections/0/p95/Results.h5", 'r')
    total_95_detected = np.round(f_p95_detected['0']['Total'][()][:, 0][::24][0:90].flatten())

    # we do the same for the new infecitons (same as above but folder   new_detected_infections)
    f_p50_diff_detected = h5py.File(
        path+"/new_detected_infections/0/p50/Results.h5", 'r')
    total_50_diff_detected = np.round(f_p50_diff_detected['0']['Total'][()].flatten())

    f_p95_diff_detected = h5py.File(
        path+"/new_detected_infections/0/p95/Results.h5", 'r')
    total_95_diff_detected = np.round(f_p95_diff_detected['0']['Total'][()].flatten())

    f_p05_diff_detected = h5py.File(
        path+"/new_detected_infections/0/p05/Results.h5", 'r')
    total_05_diff_detected = np.round(f_p05_diff_detected['0']['Total'][()].flatten())

    f_p25_diff_detected = h5py.File(
        path+"/new_detected_infections/0/p25/Results.h5", 'r')
    total_25_diff_detected = np.round(f_p25_diff_detected['0']['Total'][()].flatten())
    
    f_p75_diff_detected = h5py.File(
        path+"/new_detected_infections/0/p75/Results.h5", 'r')
    total_75_diff_detected = np.round(f_p75_diff_detected['0']['Total'][()].flatten())

    # we need to sum every 24 entries to get the daily amount
    total_50_detected = np.cumsum(total_50_detected, axis=0)
    total_50_detected = np.diff(total_50_detected)[0:90]

    total_95_detected = np.cumsum(total_95_detected, axis=0)
    total_95_detected = np.diff(total_95_detected)[0:90]

    total_05_detected = np.cumsum(total_05_detected, axis=0)
    total_05_detected = np.diff(total_05_detected)[0:90]

    total_75_detected = np.cumsum(total_75_detected, axis=0)
    total_75_detected = np.diff(total_75_detected)[0:90]

    total_25_detected = np.cumsum(total_25_detected, axis=0)
    total_25_detected = np.diff(total_25_detected)[0:90]

    time = f_p50_detected['0']['Time'][()]
    time = time[::24]
    time = time[0:90]

    # we plot the cumulative detected infections
    axs[0,1].plot(df_abb_icu['Date'][0:89], total_50_detected, color='tab:blue', label='Simulation')
    axs[0,1].plot(df_abb_icu['Date'][0:89], df_total_detected[0:89], 'v', color='tab:red', label='Cumulative reported detected infections', linewidth=4)
    axs[0,1].fill_between(df_abb_icu['Date'][0:89], total_75_detected, total_25_detected,
                            alpha=0.5, color='tab:blue', label='50% Confidence interval')
    axs[0,1].fill_between(df_abb_icu['Date'][0:89], total_95_detected, total_05_detected,
                            alpha=0.25, color='tab:blue', label='90% Confidence interval')
    axs[0,1].set_ylabel('Cumulative detected infections', fontsize=fontsize)
    axs[0,1].set_title('Cumulative detected infections', fontsize=fontsize+4)
    axs[0,1].tick_params(axis='both', which='major', labelsize=fontsize-4) 
    axs[0,1].tick_params(axis='both', which='minor', labelsize=fontsize-4)
    axs[0,1].legend(fontsize=fontsize-4)
    

    # we want to plot the new detected infections
    axs[1,1].plot(df_abb_icu['Date'][0:90], total_50_diff_detected, color='tab:blue', label='Simulation')
    axs[1,1].plot(df_abb_icu['Date'][0:90], df_detected_diff[0:90], 'v', color='tab:red', label='Daily reported detected infections', linewidth=4)
    axs[1,1].plot(df_abb_icu['Date'][0:90], df_detected_rolling_diff, '--', color='tab:red', label='Reported detected infections 7 day rolling average', linewidth=4)
    axs[1,1].fill_between(df_abb_icu['Date'][0:90], total_75_diff_detected, total_25_diff_detected,
                            alpha=0.5, color='tab:blue', label='50% Confidence interval')
    axs[1,1].fill_between(df_abb_icu['Date'][0:90], total_95_diff_detected, total_05_diff_detected,
                            alpha=0.25, color='tab:blue', label='90% Confidence interval')
    axs[1,1].set_ylabel('New detected infections', fontsize=fontsize)
    axs[1,1].set_title('New detected infections', fontsize=fontsize+4)
    axs[1,1].tick_params(axis='both', which='major', labelsize=fontsize-4)
    axs[1,1].tick_params(axis='both', which='minor', labelsize=fontsize-4)
    axs[1,1].legend(fontsize=fontsize-4)

    # we also need to calculate the MSE for the fitting
    mse_deaths = ((mse_death - mse_sim_dead)**2).mean()
    mse_icu = ((ICU_real[0:90] - ICU_Simulation[0:90])**2).mean()
    mse_detected = ((df_total_detected[0:89] - total_50_detected)**2).mean()

    mse_final = mse_deaths + mse_icu*0.1 + mse_detected*0.01*0.01*3
    print('MSE Deaths: ', mse_deaths)
    print('MSE ICU: ', mse_icu) 
    print('MSE Detected: ', mse_detected)
    print('MSE Final: ', mse_final)


   
    plt.savefig('fitting_plots.png', dpi=300)
   
   


if __name__ == "__main__":
    # path = "/Users/david/Documents/HZI/memilio/data/results_last_run"
    # path = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/results_last_run"
    # path = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/cluster_results/vorlaufige_ergebnisse/results_2024-08-28113051"
    path = "/Users/saschakorf/Nosynch/Arbeit/memilio/memilio/cpp/examples/results/results_last_run"

    if (len(sys.argv) > 1):
        n_runs = sys.argv[1]
    else:
        n_runs = len([entry for entry in os.listdir(path)
                     if os.path.isfile(os.path.join(path, entry))])
    plot_infection_states_results(path)
    plot_infections_loc_types_avarage(path)
    # plot_icu(path)
    # plot_dead(path)
    # plot_cumulative_detected_infections(path)
    # plot_positive_and_done_test(path)

    # plot_estimated_reproduction_number(path)
    # plot_fitting_plots(path)