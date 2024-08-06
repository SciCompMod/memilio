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


def plot_infections_loc_types_avarage(path):
    # 50-percentile
    f_p50 = h5py.File(
        path+"/infection_per_location_type_per_age_group/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]

    # 25-percentile
    f_p25 = h5py.File(
        path+"/infection_per_location_type_per_age_group/p25/Results.h5", 'r')
    p25_bs = f_p25['0']
    total_25 = p25_bs['Total'][()]

    # 75-percentile
    f_p75 = h5py.File(
        path + "/infection_per_location_type_per_age_group/p75/Results.h5", 'r')
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

def plot_dead(path):
    # we will have a seperate plot the cumulative infected individuals, cumulative symptomatic individuals and cumulative dead individual
    # we need to load the data
    f_p50 = h5py.File(
        path+"/infection_state_per_age_group/p50/Results.h5", 'r')
    p50_bs = f_p50['0']

    age_group_access = ['Group1', 'Group2', 'Group3',
                        'Group4', 'Group5', 'Group6', 'Total']

    # we need the real data json file cases_all_state_repdate_ma7
    df_abb = pd.read_json(
        path+"/../pydata/Germany/cases_all_county_age_ma7.json")

    # we just need the columns cases and date
    # we need to offset the dates by 19 day
    df_abb['Date'] = df_abb['Date'] + pd.DateOffset(days=19)
    # we need just the dates bewteen 2021-03-01 and 2021-06-01
    df_abb = df_abb[(df_abb['Date'] >= '2021-03-01') &
                    (df_abb['Date'] <= '2021-06-01')]
    # we just need the cases with id 3101
    df_abb = df_abb[df_abb['ID_County'] == 3101]
    df_abb['Deaths'] = np.floor(df_abb[['Deaths']].to_numpy())

    # we need the amount of dead persons for each age group: These are A00-A04, A05-A14, A15-A34, A35-A59, A60-A79, A80+
    age_groups = ['A00-A04', 'A05-A14', 'A15-A34', 'A35-A59', 'A60-A79', 'A80+']
    # we need to sum up the amount of dead persons for each age group

    # we want the deaths for the age groups
    df_abb = df_abb[['Date', 'Deaths', 'Age_RKI']]
    # we want a plot with 2 rows. Second row has a plot with each age group and the simulated and real dead persons
    # First row has the cumulative dead persons
    fig = plt.figure('Dead')
    fig.set_figwidth(20)
    fig.set_figheight(9)
    gs = fig.add_gridspec(2,6)

    # we need the cumulative dead persons
    ax = fig.add_subplot(gs[0, :])
    df_total_dead = df_abb.groupby('Date').sum()[0:90]
    y_real = df_total_dead['Deaths'].to_numpy()

    # we calculate the RMSE
    rmse_dead = np.sqrt(((y_real- p50_bs['Total'][()][:, 7][::24][0:90])**2).mean())
    # we need to plot the cumulative dead persons from the real world and from the simulation
    ax.plot(df_total_dead.index, y_real, color='tab:blue')
    ax.plot(df_total_dead.index, p50_bs['Total'][()][:, 7][::24][0:90], color='tab:red')
    ax.text(0.25, 0.8, 'RMSE: '+str(float("{:.2f}".format(rmse_dead))), horizontalalignment='center',
            verticalalignment='center', transform=plt.gca().transAxes, color='pink', fontsize=15)
    ax.set_xlabel('Date')
    ax.set_label('Number of dead')
    ax.set_title('Cumulative Dead')
    ax.legend(['Real dead', 'Simulated dead'])

    # now for each age group
    for i, age_group in zip(range(6), age_group_access):
        ax = fig.add_subplot(gs[1, i])
        # we need the amount of dead persons for each age group 
        df_abb_age_group = df_abb[df_abb['Age_RKI'] == age_groups[i]][0:90]
        y_real =  np.floor(df_abb_age_group['Deaths'].to_numpy())
        # we need to plot the dead persons from the real world and from the simulation
        ax.plot(df_abb_age_group['Date'], y_real, color='tab:blue')
        ax.plot(df_abb_age_group['Date'], p50_bs[age_group_access[i]][()][:, 7][::24][0:90], color='tab:red')
        ax.set_title('Dead, Age{}'.format(i))
        ax.set_xlabel('Date')
        ax.set_xticks(df_abb_age_group['Date'][::50])
        ax.legend(['Real dead, Age{}'.format(i), 'Simulated dead, Age{}'.format(i)])
    
    plt.show()
   
def plot_icu(path):
    df_abb = pd.read_json(path+"/../pydata/Germany/county_divi_ma7.json")
    perc_of_critical_in_icu = 0.45

    # we just need the columns ICU_low and ICU_hig
    df_abb = df_abb[['ID_County', 'ICU', 'Date']]

    df_abb = df_abb[df_abb['ID_County'] == 3101]
    # we need just the dates bewteen 2021-03-01 and 2021-06-01
    df_abb = df_abb[(df_abb['Date'] >= '2021-03-01') &
                    (df_abb['Date'] <= '2021-06-01')]

    # we plot this against this the Amount of persons in the ICU from our model
    f_p50 = h5py.File(
        path+"/infection_state_per_age_group/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]
    # we need just every 24th value
    total_50 = total_50[::24]
    # we just take the first 90 days
    total_50 = total_50[0:90]

    ICU_Simulation = np.floor(total_50[:, 5]*perc_of_critical_in_icu)
    ICU_Real = np.floor(df_abb['ICU'][0:90])

    #smooth the data
    # ICU_Real = gaussian_filter1d(ICU_Real, sigma=1, mode='nearest')
    # ICU_Simulation = gaussian_filter1d(ICU_Simulation, sigma=1, mode='nearest')



    # we calculate the RMSE
    rmse_ICU = np.sqrt(((ICU_Real - ICU_Simulation)**2).mean())

    # plot the ICU beds and the ICU beds taken
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    fig.set_figwidth(20)
    fig.set_figheight(9)
    # we plot the ICU_low and the ICU_high
    ax.plot(df_abb['Date'][0:90], ICU_Real, color='tab:blue')
    ax.plot(df_abb['Date'][0:90], ICU_Simulation, color='tab:red')
    # we also write the rmse
    ax.text(0.25, 0.8, 'RMSE: '+str(float("{:.2f}".format(rmse_ICU))), horizontalalignment='center',
            verticalalignment='center', transform=plt.gca().transAxes, color='pink', fontsize=15)
    ax.set_xlabel('Date')
    ax.set_ylabel('Number of persons in an ICU')
    ax.title.set_text('Simulated and real ICU beds occupied')
    ax.legend(['Real in ICU, smoothed','Simulated in ICU, Perc of Critical:'+str(perc_of_critical_in_icu)])
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
        path+"/results_last_run/infection_state_per_age_group/p50/Results.h5", 'r')
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

def infer_positive_tests(path):
    # First way: we just take x amount of eahc compartment and fit this to the positive tested category.
    # A few assumptions:
    # 1. The specificities of the tests are 99.9% this means that we have 0.1% false positives. We just assme, that nonsymptomatics test themselves very rarely, e.g. we just take a fraction of symptomatic persons.
    # 2. We just have a around 60% sensitivity, this means that we have 40% false negatives.
    # 3. We assume that it is way more likely to test yourself if you are symptomatic than if you are asymptomatic

    # we need every compartment of the model
    f_p50 = h5py.File(
        path+"/infection_state_per_age_group/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]
    total_50 = total_50[::24]
    total_50 = total_50[0:90]

    time = p50_bs['Time'][()]
    time = time[::24]
    time = time[0:90]


    # we need the real data from the json file cases_all_county_age_repdate_ma7.json
    df_abb = pd.read_json(
        path+"/../pydata/Germany/cases_infected_county_repdate_ma7.json")
    # we just need the columns cases and date
    df_abb = df_abb[['Date', 'Confirmed', 'ID_County']]
    # we need just the dates bewteen 2021-03-01 and 2021-06-01
    df_abb = df_abb[(df_abb['Date'] >= '2021-03-01') &
                    (df_abb['Date'] <= '2021-06-01')]
    # we just take the first 90 days
    df_abb = df_abb[0:90]
    # we need the amount of new positive tests each day insetad of cumulative
    df_abb['Confirmed'] = df_abb['Confirmed'].diff()
  
    sensitivity = 0.69
    specificity = 0.95
    # we need to derive the lowest rmse for the real positive tests and the inferred positive tests we use a grid search
    # we need to find the best r_sns and lt_sympt
    best_rmse = 1000000000
    best_r_sns = 0
    best_lt_sympt = 0
    for r_sns in np.linspace(1, 100, 100):
        for lt_sympt in np.linspace(0.005, 0.3, 500):
            total_positive_tests, inferred_positive_tests_sympt, inferred_positive_tests_asympt = calc_positive_tests_overall(total_50, sensitivity, specificity, r_sns, lt_sympt)
            rmse = np.sqrt(((df_abb['Confirmed'] - total_positive_tests)**2).mean())
            if rmse < best_rmse:
                best_rmse = rmse
                best_r_sns = r_sns
                best_lt_sympt = lt_sympt
    
    print("Best RMSE: ", best_rmse)
    print("Best r_sns: ", best_r_sns)
    print("Best lt_sympt: ", best_lt_sympt)


    total_positive_tests, inferred_positive_tests_sympt, inferred_positive_tests_asympt = calc_positive_tests_overall(total_50, sensitivity, specificity, best_r_sns, best_lt_sympt)


    # we save the assumed tests done
    # assumed__amount_of_test = (total_50[:, 3]*lt_sympt+total_50[:, 4]*lt_sympt+total_50[:, 5]*lt_sympt)+(
    #     total_50[:, 0]*lt_asympt)+(total_50[:, 1]+total_50[:, 2])*lt_asympt

    

    # we plot this
    # we plot the tests positive and the real cases
    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in time]
    xx = [xx[i].strftime('%Y-%m-%d') for i in range(len(xx))]

    plt.gca().set_xticks(time[::5])
    plt.gca().set_xticklabels(xx[::5])
    plt.gcf().autofmt_xdate()



    plt.plot(xx, inferred_positive_tests_sympt, color='tab:red')
    plt.plot(xx, inferred_positive_tests_asympt, color='tab:blue')
    plt.plot(xx, total_positive_tests, color='tab:green')
    plt.plot(xx, df_abb['Confirmed'], color='tab:orange')
    plt.xlabel('time (days)')
    plt.ylabel('Number of positive tests')
    plt.legend(['Assumed positive from Symptomatic', 'Assumed positive from Asymptomatic',
              'Assumed positive from all Persons', 'Real positive tests'])
 
    # we also write calculated best rmse, r_sns and lt_sympt into the title
    plt.title('Positive tests inferred from model and real positive tests, RMSE: '+str(float("{:.2f}".format(best_rmse)))+' r_sns: '+str(float("{:.2f}".format(best_r_sns)))+' lt_sympt: '+str(float("{:.3f}".format(best_lt_sympt))))
    
    

    plt.show()

def plot_estimated_reproduction_number(path):
    f_p50 = h5py.File(
        path+"/estimated_reproduction_number/p50/Results.h5", 'r')
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
    f_p50 = h5py.File(
        path+"/cumulative_detected_infections/p50/Results.h5", 'r')
    p50_bs = f_p50['0']
    total_50 = p50_bs['Total'][()]
    total_50 = total_50[::24]
    total_50 = total_50[0:90].flatten()
    # we smooth this with a gaussian filter
    total_50 = gaussian_filter1d(total_50, sigma=1, mode='nearest')
    time = p50_bs['Time'][()]
    time = time[::24]
    time = time[0:90]

    # we also plot the amount of new detected infections with diff
    total_50_diff = np.diff(total_50)


    # we plot this
    # we plot the tests positive and the real cases
    plt.plot(time, total_50, color='tab:red')
    plt.xlabel('time (days)')
    plt.ylabel('Cumulative Amount of detected infections')
    plt.title('Cumulative detected infections')
    # also plot the new detected infections with the same color but dashed
    plt.plot(time[1:], total_50_diff, color='tab:red', linestyle='dashed')
    plt.show()

def plot_positive_and_done_test(path):
    f_p50_positive = h5py.File(
        path+"/positive_test_per_location_type_per_age_group/p50/Results.h5", 'r')
    p50_bs_positive = f_p50_positive['0']
    total_50_positive = p50_bs_positive['Total'][()]

    f_p50_done = h5py.File(
        path+"/test_per_location_type_per_age_group/p50/Results.h5", 'r')
    p50_bs_done = f_p50_done['0']
    total_50_done = p50_bs_done['Total'][()]

    time = p50_bs_positive['Time'][()]
    time = time[::24]
    time = time[0:90]


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
    total_50_done = gaussian_filter1d(total_50_done, sigma=1, mode='nearest')

    # we plot this
    # we plot the tests positive and the real cases
    start_date = datetime.strptime('2021-03-01', '%Y-%m-%d')
    xx = [start_date + pd.Timedelta(days=int(i)) for i in range(90)]
    xx = [xx[i].strftime('%Y-%m-%d') for i in range(len(xx))]
    plt.gcf().autofmt_xdate()

    plt.plot(xx, total_50_positive, color='tab:red')
    plt.plot(xx, total_50_done, color='tab:blue')
    plt.xlabel('time (days)')
    plt.ylabel('Number of tests')
    plt.legend(['Positive tests', 'Done tests'])
    plt.title('Positive and done tests')
    plt.show()




if __name__ == "__main__":
    # path = "/Users/david/Documents/HZI/memilio/data/results_last_run_last_run"
    path = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/results_last_run"
    # path = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/cluster_results/1/results_last_run"
    # path = r"C:\Users\korf_sa\Documents\rep\data\results_last_run"

    if (len(sys.argv) > 1):
        n_runs = sys.argv[1]
    else:
        n_runs = len([entry for entry in os.listdir(path)
                     if os.path.isfile(os.path.join(path, entry))])
    # plot_infection_states_results(path)
    # plot_infections_loc_types_avarage(path)
    # plot_icu(path)
    # plot_dead(path)
    # infer_positive_tests(path)
    # plot_estimated_reproduction_number(path)
    # plot_cumulative_detected_infections(path)
    plot_positive_and_done_test(path)
