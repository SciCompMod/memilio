import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import datetime

import pycode.epidemiology.epidata.getDIVIData as getDIVIData
import pycode.epidemiology.epidata.getRKIData as getRKIData
import pycode.epidemiology.epidata.defaultDict as dd

"""function to get RKI and DIVI data on a specific given day, json files with data are downloaded in file 'Germany' if not downloaded before
additional one can plot RKI and DIVI data and age resolved incidence-values- just run this file."""

yesterday = pd.Timestamp(datetime.date.today()) - pd.DateOffset(days=1)
endday_divi = pd.Timestamp('2021.01.04')

# if True, plots will be saved in folder RKI_plots
saveplot_RKI_and_DIVI = False
# if you already downloaded RKI and DIVI-Data in Germany folder you can just read data,else it will be downloaded
read_data_global = False
path_for_data = os.path.join(dd.defaultDict['out_folder'], 'Germany/')
# population data from Statistisches Bundesamt 31.12.2019 as RKI used since 08.10.2020
# age_group_population_total=[3961376,7429883,19117865,28919134,18057318,5681135]
# population data from 31.12.2018 as RKI used from beginning until 08.10.2020
age_group_population_total = [3926397, 7364418, 19213113, 29137839, 17988340, 5389106]


def get_day_confirmed(day=pd.Timestamp('2020.02.02'), read_data=False):
    """get confirmed value of whole Germany on input day in format [confirmedvalue]"""
    if day < pd.Timestamp('2019.12.31'):
        return [0]
    if not read_data:
        print('Download RKI Data from the Internet, takes some time')
        getRKIData.get_rki_data()
    df = pd.read_json(path_for_data + "infected_rki.json")
    start_date = day
    mask = (df['Date'] == start_date)
    if df.loc[mask].empty:
        # on this day, no additional people got confirmed, so one have to take confirmed value from day before
        start_date = start_date - pd.DateOffset(days=1)
        return get_day_confirmed(day=start_date,read_data=True)
    return df.loc[mask].iloc[:, 1].values


def get_timeperiod_RKI(startday=pd.Timestamp('2020.02.02'), endday=pd.Timestamp('2020.02.20'), comp='Confirmed',
                       read_data=False):
    """get confirmed or death numbers in Germany from days startdate to enddate in dataframe format"""
    if not read_data:
        print('Download RKI Data from the Internet, takes some time')
        getRKIData.get_rki_data()
    if comp == 'Death':
        df = pd.read_json(path_for_data + 'deaths_rki.json')
    else:
        df = pd.read_json(path_for_data + "infected_rki.json")
    start_date = startday
    end_date = endday
    mask = (df['Date'] >= start_date) & (df['Date'] <= end_date)
    if len(df.loc[mask].iloc[:, 1].values) < (endday - startday).days:
        list = []
        for i in range((endday - startday).days):
            if comp == 'Death':
                [app] = get_day_death(day=start_date + pd.DateOffset(days=i),read_data=True)
            else:
                [app] = get_day_confirmed(day=start_date + pd.DateOffset(days=i),read_data=True)
            list.append(app)
        return list
    return df.loc[mask].iloc[:, 1].values


def get_day_death(day=pd.Timestamp('2020.02.02'), read_data=False):
    """get death value of whole Germany on input day in format [deathvalue]"""
    if day < pd.Timestamp('2020.01.16'):
        return [0]
    if not read_data:
        print('Download RKI Data from the Internet, takes some time')
        getRKIData.get_rki_data()
    df = pd.read_json(path_for_data + "deaths_rki.json")
    start_date = day
    mask = (df['Date'] == start_date)

    if df.loc[mask].empty:
        start_date = start_date - pd.DateOffset(days=1)
        return get_day_death(day=start_date,read_data=True)
    return df.loc[mask].iloc[:, 1].values


def plot_RKI(daystart=pd.Timestamp('2020.01.28'), simulationperiod=100, comp='Confirmed', read_data=False):
    '''prints RKI Data Beginning on daystart for simulationperiod days in given compartment
    options for comp: 'Confirmed','Death'''
    if not (daystart + pd.DateOffset(days=simulationperiod) <= yesterday):
        simulationperiod = (yesterday - daystart).days
    list = get_timeperiod_RKI(startday=daystart, endday=daystart + pd.DateOffset(days=simulationperiod), comp=comp,
                              read_data=read_data)

    list = np.asarray(list)
    t = range(len(list))
    datelist = np.array(
        pd.date_range(daystart.date(), periods=simulationperiod, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int(simulationperiod / 10) + 1) * 10)
    tick_range[-1] -= 1
    fig, ax = plt.subplots()

    ax.plot(t, list[:])
    ax.set_title(comp + " RKI-Data all ages")
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45)
    ax.set(ylabel='number of people')
    ax.grid(linestyle='dotted')
    plt.show()
    if saveplot_RKI_and_DIVI:
        fig.savefig('RKI_plots/' + 'RKI_plot_until_' + yesterday.strftime(
            "%d.%m.%Y") + '_' + comp + '.png')


def plot_RKI_last7days(daystart=pd.Timestamp('2020.01.28'), simulationperiod=100, comp='Confirmed', read_data=False):
    ''' function plotting on each day the the average Confirmed RKI-Data from last 7 days'''
    if not (daystart + pd.DateOffset(days=simulationperiod) <= yesterday):
        simulationperiod = (yesterday - daystart).days

    list = get_timeperiod_RKI(startday=daystart, endday=daystart + pd.DateOffset(days=simulationperiod), comp=comp,
                              read_data=read_data)
    list = np.asarray(list)
    listerg = np.concatenate((list, [0] * 6))
    for i in range(6):
        list = np.concatenate((np.array([0]), list))
        listerg[0:len(list)] = listerg[0:len(list)] + list
    listerg = listerg / 7
    erg = listerg[6:]
    erg = erg[0:len(erg) - 6]

    t = range(len(erg))
    datelist = np.array(
        pd.date_range((daystart + pd.DateOffset(days=6)).date(), periods=simulationperiod - 6, freq='D').strftime(
            '%y-%m-%d').tolist())
    tick_range = (np.arange(int((simulationperiod - 6) / 10) + 1) * 10)
    tick_range[-1] -= 1
    fig, ax = plt.subplots()

    ax.plot(t, erg[:])
    ax.set_title(comp + ' smoothed (average of last 7 days)')
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45)
    ax.set(ylabel='number of people')
    ax.grid(linestyle='dotted')
    plt.show()
    if saveplot_RKI_and_DIVI:
        fig.savefig('RKI_plots/' + 'RKI_smoothed_plot_until_' + yesterday.strftime(
            "%d.%m.%Y") + "_" + comp + '.png')


# ageresolved:

def get_day_age(day=pd.Timestamp('2020.02.02'), comp='Confirmed', read_data=False, interpol_unknown=True,
                search_for_last=True):
    """get age resolved values for comp=Confirmed or Death on input day in format [value] with
     'unknown' age distributed to other age groups if interpol_unknown true
     if search_for_last is True, then one will get the last change of specific age groups, to have real age data
      (for the plot this is very time intensive for every single day, so we have the timeperiod function to solve this problem)"""
    # RKI age groups: 0-4; 5-14; 15-34;35-59,60-79,80+
    if day <= pd.Timestamp('2020.01.01'):
        return np.array([0] * 6)
    if not read_data:
        print('Download RKI Data from the Internet, takes some time')
        getRKIData.get_rki_data()
    df = pd.read_json(path_for_data + "all_age_rki.json")
    mask = (df['Date'] == day)

    if comp == 'Death':
        ret = df.loc[mask].iloc[:, 3].values
    else:
        ret = df.loc[mask].iloc[:, 2].values

    if len(ret) < 6:
        # fill with zeros to always return an array length 6
        new = np.array([0] * 6)
        ages = ['A00-A04', 'A05-A14', 'A15-A34', 'A35-A59', 'A60-A79', 'A80+']
        if search_for_last:
            daybefore = day - pd.DateOffset(days=1)
            datadaybefore = get_day_age(day=daybefore, comp=comp,read_data=read_data, interpol_unknown=interpol_unknown)
        j = 0
        for i in range(6):
            if not df.loc[mask & (df['Age_RKI'] == ages[i])].empty:
                new[i] = ret[j]
                j = j + 1
            elif search_for_last:
                new[i] = datadaybefore[i]
        if not j == len(ret):
            print('Something went wrong in get_RKI_into_panda.py function: get_day_age')
        return new

    if len(ret) == 7 and interpol_unknown:
        new = np.array([0] * 6)
        new[0:6] = ret[0:6]
        help = np.array([(1 / sum(age_group_population_total))] * len(age_group_population_total))
        new = new + age_group_population_total * help * ret[6]
        ret = new

    return ret


def get_timeperiod_age_RKI(daystart=pd.Timestamp('2020.01.28'), simulationperiod=100, comp='Confirmed',
                           read_data=False):
    '''get a time period of age resolved data, to make plot_age more efficient compared to get every single day in
    plot_age and always search for the last change in specific age group'''
    list = []
    t = range(simulationperiod)
    for i in t:
        if i == 0:
            ret = get_day_age(day=daystart, comp=comp, read_data=read_data, interpol_unknown=True, search_for_last=True)
            read_data=True
        else:
            ret = get_day_age(day=daystart + pd.DateOffset(days=i), comp=comp, read_data=read_data, interpol_unknown=False,
                              search_for_last=False)
            if 0 in ret:
                iter = len(ret)
                if iter == 7:
                    iter = 6
                for j in range(iter):
                    if ret[j] == 0:
                        ret[j] = list[-1][j]
            if len(ret) == 7:
                new = np.array([0] * 6)
                new[0:6] = ret[0:6]
                help = np.array([(1 / sum(age_group_population_total))] * len(age_group_population_total))
                new = new + age_group_population_total * help * ret[6]
                ret = new
        list.append(ret.tolist())
    return list


def plot_age(daystart=pd.Timestamp('2020.01.28'), simulationperiod=100, comp='Confirmed', read_data=False):
    ages = ["[0,5)", "[5,15)", "[15,35)", "[35,60)", "[60,80)", "80+"]
    if not (daystart + pd.DateOffset(days=simulationperiod) <= yesterday):
        simulationperiod = (yesterday - daystart).days
    list = get_timeperiod_age_RKI(daystart=daystart, simulationperiod=simulationperiod, comp=comp, read_data=read_data)

    list = np.asarray(list)
    t = range(len(list))
    datelist = np.array(
        pd.date_range(daystart.date(), periods=simulationperiod, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int(simulationperiod / 10) + 1) * 10)
    tick_range[-1] -= 1
    fig, ax = plt.subplots()

    for age in range(len(ages)):
        ax.plot(t, list[:, age], label=ages[age])
    ax.set_title(comp + " age resolved RKI-Data")
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45)
    ax.set(ylabel='number of people')
    ax.legend()
    ax.grid(linestyle='dotted')
    plt.show()
    if saveplot_RKI_and_DIVI:
        fig.savefig('RKI_plots/RKI_ageresolved_plot_until_' + yesterday.strftime(
            "%d.%m.%Y") + '_' + comp + '.png')


def plot_age_incidence(daystart=pd.Timestamp('2020.01.28'), simulationweeks=28, read_data=False):
    '''one can compare these incidence plot with RKI report from 2020-09-15 '''
    simulationperiod = simulationweeks * 7

    ages = ["[0,5)", "[5,15)", "[15,35)", "[35,60)", "[60,80)", "80+"]
    if not (daystart + pd.DateOffset(days=simulationperiod) <= yesterday):
        simulationperiod = yesterday.dayofyear - daystart.dayofyear

    list = get_timeperiod_age_RKI(daystart=daystart, simulationperiod=simulationperiod, comp='Confirmed',
                                  read_data=read_data)

    list = np.asarray(list)
    t = np.arange(daystart.isocalendar()[1], daystart.isocalendar()[1] + simulationweeks)
    result = np.zeros((simulationweeks, len(ages)))
    for i in range(simulationweeks):
        # calculate new reported cases per week per age group
        result[i, :] = list[(i + 1) * 7 - 1, :] - list[i * 7, :]
    for j in range(len(ages)):
        # new reported cases per 100000 residents
        result[:, j] = result[:, j] * (100000 / age_group_population_total[j])

    fig, ax = plt.subplots()

    for age in range(len(ages)):
        ax.plot(t, result[:, age], label=ages[age])
    ax.set_title("Icidence-values")
    ax.set(ylabel='new reported cases per 100000 residents in calendar week', xlabel='calendar week')
    ax.legend()
    ax.grid(linestyle='dotted')
    plt.show()
    if saveplot_RKI_and_DIVI:
        fig.savefig('RKI_plots/RKI_ageresolved_plot_until_' + yesterday.strftime(
            "%d.%m.%Y") + '_incidence.png')


def plot_age_incidence_daily(daystart=pd.Timestamp('2020.01.28'), simulationperiod=200, read_data=False):
    ages = ["[0,5)", "[5,15)", "[15,35)", "[35,60)", "[60,80)", "80+"]
    if not (daystart + pd.DateOffset(days=simulationperiod) <= yesterday):
        simulationperiod = (yesterday - daystart).days

    list = get_timeperiod_age_RKI(daystart=daystart - pd.DateOffset(days=7), simulationperiod=simulationperiod + 7,
                                  comp='Confirmed', read_data=read_data)
    list = np.asarray(list)
    t = range(len(list) - 7)
    result = np.zeros((simulationperiod, len(ages)))
    for i in range(simulationperiod):
        # calculate new reported cases per week per age group
        result[i, :] = list[i + 6, :] - list[i, :]
    for j in range(len(ages)):
        # new reported cases per 100000 residents
        result[:, j] = result[:, j] * (100000 / age_group_population_total[j])
    datelist = np.array(
        pd.date_range(daystart.date(), periods=simulationperiod, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int(simulationperiod / 10) + 1) * 10)
    tick_range[-1] -= 1
    fig, ax = plt.subplots()

    for age in range(len(ages)):
        ax.plot(t, result[:, age], label=ages[age])
    ax.set_title("Icidence-values daily")
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45)
    ax.set(ylabel='new reported cases per 100000 residents in calendar week', xlabel='calendar week')
    ax.legend()
    ax.grid(linestyle='dotted')
    plt.show()
    if saveplot_RKI_and_DIVI:
        fig.savefig('RKI_plots/RKI_ageresolved_plot_until_' + yesterday.strftime(
            "%d.%m.%Y") + '_incidence_daily.png')


# countyresolved
def get_day_county(day=pd.Timestamp('2020.02.02'), county='SK Köln', comp='Confirmed',read_data=False):
    """get values on date 'day' from compartement 'comp' (options: 'Confirmed' , 'Death') in county 'county' ,output in format [value]"""
    if not read_data:
        print('Download RKI Data from the Internet, takes some time')
        getRKIData.get_rki_data()
    df = pd.read_json(path_for_data+"all_county_rki.json")
    start_date = day
    end_date = day + pd.DateOffset(days=1)
    mask = (df['Date'] > start_date) & (df['Date'] <= end_date) & (df['County'] == county)
    if df.loc[mask].empty:
        print("warning: There are no RKI Data for " + county + " on " + day.strftime("%d.%m.%Y"))
        return np.array([0])
    if comp == 'Death':
        return df.loc[mask].iloc[:, 4].values
    else:
        return df.loc[mask].iloc[:, 3].values


# get DIVI data
def get_day_DIVI(day=pd.Timestamp('2020.03.27'),read_data=False):
    """get Divi data ICU on input day, data starting on 24.4.2020 """
    if day < pd.Timestamp('2020.04.24'):
        print("error: DIVI dataset starts on 24.04.2020; you asked for " + day.strftime(
            "%d.%m.%Y") + " please search for another solution")
        return np.array([0])
    else:
        if not read_data:
            print('Download DIVI Data from the Internet, takes some time')
            getDIVIData.get_divi_data(end_date=endday_divi)
        df = pd.read_json(path_for_data+"germany_divi.json")
        start_date = day
        end_date = day + pd.DateOffset(days=1)
        mask = (df['Date'] > start_date) & (df['Date'] <= end_date)
        if df.loc[mask].empty:
            start_date = start_date - pd.DateOffset(days=1)
            return get_day_DIVI(day=start_date,read_data=True)
        return df.loc[mask].iloc[:, 1].values


def plot_DIVI(daystart=pd.Timestamp('2020.04.25'), simulationperiod=100,read_data=False):
    '''prints DIVI Data Beginning on daystart for simulationperiod days '''
    list = []
    if daystart + pd.DateOffset(days=simulationperiod) <= endday_divi:
        t = range(simulationperiod + 1)
    else:
        simulationperiod = (endday_divi - daystart).days
        t = range(simulationperiod + 1)
    for i in t:
        ret = get_day_DIVI(day=daystart + pd.DateOffset(days=i),read_data=read_data)
        read_data=True
        list.append(ret.tolist())

    list = np.asarray(list)
    datelist = np.array(
        pd.date_range(daystart.date(), periods=simulationperiod, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int(simulationperiod / 10) + 1) * 10)
    tick_range[-1] -= 1
    fig, ax = plt.subplots()

    ax.plot(t, list[:])
    ax.set_title('DIVI-Data: people in Intensive Care Units')
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45)
    ax.set(ylabel='number of people')
    ax.grid(linestyle='dotted')
    plt.show()
    if saveplot_RKI_and_DIVI:
        fig.savefig('RKI_plots/DIVI_plot_until_' + endday_divi.strftime("%d.%m.%Y") + '.png')


def get_day_countyDIVI(day=pd.Timestamp('2020.03.27'), county='SK Köln',read_data=False):
    """get Divi data ICU on input day ,data starting on 24.4.2020"""
    if day < pd.Timestamp('2020.04.24'):
        print("error: DIVI dataset starts on 24.04.2020; you asked for " + day.strftime(
            "%d.%m.%Y") + " please search for another solution")
    else:
        if not read_data:
            print('Download DIVI Data from the Internet, takes some time')
            getDIVIData.get_divi_data()
        df = pd.read_json(path_for_data+"county_divi.json")
        start_date = day
        end_date = day + pd.DateOffset(days=1)
        mask = (df['Date'] > start_date) & (df['Date'] <= end_date) & (df['County'] == county)
        if df.loc[mask].empty:
            return np.array([0])
        else:
            return df.loc[mask].iloc[:, 2].values


if __name__ == "__main__":
    if saveplot_RKI_and_DIVI:
        if not os.path.exists('RKI_plots'):
            os.makedirs('RKI_plots')

    simulationperiod = 400
    simulationweeks = 40
    daystart = pd.Timestamp('2020.03.02')
    read_data=read_data_global
    if True:
        plot_age_incidence(daystart=daystart, simulationweeks=simulationweeks,read_data=read_data)
        read_data=True
        plot_age_incidence_daily(daystart=daystart, simulationperiod=simulationweeks * 7,read_data=read_data)
    if True:
        plot_RKI_last7days(daystart=daystart, simulationperiod=simulationperiod,read_data=read_data)
        read_data=True
        plot_RKI_last7days(daystart=daystart, simulationperiod=simulationperiod, comp='Death',read_data=read_data)
        plot_RKI(daystart=daystart, simulationperiod=simulationperiod,read_data=read_data)
        plot_RKI(daystart=daystart, simulationperiod=simulationperiod, comp='Death',read_data=read_data)
        plot_age(daystart=daystart, simulationperiod=simulationperiod,read_data=read_data)
        plot_age(daystart=daystart, simulationperiod=simulationperiod, comp='Death',read_data=read_data)
    if True:
        read_data=read_data_global
        plot_DIVI(simulationperiod=simulationperiod,read_data=read_data)
