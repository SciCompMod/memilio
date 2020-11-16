import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pycode.epidemiology.epidata.getDIVIData as getDIVIData
import pycode.epidemiology.epidata.getRKIData as getRKIData
import datetime

"""function to get RKI and DIVI data on a specific given day, json files with data are downloaded in file 'Germany' if not downloaded before
additional one can plot RKI and DIVI data- just run this file."""

#yesterday = pd.Timestamp(datetime.date.today()) - pd.DateOffset(days=1)
yesterday = pd.Timestamp('2020.10.28')

# if True, plots will be saved in folder RKI_plots
saveplot_RKI_and_DIVI = False


def get_day_confirmed(day=pd.Timestamp('2020.02.02')):
    """get confirmed value of whole Germany on input day in format [confirmedvalue]"""
    if day < pd.Timestamp('2019.12.31'):
        return [0]
    if not os.path.exists('Germany') or not os.path.exists('Germany/infected_rki.json'):
        print('Download RKI Data from the Internet, takes some time')
        getRKIData.get_rki_data(out_folder='../examples', make_plot=False)
    df = pd.read_json("Germany/infected_rki.json")
    start_date = day
    end_date = day + pd.DateOffset(days=1)
    mask = (df['Date'] > start_date) & (df['Date'] <= end_date)
    if df.loc[mask].empty:
        #on this day, no additional people got confirmed, so one have to take confirmed value from day before
        start_date = start_date - pd.DateOffset(days=1)
        return get_day_confirmed(day=start_date)
    return df.loc[mask].iloc[:, 1].values


def get_timeperiod_RKI(startday=pd.Timestamp('2020.02.02'), endday=pd.Timestamp('2020.02.20'), comp='Confirmed'):
    """get confirmed or death numbers in Germany from days startdate to enddate in dataframe format"""
    if comp == 'Death':
        if not os.path.exists('Germany') or not os.path.exists('Germany/deaths_rki.json'):
            print('Download RKI Data from the Internet, takes some time')
            getRKIData.get_rki_data(out_folder='../examples', make_plot=False)
        df = pd.read_json('Germany/deaths_rki.json')
    else:
        if not os.path.exists('Germany') or not os.path.exists('Germany/infected_rki.json'):
            print('Download RKI Data from the Internet, takes some time')
            getRKIData.get_rki_data(out_folder='../examples', make_plot=False)
        df = pd.read_json("Germany/infected_rki.json")
    start_date = startday
    end_date = endday
    mask = (df['Date'] >= start_date) & (df['Date'] <= end_date)
    if len(df.loc[mask].iloc[:, 1].values) < (endday - startday).days:
        list = []
        for i in range((endday - startday).days):
            if comp == 'Death':
                [app] = get_day_death(day=start_date + pd.DateOffset(days=i))
            else:
                [app] = get_day_confirmed(day=start_date + pd.DateOffset(days=i))
            list.append(app)
        return list
    return df.loc[mask].iloc[:, 1].values


def get_day_death(day=pd.Timestamp('2020.02.02')):
    """get death value of whole Germany on input day in format [deathvalue]"""
    if day < pd.Timestamp('2020.01.16'):
        return [0]
    if not os.path.exists('Germany') or not os.path.exists('Germany/deaths_rki.json'):
        print('Download RKI Data from the Internet, takes some time')
        getRKIData.get_rki_data(out_folder='../examples', make_plot=False)
    df = pd.read_json("Germany/deaths_rki.json")
    start_date = day
    end_date = day + pd.DateOffset(days=1)
    mask = (df['Date'] > start_date) & (df['Date'] <= end_date)
    if df.loc[mask].empty:
        start_date = start_date - pd.DateOffset(days=1)
        return get_day_death(day=start_date)
    return df.loc[mask].iloc[:, 1].values


def plot_RKI(daystart=pd.Timestamp('2020.01.28'), simulationperiod=100, comp='Confirmed'):
    '''prints RKI Data Beginning on daystart for simulationperiod days in given compartment
    options for comp: 'Confirmed','Death'''
    if not (daystart + pd.DateOffset(days=simulationperiod) <= yesterday):
        simulationperiod = yesterday.dayofyear - daystart.dayofyear
    t = range(simulationperiod)
    list = get_timeperiod_RKI(startday=daystart, endday=daystart + pd.DateOffset(days=simulationperiod), comp=comp)

    list = np.asarray(list)
    datelist = np.array(
        pd.date_range(daystart.date(), periods=simulationperiod, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int(simulationperiod / 10) + 1) * 10)
    tick_range[-1] -= 1
    fig, ax = plt.subplots()

    ax.plot(t, list[:])
    ax.set_title(comp)
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45)
    ax.set(ylabel='number of people')
    ax.grid(linestyle='dotted')
    if saveplot_RKI_and_DIVI:
        fig.savefig('RKI_plots/' + 'RKI_plot_until_' + yesterday.strftime(
            "%d.%m.%Y") + '_' + comp + '.png')


def plot_RKI_last7days(daystart=pd.Timestamp('2020.01.28'), simulationperiod=100, comp='Confirmed'):
    '''function plotting on each day the the average Confirmed RKI-Data from last 7 days'''
    if not(daystart + pd.DateOffset(days=simulationperiod) <= yesterday):
        simulationperiod = yesterday.dayofyear - daystart.dayofyear
    t = range(simulationperiod - 6)

    list = get_timeperiod_RKI(startday=daystart, endday=daystart + pd.DateOffset(days=simulationperiod), comp=comp)
    list = np.asarray(list)
    listerg = np.concatenate((list, [0] * 6))
    for i in range(6):
        list = np.concatenate((np.array([0]), list))
        listerg[0:len(list)] = listerg[0:len(list)] + list
    listerg = listerg / 7
    erg = listerg[6:]
    erg = erg[0:len(erg) - 6]

    datelist = np.array(
        pd.date_range((daystart + pd.DateOffset(days=6)).date(), periods=simulationperiod - 6, freq='D').strftime(
            '%m-%d').tolist())
    tick_range = (np.arange(int((simulationperiod - 6 )/ 10) + 1) * 10)
    tick_range[-1] -= 1
    fig, ax = plt.subplots()

    ax.plot(t, erg[:])
    ax.set_title(comp+' smoothed (average of last 7 days)')
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45)
    ax.set(ylabel='number of people')
    ax.grid(linestyle='dotted')
    if saveplot_RKI_and_DIVI:
        fig.savefig('RKI_plots/' + 'RKI_smoothed_plot_until_' + yesterday.strftime(
            "%d.%m.%Y") +"_"+ comp+ '.png')


# ageresolved:

def get_day_age(day=pd.Timestamp('2020.02.02'), comp='Confirmed', interpol_unknown=True,search_for_last=True):
    """get age resolved values for comp=Confirmed or Death on input day in format [value] with
     'unknown' age distributed to other age groups if interpol_unknown true
     if search_for_last is True, then one will get the last change of specific age groups, to have real age data
      (for the plot this is very time intensive for every single day, so we have the timeperiod function to solve this problem)"""
    # RKI age groups: 0-4; 5-14; 15-34;35-59,60-79,80+
    if day <= pd.Timestamp('2020.01.01'):
        return np.array([0] * 6)
    if not os.path.exists('Germany') or not os.path.exists('Germany/all_age_rki.json'):
        print('Download RKI Data from the Internet, takes some time')
        getRKIData.get_rki_data(out_folder='../examples', make_plot=False)
    df = pd.read_json("Germany/all_age_rki.json")
    start_date = day
    end_date = day + pd.DateOffset(days=1)
    mask = (df['Date'] > start_date) & (df['Date'] <= end_date)

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
            datadaybefore = get_day_age(day=daybefore, comp=comp, interpol_unknown=interpol_unknown)
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
        new = new + (ret[6] / 6)
        ret = new

    return ret

def get_timeperiod_age_RKI(daystart=pd.Timestamp('2020.01.28'), simulationperiod=100, comp='Confirmed'):
    '''get a time period of age resolved data, to make plot_age more efficient compared to get every single day in
    plot_age and always search for the last change in specific age group'''
    list = []
    t = range(simulationperiod)
    for i in t:
        if i==0:
            ret = get_day_age(day=daystart + pd.DateOffset(days=i), comp=comp, interpol_unknown=True,search_for_last=True)
        else:
            ret = get_day_age(day=daystart + pd.DateOffset(days=i), comp=comp, interpol_unknown=False,search_for_last=False)
            if 0 in ret:
                iter=len(ret)
                if iter==7:
                    iter=6
                for j in range(iter):
                    if ret[j]==0:
                        ret[j]=list[-1][j]
            if len(ret) == 7:
                new = np.array([0] * 6)
                new[0:6] = ret[0:6]
                new = new + (ret[6] / 6)
                ret = new
        list.append(ret.tolist())
    return list

def plot_age(daystart=pd.Timestamp('2020.01.28'), simulationperiod=100, comp='Confirmed'):
    ages = ["[0,5)", "[5,15)", "[15,35)", "[35,60)", "[60,80)", "80+"]
    if daystart + pd.DateOffset(days=simulationperiod) <= yesterday:
        t = range(simulationperiod)
    else:
        simulationperiod = yesterday.dayofyear - daystart.dayofyear
        t = range(simulationperiod)
    list=get_timeperiod_age_RKI(daystart=daystart, simulationperiod=simulationperiod, comp=comp)

    list = np.asarray(list)
    # list=np.transpose(list)
    datelist = np.array(
        pd.date_range(daystart.date(), periods=simulationperiod, freq='D').strftime('%m-%d').tolist())
    tick_range = (np.arange(int(simulationperiod / 10) + 1) * 10)
    tick_range[-1] -= 1
    fig, ax = plt.subplots()

    for age in range(len(ages)):
        ax.plot(t, list[:, age], label=ages[age])
    ax.set_title(comp)
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45)
    ax.set(ylabel='number of people')
    ax.legend()
    ax.grid(linestyle='dotted')
    if saveplot_RKI_and_DIVI:
        fig.savefig('RKI_plots/RKI_ageresolved_plot_until_' + yesterday.strftime(
            "%d.%m.%Y") + '_' + comp + '.png')


# countyresolved
def get_day_county(day=pd.Timestamp('2020.02.02'), county='SK Köln', comp='Confirmed'):
    """get values on date 'day' from compartement 'comp' (options: 'Confirmed' , 'Death') in county 'county' ,output in format [value]"""
    if not os.path.exists('Germany') or not os.path.exists('Germany/all_county_rki.json'):
        print('Download RKI Data from the Internet, takes some time')
        getRKIData.get_rki_data(out_folder='../examples', make_plot=False)
    df = pd.read_json("Germany/all_county_rki.json")
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
def get_day_DIVI(day=pd.Timestamp('2020.03.27')):
    """get Divi data ICU on input day, data starting on 24.4.2020 """
    if day < pd.Timestamp('2020.04.24'):
        print("error: DIVI dataset starts on 24.04.2020; you asked for " + day.strftime(
            "%d.%m.%Y") + " please search for another solution")
        return np.array([0])
    else:
        if not os.path.exists('Germany') or not os.path.exists('Germany/germany_divi.json'):
            print('Download DIVI Data from the Internet, takes some time')
            getDIVIData.get_divi_data(out_folder='../examples', end_date=yesterday)
        df = pd.read_json("Germany/germany_divi.json")
        start_date = day
        end_date = day + pd.DateOffset(days=1)
        mask = (df['Date'] > start_date) & (df['Date'] <= end_date)
        if df.loc[mask].empty:
            start_date = start_date - pd.DateOffset(days=1)
            return get_day_DIVI(day=start_date)
        return df.loc[mask].iloc[:, 1].values


def plot_DIVI(daystart=pd.Timestamp('2020.04.25'), simulationperiod=100):
    '''prints DIVI Data Beginning on daystart for simulationperiod days '''
    list = []
    if daystart + pd.DateOffset(days=simulationperiod) <= yesterday:
        t = range(simulationperiod + 1)
    else:
        simulationperiod = yesterday.dayofyear - daystart.dayofyear
        t = range(simulationperiod + 1)
    for i in t:
        ret = get_day_DIVI(day=daystart + pd.DateOffset(days=i))
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
    if saveplot_RKI_and_DIVI:
        fig.savefig('RKI_plots/DIVI_plot_until_' + yesterday.strftime("%d.%m.%Y") + '.png')


def get_day_countyDIVI(day=pd.Timestamp('2020.03.27'), county='SK Köln'):
    """get Divi data ICU on input day ,data starting on 24.4.2020"""
    if day < pd.Timestamp('2020.04.24'):
        print("error: DIVI dataset starts on 24.04.2020; you asked for " + day.strftime(
            "%d.%m.%Y") + " please search for another solution")
    else:
        if not os.path.exists('Germany') or not os.path.exists('Germany/county_divi.json'):
            print('Download DIVI Data from the Internet, takes some time')
            getDIVIData.get_divi_data(out_folder='../examples', end_date=yesterday)
        df = pd.read_json("Germany/county_divi.json")
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

    simulationperiod=300
    if True:
        plot_RKI_last7days(simulationperiod=simulationperiod)
        plot_RKI_last7days(simulationperiod=simulationperiod,comp='Death')
        plot_DIVI(simulationperiod=simulationperiod)
        plot_age(simulationperiod=simulationperiod)
        plot_age(simulationperiod=simulationperiod, comp='Death')
        plot_RKI(daystart=pd.Timestamp('2019.12.30'),simulationperiod=simulationperiod)
        plot_RKI(simulationperiod=simulationperiod, comp='Death')
    plt.show()
