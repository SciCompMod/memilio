"""@plot_RKI_data.py
function to get RKI and DIVI data plots
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import datetime

import epidemiology.epidata.getDIVIData as getDIVIData
import epidemiology.epidata.getRKIData as getRKIData
import epidemiology.epidata.defaultDict as dd

yesterday = pd.Timestamp(datetime.date.today()) - pd.DateOffset(days=1)

# folder where divi and rki data can be found after downloading with default out_folder
data_folder = os.path.join(dd.defaultDict['out_folder'],"pydata", "Germany")


def get_Data(endday_divi=yesterday, moving_average=True):
    """function to get DIVI and RKI data
    @param endday_divi divi data will be loaded until this day
    @param moving_average Defines if moving average is used"""
    print('Download RKI Data from the Internet, takes some time')
    getRKIData.get_rki_data(out_folder=dd.defaultDict['out_folder'], moving_average=moving_average)
    print('Download DIVI Data from the Internet, takes some time')
    getDIVIData.get_divi_data(out_folder=dd.defaultDict['out_folder'], end_date=endday_divi)


def plot_RKI(moving_average, daystart=pd.Timestamp('2020.01.28'), simulationperiod=100, saveplot=False):
    """plot RKI data: Cumulative Cofirmed Cases, Deaths and for each 7-days moving average
        @param moving_average Defines if moving-average is used
        @param daystart Day at which should be started in timestamp format
        @param simulationperiod number in integer format of days for which data should be plotted
        @param saveplot boolean value; says if plot should be saved """
    df = pd.read_json(os.path.join(data_folder, "infected_rki.json"))
    if not (daystart + pd.DateOffset(days=simulationperiod) <= yesterday):
        simulationperiod = (yesterday - daystart).days
    mask = (df['Date'] >= daystart) & (df['Date'] <= daystart + pd.DateOffset(days=simulationperiod))

    fig_name = 'RKI_confirmed'
    plt.figure(fig_name)
    plt.plot(df.loc[mask]['Date'], df.loc[mask]["Confirmed"])
    plt.title('Cumulative Confirmed cases RKI')
    plt.xlabel('Date')
    plt.xticks(rotation=25)
    plt.ylabel('Confirmed cases')
    plt.ylim(bottom=0)
    plt.grid(True, linestyle='--')
    plt.tight_layout()
    if saveplot:
        plt.savefig(os.path.join('Plots', fig_name+".png"))

    if moving_average:
        df = pd.read_json(os.path.join(data_folder, "infected_ma_rki.json"))
        mask = (df['Date'] >= daystart) & (df['Date'] <= daystart + pd.DateOffset(days=simulationperiod))
        fig_name = 'RKI_moving_average_confirmed'
        plt.figure(fig_name)
        plt.plot(df.loc[mask]['Date'], df.loc[mask]["Confirmed"])
        plt.title('7-days moving average cumulative Confirmed cases RKI')
        plt.xlabel('Date')
        plt.xticks(rotation=25)
        plt.ylabel('Confirmed cases')
        plt.ylim(bottom=0)
        plt.grid(True, linestyle='--')
        plt.tight_layout()
        if saveplot:
            plt.savefig('Plots/plot_RKI_ma_Confirmed.png')

    df = pd.read_json(os.path.join(data_folder, "deaths_rki.json"))
    mask = (df['Date'] >= daystart) & (df['Date'] <= daystart + pd.DateOffset(days=simulationperiod))

    fig_name = 'RKI_deaths'
    plt.figure(fig_name)
    plt.plot(df.loc[mask]['Date'], df.loc[mask]["Deaths"])
    plt.title('Cumulative Deaths RKI')
    plt.xlabel('Date')
    plt.xticks(rotation=25)
    plt.ylabel('Deaths')
    plt.ylim(bottom=0)
    plt.grid(True, linestyle='--')
    plt.tight_layout()
    if saveplot:
        plt.savefig(os.path.join('Plots', fig_name+".png"))

    if moving_average:
        df = pd.read_json(os.path.join(data_folder, "deaths_ma_rki.json"))
        mask = (df['Date'] >= daystart) & (df['Date'] <= daystart + pd.DateOffset(days=simulationperiod))
        fig_name = 'RKI_moving_average_deaths'
        plt.figure(fig_name)
        plt.plot(df.loc[mask]["Date"], df.loc[mask]["Deaths"])
        plt.title('7-days moving average cumulative Deaths cases RKI')
        plt.xlabel('Date')
        plt.xticks(rotation=25)
        plt.ylabel('Deaths')
        plt.ylim(bottom=0)
        plt.grid(True, linestyle='--')
        plt.tight_layout()
        if saveplot:
            plt.savefig(os.path.join('Plots', fig_name+".png"))


def plot_RKI_age(moving_average, daystart=pd.Timestamp('2020.01.28'), simulationperiod=100, saveplot=False):
    """plot RKI data age resolved: Cumulative Cofirmed Cases, Deaths and for each 7-days moving average
    @param moving_average Defines if moving-average is used
    @param daystart Day at which should be started in timestamp format
    @param simulationperiod number in integer format of days for which data should be plotted
    @param saveplot boolean value; says if plot should be saved """
    if not (daystart + pd.DateOffset(days=simulationperiod) <= yesterday):
        simulationperiod = (yesterday - daystart).days
    df = pd.read_json(os.path.join(data_folder, "all_age_rki.json"))
    mask = (df['Date'] >= daystart) & (df['Date'] <= daystart + pd.DateOffset(days=simulationperiod))
    ages = dd.age_rki_list
    fig_name = 'RKI_confirmed_ageresolved'
    plt.figure(fig_name)
    for value in ages:
        plt.plot(df.loc[mask & (df["Age_RKI"] == value)]['Date'], df.loc[mask & (df["Age_RKI"] == value)]["Confirmed"])
    plt.legend(ages)
    plt.title('cumulative Confirmed cases RKI age resolved')
    plt.xlabel('Date')
    plt.xticks(rotation=25)
    plt.ylabel('Confirmed cases')
    plt.ylim(bottom=0)
    plt.grid(True, linestyle='--')
    plt.tight_layout()
    if saveplot:
        plt.savefig(os.path.join('Plots', fig_name+".png"))

    fig_name = 'RKI_deaths_ageresolved'
    plt.figure(fig_name)
    for value in ages:
        plt.plot(df.loc[mask & (df["Age_RKI"] == value)]['Date'], df.loc[mask & (df["Age_RKI"] == value)]["Deaths"])
    plt.legend(ages)
    plt.title('cumulative deaths RKI age resolved')
    plt.xlabel('Date')
    plt.xticks(rotation=25)
    plt.ylabel('Deaths')
    plt.ylim(bottom=0)
    plt.grid(True, linestyle='--')
    plt.tight_layout()
    if saveplot:
        plt.savefig(os.path.join('Plots', fig_name+".png"))

    if moving_average:
        df = pd.read_json(os.path.join(data_folder, "all_age_ma_rki.json"))
        mask = (df['Date'] >= daystart) & (df['Date'] <= daystart + pd.DateOffset(days=simulationperiod))
        fig_name = 'RKI_ma_confirmed_ageresolved'
        plt.figure(fig_name)
        for value in ages:
            plt.plot(df.loc[mask & (df["Age_RKI"] == value)]['Date'], df.loc[mask & (df["Age_RKI"] == value)]["Confirmed"])
        plt.legend(ages)
        plt.title('7-days moving average cumulative Confirmed cases RKI age resolved')
        plt.xlabel('Date')
        plt.xticks(rotation=25)
        plt.ylabel('Confirmed cases')
        plt.ylim(bottom=0)
        plt.grid(True, linestyle='--')
        plt.tight_layout()
        if saveplot:
            plt.savefig(os.path.join('Plots', fig_name+".png"))

        fig_name = 'RKI_ma_deaths_ageresolved'
        plt.figure(fig_name)
        for value in ages:
            plt.plot(df.loc[mask & (df["Age_RKI"] == value)]['Date'], df.loc[mask & (df["Age_RKI"] == value)]["Deaths"])
        plt.legend(ages)
        plt.title('7-days moving average cumulative deaths RKI age resolved')
        plt.xlabel('Date')
        plt.xticks(rotation=25)
        plt.ylabel('Deaths')
        plt.ylim(bottom=0)
        plt.grid(True, linestyle='--')
        plt.tight_layout()
        if saveplot:
            plt.savefig(os.path.join('Plots', fig_name+".png"))


def plot_RKI_county(moving_average, daystart=pd.Timestamp('2020.01.28'), simulationperiod=100, county='SK Köln', saveplot=False):
    """plot RKI data for one county: Cumulative Cofirmed Cases, Deaths and for each 7-days moving average
    @param moving_average Defines if moving-average is used
    @param daystart Day at which should be started in timestamp format
    @param simulationperiod number in integer format of days for which data should be plotted
    @param county county, for which data chould be plotted
    @param saveplot boolean value; says if plot should be saved """
    if not (daystart + pd.DateOffset(days=simulationperiod) <= yesterday):
        simulationperiod = (yesterday - daystart).days
    df = pd.read_json(os.path.join(data_folder, "all_county_rki.json"))
    mask = (df['Date'] >= daystart) & (df['Date'] <= daystart + pd.DateOffset(days=simulationperiod)) & (
            df["County"] == county)
    fig_name = 'RKI_confirmed_county_' + county.replace(" ", "_")
    plt.figure(fig_name)
    plt.plot(df.loc[mask]["Date"], df.loc[mask]["Confirmed"])
    plt.title('cumulative Confirmed cases ' + county)
    plt.xlabel('Date')
    plt.xticks(rotation=25)
    plt.ylabel('Confirmed cases')
    plt.ylim(bottom=0)
    plt.grid(True, linestyle='--')
    plt.tight_layout()
    if saveplot:
        plt.savefig(os.path.join('Plots', fig_name+".png"))

    fig_name = 'RKI_deaths_county_' + county.replace(" ", "_")
    plt.figure(fig_name)
    plt.plot(df.loc[mask]["Date"], df.loc[mask]["Deaths"])
    plt.title('cumulative deaths ' + county)
    plt.xlabel('Date')
    plt.xticks(rotation=25)
    plt.ylabel('Deaths')
    plt.ylim(bottom=0)
    plt.grid(True, linestyle='--')
    plt.tight_layout()
    if saveplot:
        plt.savefig(os.path.join('Plots', fig_name+".png"))

    if moving_average:
        df = pd.read_json(os.path.join(data_folder, "all_county_ma_rki.json"))
        mask = (df['Date'] >= daystart) & (df['Date'] <= daystart + pd.DateOffset(days=simulationperiod)) & (
                df["County"] == county)
        fig_name = 'RKI_ma_confirmed_county_' + county.replace(" ", "_")
        plt.figure(fig_name)
        plt.plot(df.loc[mask]["Date"], df.loc[mask]["Confirmed"])
        plt.title('7-days moving average cumulative Confirmed cases ' + county)
        plt.xlabel('Date')
        plt.xticks(rotation=25)
        plt.ylabel('Confirmed cases')
        plt.ylim(bottom=0)
        plt.grid(True, linestyle='--')
        plt.tight_layout()
        if saveplot:
            plt.savefig(os.path.join('Plots', fig_name+".png"))

        fig_name = 'RKI_ma_deaths_county_' + county.replace(" ", "_")
        plt.figure(fig_name)
        plt.plot(df.loc[mask]["Date"], df.loc[mask]["Deaths"])
        plt.title('7-days moving average cumulative deaths ' + county)
        plt.xlabel('Date')
        plt.xticks(rotation=25)
        plt.ylabel('Deaths')
        plt.ylim(bottom=0)
        plt.grid(True, linestyle='--')
        plt.tight_layout()
        if saveplot:
            plt.savefig(os.path.join('Plots', fig_name+".png"))


def plot_DIVI_data(daystart=pd.Timestamp('2020.04.26'), simulationperiod=100, endday_divi=yesterday, saveplot=False):
    """plot DIVI data
        @param daystart Day at which should be started in timestamp format
        @param simulationperiod number in integer format of days for which data should be plotted
        @param endday_divi last available day of divi data in timestamp format
        @param saveplot boolean value; says if plot should be saved """
    if daystart < pd.Timestamp('2020.04.26'):
        print("error: DIVI dataset including ventilated data starts on 26.04.2020; you asked for " + daystart.strftime(
            "%d.%m.%Y") + " please search for another solution")
    else:
        df = pd.read_json(os.path.join(data_folder, "germany_divi.json"))
        if not (daystart + pd.DateOffset(days=simulationperiod) <= endday_divi):
            simulationperiod = (endday_divi - daystart).days
        mask = (df['Date'] >= daystart) & (df['Date'] <= daystart + pd.DateOffset(days=simulationperiod))
        fig_name = 'DIVI_plot'
        plt.figure(fig_name)
        plt.plot(df.loc[mask]['Date'], df.loc[mask]["ICU"])
        plt.plot(df.loc[mask]['Date'], df.loc[mask]["ICU_ventilated"])
        plt.title('DIVI data')
        plt.legend(['ICU', 'ICU_ventilated'])
        plt.xticks(rotation=25)
        plt.xlabel('Date')
        plt.ylabel('numer of people')
        plt.ylim(bottom=0)
        plt.grid(True, linestyle='--')
        plt.tight_layout()
        if saveplot:
            plt.savefig(os.path.join('Plots', fig_name+".png"))


def main():
    """plot all available plts
        @param read_data boolean, if True data will just be read from available data downloaded before,
         else data will be downloaded
        @param saveplot boolean value; says if plot should be saved """

    read_data = True
    saveplot = False
    endday_divi = pd.Timestamp('2021.01.04')
    moving_average = False

    if not read_data:
        get_Data(endday_divi, moving_average)
    if saveplot:
        if not os.path.exists('Plots'):
            os.makedirs('Plots')

    simulationperiod = 400
    daystart = pd.Timestamp('2020.03.02')
    plot_RKI(moving_average, daystart, simulationperiod, saveplot)
    plot_RKI_age(moving_average, daystart, simulationperiod, saveplot)
    plot_RKI_county(moving_average, daystart, simulationperiod, county='SK Köln', saveplot=saveplot)
    plot_DIVI_data(simulationperiod=simulationperiod, saveplot=saveplot, endday_divi=endday_divi)
    if not saveplot:
        plt.show()


if __name__ == "__main__":
    main()
