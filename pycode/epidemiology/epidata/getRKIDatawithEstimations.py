## @file getRKIDatawithEstimations.py
#
# @brief Estimates recovered and deaths from data from RKI and john hopkins together
#

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from epidemiology.epidata import getRKIData as grd
from epidemiology.epidata import getJHData as gjd


def get_rki_data_with_estimations(read_data=dd.defaultDict['read_data'],
                                  out_form=dd.defaultDict['out_form'],
                                  out_folder=dd.defaultDict['out_folder'],
                                  make_plot=dd.defaultDict['make_plot']):

    """! Function to estimate recovered and deaths from combination of RKI and JH data

    From the John-Hopkins (JH) data the fraction revered/confirmed and deaths/confiremd are calculated
    With this fraction every existing RKI data is scaled.
    The new columns recovered_estimated and deaths_estimated are added.

    @param read_data False [Default] or True. Defines if data is read from file or downloaded.
    @param out_form json [Default]
    @param out_folder Folder where data is written to.
    @param make_plot [Optional] RKI and estimated data can be compared by plots
    """

    data_path = os.path.join(out_folder, 'Germany/')

    if not read_data:

       # get rki data
       grd.get_rki_data(read_data, out_form, out_folder, False)

       # get data from John Hopkins University
       gjd.get_jh_data(read_data, out_form, out_folder)

    # Now we now which data is generated and we can use it
    # read in jh data

    jh_data_file = os.path.join(data_path, 'whole_country_Germany_jh.json')
    df_jh = pd.read_json(jh_data_file)

    confirmed = dd.EngEng['confirmed']
    recovered = dd.EngEng['recovered']
    deaths = dd.EngEng['deaths']
    date = dd.EngEng['date']

    fr_r_c = "fraction_recovered_confirmed"
    fr_d_c = "fraction_deaths_confirmed"

    recovered_estimated = recovered + "_estimated"
    recovered_after14days = recovered + "_after14days"
    deaths_estimated = deaths+ "_estimated"
    dstr = date + "_string"
    week = "Week"

    #delta = timedelta(days=14)

    df_jh[fr_r_c] =  df_jh[recovered] / df_jh[confirmed]
    df_jh[fr_d_c] = df_jh[deaths] / df_jh[confirmed]

    # There might be some NaN due to division with 0
    # Nan is replaced by 0
    df_jh[fr_r_c] = df_jh[fr_r_c].fillna(0)
    df_jh[fr_d_c] = df_jh[fr_d_c].fillna(0)

    # get data from rki and make new data
    rki_files_to_change = ["all_germany_rki", "all_gender_rki", "all_age_rki",
                           "all_state_rki", "all_state_gender_rki", "all_state_age_rki",
                           "all_county_rki", "all_county_gender_rki", "all_county_age_rki"]

    for file_to_change in rki_files_to_change:
        # read data of rki file
        rki_data_file = os.path.join(data_path, file_to_change + ".json")
        df_rki = pd.read_json(rki_data_file)

        # generate new columns to store estimated values
        # TODO Add also column infected and calculate it in the end
        df_rki[recovered_estimated] = np.nan
        #df_rki[recovered_after14days] = np.nan
        df_rki[deaths_estimated] = np.nan

        # convert time stamp of rki data
        df_rki[dstr] = df_rki[date].dt.strftime('%Y-%m-%d')

        for i in range(len(df_jh)):

            date_jh = df_jh.loc[i, date].strftime('%Y-%m-%d')
            fraction_rec_conf = df_jh.loc[i, fr_r_c]
            fraction_deaths_conf = df_jh.loc[i, fr_d_c]

            # go to date and calculate estimation
            df_rki.loc[(df_rki[dstr] == date_jh), recovered_estimated] \
                = np.round(fraction_rec_conf*df_rki.loc[(df_rki[dstr] == date_jh), confirmed])

            # TODO: Check how recovered would be, if everyone would be recovered in 14 days
            #now_confirmed = df_rki.loc[(df_rki[dstr] == date_jh), confirmed]
            # shift confirmed to recovered 14 days later
            #date_after14days = (datetime.strptime(date_jh, '%Y-%m-%d')+delta).strftime('%Y-%m-%d')
            #df_rki.loc[(df_rki[dstr] == date_after14days), recovered_after14days] = now_confirmed

            df_rki.loc[(df_rki[dstr] == date_jh), deaths_estimated] \
                = np.round(fraction_deaths_conf * df_rki.loc[(df_rki[dstr] == date_jh), confirmed])

        df_rki = df_rki.drop([dstr], 1)
        gd.write_dataframe(df_rki, data_path, file_to_change + "_estimated", out_form)

        # check if calculation is meaningful
        # TODO Add jh data to whole germany plot

        if(make_plot == True):
           df_rki.plot(x=date, y = [recovered, recovered_estimated],
                       title = 'COVID-19 check recovered for '+ file_to_change,
                       grid = True, style = '-o')
           plt.tight_layout()
           plt.show()

           df_rki.plot(x=date, y=[deaths, deaths_estimated],
                       title='COVID-19 check deaths '+ file_to_change,
                       grid=True, style='-o')
           plt.tight_layout()
           plt.show()

           df_rki[week] = df_rki[date].dt.isocalendar().week

           # TODO download and plot the rki file where there are the real number of deaths dependent on week number.
           # df_rki_week = df_rki.groupby(week).agg({deaths: sum, deaths_estimated: sum}).reset_index()

           #url = "https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Projekte_RKI/COVID-19_Todesfaelle.xlsx?__blob=publicationFile"

           #df_real_deaths_per_week = pd.read_excel(url)

           #df_rki_week.plot(x=week, y=[deaths, deaths_estimated], title='COVID-19 check deaths dependent on wek number', grid=True,
           #            style='-o')
           #plt.tight_layout()
           #plt.show()

        # TODO: think about to add age dependent weight function


def main():
    """! Main program entry."""

    [read_data, out_form, out_folder, make_plot] = gd.cli("rkiest")
    get_rki_data_with_estimations(read_data, out_form, out_folder, make_plot)


if __name__ == "__main__":
    main()
