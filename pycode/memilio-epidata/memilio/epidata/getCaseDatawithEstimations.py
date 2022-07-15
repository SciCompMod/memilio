#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Kathrin Rack, Wadim Koslow
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
"""
@file getCaseDatawithEstimations.py

@brief Estimates recovered and deaths from case data from RKI and Johns Hopkins (JH) University together
"""

import os
from datetime import datetime
#from datetime import timedelta
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import requests

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from memilio.epidata import getCaseData as gcd
from memilio.epidata import getJHData as gjd


def get_case_data_with_estimations(data_folder,
                                   read_data=dd.defaultDict['read_data'],
                                   file_format=dd.defaultDict['file_format'],
                                   no_raw=dd.defaultDict['no_raw'],
                                   make_plot=dd.defaultDict['make_plot']):
    """! Function to estimate recovered and deaths from combination of case data from RKI and JH data
WARNING: This file is experimental and has not been tested.
    From the John-Hopkins (JH) data the fraction recovered/confirmed and deaths/confirmed
    are calculated for whole Germany
    With this fraction every existing case data from RKI is scaled.
    The new columns recovered_estimated and deaths_estimated are added.

    @param data_folder Folder where data is written to.
    @param read_data False [Default] or True. Defines if data is read from file or downloaded.
    @param file_format json [Default]
    @param no_raw True or False [Default]. Defines if unchanged raw data is saved or not.
    @param make_plot [Optional] case data from RKI and estimated data can be compared by plots
    """

    data_path = os.path.join(data_folder, 'Germany/')

    if not read_data:
        impute_dates = False
        make_plot_cases = False
        moving_average = False
        split_berlin = False

        # get case data
        gcd.get_case_data(
            data_path, read_data, file_format, no_raw, impute_dates,
            make_plot_cases, moving_average, no_raw, split_berlin)

        # get data from John Hopkins University
        gjd.get_jh_data(data_path, read_data, file_format, no_raw)

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
    # recovered_after14days = recovered + "_after14days"
    deaths_estimated = deaths + "_estimated"
    dstr = date + "_string"
    # week = "Week"
    # delta = timedelta(days=14)

    df_jh[fr_r_c] = df_jh[recovered] / df_jh[confirmed]
    df_jh[fr_d_c] = df_jh[deaths] / df_jh[confirmed]

    # There might be some NaN due to division with 0
    # Nan is replaced by 0
    df_jh[fr_r_c] = df_jh[fr_r_c].fillna(0)
    df_jh[fr_d_c] = df_jh[fr_d_c].fillna(0)

    # get case data from rki and make new data
    case_files_to_change = [
        "cases_all_germany", "cases_all_gender", "cases_all_age",
        "cases_all_state", "cases_all_state_gender", "cases_all_state_age",
        "cases_all_county", "cases_all_county_gender", "cases_all_county_age"]

    for file_to_change in case_files_to_change:
        # read data of case file
        case_data_file = os.path.join(data_path, file_to_change + ".json")
        try:
            df_cases = pd.read_json(case_data_file)
        except ValueError as e:
            print("WARNING: The file ", file_to_change + ".json does not exist.")
            continue

        # generate new columns to store estimated values
        # TODO Add also column infected and calculate it in the end
        df_cases[recovered_estimated] = np.nan
        # df_cases[recovered_after14days] = np.nan
        df_cases[deaths_estimated] = np.nan

        # convert time stamp of case data
        df_cases[dstr] = df_cases[date].dt.strftime('%Y-%m-%d')

        for i in range(len(df_jh)):
            date_jh = df_jh.loc[i, date].strftime('%Y-%m-%d')
            fraction_rec_conf = df_jh.loc[i, fr_r_c]
            fraction_deaths_conf = df_jh.loc[i, fr_d_c]

            # go to date and calculate estimation
            df_cases.loc[(df_cases[dstr] == date_jh), recovered_estimated] = np.round(
                fraction_rec_conf * df_cases.loc[(df_cases[dstr] == date_jh), confirmed])

            # TODO: Check how recovered would be, if everyone would be recovered in 14 days
            # now_confirmed = df_cases.loc[(df_cases[dstr] == date_jh), confirmed]
            # shift confirmed to recovered 14 days later
            # date_after14days = (datetime.strptime(date_jh, '%Y-%m-%d')+delta).strftime('%Y-%m-%d')
            # df_cases.loc[(df_cases[dstr] == date_after14days), recovered_after14days] = now_confirmed

            df_cases.loc[(df_cases[dstr] == date_jh), deaths_estimated] = np.round(
                fraction_deaths_conf * df_cases.loc[(df_cases[dstr] == date_jh), confirmed])

        df_cases = df_cases.drop([dstr], axis=1)
        gd.write_dataframe(df_cases, data_path,
                           file_to_change + "_estimated", file_format)

        # check if calculation is meaningful
        # TODO Add jh data to whole germany plot
        if make_plot:
            df_cases.plot(
                x=date, y=[recovered, recovered_estimated],
                title='COVID-19 check recovered for ' + file_to_change,
                grid=True, style='-o')
            plt.tight_layout()
            plt.show()

            df_cases.plot(x=date, y=[deaths, deaths_estimated],
                          title='COVID-19 check deaths ' + file_to_change,
                          grid=True, style='-o')
            plt.tight_layout()

            plt.show()

        if file_to_change == "cases_all_germany":
            compare_estimated_and_rki_deathsnumbers(
                df_jh, df_cases, data_path, read_data, make_plot)
            get_weekly_deaths_data_age_gender_resolved(
                data_path, read_data=True)
            # df_cases[week] = df_cases[date].dt.isocalendar().week
    #if make_plot:
    #    plt.show()


def compare_estimated_and_rki_deathsnumbers(
        df_jh, df_cases, data_path, read_data, make_plot):
    """! Comparison of estimated values and monthly case values from rki

    From the daily number of deaths the value before is subtracted to have the actual value of the day not
    an accumulation.
    From the daily values, the weekly values are calculated by calculating the sum for the same week.
    This is done with the values from the original case data from RKI and for the estimated values
    (RKI-confirmed times JH fraction).
    Furthermore, we read the excel table with weekly case data given by the RKI additionally.

    From this comparison we can see how good the estimated values are.

    @param df_cases Pandas dataframe with case data
    @param data_path Path where to store the file.
    @param read_data False or True. Defines if data is read from file or downloaded.
    @param make_plot Defines if plots are generated

    """
    df_cases['Date'] = pd.to_datetime(df_cases['Date'], format="%Y-%m-%d")
    # we set january 2020 to week 1
    # 2020 had 53 weeks
    # meaning weak 45 is first week in 2021
    df_cases["week"] = df_cases['Date'].dt.isocalendar(
    ).week + (df_cases['Date'].dt.isocalendar().year-2020)*53
    df_jh["week"] = df_jh['Date'].dt.isocalendar(
    ).week + (df_jh['Date'].dt.isocalendar().year - 2020)*53
    # want to have daily deaths numbers, not accumulated
    df_cases["deaths_daily"] = df_cases['Deaths'] - \
        df_cases['Deaths'].shift(periods=1, fill_value=0)
    df_jh["deaths_daily"] = df_jh['Deaths'] - \
        df_jh['Deaths'].shift(periods=1, fill_value=0)
    df_cases["deaths_estimated_daily"] = df_cases['Deaths_estimated'] - \
        df_cases['Deaths_estimated'].shift(periods=1, fill_value=0)
    df_cases_week = df_cases.groupby("week").agg(
        {"deaths_daily": sum, "deaths_estimated_daily": sum}).reset_index()
    df_jh_week = df_jh.groupby("week").agg({"deaths_daily": sum}).reset_index()
    df_cases_week.rename(
        columns={'deaths_daily': 'Deaths_weekly',
                 'deaths_estimated_daily': 'Deaths_estimated_weekly'},
        inplace=True)
    df_jh_week.rename(
        columns={'deaths_daily': 'Deaths_weekly'}, inplace=True)

    # download weekly deaths numbers from rki
    if not read_data:
        download_weekly_deaths_numbers(data_path)

    df_real_deaths_per_week = gd.loadExcel(
        "Cases_deaths_weekly", apiUrl=data_path, extension=".xlsx",
        param_dict={"sheet_name": 'COVID_Todesfälle', "header": 0,
                    "engine": 'openpyxl'})
    df_real_deaths_per_week.rename(
        columns={'Sterbejahr': 'year', 'Sterbewoche': 'week',
                 'Anzahl verstorbene COVID-19 Fälle': 'confirmed_deaths_weekly'},
        inplace=True)
    df_real_deaths_per_week['confirmed_deaths_weekly'] = pd.to_numeric(
        df_real_deaths_per_week['confirmed_deaths_weekly'], errors='coerce')
    # values which can't be transformed to numeric number are a String called '<4' in dataframe
    # (probably for data protection reasons)
    # So get random numbers between 1 and 3 to fill these missing values
    help_random = np.random.randint(
        1, 4, (df_real_deaths_per_week['confirmed_deaths_weekly'].shape[0]))
    random_fill = pd.Series(help_random)
    df_real_deaths_per_week['confirmed_deaths_weekly'] =\
        df_real_deaths_per_week['confirmed_deaths_weekly'].fillna(random_fill)
    # we set january 2020 to week 1
    # 2020 had 53 weeks
    # meaning weak 54 is first week in 2021
    df_real_deaths_per_week.loc[df_real_deaths_per_week.year ==
                                2021, 'week'] += 53

    #combine both dataframes to one dataframe
    df_cases_week = df_cases_week.merge(
        df_real_deaths_per_week, how='outer', on="week")
    del df_cases_week['year']
    # confirmed_deaths_weekly has lass rows than other data so fill NaN's with zeros
    df_cases_week['confirmed_deaths_weekly'] = df_cases_week['confirmed_deaths_weekly'].fillna(
        0)

    # safe information in dataframe in json-file
    gd.write_dataframe(df_cases_week, data_path, 'cases_weekly_deaths', 'json')

    if make_plot:
        df_cases_week['confirmed_deaths_accumulated'] = df_cases_week['confirmed_deaths_weekly'].cumsum()
        df_cases_week['Deaths_accumulated'] = df_cases_week['Deaths_weekly'].cumsum()
        df_cases_week['Deaths_estimated_accumulated'] = df_cases_week['Deaths_estimated_weekly'].cumsum()

        df_jh_week['Deaths_weekly_accumulated'] = df_jh_week['Deaths_weekly'].cumsum()

        #plot
        df_cases_week.plot(
            x="week",
            y=["Deaths_weekly", "Deaths_estimated_weekly",
               "confirmed_deaths_weekly"],
            title='COVID-19 check deaths dependent on week number', grid=True,
            style='-o')
        plt.plot(df_jh_week["week"], df_jh_week["Deaths_weekly"])
        plt.legend(["cases daily", "estimated with JH",
                    "weekly cases from RKI", "JH"])
        plt.tight_layout()
        df_cases_week.plot(
            x="week",
            y=["Deaths_accumulated", "Deaths_estimated_accumulated",
               "confirmed_deaths_accumulated"],
            title='COVID-19 check deaths accumulated dependent on week number',
            grid=True, style='-o')
        plt.plot(df_jh_week["week"], df_jh_week["Deaths_weekly_accumulated"])
        plt.legend(["cases daily", "estimated with JH",
                    "weekly cases from RKI", "JH"])
        plt.tight_layout()
        plt.show()

    # TODO: think about to add age dependent weight function


def get_weekly_deaths_data_age_gender_resolved(data_path, read_data):
    """! Read case data from rki excel file

    @param data_path Path where to store the file.
    @param read_data False or True. Defines if data is read from file or downloaded.
    """

    if not read_data:
        download_weekly_deaths_numbers(data_path)

    df_real_deaths_per_week_age = gd.loadExcel(
        'Cases_deaths_weekly', apiUrl=data_path, extension='.xlsx',
        param_dict={'sheet_name': 'COVID_Todesfälle_KW_AG10', "header": 0,
                    "engine": 'openpyxl'})
    df_real_deaths_per_week_gender = gd.loadExcel(
        'Cases_deaths_weekly', apiUrl=data_path, extension='.xlsx',
        param_dict={'sheet_name': 'COVID_Todesfälle_KW_AG20_G', "header": 0,
                    "engine": 'openpyxl'})
    df_real_deaths_per_week_age.rename(
        columns={'Sterbejahr': 'year', 'Sterbewoche': 'week',
                 'AG 0-9 Jahre': 'age 0-9 years',
                 'AG 10-19 Jahre': 'age 10-19 years',
                 'AG 20-29 Jahre': 'age 20-29 years',
                 'AG 30-39 Jahre': 'age 30-39 years',
                 'AG 40-49 Jahre': 'age 40-49 years',
                 'AG 50-59 Jahre': 'age 50-59 years',
                 'AG 60-69 Jahre': 'age 60-69 years',
                 'AG 70-79 Jahre': 'age 70-79 years',
                 'AG 80-89 Jahre': 'age 80-89 years',
                 'AG 90+ Jahre': 'age 90+ years'},
        inplace=True)
    df_real_deaths_per_week_gender.rename(
        columns={'Sterbejahr': 'year', 'Sterbewoche': 'week',
                 'Männer, AG 0-19 Jahre': 'male, age 0-19 years',
                 'Männer, AG 20-39 Jahre': 'male, age 20-39 years',
                 'Männer, AG 40-59 Jahre': 'male, age 40-59 years',
                 'Männer, AG 60-79 Jahre': 'male, age 60-79 years',
                 'Männer, AG 80+ Jahre': 'male, age 80+ years',
                 'Frauen, AG 0-19 Jahre': 'female, age 0-19 years',
                 'Frauen, AG 20-39 Jahre': 'female, age 20-39 years',
                 'Frauen, AG 40-59 Jahre': 'female, age 40-59 years',
                 'Frauen, AG 60-79 Jahre': 'female, age 60-79 years',
                 'Frauen, AG 80+ Jahre': 'female, age 80+ years'},
        inplace=True)

    for df_real_deaths_per_week in [
            df_real_deaths_per_week_age, df_real_deaths_per_week_gender]:
        for column in df_real_deaths_per_week:
            df_real_deaths_per_week[column] = pd.to_numeric(
                df_real_deaths_per_week[column], errors='coerce')

            # values which can't be transformed to numeric number are a String called '<4' in dataframe
            # (probably for data protection reasons)
            # So get random numbers between 1 and 3 to fill these missing values
            help_random = np.random.randint(
                1, 4, (df_real_deaths_per_week[column].shape[0]))
            random_fill = pd.Series(help_random)
            df_real_deaths_per_week[column] = df_real_deaths_per_week[column].fillna(
                random_fill)

    # safe information in dataframe in json-file
    gd.write_dataframe(
        df_real_deaths_per_week_age, data_path,
        'cases_weekly_deaths_age_resolved', 'json')
    gd.write_dataframe(
        df_real_deaths_per_week_gender, data_path,
        'cases_weekly_deaths_gender_resolved', 'json')


def download_weekly_deaths_numbers(data_path):
    """!Downloads excel file from RKI webpage

    @param data_path Path where to store the file.
    """

    name_file = "Cases_deaths_weekly.xlsx"
    url = "https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Projekte_RKI/" \
          "COVID-19_Todesfaelle.xlsx?__blob=publicationFile"

    # data_path: path where to safe Excel-file
    r = requests.get(url)
    filename = os.path.join(data_path, name_file)
    with open(filename, 'wb') as output_file:
        output_file.write(r.content)


def main():
    """! Main program entry."""

    path = os.path.join(os.getcwd(), 'data', 'pydata')
    arg_dict = gd.cli("cases_est")
    get_case_data_with_estimations(path, **arg_dict)


if __name__ == "__main__":
    main()
