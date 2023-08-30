#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Patrick Lenz
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
@file getHospitalizationData.py
@brief Downloads the hospitalization data of the Robert Koch-Institute (RKI) and provides it in different ways.

The raw hospitalization data can be found at
https://github.com/robert-koch-institut/COVID-19-Hospitalisierungen_in_Deutschland
"""

# Imports
import os
from datetime import timedelta

import numpy as np
import pandas as pd

from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import modifyDataframeSeries as mdfs


def hospit_sanity_checks(df):
    """! Checks the sanity of the hospitalization_data dataframe

    Checks if type of the given data is a dataframe
    Checks if the headers of the dataframe are those which are needed

    @param df The dataframe which has to be checked
    """
    # test if dataframe is empty
    if df.empty:
        raise gd.DataError(
            "Download of Hospitalization Data failed. File is empty.")

    actual_strings_list = df.columns.tolist()
    # check number of data categories
    if len(actual_strings_list) != 6:
        print("Warning: Number of data categories changed.")

    # these strings need to be in the header
    test_strings = {
        'Datum', 'Bundesland', 'Bundesland_Id', 'Altersgruppe',
        '7T_Hospitalisierung_Faelle', '7T_Hospitalisierung_Inzidenz'}
    # check if headers are those we want
    for name in test_strings:
        if (name not in actual_strings_list):
            raise gd.DataError("Error: Data categories have changed.")


def get_hospitailzations_per_day(seven_days_values):
    """! Gets the daily cases of hospitalizations from the seven day sum.

    A zero filled array is created where the one day data is stored.
    For each calculated daily case the copied array is adjusted to include only the remaining values.
    Whenever the seven day sum changes to the following day, the difference is the same as the daily case difference from last week.
    The computation can be done forward and backward. Both is done until the array is flattened.
    If there are constant cases left on each day, they are divided on each date by one seventh.
    After that, a few tests are done to check if all cases were distributed correctly.

    @param seven_days_values Array. Total hospitalizations over the last seven days

    @return daily_values Hospitalizations per day.
    """

    daily_values = np.zeros(len(seven_days_values), dtype=float)
    to_split = seven_days_values.copy()
    run = 1
    while sum(to_split[6:]) != 0 and run < 5:
        # backward computation
        if max(to_split[6:]) != min(to_split[6:]):
            backward = np.zeros(len(daily_values), dtype=float)
            for i in range(1, len(daily_values)-6):
                if to_split[-i-1] > to_split[-i]:
                    backward[-i-7] = to_split[-i-1]-to_split[-i]
                    for day in range(7):
                        try:
                            to_split[-i-1-day] -= backward[-7-i]
                        except IndexError:
                            pass
            daily_values += backward
        # start at first known value
        if max(to_split[6:]) > 0:
            forward = np.zeros(len(daily_values), dtype=float)
            for i in range(7, len(to_split)):
                if to_split[i-1] < to_split[i]:
                    forward[i] += to_split[i]-to_split[i-1]
                    for day in range(7):
                        try:
                            to_split[i+day] -= forward[i]
                        except IndexError:
                            pass
            daily_values += forward

        if max(to_split[6:]) == min(to_split[6:]):
            daily_values += max(to_split[6:])/7
            to_split[6:] -= max(to_split[6:])
        run += 1

    # break after 5 runs to prevent endless loop
    if run == 5:
        print("Can't get hospitalizations per day from incidence.")
    if len(daily_values[daily_values < 0]) > 0:
        raise gd.DataError('Negative hospitalizations found.')
    # check that daily values are calculated correctly
    check = np.zeros(len(daily_values)+7, dtype=float)
    for i in range(len(daily_values)):
        for day in range(7):
            check[i+day] += daily_values[i]
    if sum(check[6:-7]-seven_days_values[6:]) > 10**-7:
        raise gd.DataError("Check failed.")

    return daily_values


def get_hospitalization_data(read_data=dd.defaultDict['read_data'],
                             file_format=dd.defaultDict['file_format'],
                             out_folder=dd.defaultDict['out_folder'],
                             no_raw=dd.defaultDict['no_raw'],
                             start_date=dd.defaultDict['start_date'],
                             end_date=dd.defaultDict['end_date'],
                             impute_dates=dd.defaultDict['impute_dates'],
                             moving_average=dd.defaultDict['moving_average'],
                             make_plot=dd.defaultDict['make_plot']
                             ):
    """! Downloads or reads the RKI hospitalization data and writes them in different files.

    Available data starts from 2020-03-01.
    If it does not already exist, the folder Germany is generated in the given out_folder.
    If read_data == True and the file "RKIHospitFull.json" exists, the data is read from this file
    and stored in a pandas dataframe. If read_data = True and the file does not exist the program is stopped.

    The downloaded dataframe is written to the file "RKIHospitFull".
    After that, the columns are renamed to english.
    From the sum of the cases of the last seven days the daily cases are calculated.
    Afterwards, the data is stored in four different files:
    "hospit_state_age", "hospit_state", "hospit_germany_age" and "hospit_germany"
    for states or germany and age groups.

    @param read_data True or False. Defines if data is read from file or downloaded.  Default defined in defaultDict.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    @param out_folder Folder where data is written to. Default defined in defaultDict.
    @param no_raw True or False. Defines if unchanged raw data is saved or not. Default defined in defaultDict.
    @param start_date Date of first date in dataframe. Default defined in defaultDict.
    @param end_date Date of last date in dataframe. Default defined in defaultDict.
    @param impute_dates True or False. Defines if values for dates without new information are imputed. Default defined in defaultDict.
        Here Dates are always imputed so False changes nothing.
    @param moving_average [Currently not used] Integers >=0. Applies an 'moving_average'-days moving average on all time series
        to smooth out weekend effects.  Default defined in defaultDict.
    @param make_plot [currently not used] True or False. Defines if plots are generated with matplotlib. Default defined in defaultDict.
    """
    impute_dates = True
    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    # get raw dataframe
    filename = "RKIHospitFull"
    url = "https://raw.githubusercontent.com/robert-koch-institut/COVID-19-Hospitalisierungen_in_Deutschland/master/Aktuell_Deutschland_COVID-19-Hospitalisierungen.csv"
    path = os.path.join(directory + filename + ".json")
    df_raw = gd.get_file(path, url, read_data, param_dict={}, interactive=True)

    hospit_sanity_checks(df_raw)

    if not no_raw:
        gd.write_dataframe(df_raw, directory, filename, file_format)
    df_data = df_raw.copy()
    # drop unwanted columns and rows, rename and sort dataframe
    df_data.rename(dd.GerEng, axis=1, inplace=True)
    df_data.rename(columns={'Datum': dd.EngEng['date']}, inplace=True)
    df_data = df_data.drop(columns=['State', '7T_Hospitalisierung_Inzidenz'])
    df_data = df_data.sort_values(
        by=[dd.EngEng['date'],
            dd.EngEng['idState'],
            dd.EngEng['ageRKI']]).reset_index(
        drop=True)
    # impute 6 days before min_date to split up the seven day cases
    df_data = mdfs.impute_and_reduce_df(
        df_old=df_data,
        group_by_cols={dd.EngEng['idState']:
                       df_data[dd.EngEng['idState']].unique(),
                       dd.EngEng['ageRKI']:
                       df_data[dd.EngEng['ageRKI']].unique()},
        mod_cols=['7T_Hospitalisierung_Faelle'],
        impute='zeros', moving_average=0, min_date=pd.to_datetime(
            min(df_data.Date)).date() - timedelta(6),
        max_date='', start_w_firstval=False)

    # get data for each day
    # for each state and age group seperately
    df_daily = pd.DataFrame()
    for age in df_data[dd.EngEng['ageRKI']].unique():
        df_age = df_data[df_data[dd.EngEng['ageRKI']] == age]
        for stateid in df_data[dd.EngEng['idState']].unique():
            df_age_stateid = df_age[df_age[dd.EngEng['idState']]
                                    == stateid].copy()
            # get hospitalizations per day from incidence
            seven_days_values = df_age_stateid['7T_Hospitalisierung_Faelle'].values
            daily_values = get_hospitailzations_per_day(seven_days_values)
            # save data in dataframe
            df_age_stateid['hospitalized'] = daily_values
            df_age_stateid = df_age_stateid.drop(
                ['7T_Hospitalisierung_Faelle'], axis=1)
            df_daily = pd.concat(
                [df_daily.reset_index(drop=True),
                 df_age_stateid.reset_index(drop=True)],
                join='outer')

    df_daily = mdfs.extract_subframe_based_on_dates(
        df_daily, start_date, end_date)
    # write dataframe with all states
    # drop columns for states and all age groups
    df_state_age = df_daily[df_daily[dd.EngEng['idState']] != 0]
    df_state_age = df_state_age[df_state_age[dd.EngEng['ageRKI']] != '00+']
    filename = gd.append_filename(
        'hospit_state_age', impute_dates, moving_average=0)
    gd.write_dataframe(df_state_age, directory, filename, file_format)
    # write dataframe for germany and all age groups
    df_germany_age = df_daily[df_daily[dd.EngEng['idState']] == 0]
    df_germany_age = df_germany_age[
        df_germany_age[dd.EngEng['ageRKI']] != '00+'].drop(
        columns=[dd.EngEng['idState']])
    filename = gd.append_filename(
        'hospit_germany_age', impute_dates, moving_average=0)
    gd.write_dataframe(df_germany_age, directory, filename, file_format)
    # write dataframe for states and combined age groups
    df_state = df_daily[df_daily[dd.EngEng['idState']] != 0]
    df_state = df_state[df_state[dd.EngEng['ageRKI']] == '00+'].drop(
        columns=[dd.EngEng['ageRKI']])
    filename = gd.append_filename(
        'hospit_state', impute_dates, moving_average=0)
    gd.write_dataframe(df_state, directory, filename, file_format)
    # write dataframe for germany and combined age groups
    df_germany = df_daily[df_daily[dd.EngEng['idState']] == 0]
    df_germany = df_germany[df_germany[dd.EngEng['ageRKI']] == '00+'].drop(
        columns=[dd.EngEng['idState'],
                 dd.EngEng['ageRKI']])
    filename = gd.append_filename(
        'hospit_germany', impute_dates, moving_average=0)
    gd.write_dataframe(df_germany, directory, filename, file_format)


def main():
    """! Main program entry."""
    arg_dict = gd.cli('hospitalization')
    get_hospitalization_data(**arg_dict)


if __name__ == "__main__":

    main()
