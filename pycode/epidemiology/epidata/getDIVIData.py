#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Kathrin Rack, Lena Ploetzke, Martin J. Kuehn
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
@file getDIVIData.py

@brief Data of the DIVI
about Sars-CoV2 is downloaded.
This data contains the number of Covid19 patients in intensive care
and the number of those that are additionally ventilated.

DIVI - Deutsche interdisziplinäre Vereinigung für Intensiv- und Notfallmedizin

data explanation:
- reporting_hospitals is the number of reporting hospitals
- ICU is the number of covid patients in reporting hospitals
- ICU_ventilated is the number of ventilated covid patients in reporting hospitals
- free_ICU is the number of free ICUs in reporting hospitals
- occupied_ICU is the number of occupied ICUs in in reporting hospitals

ID_County and ID_State is defined by the "Amtlicher Gemeindeschlüssel (AGS)"
which is also used in the RKI data as ID_County and ID_State
https://de.wikipedia.org/wiki/Liste_der_Landkreise_in_Deutschland.

Specific features about the data:
The column "faelle_covid_im_bundesland" exits only in the data from the first day (24.4)
The column ICU does not exist for the 24.4.
ICU_ventilated does not exist for the 24.4. and 25.4.
"""

import os
import sys
import bisect
from datetime import timedelta, date, datetime
import pandas as pd

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from epidemiology.epidata import geoModificationGermany as geoger
from epidemiology.epidata import modifyDataframeSeries


def cut_of_dates(df, start_date, end_date):

    startdelta = (start_date - dd.defaultDict['start_date']).days
    enddelta = (dd.defaultDict['end_date'] - end_date).days
    for i in range(startdelta):
        lowest_date = dd.defaultDict['start_date'] + timedelta(i)
        dt = lowest_date
        df_new = df[df.Date != datetime.strftime(dt, '%Y-%m-%d')]
        df = df_new
    for i in range(enddelta):
        highest_date = dd.defaultDict['end_date'] - timedelta(i)
        dt = highest_date
        df_new = df[df.Date != datetime.strftime(dt, '%Y-%m-%d')]
        df = df_new
    df.reset_index(drop=True, inplace=True)

    return df


def get_divi_data(read_data=dd.defaultDict['read_data'],
                  file_format=dd.defaultDict['file_format'],
                  out_folder=dd.defaultDict['out_folder'],
                  no_raw=dd.defaultDict['no_raw'],
                  end_date=dd.defaultDict['end_date'],
                  start_date=dd.defaultDict['start_date'],
                  impute_dates=dd.defaultDict['impute_dates'],
                  moving_average=dd.defaultDict['moving_average']
                  ):
    """! Downloads or reads the DIVI ICU data and writes them in different files.

    Available data starts from 2020-04-24.
    If the given start_date is earlier, it is changed to this date and a warning is printed.
    If it does not already exist, the folder Germany is generated in the given out_folder.
    If read_data == True and the file "FullData_DIVI.json" exists, the data is read form this file
    and stored in a pandas dataframe.
    Otherwise the program is stopped.

    The dataframe is written to the file filename = "FullData_DIVI".
    Than the columns are renamed to English and the state and county names are added.
    Afterwards, three kinds of structuring of the data are done.
    We obtain the chronological sequence of ICU and ICU_ventilated
    stored in the files "county_divi".json", "state_divi.json" and "germany_divi.json"
    for counties, states and whole Germany, respectively.

    @param read_data False [Default] or True. Defines if data is read from file or downloaded.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    "False [Default]" if it is downloaded for all dates from start_date to end_date.
    @param out_folder Folder where data is written to.
    @param no_raw True or False [Default]. Defines if unchanged raw data is saved or not.
    @param start_date [Optional] Date of first data in dataframe[Default = 2020.4.24].
    @param end_date [Optional] Date to last data in dataframe [Default = today].
    @param impute_dates True or False [Default]. Defines if values for dates without new information are imputed.
    @param moving_average 0 [Default] or >0. Applies an 'moving_average'-days moving average on all time series
        to smooth out weekend effects.
    """

    # First csv data on 24-04-2020
    if start_date < date(2020, 4, 24):
        print("Warning: First data available on 2020-04-24. "
              "You asked for " + start_date.strftime("%Y-%m-%d") + ".")
        start_date = date(2020, 4, 24)

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    filename = "FullData_DIVI"

    if read_data:
        # read json file for already downloaded data
        file_in = os.path.join(directory, filename + ".json")

        try:
            df = pd.read_json(file_in)
        except ValueError:
            exit_string = "Error: The file: " + file_in + " does not exist. "\
                          "Call program without -r flag to get it."
            sys.exit(exit_string)
    else:
        try:
            df = gd.loadCsv(
                'zeitreihe-tagesdaten',
                apiUrl='https://diviexchange.blob.core.windows.net/%24web/',
                extension='.csv')
        except:
            exit_string = "Error: Download link for Divi data has changed."
            sys.exit(exit_string)

    if not df.empty:
        if not no_raw:
            gd.write_dataframe(df, directory, filename, file_format)
    else:
        exit_string = "Something went wrong, dataframe is empty."
        sys.exit(exit_string)
    df_raw = df.copy()
    gdd_sanity_checks(df_raw)
    df.rename(columns={'date': dd.EngEng['date']}, inplace=True)
    df.rename(dd.GerEng, axis=1, inplace=True)

    df['Date'] = pd.to_datetime(df['Date'], format='%Y-%m-%d %H:%M:%S')
    df = cut_of_dates(df, start_date, end_date)

    # insert names of  states
    df.insert(loc=0, column=dd.EngEng["idState"], value=df[dd.EngEng["state"]])
    for item in geoger.get_state_names_and_ids():
        df.loc[df[dd.EngEng["idState"]] == item[1],
               [dd.EngEng["state"]]] = item[0]

    # insert names of counties
    df.insert(loc=3, column=dd.EngEng["county"],
              value=df[dd.EngEng["idCounty"]])
    for item in geoger.get_county_names_and_ids(merge_eisenach=False):
        df.loc[df[dd.EngEng["idCounty"]] == item[1],
               [dd.EngEng["county"]]] = item[0]

    # remove leading zeros for ID_County (if not yet done)
    df['ID_County'] = df['ID_County'].astype(int)
    # add missing dates (and compute moving average)
    if (impute_dates == True) or (moving_average > 0):
        df = modifyDataframeSeries.impute_and_reduce_df(
            df,
            {dd.EngEng["idCounty"]: geoger.get_county_ids(
                merge_eisenach=False)},
            [dd.EngEng["ICU"], dd.EngEng["ICU_ventilated"]],
            impute='forward', moving_average=moving_average)

    # add names etc for empty frames (counties where no ICU beds are available)
    countyid_to_name = geoger.get_countyid_to_name()
    stateid_to_name = geoger.get_stateid_to_name()
    countyid_to_stateid = geoger.get_countyid_to_stateid_map()
    for id in df.loc[df.isnull().any(axis=1), dd.EngEng['idCounty']].unique():
        stateid = countyid_to_stateid[id]
        df.loc[df[dd.EngEng['idCounty']] == id,
               [dd.EngEng['idState'],
                dd.EngEng['state'],
                dd.EngEng['county']]] = [stateid, stateid_to_name[stateid],
                                         countyid_to_name[id]]

    # write data for counties to file
    df_counties = df[[dd.EngEng["idCounty"],
                      dd.EngEng["county"],
                      dd.EngEng["ICU"],
                      dd.EngEng["ICU_ventilated"],
                      dd.EngEng["date"]]].copy()
    # merge Eisenach and Wartburgkreis from DIVI data
    df_counties = geoger.merge_df_counties_all(
        df_counties, sorting=[dd.EngEng["idCounty"], dd.EngEng["date"]])
    # save
    filename = "county_divi"
    filename = gd.append_filename(filename, impute_dates, moving_average)
    gd.write_dataframe(df_counties, directory, filename, file_format)

    # write data for states to file
    df_states = df.groupby(
        [dd.EngEng["idState"],
         dd.EngEng["state"],
         dd.EngEng["date"]]).agg(
        {dd.EngEng["ICU"]: sum, dd.EngEng["ICU_ventilated"]: sum})
    df_states = df_states.reset_index()
    df_states.sort_index(axis=1, inplace=True)

    filename = "state_divi"
    filename = gd.append_filename(filename, impute_dates, moving_average)
    gd.write_dataframe(df_states, directory, filename, file_format)

    # write data for germany to file
    df_ger = df.groupby(["Date"]).agg({"ICU": sum, "ICU_ventilated": sum})
    df_ger = df_ger.reset_index()
    df_ger.sort_index(axis=1, inplace=True)

    filename = "germany_divi"
    filename = gd.append_filename(filename, impute_dates, moving_average)
    gd.write_dataframe(df_ger, directory, filename, file_format)

    return(df_raw, df_counties, df_states, df_ger)


def gdd_sanity_checks(df = pd.DataFrame()):

    # get actual headers
    try:
        actual_strings_list = df.columns.tolist()
    except:
        exit_string = "Error: no dataframe given."
        sys.exit(exit_string)

    if len(actual_strings_list) != 11:
        exit_string = "Error: Number of data categories changed."
        sys.exit(exit_string)

    # These strings need to be in the header
    test_strings = {
        "date", "bundesland", "gemeindeschluessel", "faelle_covid_aktuell",
        "faelle_covid_aktuell_invasiv_beatmet"}

    # Compare
    for name in test_strings:
        if(name not in actual_strings_list):
            exit_string = "Error: Data categories have changed."
            sys.exit(exit_string)

    if (len(df.index) < 200000) or (len(df.index) > 500000):
        exit_string = "Error: unexpected length of dataframe."
        sys.exit(exit_string)


def main():
    """ Main program entry."""

    arg_dict = gd.cli('divi',)
    get_divi_data(**arg_dict)


if __name__ == "__main__":
    main()
