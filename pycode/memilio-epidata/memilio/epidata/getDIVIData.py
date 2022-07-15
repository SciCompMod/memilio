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
"""

import os
from datetime import date
import pandas as pd

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import modifyDataframeSeries


def get_divi_data(data_folder,
                  read_data=dd.defaultDict['read_data'],
                  file_format=dd.defaultDict['file_format'],
                  no_raw=dd.defaultDict['no_raw'],
                  end_date=dd.defaultDict['end_date'],
                  start_date=dd.defaultDict['start_date'],
                  impute_dates=dd.defaultDict['impute_dates'],
                  moving_average=dd.defaultDict['moving_average']
                  ):
    """! Downloads or reads the DIVI ICU data and writes them in different files.

    Available data starts from 2020-04-24.
    If the given start_date is earlier, it is changed to this date and a warning is printed.
    If it does not already exist, the folder Germany is generated in the given data_folder.
    If read_data == True and the file "FullData_DIVI.json" exists, the data is read form this file
    and stored in a pandas dataframe. If read_data = True and the file does not exist the program is stopped.

    The downloaded dataframe is written to the file "FullData_DIVI".
    After that, the columns are renamed to English and the state and county names are added.
    Afterwards, three kinds of structuring of the data are done.
    We obtain the chronological sequence of ICU and ICU_ventilated
    stored in the files "county_divi".json", "state_divi.json" and "germany_divi.json"
    for counties, states and whole Germany, respectively.

    @param data_folder Folder where data is written to.
    @param read_data False [Default] or True. Defines if data is read from file or downloaded.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    "False [Default]" if it is downloaded for all dates from start_date to end_date.
    @param no_raw True or False [Default]. Defines if unchanged raw data is saved or not.
    @param start_date [Optional] Date of first date in dataframe. Default defined in defaultDict.
    @param end_date [Optional] Date of last date in dataframe. Default defined in defaultDict.
    @param impute_dates True or False [Default]. Defines if values for dates without new information are imputed.
    @param moving_average 0 [Default] or >0. Applies an 'moving_average'-days moving average on all time series
        to smooth out weekend effects.
    """

    # First csv data on 24-04-2020
    if start_date < date(2020, 4, 24):
        print("Warning: First data available on 2020-04-24. "
              "You asked for " + start_date.strftime("%Y-%m-%d") + ".")
        start_date = date(2020, 4, 24)

    directory = os.path.join(data_folder, 'Germany/')
    gd.check_dir(directory)

    filename = "FullData_DIVI"

    if read_data:
        # read json file for already downloaded data
        file_in = os.path.join(directory, filename + ".json")

        try:
            df_raw = pd.read_json(file_in)
        except ValueError:
            raise FileNotFoundError("Error: The file: " + file_in +
                                    " does not exist. Call program without"
                                    " -r flag to get it.")
    else:
        try:
            df_raw = gd.loadCsv(
                'zeitreihe-tagesdaten',
                apiUrl='https://diviexchange.blob.core.windows.net/%24web/',
                extension='.csv')
        except Exception as err:
            raise FileNotFoundError(
                "Error: Download link for Divi data has changed.") from err

    if not df_raw.empty:
        if not no_raw:
            gd.write_dataframe(df_raw, directory, filename, file_format)
    else:
        raise gd.DataError("Something went wrong, dataframe is empty.")
    df = df_raw.copy()
    divi_data_sanity_checks(df_raw)
    df.rename(columns={'date': dd.EngEng['date']}, inplace=True)
    df.rename(dd.GerEng, axis=1, inplace=True)

    df[dd.EngEng['date']] = pd.to_datetime(
        df[dd.EngEng['date']], format='%Y-%m-%d %H:%M:%S')
    df = modifyDataframeSeries.extract_subframe_based_on_dates(
        df, start_date, end_date)

    # remove leading zeros for ID_County (if not yet done)
    df['ID_County'] = df['ID_County'].astype(int)
    # add missing dates (and compute moving average)
    if (impute_dates == True) or (moving_average > 0):
        df = modifyDataframeSeries.impute_and_reduce_df(
            df,
            {dd.EngEng["idCounty"]: geoger.get_county_ids()},
            [dd.EngEng["ICU"], dd.EngEng["ICU_ventilated"]],
            impute='forward', moving_average=moving_average)

    # add names etc for empty frames (counties where no ICU beds are available)
    countyid_to_stateid = geoger.get_countyid_to_stateid_map()
    for id in df.loc[df.isnull().any(axis=1), dd.EngEng['idCounty']].unique():
        stateid = countyid_to_stateid[id]
        df.loc[df[dd.EngEng['idCounty']] == id, dd.EngEng['idState']] = stateid

    df = geoger.insert_names_of_states(df)
    df = geoger.insert_names_of_counties(df)

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


def divi_data_sanity_checks(df=pd.DataFrame()):
    """! Checks the sanity of the divi_data dataframe

    Checks if type of the given data is a dataframe
    Checks if the headers of the dataframe are those which are needed
    Checks if the size of the dataframe is not unusual

    @param df The dataframe which has to be checked
    """
    # get actual headers
    actual_strings_list = df.columns.tolist()
    # check number of data categories
    if len(actual_strings_list) != 11:
        raise gd.DataError("Error: Number of data categories changed.")

    # These strings need to be in the header
    test_strings = {
        "date", "bundesland", "gemeindeschluessel", "faelle_covid_aktuell",
        "faelle_covid_aktuell_invasiv_beatmet"}

    # check if headers are those we want
    for name in test_strings:
        if(name not in actual_strings_list):
            raise gd.DataError("Error: Data categories have changed.")
    # check if size of dataframe is not unusal
    # data colletion starts at 24.04.2020
    num_dates = (date.today() - date(2020, 4, 24)).days
    min_num_data = 390*num_dates  # not all 400 counties report every day
    max_num_data = 400*num_dates
    if (len(df) < min_num_data) or (len(df) > max_num_data):
        raise gd.DataError("Error: unexpected length of dataframe.")


def main():
    """ Main program entry."""

    path = os.path.join(os.getcwd(), 'data', 'pydata')
    arg_dict = gd.cli('divi',)
    get_divi_data(path, **arg_dict)


if __name__ == "__main__":
    main()
