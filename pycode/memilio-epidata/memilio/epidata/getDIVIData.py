#############################################################################
# Copyright (C) 2020-2025 MEmilio
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
:strong:`getDIVIData.py`

Data of the DIVI
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
from typing import Tuple, Dict

import pandas as pd

from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import modifyDataframeSeries as mdfs


def fetch_divi_data(
        directory: str,
        filename: str,
        conf_obj,
        read_data: bool = dd.defaultDict['read_data'],
        file_format: str = dd.defaultDict['file_format'],
) -> pd.DataFrame:
    """ Downloads or reads the DIVI ICU data and writes them in different files.

    If it does not already exist, the folder Germany is generated in the given out_folder.
    If read_data == True and the file "FullData_DIVI.json" exists, the data is read form this file
    and stored in a pandas dataframe. If read_data = True and the file does not exist the program is stopped.
    The downloaded dataframe is written to the file "FullData_DIVI".

    :param directory: str
        Path to the output directory
    :param conf_obj: configuration object
    :param filename: str
        File format which is used for writing the data. Default defined in defaultDict.
    :param read_data: bool. True or False. Defines if data is read from file or downloaded. Default defined in defaultDict. (Default value = dd.defaultDict['read_data'])
    :param file_format: str. File format which is used for writing the data. Default defined in defaultDict. (Default value = dd.defaultDict['file_format'])
    :returns: Tuple[df_raw, start_date] Tuple. Contains the fetched data as well as the adjusted starting date

    """
    no_raw = conf_obj.no_raw

    url = "https://raw.githubusercontent.com/robert-koch-institut/" \
          "Intensivkapazitaeten_und_COVID-19-Intensivbettenbelegung_in_Deutschland/" \
          "main/Intensivregister_Landkreise_Kapazitaeten.csv"
    path = os.path.join(directory + filename + ".json")
    df_raw = gd.get_file(path, url, read_data, param_dict={},
                         interactive=conf_obj.interactive)
    if not df_raw.empty:
        if not no_raw:
            gd.write_dataframe(df_raw, directory, filename, file_format)
    else:
        raise gd.DataError("Something went wrong, dataframe is empty.")
    if conf_obj.checks is True:
        divi_data_sanity_checks(df_raw)
    else:
        gd.default_print(
            "Warning", "Sanity checks for DIVI data have not been executed.")
    return df_raw


def preprocess_divi_data(df_raw: pd.DataFrame,
                         conf_obj,
                         start_date: date = date(2020, 4, 24),
                         end_date: date = dd.defaultDict['end_date'],
                         impute_dates: bool = dd.defaultDict['impute_dates'],
                         moving_average: int = dd.defaultDict['moving_average'],
                         ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """ Processing of the downloaded data
        * the columns are renamed to English and the state and county names are added.

    :param df_raw: pd.DataFrame
    :param conf_obj: configuration object
    :param start_date: date  The first date in dataframe. Default value = date(2020, 4, 24).
    :param end_date: date The last date in dataframe. Default defined in defaultDict. (Default value = dd.defaultDict['end_date'])
    :param impute_dates: bool  Defines if values for dates without new information are imputed. Default defined in defaultDict. (Default value = dd.defaultDict['impute_dates'])
    :param moving_average: int  Integers >=0.Applies an 'moving_average'-days moving average on all time seriesto smooth out effects of irregular reporting. Default defined in defaultDict. (Default value = dd.defaultDict['moving_average'])
    :returns: df pd.DataFrame  processed divi data

    """
    # First csv data on 24-04-2020
    if start_date < date(2020, 4, 24):
        gd.default_print(
            'Warning',
            "First data available on 2020-04-24. "
            "You asked for " +
            start_date.strftime("%Y-%m-%d") +
            ". Changed it to 2020-04-24.")
        start_date = date(2020, 4, 24)
     # Dataset will no longer be updated from CW29 2024 on.
    if end_date > date(2024, 7, 21):
        gd.default_print(
            'Warning',
            "Dataset will no longer be updated after 2024-07-21. "
            "You asked for " +
            start_date.strftime("%Y-%m-%d") +
            ". Use the provided data carefully.")

    if conf_obj.checks is True:
        divi_data_sanity_checks(df_raw)
    else:
        gd.default_print(
            "Warning", "Sanity checks for DIVI data have not been executed.")
    df = df_raw.rename(dd.GerEng, axis=1, inplace=False)

    try:
        df[dd.EngEng['date']] = pd.to_datetime(
            df[dd.EngEng['date']], format="ISO8601")
    except ValueError:
        try:
            df[dd.EngEng['date']] = pd.to_datetime(
                df[dd.EngEng['date']], format="%Y-%m-%d %H:%M:%S")
        except BaseException:
            raise gd.DataError(
                "Time data can't be transformed to intended format")

    # remove leading zeros for ID_County (if not yet done)
    df['ID_County'] = df['ID_County'].astype(int)
    # add missing dates (and compute moving average)
    if (impute_dates is True) or (moving_average > 0):
        df = mdfs.impute_and_reduce_df(
            df, {dd.EngEng["idCounty"]: df[dd.EngEng["idCounty"]].unique()},
            [dd.EngEng["ICU"],
             dd.EngEng["ICU_ventilated"]],
            impute='forward', moving_average=moving_average,
            min_date=start_date, max_date=end_date)

    # add names etc for empty frames (counties where no ICU beds are available)
    countyid_to_stateid = geoger.get_countyid_to_stateid_map()
    for county_id in df.loc[df.isna().any(axis=1), dd.EngEng['idCounty']].unique():
        state_id = countyid_to_stateid[county_id]
        df.loc[df[dd.EngEng['idCounty']] == county_id,
               dd.EngEng['idState']] = state_id

    # extract subframe of dates
    df = mdfs.extract_subframe_based_on_dates(df, start_date, end_date)

    return df, df_raw


def write_divi_data(df: pd.DataFrame,
                    directory: str,
                    conf_obj,
                    file_format: str = dd.defaultDict['file_format'],
                    impute_dates: bool = dd.defaultDict['impute_dates'],
                    moving_average: int = dd.defaultDict['moving_average'],
                    ) -> Dict:
    """ Write the divi data into json files

    Three kinds of structuring of the data are done.
    We obtain the chronological sequence of ICU and ICU_ventilated
    stored in the files "county_divi".json", "state_divi.json" and "germany_divi.json"
    for counties, states and whole Germany, respectively.

    :param df: pd.DataFrame. Dataframe containing processed divi data
    :param directory: str
        Path to the output directory
    :param conf_obj: configuration object
    :param file_format: str. File format which is used for writing the data. Default defined in defaultDict. (Default value = dd.defaultDict['file_format'])
    :param impute_dates: bool True or False. Defines if values for dates without new information are imputed. Default defined in defaultDict. (Default value = dd.defaultDict['impute_dates'])
    :param moving_average: int Integers >=0. Applies an 'moving_average'-days moving average on all time series to smooth out effects of irregular reporting. Default defined in defaultDict. (Default value = dd.defaultDict['moving_average'])
    :returns: data_dict Dict Dictionary containing datasets at county, state and national level

    """
    # write data for counties to file
    df_counties = df[[dd.EngEng["idCounty"],
                      dd.EngEng["county"],
                      dd.EngEng["ICU"],
                      dd.EngEng["ICU_ventilated"],
                      dd.EngEng["date"]]]
    # merge Eisenach and Wartburgkreis from DIVI data
    df_counties = geoger.merge_df_counties_all(
        df_counties, sorting=[dd.EngEng["idCounty"], dd.EngEng["date"]])
    # save

    # write data for states to file
    df_states = df.groupby(
        [dd.EngEng["idState"],
         dd.EngEng["state"],
         dd.EngEng["date"]]).agg(
        {dd.EngEng["ICU"]: "sum", dd.EngEng["ICU_ventilated"]: "sum"})
    df_states.reset_index(inplace=True)
    df_states.sort_index(axis=1, inplace=True)

    # write data for germany to file
    df_ger = df.groupby(["Date"]).agg({"ICU": "sum", "ICU_ventilated": "sum"})
    df_ger.reset_index(inplace=True)
    df_ger.sort_index(axis=1, inplace=True)

    if not conf_obj.to_dataset:
        filename = "county_divi"
        filename = gd.append_filename(filename, impute_dates, moving_average)
        gd.write_dataframe(df_counties, directory, filename, file_format)

        filename = "state_divi"
        filename = gd.append_filename(filename, impute_dates, moving_average)
        gd.write_dataframe(df_states, directory, filename, file_format)

        filename = "germany_divi"
        filename = gd.append_filename(filename, impute_dates, moving_average)
        gd.write_dataframe(df_ger, directory, filename, file_format)

    data_dict = {
        "counties": df_counties,
        "states": df_states,
        "Germany": df_ger
    }
    return data_dict


def get_divi_data(read_data: bool = dd.defaultDict['read_data'],
                  file_format: str = dd.defaultDict['file_format'],
                  out_folder: str = dd.defaultDict['out_folder'],
                  start_date: date = date(2020, 4, 24),
                  end_date: date = dd.defaultDict['end_date'],
                  impute_dates: bool = dd.defaultDict['impute_dates'],
                  moving_average: int = dd.defaultDict['moving_average'],
                  **kwargs
                  ):
    """ Downloads or reads the DIVI ICU data and writes them in different files.

    Available data starts from 2020-04-24.
    If the given start_date is earlier, it is changed to this date and a warning is printed.
    It has been announced that the dataset will no longer be updated from 2024-07-21 (CW29).
    If end_date is later, a warning is displayed.
    If it does not already exist, the folder Germany is generated in the given out_folder.
    If read_data == True and the file "FullData_DIVI.json" exists, the data is read form this file
    and stored in a pandas dataframe. If read_data = True and the file does not exist the program is stopped.

    The downloaded dataframe is written to the file "FullData_DIVI".
    After that, the columns are renamed to English and the state and county names are added.
    Afterwards, three kinds of structuring of the data are done.
    We obtain the chronological sequence of ICU and ICU_ventilated
    stored in the files "county_divi".json", "state_divi.json" and "germany_divi.json"
    for counties, states and whole Germany, respectively.

    :param read_data: True or False. Defines if data is read from file or downloaded. Default defined in defaultDict. (Default value = dd.defaultDict['read_data'])
    :param file_format: File format which is used for writing the data. Default defined in defaultDict. (Default value = dd.defaultDict['file_format'])
    :param out_folder: Folder where data is written to. Default defined in defaultDict. (Default value = dd.defaultDict['out_folder'])
    :param start_date: Date of first date in dataframe. Default value = ``date(2020, 4, 24)``.
    :param end_date: Date of last date in dataframe. Default defined in defaultDict. (Default value = dd.defaultDict['end_date'])
    :param impute_dates: True or False. Defines if values for dates without new information are imputed. Default defined in defaultDict. (Default value = dd.defaultDict['impute_dates'])
    :param moving_average: Integers >=0. Applies an 'moving_average'-days moving average on all time series
        to smooth out effects of irregular reporting. Default defined in defaultDict. (Default value = dd.defaultDict['moving_average'])
    :param **kwargs: 

    """
    conf = gd.Conf(out_folder, **kwargs)
    out_folder = conf.path_to_use

    directory = os.path.join(out_folder, 'Germany', 'pydata')
    gd.check_dir(directory)

    filename = "FullData_DIVI"

    downloaded_data_df = fetch_divi_data(
        directory=directory,
        conf_obj=conf,
        filename=filename,
        read_data=read_data,
        file_format=file_format,
    )

    preprocess_df, df_raw = preprocess_divi_data(
        conf_obj=conf,
        df_raw=downloaded_data_df,
        start_date=start_date,
        end_date=end_date,
        impute_dates=impute_dates,
        moving_average=moving_average,
    )
    datasets = write_divi_data(
        df=preprocess_df,
        directory=directory,
        file_format=file_format,
        impute_dates=impute_dates,
        moving_average=moving_average,
        conf_obj=conf,
    )
    datasets['raw_data'] = df_raw
    return datasets


def divi_data_sanity_checks(df: pd.DataFrame) -> None:
    """ Checks the sanity of the divi_data dataframe

    Checks if type of the given data is a dataframe
    Checks if the headers of the dataframe are those which are needed
    Checks if the size of the dataframe is not unusual

    :param df: The dataframe which has to be checked pd.DataFrame

    """
    # there should be 13 headers
    num_headers = 13
    # get actual headers
    actual_strings_list = df.columns.tolist()
    # check number of data categories
    if len(actual_strings_list) != num_headers:
        raise gd.DataError("Error: Number of data categories changed.")

    # These strings need to be in the header
    test_strings = {
        "datum", "bundesland_id", "landkreis_id", "bundesland_name",
        "landkreis_name", "faelle_covid_aktuell",
        "faelle_covid_aktuell_invasiv_beatmet"}

    # check if headers are those we want
    for name in test_strings:
        if name not in actual_strings_list:
            raise gd.DataError("Error: Data categories have changed.")
    # check if size of dataframe is not unusal
    # data colletion starts at 24.04.2020
    # TODO: Number of reporting counties get less with time.
    # Maybe we should look for a new method to sanitize the size of the
    # DataFrame.
    num_dates = (date.today() - date(2020, 4, 24)).days
    min_num_data = 380 * num_dates  # not all 400 counties report every day
    max_num_data = 400 * num_dates
    if (len(df) < min_num_data) or (len(df) > max_num_data):
        raise gd.DataError("Error: unexpected length of dataframe.")


def main():
    """Main program entry."""

    arg_dict = gd.cli('divi',)
    get_divi_data(**arg_dict)


if __name__ == "__main__":
    main()
