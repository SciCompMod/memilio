#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Kathrin Rack, Lena Ploetzke
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
from datetime import timedelta, date
import pandas

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd


def adjust_data(df, date_of_data):
    """!Adjusts data such that the data structure looks the same for all dates

    For the first available dates of year 2020, there are some differences and thus special treatment is needed:
    - 24.4.: Rename column 'kreis' to 'gemeindeschluessel'
    - 25.4.: Add id bundesland, which is extracted from the given column "gemeindeschluessel"
    - 24.4.-27.4.: Add date as a column
    - 29.9.: Empty column has to be removed

    @param df A pandas data frame
    @param date_of_data The date for the data stored in df
    @return changed pandas dataframe
    """

    # rename column 'kreis' of first date to match data of following days
    if date_of_data == date(2020, 4, 24):
        df.rename({'kreis': 'gemeindeschluessel'}, axis=1, inplace=True)
        # Add empty column ICU_ventilated
        df.insert(loc=len(df.columns), column='faelle_covid_aktuell_beatmet', value='')

    # remove empty column
    if date_of_data == date(2020, 4, 29):
        df.drop(columns='Unnamed: 0', inplace=True)

    # add dates for data until 27.4.
    if date_of_data <= date(2020, 4, 27):
        date_str = date_of_data.strftime("%Y-%m-%d") + " 09:15:00"
        df.insert(loc=len(df.columns), column='daten_stand', value=date_str)

    # add 'bundesland' for data from 25.4.
    if date_of_data == date(2020, 4, 25):
        df.insert(loc=0, column='bundesland', value=df['gemeindeschluessel'].floordiv(1000))
        df.insert(loc=len(df.columns), column='faelle_covid_aktuell_beatmet', value='')

    return df


def call_call_url(url_prefix, call_number):
    """!Calls, if possible a url and breaks if it is not possible to call the url,
    e.g., due to a not working internet connection
    or returns an empty dataframe if url is not the correct one.

    @param url_prefix Date specific part of the url
    @param call_number Irregular number which is needed for the url

    @return pandas dataframe which is either empty or contains requested data
    """

    call_url = "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv" + "/viewdocument/" \
               + str(call_number) + "/divi-intensivregister-" + url_prefix

    # empty data frame
    df = pandas.DataFrame()

    try:
        df = pandas.read_csv(call_url)
    except:
        pass

    return df


def nearest_earlier_date(date_list, date_given):
    """!Finds nearest but earlier date to a given date from a list of dates.

    @param date_list Iterable object of dates
    @param date_given Single date

    @return date from date_list which is closest but earlier to date_given.
    """

    # should be sorted but for case
    date_list.sort()
    # find all dates before date_given
    index = bisect.bisect(date_list, date_given)
    # get date which is closest
    return min(date_list[:index], key=lambda x: abs(x - date_given))


def download_data_for_one_day(last_number, download_date):
    """!Downloads data for given date

    This function checks if the given date is a key of a dictionary where dates
    and their corresponding call_numbers are stored.
    This is stored there if the difference of the call_number to the call_number of the date before is larger than
    1 or 2.
    If the date is not part of the call_number_dict the new call_number is calculated.
    First, the last_number is successively increased by one until a difference of 300 is reached.
    Then, the last_number is successively decreased by one until a difference of 300 is reached.
    At last the last_number without a change is tried.
    If data could be downloaded by the function call_call_url and the difference was 1 or 2,
    the data is simply given back.
    If the difference is different an additional message is printed in the end of the program,
    which can be used to copy the date and the corresponding call_number directly into the call_number_dict
    to decrease runtime of the program.
    Furthermore, in this function another specific part of the url, the call_time, is estimated from the given date.
    If the date is before 2020-6-5 the time "-09-15" and afterwards "-12-15" has to be added to the date.
    Moreover for dates in the range [2020-6-12, 2020-6-5], an additional "-2" has to be added after the call_time.
    Thus, the url_prefix = call_date + call_time + ext and is then given to the call_call_url fct.

    @param last_number This is the call_number which is needed for the download url
    of the date 1 day before this one
    @param download_date The date for which the data should be downloaded
    @return List of call_number of the download_data, the pandas dataframe, and a string which is either empty or
    contains what should be added to "call_number_dict"
    """

    # define call numbers for dates where call number doesn't increase by 1
    call_number_dict = {date(2020, 4, 24): 3974,
                        date(2020, 5, 6): 3691,
                        date(2020, 6, 5): 3842,
                        date(2020, 6, 12): 3906,
                        date(2020, 6, 26): 3774,
                        date(2020, 6, 28): 3777,
                        date(2020, 6, 29): 3830,
                        date(2020, 6, 30): 3840,
                        date(2020, 7, 1): 3928,
                        date(2020, 7, 2): 3951,
                        date(2020, 7, 6): 3958,
                        date(2020, 7, 7): 3961,
                        date(2020, 7, 8): 3964,
                        date(2020, 7, 10): 3987,
                        date(2020, 7, 11): 4065,
                        date(2020, 7, 12): 4069,
                        date(2020, 7, 13): 4073,
                        date(2020, 7, 14): 4081,
                        date(2020, 7, 15): 4108,
                        date(2020, 7, 17): 4127,
                        date(2020, 7, 18): 4143,
                        date(2020, 7, 20): 4156,
                        date(2020, 7, 22): 4162,
                        date(2020, 7, 27): 4173,
                        date(2020, 7, 28): 4177,
                        date(2020, 7, 29): 4181,
                        date(2020, 7, 30): 4184,
                        date(2020, 8, 3): 4204,
                        date(2020, 8, 4): 4316,
                        date(2020, 8, 5): 4319,
                        date(2020, 8, 6): 4324,
                        date(2020, 8, 8): 4327,
                        date(2020, 8, 9): 4330,
                        date(2020, 8, 11): 4412,
                        date(2020, 8, 13): 4430,
                        date(2020, 8, 14): 4680,
                        date(2020, 8, 15): 4909,
                        date(2020, 8, 17): 4977,
                        date(2020, 8, 18): 4996,
                        date(2020, 8, 24): 5008,
                        date(2020, 8, 25): 5011,
                        date(2020, 8, 28): 5020,
                        date(2020, 8, 31): 5026,
                        date(2020, 9, 5): 5036,
                        date(2020, 9, 10): 5048,
                        date(2020, 9, 12): 5052,
                        date(2020, 9, 15): 5058,
                        date(2020, 9, 16): 5073,
                        date(2020, 9, 18): 5082,
                        date(2020, 9, 20): 5087,
                        date(2020, 9, 23): 5093,
                        date(2020, 9, 28): 5103,
                        date(2020, 10, 4): 5115,
                        date(2020, 10, 6): 5120,
                        date(2020, 10, 14): 5136,
                        date(2020, 10, 17): 5142,
                        date(2020, 10, 20): 5149,
                        date(2020, 10, 21): 5152,
                        date(2020, 10, 23): 5158,
                        date(2020, 10, 25): 5163,
                        date(2020, 10, 28): 5174,
                        date(2020, 10, 29): 5180,
                        date(2020, 11, 6): 5199,
                        date(2020, 11, 12): 5212,
                        date(2020, 11, 14): 5219,
                        date(2020, 11, 16): 5224,
                        date(2020, 11, 17): 5229,
                        date(2020, 11, 18): 5232,
                        date(2020, 11, 20): 5238,
                        date(2020, 11, 21): 5253,
                        date(2020, 11, 24): 5260,
                        date(2020, 11, 25): 5263,
                        date(2020, 12, 2): 5286,
                        date(2020, 12, 3): 5289,
                        date(2020, 12, 4): 5292,
                        date(2020, 12, 8): 5300,
                        date(2020, 12, 11): 5308,
                        date(2020, 12, 17): 5321,
                        date(2020, 12, 19): 5326,
                        date(2020, 12, 22): 5333,
                        date(2020, 12, 30): 5350,
                        }
    start_date_differs = False

    call_date = download_date.strftime("%Y-%m-%d")

    # first links have different upload time
    time_shift_date = date(2020, 6, 5)
    if download_date < time_shift_date:
        call_time = "-09-15"
    elif download_date == date(2020, 9, 14):
        call_time = "-14-15"
    else:
        call_time = "-12-15"

    # need extension "-2" for dates between 12.6. and 25.6.
    ext = ""
    if date(2020, 6, 12) <= download_date <= date(2020, 6, 25):
        ext = "-2"

    url_prefix = call_date + call_time + ext

    # construction of link is different for different time periods
    # for the latest dates no call number is needed
    # number in url starts at values in call_number_dict and increases every day by 1

    sign_dict = {0: 1,
                 1: -1}

    call_string = ""

    if download_date in call_number_dict.keys():

        call_number = call_number_dict[download_date]

        df = call_call_url(url_prefix, call_number)

        # if no data has been added, but numbers are in dict something else is wrong, e.g. internet connection
        if df.empty:
            exit_string = "Something went wrong with download of data for date " + str(download_date) \
                          + ", although it is part of the call_number_dict."
            sys.exit(exit_string)
    else:
        # case where start_date is not 24-04-2020 and is not in dict
        # then we search to date which is closest to given date and smaller
        if last_number == 0:
            nd = nearest_earlier_date(list(call_number_dict.keys()), download_date)
            last_number = call_number_dict[nd]
            # make sure no output about new drifting number is added for this case
            start_date_differs = True

        # It is most likely that the difference is between 1 and 2
        for sign in range(2):
            for delta in range(1, 300):
                call_number = last_number + sign_dict[sign]*delta

                # for delta 1 and 2 the number is not saved in dict,
                # because it does not take so long to get those numbers
                if sign == 0 and delta != 1 and delta != 2 and not start_date_differs:
                    call_string = "date(" + download_date.strftime("%Y, %-m, %-d") + "): " + str(call_number) + "," \
                                  + "\n"

                df = call_call_url(url_prefix, call_number)

                if not df.empty:
                    return [call_number, df, call_string]

        # case with same call_number, which is very unlikely
        call_number = last_number
        call_string = "date(" + download_date.strftime("%Y, %-m, %-d") + "): " + str(call_number) + "," + "\n"
        df = call_call_url(url_prefix, call_number)

    return [call_number, df, call_string]


def get_divi_data(read_data=dd.defaultDict['read_data'],
                  file_format=dd.defaultDict['file_format'],
                  out_folder=dd.defaultDict['out_folder'],
                  no_raw=dd.defaultDict['no_raw'],
                  end_date=dd.defaultDict['end_date'],
                  start_date=dd.defaultDict['start_date'],
                  update_data=dd.defaultDict['update_data'],
                  ):
    """!Function to get the divi data.

    Available data start from 2020-4-24.
    If the given start_date is earlier, it is changed to this date and a warning is printed.
    If it does not already exist, the folder Germany is generated in the given out_folder.
    If read_data == True and the file "FullData_DIVI.json" exists, the data is read form this file
    and stored in a pandas dataframe.
    Otherwise the program is stopped.

    If update_data == True  the same happens as for read_data == True.
    Furthermore, if the last date in read data is yesterday, the data of today is downloaded.
    For this download an easier link can be used than the one where we have to find the mysterious call_number.
    If the data has not yet been uploaded, the program is stopped with a message to try again later.
    If there is more data missing than today, the parameter start_date is changed to the first missing data,
    and the data is normally downloaded (see below).

    If data should normally be downloaded between start_date and end_date, we start with an empty pandas dataframe.
    Afterwards, for every day between start_date and end_date, both included,
    the function download_data_for_one_day is called.
    If a given back dataframe is empty, a warning is printed that this date is missing, but the program is not stopped.
    If start_date is earlier or equal 2020-04-29, the function adjust_data has to be called.
    If data has been downloaded, the dataframe of this one date is added to the dataframe with all dates.

    If the dataframe which should contain all data is empty after going through all dates, the program is stopped.
    Otherwise the dataframe is written to the file filename = "FullData_DIVI".
    Than the columns are renamed to English and the state and county names are added.
    Afterwards, three kinds of structuring of the data are done.
    We obtain the chronological sequence of ICU and ICU_ventilated
    stored in the files "county_divi".json", "state_divi.json"and "germany_divi.json"
    for counties, states and whole Germany, respectively.

    @param read_data False [Default] or True. Defines if data is read from file or downloaded.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    @param update_data "True" if existing data is updated or
    "False [Default]" if it is downloaded for all dates from start_date to end_date.
    @param out_folder Folder where data is written to.
    @param no_raw True or False [Default]. Defines if unchanged raw data is saved or not.
    @param start_date [Optional] Date to start to download data [Default = 2020.4.24].
    @param end_date [Optional] Date to stop to download data [Default = today].
    """

    delta = timedelta(days=1)
    today = date.today()

    # First csv data on 24-04-2020
    if start_date < date(2020, 4, 24):
        print("Warning: First data available on 2020-04-24. "
              "You asked for " + start_date.strftime("%Y-%m-%d") + ".")
        start_date = date(2020, 4, 24)

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    filename = "FullData_DIVI"

    if update_data or read_data:
        # read json file for already downloaded data
        file_in = os.path.join(directory, filename + ".json")

        try:
            df = pandas.read_json(file_in)
        except ValueError:
            exit_string = "Error: The file: " + file_in + " does not exist. "\
                          "Call program without -r or -u flag to get it."
            sys.exit(exit_string)

        # for read_data no data download is needed;
        # while-loop will be skipped by setting start_date to larger than end_date
        start_date = end_date + delta

        if update_data:
            if not df.empty:
                newest_date = pandas.to_datetime(df['daten_stand']).max().date()

                if (today - delta) == newest_date:
                    # just today's data is missing
                    # download data from today
                    # while loop to download further data will be skipped
                    df2 = gd.loadCsv('DIVI-Intensivregister-Tagesreport', apiUrl='https://www.divi.de/')

                    if not df2.empty:
                        # test if online data is already the one of today
                        download_date = pandas.to_datetime(df2['daten_stand']).max().date()

                        if download_date == today:
                            df = df.append(df2, ignore_index=True)
                            print("Success: Data of date " + today.strftime("%Y-%m-%d")
                                  + " has been added to dataframe")
                        else:
                            exit_string = "Data of today = " + today.strftime("%Y-%m-%d") \
                                          + " has not yet uploaded. Please, try again later."
                            sys.exit(exit_string)

                elif (today - newest_date).days > delta.days:
                    # more than today's data is missing
                    # start with the oldest missing data
                    start_date = newest_date + delta
    else:
        # Get all data:
        # start with empty dataframe
        df = pandas.DataFrame()

    last_number = 0

    new_dict_string = ""
    while start_date <= end_date:

        [last_number, df2, new_string] = download_data_for_one_day(last_number, start_date)

        if not df2.empty:

            new_dict_string = new_dict_string + new_string

            # data of first days needs adjustment
            if start_date <= date(2020, 4, 29):
                df2 = adjust_data(df2, start_date)
            df = df.append(df2, ignore_index=True)

            print("Success: Data of date " + start_date.strftime("%Y-%m-%d") + " has been added to dataframe")
        else:
            print("Warning: Data of date " + start_date.strftime("%Y-%m-%d") + " is not added to dataframe")

        start_date += delta

    # output data before renaming
    if not df.empty:
        if not new_dict_string == "":
            print("New drifting number in link found. "
                  "To decrease runtime, please copy the following "
                  "to the dictionary \"call_number_dict\" in the function \"download_data_for_one_day\": ")
            print(new_dict_string)

        if not no_raw:
            gd.write_dataframe(df, directory, filename, file_format)
    else:
        exit_string = "Something went wrong, dataframe is empty."
        sys.exit(exit_string)

    # change column names
    df.rename(dd.GerEng, axis=1, inplace=True)

    df.Date = pandas.to_datetime(df.Date, format='%Y-%m-%d %H:%M')

    # insert names of  states
    df.insert(loc=0, column="State", value=df.ID_State)
    for key, value in dd.State.items():
        df.loc[df["State"] == key, ["State"]] = value

    # insert names of counties
    df.insert(loc=0, column="County", value=df.ID_County)
    for key, value in dd.County.items():
        df.loc[df["County"] == key, ["County"]] = value

    # write data for counties to file
    df_counties = df[["County", "ID_County", "ICU", "ICU_ventilated", "Date"]].copy()
    filename = "county_divi"
    gd.write_dataframe(df_counties, directory, filename, file_format)

    # write data for states to file
    df_states = df.groupby(["ID_State", "State", "Date"]).agg({"ICU": sum, "ICU_ventilated": sum})
    df_states = df_states.reset_index()
    df_states.sort_index(axis=1, inplace=True)
    # For the sum calculation Nan is used as a 0, thus some zeros have to be changed back to NaN
    # df_states.loc[df_states["Date"] <= "2020-04-25 09:15:00", "ICU_ventilated"] = np.nan
    # TODO: Fill "faelle_covid_aktuell_im_bundesland" into ICU data
    # print(df_states.loc[df_states["Date"] == "2020-04-24 09:15:00", "ICU"])
    # print(df.groupby(["ID_State", "State"]).agg({"faelle_covid_aktuell_im_bundesland": max}).reset_index()
    # [["faelle_covid_aktuell_im_bundesland"]])
    # df_states.loc[df_states["Date"] == "2020-04-24 09:15:00", "ICU"] = df_states.loc[df_states["Date"]
    # == "2020-04-24 09:15:00", "ICU"].replace(df.groupby(["ID_State", "State"]).agg(
    # {"faelle_covid_aktuell_im_bundesland": max})[["faelle_covid_aktuell_im_bundesland"]])

    filename = "state_divi"
    gd.write_dataframe(df_states, directory, filename, file_format)

    # write data for germany to file
    df_ger = df.groupby(["Date"]).agg({"ICU": sum, "ICU_ventilated": sum})
    df_ger = df_ger.reset_index()
    df_ger.sort_index(axis=1, inplace=True)
    # for performing sum nan is treated like an zero
    # However, afterwards ist should be filled with NaN
    # TODO: should we have this? See issue #205
    # df_ger.loc[df_states["Date"] <= "2020-04-25 09:15:00", "ICU_ventilated"] = np.nan
    # TODO: Use also "faelle_covid_aktuell_im_bundesland" from 25.9.
    filename = "germany_divi"
    gd.write_dataframe(df_ger, directory, filename, file_format)

def main():
    """ Main program entry."""

    arg_dict = gd.cli('divi',)
    get_divi_data(**arg_dict)


if __name__ == "__main__":
    main()
