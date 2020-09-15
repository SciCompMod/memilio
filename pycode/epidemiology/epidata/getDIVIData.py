import os
import sys
import pandas
import numpy as np
from datetime import timedelta, date

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd

## @package getDIVIData
# Data of the DIVI -
# Deutsche interdisziplinäre Vereinigung für Intensiv- und Notfallmedizin
# about Sars-CoV2 is downloaded.
# This data contains the number of Covid19 patients in intensive care
# and the number of those that are additionally ventilated.
#
# data explanation:
# reporting_hospitals number of reporting hospitals
# ICU is the number of covid patients in reporting hospitals
# ICU_ventilated is the number of ventilated covid patients in reporting hospitals
# free_ICU is the number of free ICUs in reporting hospitals
# occupied_ICU is the number of occupied ICUs in in reporting hospitals
#
# ID_County and ID_State is as defined by the "Amtlicher Gemeindeschlüssel (AGS)"
# which is also used in the RKI data as ID_County and ID_State
# https://de.wikipedia.org/wiki/Liste_der_Landkreise_in_Deutschland
#
# Furthermore, what might be interesting about the data:
# The column "faelle_covid_im_bundesland" exits only in the data from the first day (24.4)
# The column ICU does not exist for the 24.4.
# and ICU_ventilated does not exist for the 24.4. and 25.4.


## Adjusts data such that the data structure looks the same for all dates
#
# For the first available dates of year 2020, there are some differences and thus special treatment is needed:
# - 24.4.: Rename column 'kreis' to 'gemeindeschluessel'
# - 25.4.: Add id bundesland, which is extracted from the given column "gemeindeschluessel"
# - 24.4.-27.4.: Add date as a column
# - 29.9. empty column has to be removed
#
# @param df A pandas data frame
# @param date_of_data The date for the data stored in df
# @return the changed pandas dataframe
#
def adjust_data(df, date_of_data):

    # rename column 'kreis' of first date to match data of following days
    if date_of_data == date(2020, 4, 24):
        df.rename({'kreis': 'gemeindeschluessel'}, axis=1, inplace=True)
        # Add column ICU_ventilated and fill it with Nan
        df.insert(loc=len(df.columns), column='faelle_covid_aktuell_beatmet', value=np.nan)

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
        df.insert(loc=len(df.columns), column='faelle_covid_aktuell_beatmet', value=np.nan)

    return df

## Calls, if possible a url and breaks if it is not possible to call the url,
# e.g., due to a not working internet connection
# or returns an empty dataframe if url is not the correct one.
#
# @param url_prefix Date specific part of the url
# @param call_number Irregular number which is needed for the url
#
# @return pandas dataframe which is either empty or contains requested data
#
def call_call_url(url_prefix, call_number):

    call_url = "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" \
               + url_prefix + "/viewdocument/" + str(call_number)

    # empty data frame
    df = pandas.DataFrame()

    try:
        df = pandas.read_csv(call_url)
    except OSError as e:
        exit_string = "ERROR: URL " + call_url + " could not be opened. " \
                      + "Hint: check your internet connection."

        sys.exit(exit_string)
    except:
        pass

    return df


## Downloads data for given date
#
# This function checks if the given date is a key of a dictionary where dates
# and their corresponding call_numbers are stored.
# This is stored there if the difference of the call_number to the call_number of the date before is larger than 1 or 2.
# If the date is not part of the call_number_dict the new call_number is calculated.
# First the last_number is successive increased by one until a difference of 300 is reached.
# Than the last_number is successive decreased by one until a difference of 300 is reached.
# At last the last_number without a change is tried.
# If data could be downloaded by the function call_call_url and the difference was 1 or 2.
# The data is simply given back.
# If the difference is different an additional message is printed,
# which can be used to copy directly to the call_number_dict to decrease runtime of the program.
# Furthermore, in this function another specific part of the url, the call_time, is estimated from the given date.
# If the date is before 2020-6-5 the time "-09-15" and afterwards "-12-15" has to be added to te date.
# Moreover got dates in the range [2020-6-12),2020-6-5] an additional "-2" has to be added after the call_time.
# Thus  the url_prefix = call_date + call_time + ext, which is than given to the call_call_url fct.
#
# @param last_number This is the call_number which is needed for the download url
# of the date 1 day before this one
# @param download_data The date for which the data should be downloaded
# @return List of call_number of the download_data, the pandas dataframe, and a string which is either empty or contains
# what should be added to "call_number_dict"
def download_data_for_one_day(last_number, download_date):
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
                        }

    call_date = download_date.strftime("%Y-%m-%d")

    # first links have different upload time
    time_shift_date = date(2020, 6, 5)
    if download_date < time_shift_date:
        call_time = "-09-15"
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

    else:
        # It is most likely that the difference is between 1 and 2
        for sign in range(2):
            for delta in range(1,300):

               call_number = last_number + sign_dict[sign]*delta

               # for delta 1 and 2 the number is not saved in dict,
               # because it does not take so long to get those numbers
               if sign == 0 and delta != 1 and delta != 2:
                  call_string = "date(" + download_date.strftime("%Y, %-m, %-d") + "): " + str(call_number) + "," + "\n"

               df = call_call_url(url_prefix, call_number)

               if not df.empty:
                   return [call_number, df, call_string]

        # case with same call_number, which is very unlikely
        call_number = last_number
        call_string = "date(" + download_date.strftime("%Y, %-m, %-d") + "): " + str(call_number) + "," + "\n"
        df = call_call_url(url_prefix, call_number)

    return [call_number, df, call_string]


## Function to get the divi data.
#
# Available data start from 2020-4-24.
# If the given start_date is earlier it is changed to this date and a warning is printed.
# If it does not already exist the folder Germany is generated in the given out_folder.
# If read_data == True and the file "FullData_DIVI.json" exists the data is read form this file
# and stored in a pandas dataframe.
# Otherwise the program is stopped.
#
# If update_data == True  the same happens as for read_data == True.
# Furthermore, if the last date in read data is yesterday the data of today is downloaded.
# For this download an easier link can be used than the one where we have to find the mysterious call_number
# If the data has not yet been uploaded the program is stopped with message to try again later.
# If there is more data missing than today the parameter start_date is changed to the first missing data
# and the data is normally downloaded (see below).
#
# If data should normally  be downloaded between start_date and end_date we start with an empty pandas dataframe.
# Afterwards, for everyday between start_date and end_date, both included,
# the function download_data_for_one_day is called.
# If a given back dataframe is empty, a warning is printed that this date is missing, but the program is not stopped.
# If start_date is earlier or equal 2020-04-29, the function adjust_data has to be called.
# If data has been downloaded, the dataframe of this one date is added to the dataframe with all dates.
#
# If the dataframe which should contain all data is empty after going through all dates, the program is stopped.
# Otherwise the dataframe is written to the file filename = "FullData_DIVI".
# Following the columns are renamed to English and the state and county names are added.
# Afterwards, three kind of structuring of the data are done.
# We obtain the chronological sequence of ICU and ICU_ventilated
# stored in the files "county_divi".json", "state_divi.json"and "germany_divi.json"
# for counties, states and whole Germany, respectively.
#
# @param read_data False [Default] or True. Defines if data is read from file or downloaded.
# @param update_date "True" if existing data is updated or "False [Default]" if it is downloaded.
# @param out_folder Folder where data is written to.
# @param start_date [Optional] Date to start to download data [Default = 2020.4.24].
# @param end_date [Optional] Date to stop to download data [Default = today].
#
def get_divi_data(read_data=dd.defaultDict['read_data'],
                  update_data=dd.defaultDict['update_data'],
                  make_plot=dd.defaultDict['make_plot'],
                  out_form=dd.defaultDict['out_form'],
                  out_folder=dd.defaultDict['out_folder'],
                  start_date=date(2020, 4, 24),
                  end_date=date(2020, 9, 5),
                  ):

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
                  "to the dcitionary \"call_number_dict\" in the function \"download_data_for_one_day\": ")
            print(new_dict_string)

        gd.write_dataframe(df, directory, filename, "json")
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
    df_counties = df[["County", "ID_County", "ICU", "ICU_ventilated", "Date"]]
    filename = "county_divi"
    gd.write_dataframe(df_counties, directory, filename, "json")

    # write data for states to file
    df_states = df.groupby(["ID_State", "State", "Date"]).agg({"ICU": sum, "ICU_ventilated": sum})
    df_states = df_states.reset_index()
    df_states.sort_index(axis=1, inplace=True)
    # For the sum calculation Nan is used as a 0, thus some zeros have to be changed back to NaN
    df_states.loc[df_states["Date"] <= "2020-04-25 09:15:00", "ICU_ventilated"] = np.nan
    # TODO: Fill "faelle_covid_aktuell_im_bundesland" into ICU data
    #print(df_states.loc[df_states["Date"] == "2020-04-24 09:15:00", "ICU"])
    #print(df.groupby(["ID_State", "State"]).agg({"faelle_covid_aktuell_im_bundesland": max}).reset_index()[["faelle_covid_aktuell_im_bundesland"]])
    #df_states.loc[df_states["Date"] == "2020-04-24 09:15:00", "ICU"] = df_states.loc[df_states["Date"] == "2020-04-24 09:15:00", "ICU"].replace(df.groupby(["ID_State", "State"]).agg({"faelle_covid_aktuell_im_bundesland": max})[["faelle_covid_aktuell_im_bundesland"]])

    filename = "state_divi"
    gd.write_dataframe(df_states, directory, filename, "json")

    # write data for germany to file
    df_ger = df.groupby(["Date"]).agg({"ICU": sum, "ICU_ventilated": sum})
    df_ger = df_ger.reset_index()
    df_ger.sort_index(axis=1, inplace=True)
    df_ger.loc[df_states["Date"] <= "2020-04-25 09:15:00", "ICU_ventilated"] = np.nan
    # TODO: Use also "faelle_covid_aktuell_im_bundesland" from 25.9.
    filename = "germany_divi"
    gd.write_dataframe(df_ger, directory, filename, "json")

    print("Information: DIVI data has been written to", directory)


def main():
    [read_data, update_data, make_plot, out_form, out_folder] = gd.cli('Downloads data from DIVI')
    get_divi_data(read_data, update_data, make_plot, out_form, out_folder)


if __name__ == "__main__":
    main()
