#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Martin J. Kuehn
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
import io
import os
from datetime import datetime

import pandas as pd
import requests

from memilio.epidata import customPlot
from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import modifyDataframeSeries as mdfs

# activate CoW for more predictable behaviour of pandas DataFrames
pd.options.mode.copy_on_write = True

# Downloads testing data from RKI


def download_testing_data():
    """ Downloads the Sars-CoV-2 test data sets from RKI on country
    and federal state level. Information on federal state level do not sum
    up to country-wide information since less laboratories are participating.


    :returns: dataframe array with country level information first and
        federal state level second

    """
    df_test = [[], []]

    # get country-wide testing data without resolution per federal state
    # but from much more laboratories
    url = "https://github.com/robert-koch-institut/SARS-CoV-2-PCR-Testungen_in_Deutschland/raw/main/SARS-CoV-2-PCR-Testungen_in_Deutschland.csv"
    df_test[0] = gd.get_file(url=url, read_data=False)
    # start on calender week 12/2020 as in federal states sheet, below
    df_test[0] = df_test[0].iloc[2:, :].reset_index(drop=True)

    # get testing data on federal state level (from only a subset of
    # laboratories)
    url = 'https://ars.rki.de/Docs/SARS_CoV2/Daten/data_wochenbericht.xlsx'
    header = {'User-Agent': 'Mozilla/5.0'}
    r = requests.get(url, headers=header)
    if r.status_code != 200:  # e.g. 404
        raise requests.exceptions.HTTPError("HTTPError: "+str(r.status_code))
    with io.BytesIO(r.content) as fh:
        df = pd.io.excel.ExcelFile(fh, engine=gd.Conf.excel_engine)
        sheet_names = df.sheet_names
        df_test[1] = pd.read_excel(df, sheet_name=sheet_names[3], header=[4],
                                   dtype={'Anteil positiv': float})

    return df_test

# transform calender weeks of data frames to dates using Thursday
# as the representation of each week


def transform_weeks_to_dates(df_test):
    """ Transforms the calender weeks of the two data frames obtained from
        RKI sources to dates in the middle of the corresponding week
        (i.e., Thursdays).

    :param df_test: 
    :returns: test data data frames with calender weeks replaced by dates.

    """
    # country-wide data
    df_test[0].rename(
        columns={df_test[0].columns[0]: dd.EngEng['date']}, inplace=True)
    for i in range(len(df_test[0])):
        # use %G insteaf of %Y (for year) and %V instead of %W (for month)
        # to get ISO week definition
        df_test[0].loc[i, dd.EngEng['date']] = datetime.strftime(datetime.strptime(
            str(df_test[0].loc[i, dd.EngEng['date']]).replace('W', '') + '-4', "%G-%V-%u"), "%Y-%m-%d")

    # federal state-based data
    df_test[1].rename(columns={df_test[1].columns[1]: dd.EngEng['date']}, inplace=True)
    for i in range(len(df_test[1])):
        datestr = str(df_test[1].loc[i, df_test[1].columns[2]]) + \
            "/" + str(df_test[1].loc[i, dd.EngEng['date']]) + '-4'
        # use %G insteaf of %Y (for year) and %V instead of %W (for month)
        # to get ISO week definition
        df_test[1].loc[i, dd.EngEng['date']] = datetime.strftime(
            datetime.strptime(datestr, "%V/%G-%w"), "%Y-%m-%d")
    # drop specific column on week after merge in year/date column
    df_test[1].drop(columns=df_test[1].columns[2], inplace=True)

    return df_test

# gets rki testing monitoring data resolved by federal states (which only
# is a subset of the total conducted tests)
# extrapolates the values for counties according to their population


def get_testing_data(read_data=dd.defaultDict['read_data'],
                     file_format=dd.defaultDict['file_format'],
                     out_folder=dd.defaultDict['out_folder'],
                     start_date=dd.defaultDict['start_date'],
                     end_date=dd.defaultDict['end_date'],
                     impute_dates=dd.defaultDict['impute_dates'],
                     moving_average=dd.defaultDict['moving_average'],
                     **kwargs):
    """ Downloads the RKI testing data and provides positive rates of
    testing data in different ways. Since positive rates also implicitly
    provide information on testing numbers while the opposite is
    not necessarily true without having additional information,
    only positive rates are provided.

    The data is read from the internet.
    The file is read in or stored at the folder "out_folder"/Germany/pydata.
    To store and change the data we use pandas.

    While working with the data
    - the column names are changed to English depending on defaultDict
    - The column "Date" provides information on the date of each data
        point given in the corresponding columns.

    - The data is exported in three different ways:
        - germany_testpos: Positive rates of testing for whole Germany
        - germany_states_testpos: Positive rates of testing for all
            federal states of Germany
        - germany_counties_from_states_testpos: Positive rates of testing
            for all counties of Germany, only taken from the
            values of the federal states. No extrapolation applied.

    - Missing dates are imputed for all data frames ('impute_dates' is
        not optional but always executed).
    - A central moving average of N days is optional.

    - Start and end dates can be provided to define the length of the
        returned data frames.

    :param read_data: True or False. Defines if data is read from file or downloaded. (Default value = dd.defaultDict['read_data'])
    :param file_format: File format which is used for writing the data. Default defined in defaultDict.
    :param out_folder: Folder where data is written to. Default defined in defaultDict.
    :param start_date: Date of first date in dataframe. Default defined in defaultDict.
    :param end_date: Date of last date in dataframe. Default defined in defaultDict.
    :param impute_dates: True or False. Defines if values for dates without new information are imputed. Default defined in defaultDict.
        At the moment they are always imputed.
    :param moving_average: Integers >=0. Applies an 'moving_average'-days moving average on all time series
        to smooth out effects of irregular reporting. Default defined in defaultDict.
    :param **kwargs: 

    """
    conf = gd.Conf(out_folder, **kwargs)
    out_folder = conf.path_to_use
    no_raw = conf.no_raw

    # data for all dates is automatically added
    impute_dates = True

    directory = out_folder
    directory = os.path.join(directory, 'Germany', 'pydata')
    gd.check_dir(directory)

    filename_county = "RKITestFull_Country"
    filename_state = "RKITestFull_FederalStates"

    if read_data:

        df_test = [[], []]

        county_file_in = os.path.join(directory, filename_county + ".json")
        try:
            df_test[0] = pd.read_json(county_file_in)
        # pandas>1.5 raise FileNotFoundError instead of ValueError
        except (ValueError, FileNotFoundError):
            raise FileNotFoundError("Error: The file: " + county_file_in +
                                    " does not exist. Call program without"
                                    " -r flag to get it.")

        state_file_in = os.path.join(directory, filename_state + ".json")
        try:
            df_test[1] = pd.read_json(state_file_in)
        # pandas>1.5 raise FileNotFoundError instead of ValueError
        except (ValueError, FileNotFoundError):
            raise FileNotFoundError("Error: The file: " + state_file_in +
                                    " does not exist. Call program without"
                                    " -r flag to get it.")
    else:
        df_test = download_testing_data()

    if not no_raw:
        gd.write_dataframe(
            df_test[0],
            directory, filename_county, "json")
        gd.write_dataframe(
            df_test[1],
            directory, filename_state, "json")

    # transform calender week to date by using Thursday ('-4')
    # of the corresponding calender week
    df_test = transform_weeks_to_dates(df_test)

    # rename columns
    df_test[0].rename(dd.GerEng, axis=1, inplace=True)
    df_test[1].rename(dd.GerEng, axis=1, inplace=True)

    try:
        df_test[0][dd.EngEng['date']] = pd.to_datetime(
            df_test[0][dd.EngEng['date']], format="ISO8601")
        df_test[1][dd.EngEng['date']] = pd.to_datetime(
            df_test[1][dd.EngEng['date']], format="ISO8601")
    except ValueError:
        try:
            df_test[0][dd.EngEng['date']] = pd.to_datetime(
                df_test[0][dd.EngEng['date']])
            df_test[1][dd.EngEng['date']] = pd.to_datetime(
                df_test[1][dd.EngEng['date']])
        except:
            raise gd.DataError(
                "Time data can't be transformed to intended format")

    # drop columns
    df_test[0].drop(
        columns=['tests_total', 'tests_total_accumulated', 'tests_positive',
                 'tests_positive_accumulated',
                 'laboratories_tests', 'capacities_daily',
                 'capacities_weekly_theoretically', 'capacities_weeklyweek_actually',
                 'laboratories_capacities', 'laboratories_samplebacklog',
                 'samplebacklog'], inplace=True)
    df_test[1].drop(columns='Anzahl Gesamt', inplace=True)

    # remove unknown locations
    df_test[1] = df_test[1][df_test[1].State != 'unbekannt']
    df_test[1].reset_index(drop=True, inplace=True)

    # replace state names with IDs
    df_test[1].rename(
        columns={dd.EngEng['state']: dd.EngEng['idState']}, inplace=True)
    for stateName, stateID in geoger.get_state_names_and_ids():
        df_test[1].loc[df_test[1][dd.EngEng['idState']]
                       == stateName, dd.EngEng['idState']] = stateID

    # set last values for missing dates via forward imputation
    df_test[0] = mdfs.impute_and_reduce_df(
        df_test[0],
        {},
        [dd.EngEng['testPositiveRatio']],
        impute='forward', moving_average=moving_average,
        min_date=start_date, max_date=end_date)

    # store positive rates for the whole country
    filename = 'germany_testpos'
    filename = gd.append_filename(filename, impute_dates, moving_average)
    gd.write_dataframe(df_test[0], directory, filename, file_format)

    # plot country-wide positive rates
    if conf.plot:
        # make plot
        customPlot.plot_multiple_series(
            df_test[0][dd.EngEng['date']],
            [df_test[0][dd.EngEng['testPositiveRatio']]],
            ["Germany"],
            title='Positive rate for Sars-CoV-2 testing', xlabel='Date', ylabel='Positive rate',
            fig_name="Germany_Testing_positive_rate")

    # set last values for missing dates via forward imputation
    df_test[1] = mdfs.impute_and_reduce_df(
        df_test[1],
        {dd.EngEng["idState"]: [k for k in geoger.get_state_ids()]},
        [dd.EngEng['testPositiveRatio']],
        impute='forward', moving_average=moving_average,
        min_date=start_date, max_date=end_date)
    # store positive rates for the all federal states
    filename = 'germany_states_testpos'
    filename = gd.append_filename(filename, impute_dates, moving_average)
    gd.write_dataframe(df_test[1], directory, filename, file_format)

    # plot positive rates of federal states
    if conf.plot:
        # make plot
        customPlot.plot_multiple_series(
            df_test[0][dd.EngEng['date']],
            [df_test[1].loc
             [df_test[1][dd.EngEng['idState']] == stateID,
              [dd.EngEng['testPositiveRatio']]] for stateID in geoger.get_state_ids()],
            [stateName for stateName in geoger.get_state_names()],
            title='Positive rate for Sars-CoV-2 testing', xlabel='Date', ylabel='Positive rate',
            fig_name='FederalStates_Testing_positive_rate')

    # store positive rates of federal states on county level
    # get county ids
    unique_geo_entities = geoger.get_county_ids()

    # create list of DataFrames, later to be merged
    df_test_counties = []
    states_str = dict(zip((geoger.get_state_ids(zfill=True)), range(1,
                                                                    1+len(geoger.get_state_ids()))))
    for county in unique_geo_entities:
        county_str = str(county).zfill(5)
        state_index = states_str[county_str[0:2]]
        df_local = pd.DataFrame(
            df_test[1].loc[df_test[1][dd.EngEng['idState']]
                           == state_index])
        df_local.rename(
            columns=({dd.EngEng['idState']: dd.EngEng['idCounty']}),
            inplace=True)
        df_local[dd.EngEng['idCounty']] = county
        df_test_counties.append(df_local)

    df_test_counties = pd.concat(df_test_counties)

    # store positive rates for the all federal states
    filename = 'germany_counties_from_states_testpos'
    filename = gd.append_filename(filename, impute_dates, moving_average)
    gd.write_dataframe(df_test_counties, directory, filename, file_format)


def main():
    """ Main program entry."""

    arg_dict = gd.cli("testing")
    get_testing_data(**arg_dict)


if __name__ == "__main__":

    main()
