#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
from datetime import datetime
import os
import requests
import io
import pandas as pd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from memilio.epidata import modifyDataframeSeries
from memilio.epidata import customPlot
from memilio.epidata import geoModificationGermany as geoger

# Downloads testing data from RKI

def download_testing_data():
    """! Downloads the Sars-CoV-2 test data sets from RKI on country 
    and federal state level. Information on federal state level do not sum
    up to country-wide information since less laboratories are participating.

    @return dataframe array with country level information first and 
        federal state level second
    """
    df_test = [[], []]

    # get country-wide testing data without resolution per federal state
    # but from much more laboratories
    url = 'https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Daten/Testzahlen-gesamt.xlsx?__blob=publicationFile'
    header = {'User-Agent': 'Mozilla/5.0'}
    r = requests.get(url, headers=header)
    with io.BytesIO(r.content) as fh:
        df = pd.io.excel.ExcelFile(fh, engine='openpyxl')
        sheet_names = df.sheet_names
        df_test[0] = pd.read_excel(
            df, sheet_name=sheet_names[1],
            dtype={'Positivenanteil (%)': float})
        # start on calender week 12/2020 as in federal states sheet,
        # below and remove sum at bottom
        df_test[0] = df_test[0][2:-1].reset_index()
        df_test[0] = df_test[0].drop(columns='index')

    # get testing data on federal state level (from only a subset of
    # laboratories)
    url = 'https://ars.rki.de/Docs/SARS_CoV2/Daten/data_wochenbericht.xlsx'
    header = {'User-Agent': 'Mozilla/5.0'}
    r = requests.get(url, headers=header)
    with io.BytesIO(r.content) as fh:
        df = pd.io.excel.ExcelFile(fh, engine='openpyxl')
        sheet_names = df.sheet_names
        df_test[1] = pd.read_excel(df, sheet_name=sheet_names[3], header=[4],
            dtype={'Anteil positiv': float})

    return df_test

# transform calender weeks of data frames to dates using Thursday
# as the representation of each week
def transform_weeks_to_dates(df_test):
    """! Transforms the calender weeks of the two data frames obtained from 
        RKI sources to dates in the middle of the corresponding week
        (i.e., Thursdays).

    @return test data data frames with calender weeks replaced by dates.
    """
    # country-wide data
    df_test[0] = df_test[0].rename(
        columns={df_test[0].columns[0]: dd.EngEng['date']})
    for i in range(len(df_test[0])):
        # use %G insteaf of %Y (for year) and %V instead of %W (for month)
        # to get ISO week definition
        df_test[0].loc[i, dd.EngEng['date']] = datetime.strftime(datetime.strptime(
            df_test[0].loc[i, dd.EngEng['date']] + '-4', "%V/%G-%w"), "%Y-%m-%d")

    # federal state-based data
    df_test[1] = df_test[1].rename(
        columns={df_test[1].columns[1]: dd.EngEng['date']})
    for i in range(len(df_test[1])):
        datestr = str(df_test[1].loc[i, df_test[1].columns[2]]) + \
            "/" + str(df_test[1].loc[i, dd.EngEng['date']]) + '-4'
        # use %G insteaf of %Y (for year) and %V instead of %W (for month)
        # to get ISO week definition
        df_test[1].loc[i, dd.EngEng['date']] = datetime.strftime(
            datetime.strptime(datestr, "%V/%G-%w"), "%Y-%m-%d")
    # drop specific column on week after merge in year/date column
    df_test[1] = df_test[1].drop(columns=df_test[1].columns[2])

    return df_test

# gets rki testing monitoring data resolved by federal states (which only
# is a subset of the total conducted tests)
# extrapolates the values for counties according to their population
def get_testing_data(data_folder,
                     read_data=dd.defaultDict['read_data'],
                     file_format=dd.defaultDict['file_format'],
                     no_raw=dd.defaultDict['no_raw'],
                     start_date=dd.defaultDict['start_date'],
                     end_date=dd.defaultDict['end_date'],
                     make_plot=dd.defaultDict['make_plot'],
                     moving_average=dd.defaultDict['moving_average']):
    """! Downloads the RKI testing data and provides positive rates of 
    testing data in different ways. Since positive rates also implicitly 
    provide information on testing numbers while the opposite is
    not necessarily true without having additional information, 
    only positive rates are provided.

    The data is read from the internet.
    The file is read in or stored at the folder "data_folder"/Germany/.
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

    @param data_folder Path to folder where data is written in folder 
        data_folder/Germany.
    @param read_data False [Default]. Data is always downloaded from 
        the internet.
    @param file_format File format which is used for writing the data. 
        Default defined in defaultDict.
    @param no_raw True or False [Default]. Defines if raw data is 
        saved or not.
    @param start_date [Default = '', taken from read data] Start date
        of stored data frames.
    @param end_date [Default = '', taken from read data] End date of
        stored data frames.
    @param make_plot False [Default] or True. Defines if plots are
        generated with matplotlib.
    @param moving_average 0 [Default] or Number>0. Defines the number of
        days for which a centered moving average is computed.
    """
    # data for all dates is automatically added
    impute_dates = True

    directory = os.path.join(data_folder, 'Germany/')
    gd.check_dir(directory)

    df_test = download_testing_data()

    if not no_raw:
        gd.write_dataframe(
            df_test[0],
            directory, "RKITestFull_Country", "json")
        gd.write_dataframe(
            df_test[1],
            directory, "RKITestFull_FederalStates", "json")

    # transform calender week to date by using Thursday ('-4')
    # of the corresponding calender week
    df_test = transform_weeks_to_dates(df_test)

    # rename columns
    df_test[0].rename(dd.GerEng, axis=1, inplace=True)
    df_test[1].rename(dd.GerEng, axis=1, inplace=True)

    # drop columns
    df_test[0] = df_test[0].drop(
        columns=['Anzahl Testungen', 'Positiv getestet',
                 'Anzahl übermittelnder Labore'])
    df_test[1] = df_test[1].drop(columns='Anzahl Gesamt')

    # remove unknown locations
    df_test[1] = df_test[1][df_test[1].State != 'unbekannt']
    df_test[1] = df_test[1].reset_index()
    df_test[1] = df_test[1].drop(columns='index')

    # correct positive rate to percentage
    df_test[0][dd.EngEng['positiveRate']
               ] = df_test[0][dd.EngEng['positiveRate']]/100

    # replace state names with IDs
    df_test[1] = df_test[1].rename(
        columns={dd.EngEng['state']: dd.EngEng['idState']})
    for stateName, stateID in geoger.get_state_names_and_ids():
        df_test[1].loc[df_test[1][dd.EngEng['idState']]
                       == stateName, dd.EngEng['idState']] = stateID

    # set last values for missing dates via forward imputation
    df_test[0] = modifyDataframeSeries.impute_and_reduce_df(
        df_test[0],
        {},
        [dd.EngEng['positiveRate']],
        impute='forward', moving_average=moving_average)

    # store positive rates for the whole country
    filename = 'germany_testpos'
    filename = gd.append_filename(filename, impute_dates, moving_average)
    gd.write_dataframe(df_test[0], directory, filename, file_format)

    # plot country-wide positive rates
    if make_plot:
        # make plot
        customPlot.plotList(
            df_test[0][dd.EngEng['date']],
            [df_test[0][dd.EngEng['positiveRate']]],
            ["Germany"],
            'Positive rate for Sars-CoV-2 testing', 'Date', 'Positive rate',
            "Germany_Testing_positive_rate")

    # set last values for missing dates via forward imputation
    df_test[1] = modifyDataframeSeries.impute_and_reduce_df(
        df_test[1],
        {dd.EngEng["idState"]: [k for k in geoger.get_state_ids()]},
        [dd.EngEng['positiveRate']],
        impute='forward', moving_average=moving_average)

    # store positive rates for the all federal states
    filename = 'germany_states_testpos'
    filename = gd.append_filename(filename, impute_dates, moving_average)
    gd.write_dataframe(df_test[1], directory, filename, file_format)

    # plot positive rates of federal states
    if make_plot:
        # make plot
        customPlot.plotList(
            df_test[0][dd.EngEng['date']],
            [df_test[1].loc
             [df_test[1][dd.EngEng['idState']] == stateID,
              [dd.EngEng['positiveRate']]] for stateID in geoger.get_state_ids()],
            [stateName for stateName in geoger.get_state_names()],
            'Positive rate for Sars-CoV-2 testing', 'Date', 'Positive rate',
            "FederalStates_Testing_positive_rate")

    # store positive rates of federal states on county level
    # get county ids
    unique_geo_entities = geoger.get_county_ids()

    df_test_counties = pd.DataFrame()
    states_str = dict(zip((geoger.get_state_ids(zfill=True)), range(1,
                                                      1+len(geoger.get_state_ids()))))
    for county in unique_geo_entities:
        county_str = str(county).zfill(5)
        state_index = states_str[county_str[0:2]]
        df_local = pd.DataFrame(
            df_test[1].loc[df_test[1][dd.EngEng['idState']]
                           == state_index].copy())
        df_local.rename(
            columns=({dd.EngEng['idState']: dd.EngEng['idCounty']}),
            inplace=True)
        df_local[dd.EngEng['idCounty']] = county
        df_test_counties = df_test_counties.append(df_local.copy())

     # store positive rates for the all federal states
    filename = 'germany_counties_from_states_testpos'
    filename = gd.append_filename(filename, impute_dates, moving_average)
    gd.write_dataframe(df_test_counties, directory, filename, file_format)

def main():
    """! Main program entry."""

    path = os.path.join(os.getcwd(), 'data', 'pydata')
    arg_dict = gd.cli("testing")
    get_testing_data(path, **arg_dict)

if __name__ == "__main__":

   main()
