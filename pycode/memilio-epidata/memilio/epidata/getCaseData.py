#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Kathrin Rack, Wadim Koslow, Martin J. Kuehn, Annette Lutz
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
@file getCaseData.py
@brief Downloads the case data of the Robert Koch-Institute (RKI) and provides it in different ways.

The raw case data we download can be found at
https://github.com/robert-koch-institut/SARS-CoV-2-Infektionen_in_Deutschland

Be careful: Date of deaths or recovery is not reported in original case data and will be
extrapolated in this script.
"""

# Imports
import os
from datetime import date

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import modifyDataframeSeries as mdfs
from memilio.epidata import progress_indicator


def check_for_completeness(df, merge_berlin=False, merge_eisenach=True):
    """! Checks if all counties are mentioned in the case data set

   This check had to be added due to incomplete data downloads
   It is checked if all all counties are part of the data.
   If data is incomplete the data is downloaded from another source.
   Note: There is no check if data for every day and every county is available (which can happen).

   @param df pandas dataframe to check
   @return Boolean to say if data is complete or not
   """

    if not df.empty:
        return geoger.check_for_all_counties(
            df["IdLandkreis"].unique(),
            merge_berlin, merge_eisenach)
    # if it is empty
    return False


def get_case_data(read_data=dd.defaultDict['read_data'],
                  file_format=dd.defaultDict['file_format'],
                  out_folder=dd.defaultDict['out_folder'],
                  no_raw=dd.defaultDict['no_raw'],
                  start_date=date(2020, 1, 1),
                  end_date=dd.defaultDict['end_date'],
                  impute_dates=dd.defaultDict['impute_dates'],
                  moving_average=dd.defaultDict['moving_average'],
                  make_plot=dd.defaultDict['make_plot'],
                  split_berlin=dd.defaultDict['split_berlin'],
                  rep_date=dd.defaultDict['rep_date'],
                  files='All'
                  ):
    """! Downloads the case data and provides different kind of structured data

    The data is read either from the internet or from a json file (CaseDataFull.json), stored in an earlier run.
    If the data is read from the internet, before changing anything the data is stored in CaseDataFull.json.
    If data should be downloaded, it is checked if data contains all counties.
    If not a different source is tried.
    The file is read in or stored at the folder "out_folder"/Germany/.
    To store and change the data we use pandas.

    While working with the data
    - the column names are changed to english depending on defaultDict
    - a new column "Date" is defined.
    - we are only interested in the values where the parameter NeuerFall, NeuerTodesfall, NeuGenesen are larger than 0.
    The values, when these parameters are negative are just useful,
    if one would want to get the difference to the previous day.
    For details we refer to the above mentioned webpage.
    - For all different parameters and different columns the values are added up for whole germany for every date
    and the cumulative sum is calculated. Unless something else is mentioned.
    - For Berlin all districts can be merged into one [Default]. Otherwise, Berlin is divided into multiple districts and
        different file names are used.
    - Following data is generated and written to the mentioned filename
        - All infected (current and past) for whole germany are stored in "cases_infected"
        - All deaths whole germany are stored in "cases_deaths"
        - Infected, deaths and recovered for whole germany are stored in "cases_all_germany"
        - Infected split for states are stored in "cases_infected_state"
        - Infected, deaths and recovered split for states are stored in "cases_all_state"
        - Infected split for counties are stored in "cases_infected_county(_split_berlin)"
        - Infected, deaths and recovered split for county are stored in "cases_all_county(_split_berlin)"
        - Infected, deaths and recovered split for gender are stored in "cases_all_gender"
        - Infected, deaths and recovered split for state and gender are stored in "cases_all_state_gender"
        - Infected, deaths and recovered split for county and gender are stored in "cases_all_county_gender(_split_berlin)"
        - Infected, deaths and recovered split for age are stored in "cases_all_age"
        - Infected, deaths and recovered split for state and age are stored in "cases_all_state_age"
        - Infected, deaths and recovered split for county and age are stored in "cases_all_county_age(_split_berlin)"

    @param read_data True or False. Defines if data is read from file or downloaded. Default defined in defaultDict.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    @param out_folder Folder where data is written to. Default defined in defaultDict.
    @param no_raw True or False. Defines if unchanged raw data is saved or not. Default defined in defaultDict.
    @param start_date Date of first date in dataframe. Default 2020-01-01.
    @param end_date Date of last date in dataframe. Default defined in defaultDict.
    @param impute_dates True or False. Defines if values for dates without new information are imputed. Default defined in defaultDict.
    @param moving_average Integers >=0. Applies an 'moving_average'-days moving average on all time series
        to smooth out effects of irregular reporting. Default defined in defaultDict.
    @param make_plot True or False. Defines if plots are generated with matplotlib. Default defined in defaultDict.
    @param split_berlin True or False. Defines if Berlin's disctricts are kept separated or get merged. Default defined in defaultDict.
    @param rep_date True or False. Defines if reporting date or reference date is taken into dataframe. Default defined in defaultDict.
    @param files List of strings or 'All' or 'Plot'. Defnies which files should be provided (and plotted). Default 'All'.
    """

    if files == 'All':
        files = ['infected', 'deaths', 'all_germany', 'infected_state',
                 'all_state', 'infected_county', 'all_county', 'all_gender',
                 'all_state_gender', 'all_county_gender', 'all_age',
                 'all_state_age', 'all_county_age']
    if files == 'Plot':
        # only consider plotable files
        files = ['infected', 'deaths', 'all_gender', 'all_age']

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)
    filename = "CaseDataFull"

    complete = False
    path = os.path.join(directory + filename + ".json")
    try:
        url = "https://media.githubusercontent.com/media/robert-koch-institut/" + \
            "SARS-CoV-2-Infektionen_in_Deutschland/main/Aktuell_Deutschland_SarsCov2_Infektionen.csv"
        df = gd.get_file(path, url, read_data, param_dict={}, interactive=True)
        complete = check_for_completeness(df, merge_eisenach=True)
    except:
        pass
    if complete:
        if not read_data:
            # add column with state ids
            county_to_state_map = geoger.get_countyid_to_stateid_map(
                merge_berlin=False)
            df["IdBundesland"] = df["IdLandkreis"].map(county_to_state_map)
    else:
        # try another possibility if df was empty or incomplete
        print("Note: Case data is incomplete. Trying another source.")
        try:
            url = "https://opendata.arcgis.com/datasets/66876b81065340a4a48710b062319336_0.csv"
            # if this file is encoded with utf-8 German umlauts are not displayed correctly because they take two bytes
            # utf_8_sig can identify those bytes as one sign and display it correctly
            df = gd.get_file(path, url, False, param_dict={
                             "encoding": 'utf_8_sig'}, interactive=True)
            complete = check_for_completeness(df, merge_eisenach=True)
        except:
            pass
        if not complete:
            print("Note: Case data is still incomplete. Trying a thrid source.")
            try:
                # If the data on github is not available we download the case data from rki from covid-19 datahub
                url = "https://npgeo-de.maps.arcgis.com/sharing/rest/content/" +\
                    "items/f10774f1c63e40168479a1feb6c7ca74/data"
                df = gd.get_file(path, url, False, param_dict={
                                 "encoding": 'utf_8_sig'}, interactive=True)
                df.rename(columns={'FID': "OBJECTID"}, inplace=True)
                complete = check_for_completeness(df, merge_eisenach=True)
            except:
                pass
        if not complete:
            raise FileNotFoundError(
                "Something went wrong, dataframe is empty for csv and geojson!")

        # drop columns that do not exist in data from github
        df = df.drop(["Altersgruppe2", "Datenstand", "OBJECTID",
                      "Bundesland", "Landkreis"], axis=1)
    with progress_indicator.Spinner(message='Preparing DataFrame'):
        df = df.convert_dtypes()

        # output data to not always download it
        if not no_raw:
            gd.write_dataframe(df, directory, filename, "json")

        # store dict values in parameter to not always call dict itself
        Altersgruppe = dd.GerEng['Altersgruppe']
        Geschlecht = dd.GerEng['Geschlecht']
        AnzahlFall = dd.GerEng['AnzahlFall']
        AnzahlGenesen = dd.GerEng['AnzahlGenesen']
        AnzahlTodesfall = dd.GerEng['AnzahlTodesfall']
        IdBundesland = dd.GerEng['IdBundesland']
        IdLandkreis = dd.GerEng['IdLandkreis']

        # translate column gender from German to English and standardize
        df.loc[df.Geschlecht == 'unbekannt', [
            'Geschlecht']] = dd.GerEng['unbekannt']
        df.loc[df.Geschlecht == 'W', ['Geschlecht']] = dd.GerEng['W']
        df.loc[df.Geschlecht == 'M', ['Geschlecht']] = dd.GerEng['M']
        df.loc[df.Altersgruppe == 'unbekannt', [
            'Altersgruppe']] = dd.GerEng['unbekannt']

        # change names of columns
        df.rename(dd.GerEng, axis=1, inplace=True)

        # Add column 'Date' with Date
        # = reporting date if rep_date is set
        # = reference date (date of disease onset) if IstErkrankungsbeginn = 1 else
        #       take Meldedatum (reporting date)
        if rep_date:
            df[dd.EngEng['date']] = df['Meldedatum'].astype('object')
        else:
            df[dd.EngEng['date']] = np.where(
                df['IstErkrankungsbeginn'] == 1, df['Refdatum'],
                df['Meldedatum'])

        try:
            df[dd.EngEng['date']] = pd.to_datetime(
                df[dd.EngEng['date']], format="ISO8601")
        except ValueError:
            try:
                df[dd.EngEng['date']] = pd.to_datetime(
                    df[dd.EngEng['date']], format="%Y-%m-%d")
            except:
                raise gd.DataError(
                    "Time data can't be transformed to intended format")

        # Date is either Refdatum or Meldedatum after column
        # 'IstErkrankungsbeginn' has been added. See also rep_date option.
        dateToUse = dd.EngEng['date']
        df.sort_values(dateToUse, inplace=True)

        # Manipulate data to get rid of conditions: df.NeuerFall >= 0, df.NeuerTodesfall >= 0, df.NeuGenesen >=0
        # There might be a better way
        df.loc[df.NeuerFall < 0, [AnzahlFall]] = 0
        df.loc[df.NeuerTodesfall < 0, [AnzahlTodesfall]] = 0
        df.loc[df.NeuGenesen < 0, [AnzahlGenesen]] = 0

        # get rid of unnecessary columns
        df = df.drop(['NeuerFall', 'NeuerTodesfall', 'NeuGenesen',
                      "IstErkrankungsbeginn", "Meldedatum", "Refdatum"], axis=1)

        # merge Berlin counties
        if not split_berlin:
            df = geoger.merge_df_counties(
                df, 11000, geoger.CountyMerging[11000],
                sorting=[dd.EngEng['date']],
                columns=[dd.EngEng['date'],
                         dd.EngEng['gender'],
                         dd.EngEng['idState'],
                         dd.EngEng['ageRKI']])

    # dict for all files
    # filename -> [groupby_list, .agg({}), groupby_index, groupby_cols, mod_cols]
    dict_files = {
        'infected': [dateToUse, {AnzahlFall: sum}, None, {}, ['Confirmed']],
        'deaths': [dateToUse, {AnzahlTodesfall: sum}, None, {}, ['Deaths']],
        'all_germany': [dateToUse, {AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum},
                        None, {}, ['Confirmed', 'Deaths', 'Recovered']],
        'infected_state': [[dateToUse, IdBundesland], {AnzahlFall: sum}, [IdBundesland],
                           {dd.EngEng["idState"]: [k for k, v in dd.State.items()]}, ['Confirmed']],
        'all_state': [[dateToUse, IdBundesland], {AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum},
                      [IdBundesland], {dd.EngEng["idState"]
                          : [k for k, v in dd.State.items()]},
                      ['Confirmed', 'Deaths', 'Recovered']],
        'infected_county': [[dateToUse, IdLandkreis], {AnzahlFall: sum}, [IdLandkreis],
                            {dd.EngEng["idCounty"]: sorted(set(df[dd.EngEng["idCounty"]].unique()))}, ['Confirmed']],
        'all_county': [[dateToUse, IdLandkreis], {AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum},
                       [IdLandkreis], {dd.EngEng["idCounty"]: sorted(
                           set(df[dd.EngEng["idCounty"]].unique()))},
                       ['Confirmed', 'Deaths', 'Recovered']],
        'all_gender': [[dateToUse, Geschlecht], {AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum},
                       [Geschlecht], {dd.EngEng["gender"]: list(
                           df[dd.EngEng["gender"]].unique())},
                       ['Confirmed', 'Deaths', 'Recovered']],
        'all_state_gender': [[dateToUse, IdBundesland, Geschlecht],
                             {AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}, [
                                 IdBundesland, Geschlecht],
                             {dd.EngEng["idState"]: geoger.get_state_ids(), dd.EngEng["gender"]: list(
                                 df[dd.EngEng["gender"]].unique())},
                             ['Confirmed', 'Deaths', 'Recovered']],
        'all_county_gender': [[dateToUse, IdLandkreis, Geschlecht],
                              {AnzahlFall: sum, AnzahlTodesfall: sum,
                                  AnzahlGenesen: sum}, [IdLandkreis, Geschlecht],
                              {dd.EngEng["idCounty"]: sorted(set(df[dd.EngEng["idCounty"]].unique(
                              ))), dd.EngEng["gender"]: list(df[dd.EngEng["gender"]].unique())},
                              ['Confirmed', 'Deaths', 'Recovered']],
        'all_age': [[dateToUse, Altersgruppe], {AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum},
                    [Altersgruppe], {dd.EngEng["ageRKI"]: sorted(
                        set(df[dd.EngEng["ageRKI"]].unique()))},
                    ['Confirmed', 'Deaths', 'Recovered']],
        'all_state_age': [[dateToUse, IdBundesland, Altersgruppe],
                          {AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}, [
                              IdBundesland, Altersgruppe],
                          {dd.EngEng["idState"]: geoger.get_state_ids(), dd.EngEng["ageRKI"]: sorted(
                              set(df[dd.EngEng["ageRKI"]].unique()))},
                          ['Confirmed', 'Deaths', 'Recovered']],
        'all_county_age': [[dateToUse, IdLandkreis, Altersgruppe],
                           {AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}, [
                               IdLandkreis, Altersgruppe],
                           {dd.EngEng["idCounty"]: sorted(set(df[dd.EngEng["idCounty"]].unique(
                           ))), dd.EngEng["ageRKI"]: sorted(set(df[dd.EngEng["ageRKI"]].unique()))},
                           ['Confirmed', 'Deaths', 'Recovered']]
    }
    with progress_indicator.Spinner():
        for file in files:
            if file not in dict_files.keys():
                raise gd.DataError('Error: File '+file+' cannot be written.')
            # split berlin is only relevant for county level
            if 'county' in file and split_berlin == True:
                split_berlin_local = True
            else:
                # dont append _split_berlin to filename on germany/state level
                split_berlin_local = False
            filename = 'cases_' + \
                gd.append_filename(file, impute_dates,
                                   moving_average, split_berlin_local, rep_date)
            # sum over all columns defined in dict_files
            df_local = df.groupby(dict_files[file][0]).agg(dict_files[file][1])

            if file == 'deaths':
                # only consider where deaths > 0
                df_local = df_local[df_local[AnzahlTodesfall] != 0]

            # cumulative sum over columns defined in dict_files
            if dict_files[file][2] == None:
                df_local_cs = df_local.cumsum().reset_index(drop=False)
            else:
                df_local_cs = df_local.groupby(level=[dict_files[file][0].index(
                    level_index) for level_index in dict_files[file][2]]).cumsum().reset_index()

            if impute_dates or moving_average > 0:
                df_local_cs = mdfs.impute_and_reduce_df(
                    df_local_cs, group_by_cols=dict_files[file][3],
                    mod_cols=dict_files[file][4],
                    impute='forward', moving_average=moving_average,
                    min_date=start_date, max_date=end_date)

            df_local_cs = mdfs.extract_subframe_based_on_dates(
                df_local_cs, start_date, end_date)
            gd.write_dataframe(df_local_cs, directory, filename, file_format)

            if make_plot:
                if file == 'infected':
                    # make plot
                    df_local_cs.plot(title='COVID-19 infections', grid=True,
                                     style='-o')
                    plt.tight_layout()
                    plt.show()

                if file == 'deaths':
                    df_local_cs.plot(title='COVID-19 deaths', grid=True,
                                     style='-o')
                    plt.tight_layout()
                    plt.show()

                    df.agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}) \
                        .plot(title='COVID-19 infections, deaths, recovered', grid=True,
                              kind='bar')
                    plt.tight_layout()
                    plt.show()

                if file == 'all_gender':
                    df.groupby(Geschlecht).agg(
                        {AnzahlFall: sum, AnzahlTodesfall: sum,
                         AnzahlGenesen: sum}).plot(
                        title='COVID-19 infections, deaths, recovered',
                        grid=True, kind='bar')
                    plt.tight_layout()
                    plt.show()

                if file == 'all_age':
                    df.groupby(Altersgruppe).agg(
                        {AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}).plot(
                        title='COVID-19 infections, deaths, recovered for diff ages',
                        grid=True, kind='bar')
                    plt.tight_layout()
                    plt.show()

                    # Dead by "Altersgruppe":
                    df_local = df.groupby(Altersgruppe).agg(
                        {AnzahlTodesfall: sum})

                    df_local.plot(title='COVID-19 deaths', grid=True,
                                  kind='bar')
                    plt.tight_layout()
                    plt.show()


def main():
    """! Main program entry."""

    arg_dict = gd.cli("cases")
    get_case_data(**arg_dict)


if __name__ == "__main__":
    main()
