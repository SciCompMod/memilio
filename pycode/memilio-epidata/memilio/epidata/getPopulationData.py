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
@file getPopulationData.py

@brief Downloads data about population statistic

"""
import warnings
import requests
import os
import twill
import time
import io

import numpy as np
import pandas as pd

from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd


def read_population_data(username, password, read_data, directory):
    '''! Reads Population data either from regionalstatistik.de or from directory

    A request is made using the twill package. Username and Password are required to
    sign in on regionalstatistik.de. After the sign twill navigates to the file to download.

    @param username Username to sign in at regionalstatistik.de. 
    @param password Password to sign in at regionalstatistik.de.
    @param read_data False or True. Defines if data is read from file or downloaded.
    @param directory Path to folder where data is read from.
    @return DataFrame
    '''

    filename = '12411-02-03-4'
    if not read_data:
        sign_in_url = 'https://www.regionalstatistik.de/genesis/online?Menu=Anmeldung'

        # sign in to regionalstatistik.de with given username and password
        twill.browser.user_agent = requests.utils.default_headers()[
            'User-Agent']
        twill.commands.go(sign_in_url)
        twill.commands.fv('3', 'KENNUNG', username)
        twill.commands.fv('3', 'PASSWORT', password)
        twill.commands.submit('login', '3')
        # navigate to file as in documentation
        twill.commands.follow('Themen')
        twill.commands.follow(filename[:2])
        # wait 1 second to prevent error
        time.sleep(1)
        twill.commands.follow(filename.split('-')[0])
        twill.commands.follow(filename)
        # start 'Werteabruf'
        twill.commands.submit('45', '3')
        # read csv file (1,4 for xlsx)
        twill.commands.submit('1', '5')

        df_pop_raw = pd.read_csv(io.StringIO(
            twill.browser.html), sep=';', header=6)

    else:
        data_file = os.path.join(directory, filename)
        if os.path.isfile(data_file+'.xlsx'):
            df_pop_raw = pd.read_excel(
                data_file+'.xlsx', engine='openpyxl', sheet_name=filename, header=4)
        elif os.path.isfile(data_file+'.csv'):
            df_pop_raw = pd.read_excel(data_file+'.csv', sep=';', header=6)
        else:
            raise FileNotFoundError(
                'Data file '+filename+' was not found in out_folder/Germany')

    return df_pop_raw


def export_population_dataframe(df_pop, directory, file_format, merge_eisenach):
    '''! Writes population dataframe into directory with new column names and age groups

    @param df_pop Population data DataFrame to be exported
    @param directory Directory where data is written to.
    @param file_format File format which is used for writing the data.
    @param merge_eisenach Defines whether the counties 'Wartburgkreis'
        and 'Eisenach' are listed separately or
        combined as one entity 'Wartburgkreis'.
    @return exported DataFrame
    '''

    new_cols = [
        dd.EngEng['idCounty'],
        dd.EngEng['population'],
        '<3 years', '3-5 years', '6-14 years', '15-17 years',
        '18-24 years', '25-29 years', '30-39 years', '40-49 years',
        '50-64 years', '65-74 years', '>74 years']
    df_pop_export = pd.DataFrame(columns=new_cols)
    df_pop_export[new_cols[0]] = df_pop[dd.EngEng['idCounty']]
    # <3 and 3-5
    df_pop_export[new_cols[2:4]] = df_pop[df_pop.columns[2:4]]
    # 6-14
    df_pop_export[new_cols[4]] = \
        df_pop[df_pop.columns[4:6]].sum(axis=1)
    # 15-17
    df_pop_export[new_cols[5]] = df_pop[df_pop.columns[6]]
    # 18-24
    df_pop_export[new_cols[6]] = \
        df_pop[df_pop.columns[7:9]].sum(axis=1)
    # 25-29
    df_pop_export[new_cols[7]] = df_pop[df_pop.columns[9]]
    # 30-39
    df_pop_export[new_cols[8]] = \
        df_pop[df_pop.columns[10:12]].sum(axis=1)
    # 40-49
    df_pop_export[new_cols[9]] = \
        df_pop[df_pop.columns[12:14]].sum(axis=1)
    # 50-64
    df_pop_export[new_cols[10]] = \
        df_pop[df_pop.columns[14:17]].sum(axis=1)
    # 65-74
    df_pop_export[new_cols[11]] = df_pop[df_pop.columns[17]]
    # >74
    df_pop_export[new_cols[12]] = df_pop[df_pop.columns[18]]

    df_pop_export[dd.EngEng['population']
                  ] = df_pop_export.iloc[:, 2:].sum(axis=1)

    # merge eisenach if no data available
    if '16056' in df_pop_export[dd.EngEng['idCounty']].values:
        if df_pop_export[df_pop_export[dd.EngEng['idCounty']] == '16056'][dd.EngEng['population']].values[0] == 0:
            df_pop_export = geoger.merge_df_counties_all(
                df_pop_export, sorting=[dd.EngEng["idCounty"]],
                columns=dd.EngEng["idCounty"])

    gd.check_dir(directory)

    if len(df_pop_export) == 401:
        filename = 'county_current_population_dim401'
        gd.write_dataframe(df_pop_export, directory, filename, file_format)

    if len(df_pop_export) == 400 or merge_eisenach:
        filename = 'county_current_population'

        # Merge Eisenach and Wartburgkreis
        df_pop_export = geoger.merge_df_counties_all(
            df_pop_export, sorting=[dd.EngEng["idCounty"]],
            columns=dd.EngEng["idCounty"])

        gd.write_dataframe(df_pop_export, directory, filename, file_format)

    return df_pop_export


def assign_population_data(df_pop_raw, counties, age_cols, idCounty_idx):
    '''! Assigns population data of all counties of old dataframe in new created dataframe

    In df_pop_raw there might be additional information like federal states, 
    governing regions etc. which is not necessary for the dataframe.
    Also checks for incomplete data.    

    @param df_pop_raw Raw Population DataFrame read from regionalstatistik.de
    @param counties List of counties to be assigned in new DataFrame
    @param age_cols Age groups in old DataFrame
    @param idCountyidx indexes in old DataFrame where data of corresponding county starts
    @return new DataFrame
    '''

    new_cols = {dd.EngEng['idCounty']: counties[:, 1],
                dd.EngEng['county']: counties[:, 0]}

    # number of age_cols
    num_age_groups = len(age_cols)

    # add age_cols with zero initialization to new_cols
    new_cols.update({age: 0 for age in age_cols})

    df_pop = pd.DataFrame(new_cols)

    empty_data = ['.', '-']
    for start_idx in idCounty_idx:

        county_id = df_pop_raw.loc[start_idx, dd.EngEng['idCounty']]

        # check for empty rows
        if df_pop_raw.loc[start_idx: start_idx + num_age_groups - 1,
                          dd.EngEng['number']].values.any() in empty_data:
            if not df_pop_raw.loc[start_idx: start_idx + num_age_groups - 1,
                                  dd.EngEng['number']].values.all() in empty_data:
                raise gd.DataError(
                    'Error. Partially incomplete data for county ' +
                    county_id)

        # county information needed
        elif county_id in counties[:, 1]:
            # direct assignment of population data found
            df_pop.loc[df_pop[dd.EngEng['idCounty']] == df_pop_raw.loc
                       [start_idx, dd.EngEng['idCounty']],
                       age_cols] = df_pop_raw.loc[start_idx: start_idx + num_age_groups - 1, dd.EngEng
                                                  ['number']].values.astype(int)
        # Berlin and Hamburg
        elif county_id + '000' in counties[:, 1]:
            # direct assignment of population data found
            df_pop.loc[df_pop[dd.EngEng['idCounty']] == df_pop_raw.loc
                       [start_idx, dd.EngEng['idCounty']] + '000',
                       age_cols] = df_pop_raw.loc[start_idx: start_idx + num_age_groups - 1, dd.EngEng
                                                  ['number']].values.astype(int)

        # additional information for local entities not needed
        elif county_id in ['03241001', '05334002', '10041100']:
            pass
        # Germany, federal states, and governing regions
        elif len(county_id) < 5:
            pass
        else:
            print('no data for ' + df_pop_raw.loc
                  [start_idx, dd.EngEng['idCounty']])
            raise gd.DataError(
                'Error. County ID in input population data '
                'found which could not be assigned.')

    return df_pop


def test_total_population(df_pop, age_cols):
    """! Tests if total population matches expectation
    @param df_pop Population Dataframe with all counties
    @param age_cols All age groups in DataFrame"""

    total_sum_2020 = 83155031
    total_sum_2021 = 83237124

    if df_pop[age_cols].sum().sum() != total_sum_2021:
        if df_pop[age_cols].sum().sum() == total_sum_2020:
            warnings.warn('Using data of 2020. Newer data is available.')
        else:
            raise gd.DataError('Total Population does not match expectation.')


def get_population_data(read_data=dd.defaultDict['read_data'],
                        file_format=dd.defaultDict['file_format'],
                        out_folder=dd.defaultDict['out_folder'],
                        no_raw=dd.defaultDict['no_raw'],
                        merge_eisenach=True,
                        username='',
                        password=''):
    """! Download age-stratified population data for the German counties.

    The data we use is:
    Official 'Bevölkerungsfortschreibung' 12411-02-03-4:
    'Bevölkerung nach Geschlecht und Altersgruppen (17)' 
    of regionalstatistik.de. 
    ATTENTION: The raw file cannot be downloaded 
    automatically by our scripts without an Genesis Online account. In order to
    work on this dataset, please enter your username and password or manually download it from:

    https://www.regionalstatistik.de/genesis/online -> "1: Gebiet, Bevölkerung,
    Arbeitsmarkt, Wahlen" -> "12: Bevölkerung" -> "12411 Fortschreibung des
    Bevölkerungsstandes" ->  "12411-02-03-4: Bevölkerung nach Geschlecht und 
    Altersgruppen (17) - Stichtag 31.12. - regionale Tiefe: Kreise und
    krfr. Städte". 

    Download the xlsx or csv file and put it under dd.defaultDict['out_folder'], 
    this normally is Memilio/data/pydata/Germany. 
    The folders 'pydata/Germany' have to be created if they do not exist yet. 
    Then this script can be run.

    @param read_data False or True. Defines if data is read from file or
        downloaded. Default defined in defaultDict.
    @param file_format File format which is used for writing the data.
        Default defined in defaultDict.
    @param out_folder Path to folder where data is written in folder
        out_folder/Germany. Default defined in defaultDict.
    @param no_raw True or False. Defines if unchanged raw data is written or
        not. Default defined in defaultDict. Currently not used.
    @param merge_eisenach [Default: True] or False. Defines whether the
        counties 'Wartburgkreis' and 'Eisenach' are listed separately or
        combined as one entity 'Wartburgkreis'.
    @param username Username to sign in at regionalstatistik.de. 
    @param password Password to sign in at regionalstatistik.de.
    @return DataFrame with adjusted population data for all ages to current level.
    """

    directory = os.path.join(out_folder, 'Germany')
    gd.check_dir(directory)

    df_pop_raw = read_population_data(username, password, read_data, directory)

    column_names = list(df_pop_raw.columns)
    # rename columns
    rename_columns = {
        column_names[0]: dd.EngEng['idCounty'],
        column_names[1]: dd.EngEng['county'],
        column_names[2]: dd.EngEng['ageRKI'],
        column_names[3]: dd.EngEng['number'],
        column_names[4]: dd.EngEng['male'],
        column_names[5]: dd.EngEng['female']
    }
    df_pop_raw.rename(columns=rename_columns, inplace=True)
    # remove date and explanation rows at end of table
    df_pop_raw = df_pop_raw[:np.where(df_pop_raw[dd.EngEng['idCounty']].str.contains(
        '__') == True)[0][0]].reset_index(drop=True)
    # get indices of counties first lines
    idCounty_idx = df_pop_raw.groupby(
        dd.EngEng['idCounty']).head(1).index.tolist()

    # read county list and create output data frame
    counties = np.array(geoger.get_county_names_and_ids(
        merge_berlin=True, merge_eisenach=merge_eisenach, zfill=True))
    age_cols = df_pop_raw.loc[
        idCounty_idx[0]: idCounty_idx[1] - 2,
        dd.EngEng['ageRKI']].copy().values
    for i in range(len(age_cols)):
        if i == 0:
            upper_bound = str(int(age_cols[i][
                age_cols[i].index('unter ')+6:].split(' ')[0])-1)
            age_cols[i] = '0-' + upper_bound
        elif i == len(age_cols)-1:
            lower_bound = age_cols[i].split(' ')[0]
            age_cols[i] = lower_bound + '-99'
        else:
            lower_bound = age_cols[i].split(' ')[0]
            upper_bound = str(int(age_cols[i][
                age_cols[i].index('unter ')+6:].split(' ')[0])-1)
            age_cols[i] = lower_bound + '-' + upper_bound

    df_pop = assign_population_data(
        df_pop_raw, counties, age_cols, idCounty_idx)

    test_total_population(df_pop, age_cols)

    df_pop_export = export_population_dataframe(
        df_pop, directory, file_format, merge_eisenach)

    return df_pop_export


def main():
    """! Main program entry."""

    arg_dict = gd.cli("population")
    get_population_data(**arg_dict)


if __name__ == "__main__":
    main()
