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
<<<<<<< HEAD
import sys
=======
import twill
import time
import io

>>>>>>> e13e4c705c931ad490aea943f701f918bdc8803b
import numpy as np
import pandas as pd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger


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


<<<<<<< HEAD
        try:
            # if this file is encoded with utf-8 German umlauts are not displayed correctly because they take two bytes
            # utf_8_sig can identify those bytes as one sign and display it correctly
            zensus = gd.loadCsv(
                "abad92e8eead46a4b0d252ee9438eb53_1", param_dict={"encoding":'utf_8_sig'})
        except FileNotFoundError:
            error_message = "Error: The zensus file does not exist."
            raise FileNotFoundError(error_message)
=======
def assign_population_data(df_pop_raw, counties, age_cols, idCounty_idx):
    '''! Assigns population data of all counties of old dataframe in new created dataframe
>>>>>>> e13e4c705c931ad490aea943f701f918bdc8803b

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

<<<<<<< HEAD
        empty_data = ['.', '-']
        for i in idCounty_idx:
            # county information needed
            if df_pop_raw.loc[i, dd.EngEng['idCounty']] in counties[:, 1]:
                # direct assignment of population data found
                df_pop.loc[df_pop[dd.EngEng['idCounty']] == df_pop_raw.loc
                           [i, dd.EngEng['idCounty']],
                           age_cols] = df_pop_raw.loc[i: i + nage - 1, dd.EngEng
                                                      ['number']].values
            # Berlin and Hamburg
            elif df_pop_raw.loc[i, dd.EngEng['idCounty']] + '000' in counties[:, 1]:
                # direct assignment of population data found
                df_pop.loc[df_pop[dd.EngEng['idCounty']] == df_pop_raw.loc
                           [i, dd.EngEng['idCounty']] + '000',
                           age_cols] = df_pop_raw.loc[i: i + nage - 1, dd.EngEng
                                                      ['number']].values
            # empty rows
            elif df_pop_raw.loc[i: i + nage - 1,
                                dd.EngEng['number']].values.any() in empty_data:
                if not df_pop_raw.loc[i: i + nage - 1,
                                      dd.EngEng['number']].values.all() in empty_data:
                    raise gd.DataError(
                        'Error. Partially incomplete data for county ' +
                        df_pop_raw.loc[i, dd.EngEng['idCounty']])
            # additional information for local entities not needed
            elif df_pop_raw.loc[i, dd.EngEng['idCounty']] in ['03241001', '05334002', '10041100']:
                pass
            # Germany, federal states, and governing regions
            elif len(df_pop_raw.loc[i, dd.EngEng['idCounty']]) < 5:
                pass
            else:
                print('no data for ' + df_pop_raw.loc
                      [i, dd.EngEng['idCounty']])
                raise gd.DataError(
                    'Error. County ID in input population data'
                    'found which could not be assigned.')

        if df_pop[age_cols].sum().sum() != 83155031:
            raise gd.DataError('Wrong total population size')

        new_cols = [
            dd.EngEng['idCounty'],
            dd.EngEng['population'],
            '<3 years', '3-5 years', '6-14 years', '15-17 years',
            '18-24 years', '25-29 years', '30-39 years', '40-49 years',
            '50-64 years', '65-74 years', '>74 years']
        df_pop_export = pd.DataFrame(columns=new_cols)
        df_pop_export[df_pop.columns[0]] = df_pop[dd.EngEng['idCounty']]
        # <3 and 3-5
        df_pop_export[df_pop_export.columns[2:4]] = df_pop[df_pop.columns[2:4]]
        # 6-14
        df_pop_export[df_pop_export.columns[4]] = \
            df_pop[df_pop.columns[4:6]].sum(axis=1)
        # 15-17
        df_pop_export[df_pop_export.columns[5]] = df_pop[df_pop.columns[6]]
        # 18-24
        df_pop_export[df_pop_export.columns[6]] = \
            df_pop[df_pop.columns[7:9]].sum(axis=1)
        # 25-29
        df_pop_export[df_pop_export.columns[7]] = df_pop[df_pop.columns[9]]
        # 30-39
        df_pop_export[df_pop_export.columns[8]] = \
            df_pop[df_pop.columns[10:12]].sum(axis=1)
        # 40-49
        df_pop_export[df_pop_export.columns[9]] = \
            df_pop[df_pop.columns[12:14]].sum(axis=1)
        # 50-64
        df_pop_export[df_pop_export.columns[10]] = \
            df_pop[df_pop.columns[14:17]].sum(axis=1)
        # 65-74
        df_pop_export[df_pop_export.columns[11]] = df_pop[df_pop.columns[17]]
        # >74
        df_pop_export[df_pop_export.columns[12]] = df_pop[df_pop.columns[18]]

        df_pop_export[dd.EngEng['population']
                      ] = df_pop_export.iloc[:, 2:].sum(axis=1)

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        if len(df_pop_export) == 401:
            filename = 'county_current_population_dim401'
            gd.write_dataframe(df_pop_export, directory, filename, file_format)

        if len(df_pop_export) == 400 or merge_eisenach == True:
            filename = 'county_current_population'

            # Merge Eisenach and Wartburgkreis
            df_pop_export = geoger.merge_df_counties_all(
                df_pop_export, sorting=[dd.EngEng["idCounty"]],
                columns=dd.EngEng["idCounty"])

            gd.write_dataframe(df_pop_export, directory, filename, file_format)

        return df_pop_export

    else:
        counties, zensus, reg_key = load_population_data(
            out_folder, read_data=read_data, no_raw=no_raw,
            file_format=file_format)

        # find region keys for census population data
        key = np.zeros((len(zensus)))
        for i in range(len(key)):
            for j in range(len(reg_key)):
                if zensus['Name'].values[i] == reg_key['NAME'].values.astype(
                        str)[j]:
                    if zensus[dd.invert_dict(dd.GerEng)[dd.EngEng['population']]].values[i] == round(
                            reg_key['Zensus_EWZ'].values[j] * 1000):
                        key[i] = reg_key['AGS'].values[j]

        inds = np.unique(key, return_index=True)[1]
        # columns of downloaded data which should be replaced
        male = ['M_Unter_3', 'M_3_bis_5', 'M_6_bis_14', 'M_15_bis_17',
                'M_18_bis_24', 'M_25_bis_29', 'M_30_bis_39', 'M_40_bis_49',
                'M_50_bis_64', 'M_65_bis_74', 'M_75_und_aelter']
        female = ['W_Unter_3', 'W_3_bis_5', 'W_6_bis_14', 'W_15_bis_17',
                  'W_18_bis_24', 'W_25_bis_29', 'W_30_bis_39', 'W_40_bis_49',
                  'W_50_bis_64', 'W_65_bis_74', 'W_75_und_aelter']

        if not split_gender:
            # get data from zensus file and add male and female population data
            data = np.zeros((len(inds), len(male)+2))
            data[:, 0] = key[inds].astype(int)
            data[:, 1] = zensus[dd.invert_dict(
                dd.GerEng)[dd.EngEng['population']]].values[inds].astype(int)
            for i in range(len(male)):
                data[:, i+2] = zensus[male[i]].values[inds].astype(
                    int) + zensus[female[i]].values[inds].astype(int)

            # define new columns for dataframe
            columns = [
                dd.EngEng['idCounty'],
                dd.EngEng['population'],
                '<3 years', '3-5 years', '6-14 years', '15-17 years',
                '18-24 years', '25-29 years', '30-39 years', '40-49 years',
                '50-64 years', '65-74 years', '>74 years']
=======
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
>>>>>>> e13e4c705c931ad490aea943f701f918bdc8803b
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

<<<<<<< HEAD
        # compute ratio of current and 2011 population data
        ratio = np.ones(len(data[:, 0]))
        for i in range(len(ratio)):
            for j in range(len(counties)):
                if not counties[dd.EngEng['idCounty']].isnull().values[j]:
                    try:
                        if data[i, 0] == int(
                                counties[dd.EngEng['idCounty']].values[j]):
                            ratio[i] = counties[dd.EngEng['population']].values[j]/data[i, 1]

                    except ValueError:
                        pass

        # adjust population data for all ages to current level
        data_current = np.zeros(data.shape)
        data_current[:, 0] = data[:, 0].copy()
        for i in range(len(data[0, :]) - 1):
            data_current[:, i + 1] = np.multiply(data[:, i + 1], ratio)

        # create dataframe
        df_current = pd.DataFrame(
            np.round(data_current).astype(int), columns=columns)

        directory = os.path.join(out_folder, 'Germany/')
        gd.check_dir(directory)

        if len(df_current) == 400 or merge_eisenach == True:
            # Merge Eisenach and Wartburgkreis
            df_current = geoger.merge_df_counties_all(
                df_current, sorting=[dd.EngEng["idCounty"]],
                columns=dd.EngEng["idCounty"])
            df = geoger.merge_df_counties_all(
                df, sorting=[dd.EngEng["idCounty"]],
                columns=dd.EngEng["idCounty"])

            filename = 'county_current_population'
            filename_raw = 'county_population'
        else:  # Write Dataframe without merging
            if (len(df_current) != 400) and (len(df_current) != 401):
                print('Population output only contains ' +
                      str(len(df_current)) + ' counties. Is this intended?')
            filename = 'county_current_population_dim' + str(len(df_current))
            filename_raw = 'county_population_dim' + str(len(df_current))

        gd.write_dataframe(df_current, directory, filename, file_format)
        gd.write_dataframe(df, directory, filename_raw, file_format)

        return df_current
=======
    return df_pop_export
>>>>>>> e13e4c705c931ad490aea943f701f918bdc8803b


def main():
    """! Main program entry."""

    arg_dict = gd.cli("population")
    get_population_data(**arg_dict)


if __name__ == "__main__":
    main()
