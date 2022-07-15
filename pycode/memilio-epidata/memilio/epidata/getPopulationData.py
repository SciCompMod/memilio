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

import os
import sys
import numpy as np
import pandas as pd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger


def get_new_counties(data):
    """! Creates 7 new counties that were formed since 2011 and deletes old counties

   Downloaded data is from 2011.
   However, new counties have been defined.
   Thus, the new ones have to be added, the data has to be calculated and the old counties have to be deleted.

   @param data_temp Pandas dataframe
   @return Changed data
   """

    # create 7 new counties
    data_temp = np.append(data, np.zeros((7, data.shape[1])), axis=0)

    # Göttingen
    data_temp[-7, 0] = 3159

    # Mecklenburgische Seenplatte
    data_temp[-6, 0] = 13071

    # Landkreis Rostock
    data_temp[-5, 0] = 13072

    # Vorpommern Rügen
    data_temp[-4, 0] = 13073

    # Nordwestmecklenburg
    data_temp[-3, 0] = 13074

    # Vorpommern Greifswald
    data_temp[-2, 0] = 13075

    # Ludwigslust-Parchim
    data_temp[-1, 0] = 13076

    to_delete = []

    for i in range(len(data_temp[:, 0])):
        # fuse "Göttingen" and "Osterode am Harz" into Göttingen
        if data_temp[i, 0] in [3152, 3156]:
            data_temp[-7, 1:] += data_temp[i, 1:]
            to_delete.append(i)

        # fuse "Müritz", "Neubrandenburg", "Mecklenburg-Sterlitz"
        # and "Demmin" into "Mecklenburgische Seenplatte"
        if data_temp[i, 0] in [13056, 13002, 13055, 13052]:
            data_temp[-6, 1:] += data_temp[i, 1:]
            to_delete.append(i)

        # fuse "Bad Doberan and Güstrow" into "Landkreis Rostosck"
        if data_temp[i, 0] in [13051, 13053]:
            data_temp[-5, 1:] += data_temp[i, 1:]
            to_delete.append(i)

        # fuse "Rügen", "Stralsund" and "Nordvorpommern" into
        # "Vorpommern Rügen"
        if data_temp[i, 0] in [13061, 13005, 13057]:
            data_temp[-4, 1:] += data_temp[i, 1:]
            to_delete.append(i)

        # fuse "Wismar" and "Nordwestmecklenburg" into
        # "Nordwestmecklenburg"
        if data_temp[i, 0] in [13006, 13058]:
            data_temp[-3, 1:] += data_temp[i, 1:]
            to_delete.append(i)

        # fuse "Ostvorpommern", "Uecker-Randow" and "Greifswald"
        # into "Vorpommern Greifswald"
        if data_temp[i, 0] in [13059, 13062, 13001]:
            data_temp[-2, 1:] += data_temp[i, 1:]
            to_delete.append(i)

        # fuse "Ludwigslust" and "Parchim" into "Ludwigslust-Parchim"
        if data_temp[i, 0] in [13054, 13060]:
            data_temp[-1, 1:] += data_temp[i, 1:]
            to_delete.append(i)

    data_temp = np.delete(data_temp, to_delete, 0)
    sorted_inds = np.argsort(data_temp[:, 0])
    data_temp = data_temp[sorted_inds, :]
    return data_temp


def load_population_data(data_folder,
                         read_data=dd.defaultDict['read_data'],
                         no_raw=dd.defaultDict['no_raw'],
                         file_format=dd.defaultDict['file_format']):
    """! Load of counties, zensus and reg_key files

    Data is downloaded from the following sources
   - Federal Statistical Office of Germany (Destatis/Genesis-Online), current
        data for population per county [stored in "counties"]
   - Zensus2011 data with additional information on regional keys
        [stored in "reg_key"]
   - Zensus2011 data from opendata splitted for age and gender
        [stored in "zensus"]

   Data is either downloaded or read from "data_folder"/Germany/.

   @param data_folder Path to folder where data is written in folder data_folder/Germany.
   @param read_data False or True. Defines if data is read from file or downloaded. Default defined in defaultDict.
   @param no_raw True or False. Defines if unchanged raw data is written or not. Default defined in defaultDict.
   @param file_format File format which is used for writing the data. Default defined in defaultDict.
   @return 3 Dataframes of migration, reg_key and zensus
    """

    directory = os.path.join(data_folder, 'Germany/')
    gd.check_dir(directory)

    filename_counties = 'migration'
    filename_zensus = 'zensus'
    filename_reg_key = 'reg_key'

    if read_data:

        # Read counties File
        file_in = os.path.join(directory, filename_counties + ".json")
        try:
            counties = pd.read_json(file_in)
        except ValueError:
            error_message = "Error: The file: " + file_in + \
                "could not be read. Call program without -r flag to get it."
            raise FileNotFoundError(error_message)

        # Read Zensus File
        file_in = os.path.join(directory, filename_zensus + ".json")
        try:
            zensus = pd.read_json(file_in)
        except ValueError:
            error_message = "Error: The file: " + file_in + \
                "could not be read. Call program without -r flag to get it."
            raise FileNotFoundError(error_message)

        # Read reg_key File
        file_in = os.path.join(directory, filename_reg_key + ".json")
        try:
            reg_key = pd.read_json(file_in)
        except ValueError:
            error_message = "Error: The file: " + file_in + \
                "could not be read. Call program without -r flag to get it."
            raise FileNotFoundError(error_message)
    else:
        try:
            path_counties = 'https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales/Gemeindeverzeichnis/Administrativ/04-kreise.xlsx;?__blob=publicationFile'
            counties = gd.loadExcel(
                targetFileName='', apiUrl=path_counties, extension='',
                param_dict={"sheet_name": 1, "header": 3,
                            "engine": 'openpyxl'})
        except ValueError:
            error_message = "Error: The counties file does not exist."
            raise FileNotFoundError(error_message)

        # Download zensus

        try:
            # if this file is encoded with utf-8 German umlauts are not displayed correctly because they take two bytes
            # utf_8_sig can identify those bytes as one sign and display it correctly
            zensus = gd.loadCsv(
                "abad92e8eead46a4b0d252ee9438eb53_1", param_dict={"encoding":'utf_8_sig'})
        except ValueError:
            error_message = "Error: The zensus file does not exist."
            raise FileNotFoundError(error_message)

        # Download reg_key

        try:
            path_reg_key = 'https://www.zensus2011.de/SharedDocs/Downloads/DE/Pressemitteilung/DemografischeGrunddaten/' \
                           '1A_EinwohnerzahlGeschlecht.xls?__blob=publicationFile&v=5'
            # read tables
            reg_key = gd.loadExcel(path_reg_key, apiUrl='', extension='', param_dict={
                                   "engine": None, "sheet_name": 'Tabelle_1A', "header": 12})
        except ValueError:
            error_message = "Error: The reg-key file does not exist."
            raise FileNotFoundError(error_message)

        if not no_raw:
            if not counties.empty:
                gd.write_dataframe(counties, directory,
                                   filename_counties, file_format)
            if not zensus.empty:
                gd.write_dataframe(zensus, directory,
                                   filename_zensus, file_format)
            if not reg_key.empty:
                gd.write_dataframe(reg_key, directory,
                                   filename_reg_key, file_format)

    return counties, zensus, reg_key


def get_population_data(data_folder,
                        read_data=dd.defaultDict['read_data'],
                        file_format=dd.defaultDict['file_format'],
                        no_raw=dd.defaultDict['no_raw'],
                        split_gender=False,
                        merge_eisenach=True):
    """! Download age-stratified population data for the German counties.

    There are two different data sources that can be transformed in a simple
    data file with age-resolved information per county and on which our other
    tools can continue to operate on.

    1.) Official 'Bevölkerungsfortschreibung' 12411-02-03-4:
    'Bevölkerung nach Geschlecht und Altersgruppen (17)' 
    of regionalstatistik.de. 
    ATTENTION: The raw file cannot be downloaded 
    automatically by our scripts without an Genesis Online account. In order to
    work on this dataset, please manually download it from:

    https://www.regionalstatistik.de/genesis/online -> "1: Gebiet, Bevölkerung,
    Arbeitsmarkt, Wahlen" -> "12: Bevölkerung" -> "12411 Fortschreibung des
    Bevölkerungsstandes" ->  "12411-02-03-4: Bevölkerung nach Geschlecht und 
    Altersgruppen (17) - Stichtag 31.12. - regionale Tiefe: Kreise und
    krfr. Städte". 

    Download the xlsx file and put it under data_folder, 
    this normally is Memilio/data/pydata/Germany. 
    The folders 'pydata/Germany' have to be created if they do not exist yet. 
    Then this script can be run.
    TODO: The following parameters have no effect for this source.

    2.) Combination of data from Federal Statistical Office of Germany and 
    Zensus2011 to compute very simple approximations of age-stratified number
    of inhabitants per county. Population from the Zensus data of all age groups is scaled to the
    total population of our more recent migration data by a factor which
    represents the relative increase/decrease in population size
    between 2011 and 2019 for each county"
    This data can either be downloaded automatically or read from
    "data_folder"/Germany/ if it was downloaded before.

    @param data_folder Path to folder where data is written in folder
        data_folder/Germany. 
    @param read_data False or True. Defines if data is read from file or
        downloaded. Default defined in defaultDict.
    @param file_format File format which is used for writing the data.
        Default defined in defaultDict.
    @param no_raw True or False. Defines if unchanged raw data is written or
        not. Default defined in defaultDict.
    @param split_gender [Default: False] or True. Defines whether data is
        splitted by gender
    @param merge_eisenach [Default: True] or False. Defines whether the
        counties 'Wartburgkreis' and 'Eisenach' are listed separately or
        combined as one entity 'Wartburgkreis'.
    @return DataFrame with adjusted population data for all ages to current level.
    """
    directory = os.path.join(data_folder, 'Germany')
    gd.check_dir(directory)
    filename = '12411-02-03-4'  # '12411-09-01-4-B'
    new_data_file = os.path.join(directory, filename)
    new_data_avail = os.path.isfile(new_data_file + '.xlsx')

    if new_data_avail:
        print('Information: Using new population data file ' + filename)
        df_pop_raw = gd.loadExcel(
            new_data_file, apiUrl='', extension='.xlsx',
            param_dict={"engine": "openpyxl", "sheet_name": filename, "header": 4})
        column_names = list(df_pop_raw.columns)
        # rename columns
        rename_columns = {
            column_names[0]: dd.EngEng['idCounty'],
            column_names[1]: dd.EngEng['county'],
            column_names[2]: dd.EngEng['ageRKI'],
            column_names[3]: dd.EngEng['number'],
            column_names[4]: dd.EngEng['male'],
            column_names[5]: dd.EngEng['female'],
        }
        df_pop_raw.rename(columns=rename_columns, inplace=True)
        # remove date and explanation rows at end of table
        # -33 for '12411-09-01-4-B', -38 for '12411-02-03-4'
        df_pop_raw = df_pop_raw[:-38].reset_index(drop=True)
        # get indices of counties first lines
        idCounty_idx = list(df_pop_raw.index[
            df_pop_raw[dd.EngEng['idCounty']].isna() == False])

        # read county list and create output data frame
        counties = np.array(geoger.get_county_names_and_ids(
            merge_berlin=True, merge_eisenach=False, zfill=True))
        new_cols = {dd.EngEng['idCounty']: counties[:, 1],
                    dd.EngEng['county']: counties[:, 0]}
        age_cols = df_pop_raw.loc[
            idCounty_idx[0]: idCounty_idx[1] - 2,  # -1 for '12411-09-01-4-B'
            dd.EngEng['ageRKI']].copy().values
        nage = len(age_cols)
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

        # add age_cols with zero initilization to new_cols
        new_cols.update({age: 0 for age in age_cols})
        df_pop = pd.DataFrame(new_cols)

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

        filename = 'county_current_population_dim401'
        gd.write_dataframe(df_pop_export, directory, filename, file_format)

        if merge_eisenach == True:
            filename = 'county_current_population'

            # Merge Eisenach and Wartburgkreis
            df_pop_export = geoger.merge_df_counties_all(
                df_pop_export, sorting=[dd.EngEng["idCounty"]],
                columns=dd.EngEng["idCounty"])

            gd.write_dataframe(df_pop_export, directory, filename, file_format)

        return df_pop_export

    else:
        counties, zensus, reg_key = load_population_data(
            data_folder, read_data=read_data, no_raw=no_raw,
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
        else:
            # get data from zensus file
            data = np.zeros((len(inds), len(male)+len(female)+2))
            data[:, 0] = key[inds].astype(int)
            data[:, 1] = zensus[dd.invert_dict(
                dd.GerEng)[dd.EngEng['population']]].values[inds].astype(int)
            for i in range(len(male)):
                data[:, i+2] = zensus[male[i]].values[inds].astype(int)
            for i in range(len(female)):
                data[:, i+len(male)+2] = zensus[female[i]
                                                ].values[inds].astype(int)

            # define new columns for dataframe
            columns = [
                dd.EngEng['idCounty'],
                dd.EngEng['population'],
                'M <3 years', 'M 3-5 years', 'M 6-14 years', 'M 15-17 years',
                'M 18-24 years', 'M 25-29 years', 'M 30-39 years',
                'M 40-49 years', 'M 50-64 years', 'M 65-74 years',
                'M >74 years', 'F <3 years', 'F 3-5 years', 'F 6-14 years',
                'F 15-17 years', 'F 18-24 years', 'F 25-29 years',
                'F 30-39 years', 'F 40-49 years', 'F 50-64 years',
                'F 65-74 years', 'F >74 years']

        data = get_new_counties(data)

        # create Dataframe of raw data without adjusting population
        df = pd.DataFrame(data.astype(int), columns=columns)

        # compute ratio of current and 2011 population data
        ratio = np.ones(len(data[:, 0]))
        for i in range(len(ratio)):
            for j in range(len(counties)):
                if not counties['Schlüssel-nummer'].isnull().values[j]:
                    try:
                        if data[i, 0] == int(
                                counties['Schlüssel-nummer'].values[j]):
                            ratio[i] = counties['Bevölkerung2)'].values[j]/data[i, 1]

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

        if merge_eisenach == True:
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
            filename = 'county_current_population_dim401'
            filename_raw = 'county_population_dim401'

        gd.write_dataframe(df_current, directory, filename, file_format)
        gd.write_dataframe(df, directory, filename_raw, file_format)

        return df_current


def main():
    """! Main program entry."""

    path = os.path.join(os.getcwd(), 'data', 'pydata')
    arg_dict = gd.cli("population")
    get_population_data(path, **arg_dict)


if __name__ == "__main__":
    main()
