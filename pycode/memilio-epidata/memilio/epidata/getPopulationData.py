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


def load_population_data(out_folder=dd.defaultDict['out_folder'],
                         read_data=dd.defaultDict['read_data'],
                         no_raw=dd.defaultDict['no_raw'],
                         file_format=dd.defaultDict['file_format']):
    """! Load of counties, zensus and reg_key files

    Data is downloaded from the following sources
   - Federal Statistical Office of Germany (Destatis), Genesis-Online
   current data for population per county [stored in "counties"]
   - Zensus2011 data splitted by gender for whole germany, states, counties in xls
   with additional data from 30.04.2011 (just used for reg_key?) [stored in "reg_key"]
   - Zensus2011 data from opendata splitted for age and gender [stored in "zensus"]

   Data is either downloaded or read from "out_folder"/Germany/.

   @param out_folder Path to folder where data is written in folder out_folder/Germany. Default defined in defaultDict.
   @param read_data False or True. Defines if data is read from file or downloaded. Default defined in defaultDict.
   @param no_raw True or False. Defines if unchanged raw data is written or not. Default defined in defaultDict.
   @param file_format File format which is used for writing the data. Default defined in defaultDict.
   @return 3 Dataframes of migration, reg_key and zensus
    """

    directory = os.path.join(out_folder, 'Germany/')
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
            exit_string = "Error: The file: " + file_in + " does not exist. "\
                          "Call program without -r flag to get it."
            sys.exit(exit_string)

        # Read Zensus File
        file_in = os.path.join(directory, filename_zensus + ".json")
        try:
            zensus = pd.read_json(file_in)
        except ValueError:
            exit_string = "Error: The file: " + file_in + " does not exist. "\
                          "Call program without -r flag to get it."
            sys.exit(exit_string)

        # Read reg_key File
        file_in = os.path.join(directory, filename_reg_key + ".json")
        try:
            reg_key = pd.read_json(file_in)
        except ValueError:
            exit_string = "Error: The file: " + file_in + " does not exist. "\
                          "Call program without -r flag to get it."
            sys.exit(exit_string)

    else:
        try:
            print('Trying to download data from the internet')
            path_counties = 'https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales/Gemeindeverzeichnis/Administrativ/04-kreise.xlsx;?__blob=publicationFile'
            counties = gd.loadExcel(
                targetFileName='', apiUrl=path_counties, extension='',
                param_dict={"sheet_name": 1, "header": 3,
                            "engine": 'openpyxl'})
        except ValueError:
            exit_string = "Error: The counties file does not exist."
            sys.exit(exit_string)

        # Download zensus
        try:
            # if this file is encoded with utf-8 German umlauts are not displayed correctly because they take two bytes
            # utf_8_sig can identify those bytes as one sign and display it correctly
            zensus = gd.loadCsv(
                "abad92e8eead46a4b0d252ee9438eb53_1", encoding='utf_8_sig')
        except ValueError:
            exit_string = "Error: The zensus file does not exist."
            sys.exit(exit_string)

        # Download reg_key
        try:
            path_reg_key = 'https://www.zensus2011.de/SharedDocs/Downloads/DE/Pressemitteilung/DemografischeGrunddaten/' \
                           '1A_EinwohnerzahlGeschlecht.xls?__blob=publicationFile&v=5'
            # read tables
            reg_key = gd.loadExcel(path_reg_key, apiUrl='', extension='', param_dict={
                                   "engine": None, "sheet_name": 'Tabelle_1A', "header": 12})
        except ValueError:
            exit_string = "Error: The regional key file does not exist."
            sys.exit(exit_string)

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


def get_population_data(read_data=dd.defaultDict['read_data'],
                        file_format=dd.defaultDict['file_format'],
                        out_folder=dd.defaultDict['out_folder'],
                        no_raw=dd.defaultDict['no_raw'],
                        split_gender=False,
                        merge_eisenach=True):
    """! Download data with age splitting

   Data is downloaded from the following sources
   - Federal Statistical Office of Germany (Destatis), Genesis-Online
   current data for population per county [stored in "counties"]
   - Zensus2011 data splitted by gender for whole germany, states, counties in xls
   with additional data from 30.04.2011 (just used for reg_key?) [stored in "reg_key"]
   - Zensus2011 data from opendata splitted for age and gender [stored in "zensus"]

   Data is either downloaded or read from "out_folder"/Germany/.

   Working with the data includes:
   - For the Zensus data the male and female data is added to get just the age dependence
   - Population from the Zensus data of all age groups is scaled to the total population of our more recent migration
   data by a factor which represents the relative increase/decrease in population size
   between 2011 and 2019 for each county"

   @param read_data False or True. Defines if data is read from file or downloaded. Default defined in defaultDict.
   @param file_format File format which is used for writing the data. Default defined in defaultDict.
   @param out_folder Path to folder where data is written in folder out_folder/Germany. Default defined in defaultDict.
   @param no_raw True or False. Defines if unchanged raw data is written or not. Default defined in defaultDict.
   @param split_gender [Default: False] or True. Defines whether data is splitted by gender
   @param merge_eisenach [Default: True] or False. Defines whether the counties 'Wartburgkreis' and 'Eisenach' are listed separately or combined as one entity 'Wartburgkreis'.
   @return DataFrame with adjusted population data for all ages to current level.
    """
    counties, zensus, reg_key = load_population_data(
        out_folder, read_data=read_data, no_raw=no_raw,
        file_format=file_format)

    # find region keys for census population data
    key = np.zeros((len(zensus)))
    for i in range(len(key)):
        for j in range(len(reg_key)):
            if zensus['Name'].values[i] == reg_key['NAME'].values.astype(str)[
                    j]:
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
            data[:, i+len(male)+2] = zensus[female[i]].values[inds].astype(int)

        # define new columns for dataframe
        columns = [dd.EngEng['idCounty'],
                   dd.EngEng['population'],
                   'M <3 years', 'M 3-5 years', 'M 6-14 years',
                   'M 15-17 years', 'M 18-24 years', 'M 25-29 years',
                   'M 30-39 years', 'M 40-49 years', 'M 50-64 years',
                   'M 65-74 years', 'M >74 years', 'F <3 years', 'F 3-5 years',
                   'F 6-14 years', 'F 15-17 years', 'F 18-24 years',
                   'F 25-29 years', 'F 30-39 years', 'F 40-49 years',
                   'F 50-64 years', 'F 65-74 years', 'F >74 years']

    data = get_new_counties(data)

    # compute ratio of current and 2011 population data
    ratio = np.ones(len(data[:, 0]))
    for i in range(len(ratio)):
        for j in range(len(counties)):
            if not counties['Schlüssel-nummer'].isnull().values[j]:
                try:
                    if data[i, 0] == int(counties['Schlüssel-nummer'].values[j]):
                        ratio[i] = counties['Bevölkerung2)'].values[j]/data[i, 1]
                except:
                    pass

    # adjust population data for all ages to current level
    data_current = np.zeros(data.shape)
    data_current[:, 0] = data[:, 0].copy()
    for i in range(len(data[0, :]) - 1):
        data_current[:, i + 1] = np.multiply(data[:, i + 1], ratio)

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    # create dataframe
    df_current = pd.DataFrame(
        np.round(data_current).astype(int), columns=columns)
    if merge_eisenach == True:
        # Merge Eisenach and Wartburgkreis
        df_current = geoger.merge_df_counties_all(
            df_current, sorting=[dd.EngEng["idCounty"]],
            columns=dd.EngEng["idCounty"])
        filename = 'county_current_population'

    else:
        # Write Dataframe without merging
        filename = 'county_current_population_dim401'
    gd.write_dataframe(df_current, directory, filename, file_format)

    return df_current


def main():
    """! Main program entry."""

    arg_dict = gd.cli("population")
    get_population_data(**arg_dict)


if __name__ == "__main__":
    main()
