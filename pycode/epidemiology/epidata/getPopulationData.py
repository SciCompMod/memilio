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

To download the data two functions are called
"""

import os
import sys
from collections import namedtuple
import numpy as np
import pandas
from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from epidemiology.epidata import geoModificationGermany as geoger


def get_population_data(read_data=dd.defaultDict['read_data'],
                        file_format=dd.defaultDict['file_format'],
                        out_folder=dd.defaultDict['out_folder'],
                        no_raw=dd.defaultDict['no_raw']):
    """! Downloads population data

   A list of all datasets with population data is composed.

   In a loop over this list, the different datasets are downloaded.
   The dataset consists of:
   - filename to store unchanged data
   - arcis opendata number (identifier of data)
   - Wanted columns
   - filename for results

   Data is stored and perhaps loaded from out_folder/Germany.

   @param read_data False [Default] or True. Defines if data is read from file or downloaded.
   @param file_format File format which is used for writing the data. Default defined in defaultDict.
   @param out_folder Path to folder where data is written in folder out_folder/Germany.
   @param no_raw True or False [Default]. Defines if unchanged raw data is saved or not.
   """

    print("Warning: getpopulationdata is not working correctly. A bug workaround has been applied.")

    Data = namedtuple("Data", "filename item columns_wanted filename_out")

    d1 = Data("FullDataB", '5dc2fc92850241c3be3d704aa0945d9c_2', ["LAN_ew_RS", 'LAN_ew_GEN', 'LAN_ew_EWZ'],
              "PopulStates")
    # d2 = Data("FullDataL", 'b2e6d8854d9744ca88144d30bef06a76_1', ['RS', 'GEN', 'EWZ'], "PopulCounties")

    # d = [d1, d2]
    d = [d1]

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)


    for d_i in d:
        get_one_data_set(read_data, file_format, no_raw, directory, d_i)


def get_one_data_set(read_data, file_format, no_raw, directory, d):
    """! download one dataset

   Data is either downloaded from website or loaded as specific json file from the given "directory".
   The final data is also stored in "directory".
   After load or download the columns are renamed.

   @param read_data False [Default] or True. Defines if data is read from file or downloaded.
   @param file_format File format which is used for writing the data. Default defined in defaultDict.
   @param directory Directory which wiles should be stored and loaded (out_folder/Germany).
   @param d Dataset with (filename to store unchanged data, arcis opendata number, wanted columns, filename for results)
   """

    if read_data:
        # if once dowloaded just read json file
        file = os.path.join(directory, d.filename + ".json")

        try:
            df = pandas.read_json(file)

        except ValueError:
            exit_string = "Error: The file: " + file + " does not exist. Call program without -r flag to get it."
            sys.exit(exit_string)
    else:

        # Supported data formats:
        load = {
            'csv': gd.loadCsv,
            'geojson': gd.loadGeojson
        }

        # Get data:
        df = load['csv'](d.item)

        # output data to not always download it
        if not no_raw:
            gd.write_dataframe(df, directory, d.filename, "json")

    print("Available columns:", df.columns)

    # Filter data for Bundesland/Landkreis and Einwohnerzahl (EWZ)
    dfo = df[d.columns_wanted]
    dfo = dfo.rename(columns=dd.GerEng)
    gd.write_dataframe(dfo, directory, d.filename_out, file_format)


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
    if to_delete == []:
        return data
    else:
        return data_temp


def load_age_population_data(out_folder):

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    filename_counties = 'migration'
    filename_reg_key = 'reg_key'
    filename_zensus = 'zensus'

    file_in = os.path.join(directory, filename_counties + ".json")
    try:
        counties = pandas.read_json(file_in)
    except:
        print('Local counties dataframe not found.')
        try:
            print('Trying to download from HPC server')
            path_counties = 'http://hpcagainstcorona.sc.bs.dlr.de/data/migration/'
            counties = gd.loadExcel(targetFileName = 'kreise_deu', apiUrl =  path_counties, extension = '.xlsx',
                                param_dict = {"sheet_name": 1, "header": 3})
            gd.write_dataframe(counties, directory, filename_counties, "json")
        except:
            print('No access to HPC Server.')
            try:
                print('Trying to download data from the internet')
                path_counties = 'https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales/Gemeindeverzeichnis/Administrativ/04-kreise.xlsx;?__blob=publicationFile'
                counties = gd.loadExcel(targetFileName= '', apiUrl =  path_counties, extension = '',
                                param_dict = {"sheet_name": 1, "header": 3, "engine": 'openpyxl'})
                gd.write_dataframe(counties, directory, filename_counties, "json")
            except ValueError:
                exit_string = "Error: The counties file does not exist."
                sys.exit(exit_string)

    file_in = os.path.join(directory, filename_zensus + ".json")
    try:
        zensus = pandas.read_json(file_in)
    except:
        print('Local zensus Dataframe not found.')
        try:
            print('Trying to download from the internet')
            zensus = gd.loadCsv("abad92e8eead46a4b0d252ee9438eb53_1")
            gd.write_dataframe(zensus, directory, filename_zensus, "json")
        except ValueError:
            exit_string = "Error: The zensus file does not exist."
            sys.exit(exit_string)

    file_in = os.path.join(directory, filename_reg_key + ".json")
    try:
        reg_key = pandas.read_json(file_in)
    except:
        print('Local reg_key Dataframe not found.')
        try:
            print('Trying to download from the internet')
            path_reg_key = 'https://www.zensus2011.de/SharedDocs/Downloads/DE/Pressemitteilung/DemografischeGrunddaten/' \
                           '1A_EinwohnerzahlGeschlecht.xls?__blob=publicationFile&v=5'
            # read tables
            reg_key = gd.loadExcel(path_reg_key, apiUrl='', extension = '',
                               param_dict = {"engine": None, "sheet_name": 'Tabelle_1A', "header": 12})
            gd.write_dataframe(reg_key, directory, filename_reg_key, "json")
        except ValueError:
            exit_string = "Error: The regional key file does not exist."
            sys.exit(exit_string)

    return counties, reg_key, zensus


def get_age_population_data(read_data=dd.defaultDict['read_data'],
                            file_format=dd.defaultDict['file_format'],
                            out_folder=dd.defaultDict['out_folder'],
                            no_raw=dd.defaultDict['no_raw'],
                            write_df=True,
                            merge_eisenach=True):
    """! Download data with age splitting

   Data is downloaded from the following sources
   - our own migration [stored in "counties"]
   - Zensus2011 data splitted by gender for whole germany, states, counties in xls
   with additional data from 30.04.2011 (just used for reg_key?) [stored in "reg_key"]
   - Zensus2011 data from opendata splitted for age and gender [stored in "zensus"]

   Data is either downloaded or read from "out_folder"/Germany/.

   Working with the data includes:
   - For the Zensus data the male and female data is added to get just the age dependence
   - Population from the Zensus data of all age groups is scaled to the total population of our more recent migration
   data by a factor which represents the relative increase/decrease in population size
   between 2011 and 2019 for each county"

   @param read_data False [Default] or True. Defines if data is read from file or downloaded.
   @param file_format File format which is used for writing the data. Default defined in defaultDict.
   @param out_folder Path to folder where data is written in folder out_folder/Germany.
    """
    counties, reg_key, zensus = load_age_population_data(out_folder)

    # find region keys for census population data
    key = np.zeros((len(zensus)))
    for i in range(len(key)):
        for j in range(len(reg_key)):
            if zensus.Name.values[i] == reg_key['NAME'].values.astype(str)[j]:
                if zensus.EWZ.values[i] == round(reg_key['Zensus_EWZ'].values[j]*1000):
                    key[i] = reg_key['AGS'].values[j]

    unique, inds, count = np.unique(key, return_index=True, return_counts=True)

    male = ['M_Unter_3', 'M_3_bis_5', 'M_6_bis_14', 'M_15_bis_17', 'M_18_bis_24',
            'M_25_bis_29', 'M_30_bis_39', 'M_40_bis_49', 'M_50_bis_64',
            'M_65_bis_74', 'M_75_und_aelter']
    female = ['W_Unter_3', 'W_3_bis_5', 'W_6_bis_14', 'W_15_bis_17', 'W_18_bis_24',
              'W_25_bis_29', 'W_30_bis_39', 'W_40_bis_49', 'W_50_bis_64',
              'W_65_bis_74', 'W_75_und_aelter']
    columns = [dd.EngEng['idCounty'], 'Total', '<3 years', '3-5 years', '6-14 years', '15-17 years', '18-24 years',
               '25-29 years', '30-39 years', '40-49 years', '50-64 years',
               '65-74 years', '>74 years']

    # add male and female population data
    data = np.zeros((len(inds), len(male)+2))
    data[:,0] = key[inds].astype(int)
    data[:,1] = zensus['EWZ'].values[inds].astype(int)
    for i in range(len(male)):
        data[:, i+2] = zensus[male[i]].values[inds].astype(int) + zensus[female[i]].values[inds].astype(int)

    data = get_new_counties(data)

    # compute ratio of current and 2011 population data
    ratio = np.ones(len(data[:,0]))
    for i in range(len(ratio)):
        for j in range(len(counties)):
            if not counties['Schlüssel-nummer'].isnull().values[j]:
                try:
                    if data[i,0] == int(counties['Schlüssel-nummer'].values[j]):
                        ratio[i] = counties['Bevölkerung2)'].values[j]/data[i, 1]
                except:
                    dummy = 0

    # adjust population data for all ages to current level
    data_current = np.zeros(data.shape)
    data_current[:, 0] = data[:, 0].copy()
    for i in range(len(data[0, :]) - 1):
        data_current[:, i + 1] = np.multiply(data[:, i + 1], ratio)

    #create dataframe
    df = pandas.DataFrame(data.astype(int), columns=columns)
    df_current = pandas.DataFrame(np.round(data_current).astype(int), columns=columns)
    df_401 = df_current.copy()    
    df_current_401 = df_current.copy()
    # From official county list, merge Eisenach and 
    # Wartburgkreis for now (as of oct. 2021)
    df_current = geoger.merge_df_counties_all(
        df_current, sorting=[dd.EngEng["idCounty"]],
        columns=dd.EngEng["idCounty"])

    df = geoger.merge_df_counties_all(
        df, sorting=[dd.EngEng["idCounty"]],
        columns=dd.EngEng["idCounty"])


    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    if write_df:
        gd.write_dataframe(df, directory, 'county_population', file_format)
        gd.write_dataframe(df_current, directory, 'county_current_population', file_format)
        # TODO there should be a more elegant version to write different version with Eisenach merged or not
        # or it should be prevented directly to write if Eisenach is not merged... to discuss...    
        gd.write_dataframe(df_401, directory, 'county_population', file_format)    
        gd.write_dataframe(df_current_401, directory, 'county_current_population_dim401', file_format)

    if merge_eisenach:
        return df_current
    else:
        return df_current_401


def main():
    """! Main program entry."""

    arg_dict = gd.cli("population")
    get_age_population_data(**arg_dict)
    get_population_data(**arg_dict)


if __name__ == "__main__":
    main()
