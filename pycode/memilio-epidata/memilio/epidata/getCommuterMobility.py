#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn, Lena Ploetzke
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
@file getCommuterMobility.py

@brief gets data related to county mobility from "Bundesagentur fuer Arbeit"
"""
import collections
import os
import numpy as np
import pandas as pd

from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import getPopulationData as gPd
from memilio.epidata import progress_indicator


def verify_sorted(countykey_list):
    """! verify that read countykey_list is sorted
    @param countykey_list List of county regional keys
    """
    countykey_list_is_sorted = np.all(np.array(
        countykey_list[:-1]) <= np.array(countykey_list[1:]))  # this checks if it is sorted
    if countykey_list_is_sorted:
        return True
    else:
        print('Error. Input list not sorted.')
        return False


def assign_geographical_entities(countykey_list, govkey_list):
    """! Assigns counties to governing regions based on key comparison and creates list of governing regions per state.

    Only works with sorted key lists.

    Keyword arguments:
    @param setup_dict dictionary with necessary values
    @param countykey_list List of county regional keys.
    @param govkey_list List of governing regions regional keys.

    @return countykey2govkey Hash map from county regional keys to governing region regional keys.
    @return countykey2localnumlist Hash map from county regional keys to local numbered list (per governing region).
    @return gov_county_table Table of county regional keys per governing region.
    @return state_gov_table Table of governing region regional keys per federal state.
    """

    if verify_sorted(countykey_list) == False:
        raise gd.DataError("Error. Input list not sorted.")

    # Create list of government regions with lists of counties that belong to them and list of states with government regions that belong to them; only works with sorted lists of keys.
    gov_county_table = []

    gov_index = 0
    col_index = 0
    col_list = []

    for county_key_id in range(len(countykey_list)):

        # check for belonging to currently considered government region
        if str(countykey_list[county_key_id]).startswith(str(govkey_list[gov_index])):
            # add county to current government region
            col_list.append(countykey_list[county_key_id])
            col_index += 1
        # go to next government region
        if county_key_id < len(countykey_list) - 1 and (
            not str(countykey_list[county_key_id + 1]).startswith(
                str(govkey_list[gov_index]))):
            # add government region to full table
            gov_county_table.append(col_list)
            col_list = []
            gov_index += 1
            col_index = 0
    # add last government region
    gov_county_table.append(col_list)

    if len(gov_county_table) != len(govkey_list):
        print('Error. Number of government regions wrong.')

    # create a unique hash map from county key to its government region and
    # a global key to local (in gov region) key ordering
    countykey2govkey = collections.OrderedDict()
    countykey2localnumlist = collections.OrderedDict()
    for i in range(len(gov_county_table)):
        for j in range(len(gov_county_table[i])):
            countykey2govkey[gov_county_table[i][j]] = i
            countykey2localnumlist[gov_county_table[i][j]] = j

    # create government regions list per state
    state_gov_table = []

    state_id = 1
    state_govlist_loc = []
    for govkey_id in range(len(govkey_list)):

        if str(int(govkey_list[govkey_id])).startswith(str(state_id)):
            state_govlist_loc.append(govkey_list[govkey_id])

        if govkey_id + 1 < len(govkey_list) and not str(
                int(govkey_list[govkey_id + 1])).startswith(
                str(state_id)):
            state_id += 1
            state_gov_table.append(state_govlist_loc)
            state_govlist_loc = []
    # add last state's list
    state_gov_table.append(state_govlist_loc)

    return countykey2govkey, countykey2localnumlist, gov_county_table, state_gov_table


def get_commuter_data(read_data=dd.defaultDict['read_data'],
                      file_format=dd.defaultDict['file_format'],
                      out_folder=dd.defaultDict['out_folder'],
                      no_raw=dd.defaultDict['no_raw'],
                      make_plot=dd.defaultDict['make_plot'],
                      setup_dict='',
                      ref_year=2022):
    """! Computes DataFrame of commuter migration patterns based on the Federal
    Agency of Work data.

    Keyword arguments:
    @param read_data True or False. Defines if data is read from file or downloaded.
        Only for population data. Commuter data is always downloaded. Default defined in defaultDict.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    @param out_folder Folder where data is written to. Default defined in defaultDict.
    @param no_raw [Currently not used] True or False. Defines if unchanged raw data is saved or not. Default defined in defaultDict.
    @param make_plot [Currently not used] True or False. Defines if plots are generated with matplotlib. Default defined in defaultDict.
    @param setup_dict dictionary with necessary values:
        'path': String with datapath where migration files can be found
        'abs_tol': tolerated undetected people
        'rel_tol': relative Tolerance to undetected people
    @param ref_year Year between 2013 and 2022 that specifies where the data should be taken from.
        Default value is 2022.
    @return df_commuter_migration DataFrame of commuter migration.
        df_commuter_migration[i][j]= number of commuters from county with county-id i to county with county-id j
    In commuter migration files is a cumulative value per county for number of commuters from whole Germany given.
    The printed errors are refering to the absolute and relative errors from included numbers per county in DataFrame and
    this cumulative values.
    """
    if setup_dict == '':
        abs_tol = 100  # maximum absolute error allowed per county migration
        rel_tol = 0.01  # maximum relative error allowed per county migration
        path = 'https://statistik.arbeitsagentur.de/Statistikdaten/Detail/' + \
            str(ref_year) + '12/iiia6/beschaeftigung-sozbe-krpend/'

        setup_dict = {'abs_tol': abs_tol,
                      'rel_tol': rel_tol,
                      'path': path}

    # Using the 'Einpendler' sheet to correctly distribute summed values over counties of other gov. region
    param_dict = {"sheet_name": 3, "engine": None}

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)
    mobility_dir = os.path.join(directory, 'mobility/')
    gd.check_dir(mobility_dir)

    # states 01 - 16
    states = [str(state+1).zfill(2) for state in range(16)]

    commuter_migration_files = [[] for _ in range(len(states))]

    if ref_year < 2013 or ref_year > 2022:
        raise gd.DataError('No Data available for year ' + str(ref_year) + '.')
    for state_id_file in range(len(states)):
        if ref_year <= 2020:
            # These files are in a zip folder.
            url = setup_dict['path'] + 'krpend-' + \
                states[state_id_file] + '-0-' + str(ref_year) + '12-zip.zip'
        else:
            url = setup_dict['path'] + 'krpend-' + states[state_id_file] + '-0-' + \
                str(ref_year) + '12-xlsx.xlsx?__blob=publicationFile&v=2'
        filename = 'mobility_raw_' + \
            states[state_id_file] + '_' + str(ref_year)
        filepath = os.path.join(mobility_dir) + filename + '.json'
        commuter_migration_files[state_id_file] = gd.get_file(
            filepath, url, read_data, param_dict, interactive=True)
        if not no_raw:
            gd.write_dataframe(
                commuter_migration_files[state_id_file], mobility_dir, filename, 'json')

    countykey_list = geoger.get_county_ids(merge_eisenach=True, zfill=True)
    govkey_list = geoger.get_governing_regions()

    # get population data for all countys (TODO: better to provide a corresponding method for the following lines in getPopulationData itself)
    # This is not very nice either to have the same file with either Eisenach merged or not...
    if read_data:
        population = pd.read_json(directory+'county_current_population.json')
    else:
        population = gPd.get_population_data(
            out_folder=out_folder, merge_eisenach=True, read_data=read_data)

    countypop_list = list(population[dd.EngEng["population"]])

    countykey2numlist = collections.OrderedDict(
        zip(countykey_list, list(range(len(countykey_list)))))
    govkey2numlist = collections.OrderedDict(
        zip(govkey_list, list(range(len(govkey_list)))))

    (countykey2govkey, countykey2localnumlist, gov_county_table,
     state_gov_table) = assign_geographical_entities(countykey_list, govkey_list)

    mat_commuter_migration = np.zeros(
        [len(countykey_list), len(countykey_list)])

    # maximum errors (of people not detected)
    max_abs_err = 0
    max_rel_err = 0

    # len of string to find countyIDs eg. 01001 or 16056
    len_county_id = len(countykey_list[0])

    with progress_indicator.Spinner() as p:
        for file in range(len(commuter_migration_files)):
            p.set_message('Progress '+str(file+1)+'/16')
            commuter_migration_file = commuter_migration_files[file]

            for i in range(commuter_migration_file.shape[0]):
                if (len(str(commuter_migration_file.iloc[i, 0])) == len_county_id
                        and (commuter_migration_file.iloc[i, 0]).isdigit()):
                    checksum = 0  # sum of county migration from, to be checked against sum in document

                    # make zero-filled list of counties explicitly migrated to from county considered
                    # 'implicit' migration means 'migration to' which is summed in a larger regional entity and not given in
                    # detail per county
                    counties_migratedfrom = []
                    for gov_region in range(len(gov_county_table)):
                        counties_migratedfrom.append(
                            np.zeros(len(gov_county_table[gov_region])))

                    # merge eisenach and wartburgkreis
                    commuter_migration_file.iloc[:, 2].replace(
                        '16056', '16063', inplace=True)
                    commuter_migration_file.iloc[:, 0].replace(
                        '16056', '16063', inplace=True)

                    current_col = countykey2numlist[commuter_migration_file.iloc[i, 0]]
                    curr_county_migratedto = commuter_migration_file.iloc[i, 1]
                    current_key = commuter_migration_file.iloc[i, 0]
                    # migration to itself excluded!
                    counties_migratedfrom[countykey2govkey[current_key]][
                        countykey2localnumlist[current_key]] = 1

                if not isinstance(commuter_migration_file.iloc[i, 2], float):
                    # removal of nan's, regional keys are stored as strings

                    # check if entry is a digit
                    if (commuter_migration_file.iloc[i, 2]).isdigit():
                        # explicit migration from county to county
                        if len(str(commuter_migration_file.iloc[i, 2])) == len_county_id:
                            # check if entry refers to a specific county, then set matrix value
                            current_row = countykey2numlist[commuter_migration_file.iloc[i, 2]]
                            # TODO
                            val = commuter_migration_file.iloc[i, 4]
                            mat_commuter_migration[current_row,
                                                   current_col] = val
                            checksum += val
                            counties_migratedfrom[countykey2govkey[commuter_migration_file.iloc[i, 2]]
                                                  ][countykey2localnumlist[commuter_migration_file.iloc[i, 2]]] = 1

                        # take summed values of other REMAINING counties of government region.
                        # here, some counties of the region are stated explicitly and the rest is summed.
                        elif (str(commuter_migration_file.iloc[i, 3]) == 'Übrige Kreise (Regierungsbezirk)' and str(
                                commuter_migration_file.iloc[i, 4]).isdigit()):
                            # remove trailing zeros
                            remaining_counties = str(
                                commuter_migration_file.iloc[i, 2])
                            if len(remaining_counties) > 2 and remaining_counties[2] == '0':
                                remaining_counties = remaining_counties[0:2]

                                # sum population of all counties not explicitly migrated from
                                # of the current gov region migrated from
                            dummy_pop_sum = 0
                            for mapped_county in range(
                                len(
                                    gov_county_table
                                    [govkey2numlist[remaining_counties]])):
                                if counties_migratedfrom[govkey2numlist[remaining_counties]][mapped_county] < 1:
                                    # get identifier (0-401) for county key
                                    globindex = countykey2numlist[gov_county_table[
                                        govkey2numlist[remaining_counties]][mapped_county]]
                                    # sum up
                                    dummy_pop_sum += countypop_list[globindex]

                            # distribute emigration relatively to county population where migration comes from
                            for mapped_county in range(
                                len(
                                    gov_county_table
                                    [govkey2numlist[remaining_counties]])):
                                if counties_migratedfrom[govkey2numlist[remaining_counties]][mapped_county] < 1:
                                    # get identifier (0-401) for county key
                                    globindex = countykey2numlist[gov_county_table[
                                        govkey2numlist[remaining_counties]][mapped_county]]
                                    counties_migratedfrom[govkey2numlist[remaining_counties]
                                                          ][mapped_county] = 1

                                    # set value computed relatively to county size and effective migration
                                    current_row = globindex
                                    val = commuter_migration_file.iloc[i, 4] * \
                                        countypop_list[globindex] / \
                                        dummy_pop_sum
                                    checksum += val
                                    mat_commuter_migration[current_row,
                                                           current_col] = val

                        # take summed values of ALL counties of a government region.
                        # here, no single county of the region is stated explicitly, all counties are summed together.
                        elif (commuter_migration_file.iloc[i, 2] in govkey_list and sum(
                                counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i, 2]]]) == 0):
                            # sum population of all counties not explicitly migrated to
                            # of the current gov region migrated to
                            dummy_pop_sum = 0
                            for mapped_county in range(
                                len(
                                    gov_county_table
                                    [
                                        govkey2numlist
                                        [commuter_migration_file.iloc[i, 2]]])):
                                if counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i, 2]]][mapped_county] < 1:
                                    # get identifier (0-401) for county key
                                    globindex = countykey2numlist[gov_county_table[
                                        govkey2numlist[commuter_migration_file.iloc[i, 2]]][mapped_county]]
                                    # sum up
                                    dummy_pop_sum += countypop_list[globindex]

                            # distribute emigration relatively to county population where migration comes from
                            for mapped_county in range(
                                len(
                                    gov_county_table
                                    [
                                        govkey2numlist
                                        [commuter_migration_file.iloc[i, 2]]])):
                                if counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i, 2]]][mapped_county] < 1:
                                    # get identifier (0-401) for county key
                                    globindex = countykey2numlist[gov_county_table[
                                        govkey2numlist[commuter_migration_file.iloc[i, 2]]][mapped_county]]
                                    counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i, 2]]
                                                          ][mapped_county] = 1

                                    # set value computed relatively to county size and effective migration
                                    current_row = globindex
                                    val = commuter_migration_file.iloc[i, 4] * \
                                        countypop_list[globindex] / \
                                        dummy_pop_sum
                                    checksum += val
                                    mat_commuter_migration[current_row,
                                                           current_col] = val

                        # take summed values of other REMAINING counties of a whole Bundesland
                        # here, some counties of the Bundesland are stated explicitly and the rest is summed
                        # the first or is for the case that the right first line of the incoming people directly
                        # addresses one
                        # the latter 'or's is used if no single county nor gov region of a federal state is stated
                        # explicitly
                        # although there are existent government regions in this federal state (i.e., the state itself is
                        # not considered a governement region according to gov_list)
                        elif ((str(commuter_migration_file.iloc[i, 3]) == 'Übrige Regierungsbezirke (Bundesland)' and str(
                                commuter_migration_file.iloc[i, 4]).isdigit())
                              or ((commuter_migration_file.iloc[i, 2]).isdigit() and str(
                                commuter_migration_file.iloc[i - 1][2]).startswith('nan'))
                              or (len(str(commuter_migration_file.iloc[i, 2])) == 2 and
                                  abs(float(commuter_migration_file.iloc[i, 2]) - float(
                                      commuter_migration_file.iloc[i - 1][2])) == 1)
                              or (len(str(commuter_migration_file.iloc[i, 2])) == 2 and
                                  abs(float(commuter_migration_file.iloc[i, 2]) - float(
                                      commuter_migration_file.iloc[i - 1][2])) == 2)):

                            # auxiliary key of Bundesland (key translated to int starting at zero)
                            dummy_key = int(
                                commuter_migration_file.iloc[i, 2]) - 1

                            # sum population of all counties not explicitly migrated
                            # from the current gov region migrated from
                            dummy_pop_sum = 0
                            for j in range(len(state_gov_table[dummy_key])):
                                # over all government regions not explicitly stated
                                gov_index = govkey2numlist[state_gov_table[dummy_key][j]]
                                for mapped_county in range(len(gov_county_table[gov_index])):
                                    # over all counties of the considered gov region
                                    if counties_migratedfrom[gov_index][mapped_county] < 1:
                                        # get identifier (0-401) for county key
                                        globindex = countykey2numlist[gov_county_table[gov_index]
                                                                      [mapped_county]]
                                        # sum up
                                        dummy_pop_sum += countypop_list[globindex]

                            # distribute emigration relatively to county population where migration comes from
                            for j in range(len(
                                    state_gov_table[dummy_key])):  # over all government regions not explicitly stated
                                gov_index = govkey2numlist[state_gov_table[dummy_key][j]]
                                for mapped_county in range(len(gov_county_table[gov_index])):
                                    # over all counties of the considered gov region
                                    if counties_migratedfrom[gov_index][mapped_county] < 1:
                                        # get identifier (0-401) for county key
                                        globindex = countykey2numlist[gov_county_table[gov_index]
                                                                      [mapped_county]]
                                        counties_migratedfrom[gov_index][mapped_county] = 1

                                        # set value computed relatively to county size and effective migration
                                        current_row = globindex
                                        val = commuter_migration_file.iloc[i, 4] * \
                                            countypop_list[globindex] / \
                                            dummy_pop_sum
                                        checksum += val
                                        mat_commuter_migration[current_row,
                                                               current_col] = val

                # sum of total migration 'from'
                if str(commuter_migration_file.iloc[i, 3]) == 'Einpendler aus dem Bundesgebiet':
                    abs_err = abs(
                        checksum - commuter_migration_file.iloc[i, 4])
                    if abs_err > max_abs_err:
                        max_abs_err = abs_err
                    if abs_err / checksum > max_rel_err:
                        max_rel_err = abs_err / checksum
                    if abs_err < setup_dict['abs_tol'] and abs_err / checksum < setup_dict['rel_tol']:
                        checksum = 0
                    else:
                        print('Error in calculations for county ', curr_county_migratedto,
                              '\nAccumulated values:', checksum,
                              ', correct sum:', commuter_migration_file.iloc[i, 4])
                        print('Absolute error:', abs_err,
                              ', relative error:', abs_err / checksum)

            if np.isnan(mat_commuter_migration).any():
                raise gd.DataError(
                    'NaN encountered in mobility matrix, exiting '
                    'getCommuterMobility(). Mobility data will be incomplete.')

    print('Maximum absolute error:', max_abs_err)
    print('Maximum relative error:', max_rel_err)

    countykey_list = [int(id) for id in countykey_list]
    df_commuter_migration = pd.DataFrame(
        data=mat_commuter_migration, columns=countykey_list, index=countykey_list)
    filename = 'migration_bfa_' + \
        str(ref_year) + '_dim' + str(mat_commuter_migration.shape[0])
    gd.write_dataframe(df_commuter_migration, directory, filename, file_format)

    # this is neither a very elegant nor a very general way to merge...
    # better options to be searched for!
    if 16056 in countykey_list:
        merge_id = 16063
        new_idx = countykey_list.index(geoger.CountyMerging[merge_id][0])
        old_idx = countykey_list.index(geoger.CountyMerging[merge_id][1])

        mat_commuter_migration[new_idx, :] = mat_commuter_migration[new_idx,
                                                                    :] + mat_commuter_migration[old_idx, :]
        mat_commuter_migration[:, new_idx] = mat_commuter_migration[:,
                                                                    new_idx] + mat_commuter_migration[:, old_idx]
        mat_commuter_migration[new_idx, new_idx] = 0

        mat_commuter_migration = np.delete(
            mat_commuter_migration, old_idx, axis=0)
        mat_commuter_migration = np.delete(
            mat_commuter_migration, old_idx, axis=1)

    countykey_list = geoger.get_county_ids()
    df_commuter_migration = pd.DataFrame(
        data=mat_commuter_migration, columns=countykey_list)
    df_commuter_migration.index = countykey_list
    commuter_sanity_checks(df_commuter_migration)
    filename = 'migration_bfa_' + str(ref_year)
    gd.write_dataframe(df_commuter_migration, directory, filename, file_format)
    directory = directory.split('pydata')[0] + 'mobility/'
    gd.check_dir(directory)
    gd.write_dataframe(
        df_commuter_migration, directory,
        'commuter_migration_scaled_' + str(ref_year),
        'txt', {'sep': ' ', 'index': False, 'header': False})

    return df_commuter_migration


def commuter_sanity_checks(df):
    # Dataframe should be of squared form
    if len(df.index) != len(df.columns):
        raise gd.DataError("Error. Dataframe should be of squared form.")
    # There were 401 counties at beginning of 2021 and 400 at end of 2021.
    # Check if exactly 400 counties are in dataframe.
    if not len(df) == 400:
        raise gd.DataError("Error. Size of dataframe unexpected.")


def get_neighbors_mobility(
        countyid, direction='both', abs_tol=0, rel_tol=0, tol_comb='or',
        out_folder=dd.defaultDict['out_folder'], ref_year=2022):
    '''! Returns the neighbors of a particular county ID depening on the
    commuter mobility and given absolute and relative thresholds on the number
    of commuters.

    The parameters absolute and relative tolerance decide which connections and
    neighbors are returned. If tol_comb='or', only one of this two criteria
    has to be satisfied to count the edges. If 'and' is chosen, both criteria 
    have to be satisfied.

    @param countyid ID of the county where mobility is considered and for which
        neighbors have to be returned.
    @param direction 'both' [Default], 'in', or 'out'. Defines whether 'both' or  
        'in' or 'out' commuters only are considered.
    @param abs_tol Minimum number of commuters to count the connection.
    @param rel_tol Relative tolerance with respect to the strongest connection 
        of the county to count the connections.
    @param tol_comb Defines whether absolute and relative thresholds are
        combined such that only one criterion has to be satisfied ('or') or
        both ('and').
    @param merge_eisenach [Default: True] Defines whether the counties
        'Wartburgkreis' and 'Eisenach' are listed separately or combined
        as one entity 'Wartburgkreis'.
    @param out_folder Folder where data is written to. Default defined in defaultDict.
    @param ref_year Year between 2013 and 2022 that specifies where the data should be taken from.
        Default value is 2022.
    @return Neighbors of the county with respect to mobility and the number of 
        commuters from and to the neighbors.
    '''
    # This is not very nice either to have the same file with either Eisenach merged or not...
    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)
    try:
        commuter = gd.get_file(os.path.join(
            directory, "migration_bfa_"+str(ref_year)+"_dim400.json"), read_data=True)
    except FileNotFoundError:
        print("Commuter data was not found. Download and process it from the internet.")
        commuter = get_commuter_data(out_folder=out_folder, ref_year=ref_year)

    countykey_list = commuter.columns
    commuter.index = countykey_list

    # compute in and out density
    if direction == 'both':
        commuter_all = commuter.loc[countyid, :] + commuter.loc[:, countyid]
    elif direction == 'in':
        commuter_all = commuter.loc[:, countyid]
    elif direction == 'out':
        commuter_all = commuter.loc[countyid, :]
    if tol_comb == 'and':
        neighbor_indices = np.where((commuter_all > abs_tol) & (
            commuter_all > commuter_all.max()*rel_tol))[0]
    elif tol_comb == 'or':
        neighbor_indices = np.where((commuter_all > abs_tol) | (
            commuter_all > commuter_all.max()*rel_tol))[0]

    return countykey_list[neighbor_indices], commuter_all.values[neighbor_indices]


def get_neighbors_mobility_all(
        direction='both', abs_tol=0, rel_tol=0, tol_comb='or',
        out_folder=dd.defaultDict['out_folder'], ref_year=2022):
    '''! Returns the neighbors of all counties ID depening on the
    commuter mobility and given absolute and relative thresholds on the number
    of commuters.

    The parameters absolute and relative tolerance decide which connections and
    neighbors are returned. If tol_comb='or', only one of this two criteria
    has to be satisfied to count the edges. If 'and' is chosen, both criteria 
    have to be satisfied.

    @param direction 'both' [Default], 'in', or 'out'. Defines whether 'both' or  
        'in' or 'out' commuters only are considered.
    @param abs_tol Minimum number of commuters to count the connection.
    @param rel_tol Relative tolerance with respect to the strongest connection 
        of the county to count the connections.
    @param tol_comb Defines whether absolute and relative thresholds are
        combined such that only one criterion has to be satisfied ('or') or
        both ('and')
    @param ref_year Year between 2013 and 2022 that specifies where the data should be taken from.
        Default value is 2022.
    @return Neighbors of all counties with respect to mobility.
    '''
    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)
    countyids = geoger.get_county_ids()
    neighbors_table = []
    # TODO: performance has to be improved
    for id in countyids:
        neighbors_table.append(
            get_neighbors_mobility(
                id, direction=direction, abs_tol=abs_tol,
                rel_tol=rel_tol, tol_comb=tol_comb,
                out_folder=out_folder, ref_year=ref_year))

    return dict(zip(countyids, neighbors_table))


def main():
    """! Main program entry."""

    arg_dict = gd.cli("commuter_official")
    ref_year = 2022

    abs_tol = 100  # maximum absolute error allowed per county migration
    rel_tol = 0.01  # maximum relative error allowed per county migration
    path = 'https://statistik.arbeitsagentur.de/Statistikdaten/Detail/' + \
        str(ref_year) + '12/iiia6/beschaeftigung-sozbe-krpend/'

    setup_dict = {'abs_tol': abs_tol,
                  'rel_tol': rel_tol,
                  'path': path}

    arg_dict_commuter = {**arg_dict, "setup_dict": setup_dict}

    get_neighbors_mobility(
        1001, abs_tol=0, rel_tol=0, tol_comb='or',
        out_folder=dd.defaultDict['out_folder'])

    get_commuter_data(**arg_dict_commuter)


if __name__ == "__main__":
    main()
