#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn, Lena Ploetzke, Henrik Zunker
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
import geojson
import os
import wget
import numpy as np
import pandas as pd
from zipfile import ZipFile
from shapely import geometry
from memilio.epidata import getPopulationData as gPd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import defaultDict as dd
from os.path import dirname as up


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

    # Create list of government regions with lists of counties that belong to them and list of states with government
    # regions that belong to them; only works with sorted lists of keys.
    gov_county_table = []

    gov_index = 0
    col_index = 0
    col_list = []

    for i in range(0, len(countykey_list)):

        # check for belonging to currently considered government region
        if str(countykey_list[i]).startswith(str(govkey_list[gov_index])):
            # add county to current government region
            col_list.append(countykey_list[i])
            col_index += 1
        # go to next government region
        if i < len(countykey_list) - 1 and (
            not str(countykey_list[i + 1]).startswith(
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
    for i in range(0, len(gov_county_table)):
        for j in range(0, len(gov_county_table[i])):
            countykey2govkey[gov_county_table[i][j]] = i
            countykey2localnumlist[gov_county_table[i][j]] = j

    # create government regions list per state
    state_gov_table = []

    state_id = 1
    state_govlist_loc = []
    for i in range(0, len(govkey_list)):

        if str(int(govkey_list[i])).startswith(str(state_id)):
            state_govlist_loc.append(govkey_list[i])

        if i + 1 < len(govkey_list) and not str(
                int(govkey_list[i + 1])).startswith(
                str(state_id)):
            state_id += 1
            state_gov_table.append(state_govlist_loc)
            state_govlist_loc = []
    # add last state's list
    state_gov_table.append(state_govlist_loc)

    return countykey2govkey, countykey2localnumlist, gov_county_table, state_gov_table


def get_distance(coordinates_from, coordinates_to):
    """! Calculates the distance between two coordinates

    @param coordinates_from Vector containing two start coordinates
    @param coordinates_to Vector containing two end coordinates
    """
    dist_x = np.abs(coordinates_from[0] - coordinates_to[0])
    dist_y = np.abs(coordinates_from[1] - coordinates_to[1])
    return np.sqrt(dist_x**2 + dist_y**2)


def get_counties_center_coordinates(
        path_geojson,
        file_format=dd.defaultDict['file_format'],
        out_folder=dd.defaultDict['out_folder']):
    """! Computes centers of counties based on a geojson file of germany.

    In the case of multipolygons, we choose the largest polygon and take its center.
    The results are saved in a json file.

    Keyword arguments:
    @param path_geojson geojson file which should contain the geodata for all counties in germany
    @param out_folder Path to folder where data is saved
    """

    with open(path_geojson) as f:
        file = geojson.load(f)

    countykey_list = geoger.get_county_ids(merge_eisenach=False, zfill=True)
    len_features = len(file.features)

    # If number of features not equal to number of county keys, then we set
    # merge_eisenach as True. If it stays unequal -> output error message
    if len_features != len(countykey_list):
        countykey_list = geoger.get_county_ids(
            merge_eisenach=True, zfill=True)
        if len_features != len(countykey_list):
            raise ValueError(
                'Dimension of countykeylist and given geojson is different')

    # Create list of center points coordinates.
    # List is sorted in ascending order of regional county keys.
    centers_counties = np.zeros((len_features, 2))

    # iterate over all regional counties and calculate their center.
    for i in range(0, len_features):
        features = file.features[i]

        # check if countykeys are the same in geojson and countykey_list
        if countykey_list[i] != features.properties["RS"]:
            raise ValueError(
                'Regional county keys are not the same or sorted differently')

        # Special cases are multipolygons. Here we calculate the largest polygon and take his center.
        if features['geometry']['type'] == 'MultiPolygon':
            area_polygon = 0.
            for j in features['geometry']['coordinates']:
                features['geometry']['coordinates']
                # cast coordinates to polygons and calculate center points
                poly = geometry.Polygon(j[0])
                if area_polygon < poly.area:
                    area_polygon = poly.area
                    centers_counties[i][0] = poly.centroid.x
                    centers_counties[i][1] = poly.centroid.y
        else:
            # In case of regular Polygons, cast coordinates to polygons and calculate center points
            poly = geometry.Polygon(features['geometry']['coordinates'][0])
            centers_counties[i][0] = poly.centroid.x
            centers_counties[i][1] = poly.centroid.y

    directory = os.path.join(up(out_folder), 'mobility/')
    gd.check_dir(directory)
    df_centers_counties = pd.DataFrame(
        data=centers_counties)
    filename = 'county_centers_dim' + str(centers_counties.shape[0])
    gd.write_dataframe(df_centers_counties, directory, filename, file_format)

    return None


def get_distances_from_centers(center_coordinates):
    """! Computes the distances between all counties based on the center points.

    By using center points, distances between counties can be significantly larger than they actually are.

    Keyword arguments:
    @param center_coordinates Json file which contain the center coordinates for all counties in germany
    """

    num_counties = center_coordinates.shape[0]

    if num_counties == 400:
        countykey_list = geoger.get_county_ids(merge_eisenach=True, zfill=True)
    else:
        countykey_list = geoger.get_county_ids(merge_eisenach=False, zfill=True)

    # calculate distances between counties
    distances = np.zeros((num_counties, num_counties))
    for lk1 in range(0, num_counties):
        id_lk1 = countykey_list[lk1]
        lk1_coor = [center_coordinates[0][lk1], center_coordinates[1][lk1]]
        for lk2 in range(0, num_counties):
            id_lk2 = countykey_list[lk2]

            # distance to itself
            if id_lk1 == id_lk2:
                distances[lk1][lk2] = 0.
            else:
                lk2_coor = [center_coordinates[0][lk2], center_coordinates[1][lk2]]
                distances[lk1][lk2] = get_distance(lk1_coor, lk2_coor)

    return distances


def scale_commuter_mobility(commuter_mobility, center_coordinates):
    """! Computes scaled commuter migration patterns based on the Federal
    Agency of Work data. The Scaling is as presented in https://www.medrxiv.org/content/10.1101/2020.12.18.20248509v1.full.pdf

    The Scaling is still a very heuristic approach and based on the distances between the counties.
    
    Keyword arguments:
    @param commuter_mobility DataFrame of commuter migration patterns based on the Federal Agency of Work data without any scaling
    @param center_coordinates Json file which contain the center coordinates for all counties in germany

    @return commuter_mobility Array of scaled commuter migration.
        commuter_mobility[i][j] = scaled number of commuters from county with county-id i to county with county-id j
    """

    dim_data = commuter_mobility.shape[0]
    distances_counties = get_distances_from_centers(center_coordinates)

    # Distances and commuter_mobility should have the same size
    if dim_data != distances_counties.shape[0]:
        raise ValueError(
            'Dimensions of distances and commuter mobility are different')

    for i in range(0, dim_data):
        for j in range(0, dim_data):

            # Commuting between 100-200km happens two/three times a week -> Scale with 1/2
            if distances_counties[i][j] > 1 and distances_counties[i][j] < 2:
                commuter_mobility[i][j] *= 0.5

            # Commuting +200km happens only once a week -> Scale with 1/5
            if distances_counties[i][j] > 2:
                commuter_mobility[i][j] *= 0.2

    return commuter_mobility


def get_commuter_data(setup_dict='',
                      center_coordinates="",
                      read_data=dd.defaultDict['read_data'],
                      file_format=dd.defaultDict['file_format'],
                      out_folder=dd.defaultDict['out_folder'],
                      make_plot=dd.defaultDict['make_plot'],
                      no_raw=dd.defaultDict['no_raw']):
    """! Computes DataFrame of commuter migration patterns based on the Federal
    Agency of Work data.

    Keyword arguments:
    @param setup_dict dictionary with necessary values:
        'path': String with datapath where migration files can be found
        'abs_tol': tolerated undetected people
        'rel_tol': relative Tolerance to undetected people

    @return df_commuter_migration_scaled DataFrame of commuter migration.
        df_commuter_migration_scaled[i][j]= number of commuters from county with county-id i to county with county-id j
    In commuter migration files is a cumulative value per county for number of commuters from whole Germany given.
    The printed errors are refering to the absolute and relative errors from included numbers per county in DataFrame and
    this cumulative values.
    """
    if setup_dict == '':
        ref_year = 2020
        abs_tol = 100  # maximum absolute error allowed per county migration
        rel_tol = 0.01  # maximum relative error allowed per county migration
        path = 'https://statistik.arbeitsagentur.de/Statistikdaten/Detail/' + \
            str(ref_year) + '12/iiia6/beschaeftigung-sozbe-krpend/'

        setup_dict = {'abs_tol': abs_tol,
                      'rel_tol': rel_tol,
                      'path': path}

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    countykey_list = geoger.get_county_ids(merge_eisenach=False, zfill=True)
    govkey_list = geoger.get_governing_regions()

    # get population data for all countys (TODO: better to provide a corresponding method for the following lines in getPopulationData itself)
    # This is not very nice either to have the same file with either Eisenach merged or not...

    population = gPd.get_population_data(
        out_folder=out_folder, merge_eisenach=False, read_data=read_data)

    countypop_list = list(population[dd.EngEng["population"]])

    countykey2numlist = collections.OrderedDict(
        zip(countykey_list, list(range(0, len(countykey_list)))))
    govkey2numlist = collections.OrderedDict(
        zip(govkey_list, list(range(0, len(govkey_list)))))

    (countykey2govkey, countykey2localnumlist, gov_county_table, state_gov_table) = assign_geographical_entities(
        countykey_list, govkey_list)

    mat_commuter_migration = np.zeros(
        [len(countykey_list), len(countykey_list)])

    # maxium errors (of people not detected)
    max_abs_err = 0
    max_rel_err = 0

    files = []
    for n in range(1, 17):
        # files.append('krpend_' + str(n).zfill(2) + "_0.xlsx")
        files.append(
            'krpend-' + str(n).zfill(2) + "-0-20" +
            setup_dict['path'].split('/20')[1][0:4] + "-zip.zip")

    n = 0

    for item in files:
        # Using the 'Einpendler' sheet to correctly distribute summed values over counties of other gov. region
        # This File is in a zip folder so it has to be unzipped first before it can be read.
        param_dict={"sheet_name": 3, "engine": "pyxlsb"}
        filepath = os.path.join(out_folder, 'Germany/')
        url = setup_dict['path'] + item.split('.')[0] + '.zip'
        # Unzip it
        zipfile = wget.download(url, filepath)
        with ZipFile(zipfile, 'r') as zipObj:
            zipObj.extractall(path = filepath)
        # Read the file
        filename = item.split('-20')[0] + '.xlsb'
        file = filename.replace('-','_')
        commuter_migration_file = pd.read_excel(filepath + file, **param_dict)
        # pd.read_excel(os.path.join(setup_dict['path'], item), sheet_name=3)

        # delete zip folder after extracting
        os.remove(os.path.join(filepath, item))
        # delete file after reading
        os.remove(os.path.join(filepath, file))

        counties_done = []  # counties considered as 'migration from'
        # current_row = -1  # row of matrix that belongs to county migrated from
        current_col = -1  # column of matrix that belongs to county migrated to
        checksum = 0  # sum of county migration from, to be checked against sum in document

        for i in range(0, commuter_migration_file.shape[0]):
            if (len(str(commuter_migration_file.iloc[i][0])) == 5
                    and (commuter_migration_file.iloc[i][0]).isdigit()):
                checksum = 0
                # make zero'd list of counties explicitly migrated to from county considered
                # 'implicit' migration means 'migration to' which is summed in a larger regional entity and not given in
                # detail per county
                counties_migratedfrom = []
                for j in range(0, len(gov_county_table)):
                    counties_migratedfrom.append(
                        np.zeros(len(gov_county_table[j])))

                counties_done.append(commuter_migration_file.iloc[i][0])
                current_col = countykey2numlist[commuter_migration_file.iloc[i][0]]
                curr_county_migratedto = commuter_migration_file.iloc[i][1]
                current_key = commuter_migration_file.iloc[i][0]
                # migration to itself excluded!
                counties_migratedfrom[countykey2govkey[current_key]
                                      ][countykey2localnumlist[current_key]] = 1

            if not isinstance(commuter_migration_file.iloc[i][2], float):
                # removal of nan's, regional keys are stored as strings

                # check if entry is a digit
                if (commuter_migration_file.iloc[i][2]).isdigit():
                    # explicit migration from county to county
                    if len(str(commuter_migration_file.iloc[i][2])) == 5:
                        # check if entry refers to a specific county, then set matrix value
                        current_row = countykey2numlist[commuter_migration_file.iloc[i][2]]
                        # TODO
                        val = commuter_migration_file.iloc[i][4]
                        mat_commuter_migration[current_row, current_col] = val
                        checksum += val
                        counties_migratedfrom[countykey2govkey[commuter_migration_file.iloc[i][2]]][
                            countykey2localnumlist[commuter_migration_file.iloc[i][2]]] = 1

                    # take summed values of other REMAINING counties of government region
                    # here, some counties of the region are stated explicitly and the rest is summed
                    elif (str(commuter_migration_file.iloc[i][3]) == 'Übrige Kreise (Regierungsbezirk)' and str(
                            commuter_migration_file.iloc[i][4]).isdigit()):
                        # remove trailing zeros (dummy key w/o zeros: dummy_key_wozeros)
                        dummy_key_wozeros = str(
                            commuter_migration_file.iloc[i][2])
                        if len(dummy_key_wozeros) > 2 and dummy_key_wozeros[2] == '0':
                            dummy_key_wozeros = dummy_key_wozeros[0:2]

                            # sum population of all counties not explicitly migrated from
                            # of the current gov region migrated from
                        dummy_pop_sum = 0
                        for k in range(0, len(gov_county_table[govkey2numlist[dummy_key_wozeros]])):
                            if counties_migratedfrom[govkey2numlist[dummy_key_wozeros]][k] < 1:
                                # get identifier (0-401) for county key
                                globindex = countykey2numlist[gov_county_table[govkey2numlist[dummy_key_wozeros]][k]]
                                # sum up
                                dummy_pop_sum += countypop_list[globindex]

                        # distribute emigration relatively to county population where migration comes from
                        for k in range(0, len(gov_county_table[govkey2numlist[dummy_key_wozeros]])):
                            if counties_migratedfrom[govkey2numlist[dummy_key_wozeros]][k] < 1:
                                # get identifier (0-401) for county key
                                globindex = countykey2numlist[gov_county_table[govkey2numlist[dummy_key_wozeros]][k]]
                                counties_migratedfrom[govkey2numlist[dummy_key_wozeros]][k] = 1

                                # set value computed relatively to county size and effective migration
                                current_row = globindex
                                val = commuter_migration_file.iloc[i][4] * \
                                    countypop_list[globindex] / dummy_pop_sum
                                checksum += val
                                mat_commuter_migration[current_row,
                                                       current_col] = val

                    # take summed values of ALL counties of a government region
                    # here, no single county of the region is stated explicitly, all counties are summed together
                    elif (commuter_migration_file.iloc[i][2] in govkey_list and sum(
                            counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i][2]]]) == 0):
                        # sum population of all counties not explicitly migrated to
                        # of the current gov region migrated to
                        dummy_pop_sum = 0
                        for k in range(0, len(gov_county_table[govkey2numlist[commuter_migration_file.iloc[i][2]]])):
                            if counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i][2]]][k] < 1:
                                # get identifier (0-401) for county key
                                globindex = countykey2numlist[gov_county_table[govkey2numlist[
                                    commuter_migration_file.iloc[i][2]]][k]]
                                # sum up
                                dummy_pop_sum += countypop_list[globindex]

                        # distribute emigration relatively to county population where migration comes from
                        for k in range(0, len(gov_county_table[govkey2numlist[commuter_migration_file.iloc[i][2]]])):
                            if counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i][2]]][k] < 1:
                                # get identifier (0-401) for county key
                                globindex = countykey2numlist[gov_county_table[govkey2numlist[
                                    commuter_migration_file.iloc[i][2]]][k]]
                                counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i][2]]][k] = 1

                                # set value computed relatively to county size and effective migration
                                current_row = globindex
                                val = commuter_migration_file.iloc[i][4] * \
                                    countypop_list[globindex] / dummy_pop_sum
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
                    elif ((str(commuter_migration_file.iloc[i][3]) == 'Übrige Regierungsbezirke (Bundesland)' and str(
                            commuter_migration_file.iloc[i][4]).isdigit())
                          or ((commuter_migration_file.iloc[i][2]).isdigit() and str(
                              commuter_migration_file.iloc[i - 1][2]).startswith('nan'))
                          or (len(str(commuter_migration_file.iloc[i][2])) == 2 and
                              abs(float(commuter_migration_file.iloc[i][2]) - float(
                                  commuter_migration_file.iloc[i - 1][2])) == 1)
                          or (len(str(commuter_migration_file.iloc[i][2])) == 2 and
                              abs(float(commuter_migration_file.iloc[i][2]) - float(
                                  commuter_migration_file.iloc[i - 1][2])) == 2)):

                        # auxiliary key of Bundesland (key translated to int starting at zero)
                        dummy_key = int(commuter_migration_file.iloc[i][2]) - 1

                        # sum population of all counties not explicitly migrated
                        # from the current gov region migrated from
                        dummy_pop_sum = 0
                        for j in range(0, len(state_gov_table[dummy_key])):
                            # over all government regions not explicitly stated
                            gov_index = govkey2numlist[state_gov_table[dummy_key][j]]
                            for k in range(0,
                                           len(gov_county_table[gov_index])):
                                # over all counties of the considered gov region
                                if counties_migratedfrom[gov_index][k] < 1:
                                    # get identifier (0-401) for county key
                                    globindex = countykey2numlist[gov_county_table[gov_index][k]]
                                    # sum up
                                    dummy_pop_sum += countypop_list[globindex]

                        # distribute emigration relatively to county population where migration comes from
                        for j in range(0, len(
                                state_gov_table[dummy_key])):  # over all government regions not explicitly stated
                            gov_index = govkey2numlist[state_gov_table[dummy_key][j]]
                            for k in range(0,
                                           len(gov_county_table[gov_index])):
                                # over all counties of the considered gov region
                                if counties_migratedfrom[gov_index][k] < 1:
                                    # get identifier (0-401) for county key
                                    globindex = countykey2numlist[gov_county_table[gov_index][k]]
                                    counties_migratedfrom[gov_index][k] = 1

                                    # set value computed relatively to county size and effective migration
                                    current_row = globindex
                                    val = commuter_migration_file.iloc[i][4] * \
                                        countypop_list[globindex] / \
                                        dummy_pop_sum
                                    checksum += val
                                    mat_commuter_migration[current_row,
                                                           current_col] = val

            # sum of total migration 'from'
            if str(commuter_migration_file.iloc[i][3]) == 'Einpendler aus dem Bundesgebiet':
                abs_err = abs(checksum - commuter_migration_file.iloc[i][4])
                if abs_err > max_abs_err:
                    max_abs_err = abs_err
                if abs_err / checksum > max_rel_err:
                    max_rel_err = abs_err / checksum
                if abs_err < setup_dict['abs_tol'] and abs_err / checksum < setup_dict['rel_tol']:
                    checksum = 0
                else:
                    print('Error in calculations for county ', curr_county_migratedto,
                          '\nAccumulated values:', checksum,
                          ', correct sum:', commuter_migration_file.iloc[i][4])
                    print('Absolute error:', abs_err,
                          ', relative error:', abs_err / checksum)

        n += 1
        print(' Federal state read. Progress ', n, '/ 16')
        if np.isnan(mat_commuter_migration).any():
            raise gd.DataError(
                'NaN encountered in mobility matrix, exiting '
                'getCommuterMobility(). Mobility data will be incomplete.')
    if n != 16:
        print('Error. Files missing.')

    print('Maximum absolute error:', max_abs_err)
    print('Maximum relative error:', max_rel_err)

    countykey_list = [int(id) for id in countykey_list]
    df_commuter_migration = pd.DataFrame(
        data=mat_commuter_migration, columns=countykey_list)
    df_commuter_migration.index = countykey_list
    filename = 'migration_bfa_20' + files[0].split(
        '-20')[1][0:2] + '_dim' + str(mat_commuter_migration.shape[0])
    gd.write_dataframe(df_commuter_migration, directory, filename, file_format)

    # this is neither a very elegant nor a very general way to merge...
    # better options to be searched for!
    merge_id = 16063
    new_idx = countykey_list.index(geoger.CountyMerging[merge_id][0])
    old_idx = countykey_list.index(geoger.CountyMerging[merge_id][1])

    mat_commuter_migration[new_idx, :] = mat_commuter_migration[new_idx,
                                                                :] + mat_commuter_migration[old_idx, :]
    mat_commuter_migration[:, new_idx] = mat_commuter_migration[:,
                                                                new_idx] + mat_commuter_migration[:, old_idx]
    mat_commuter_migration[new_idx, new_idx] = 0

    mat_commuter_migration = np.delete(mat_commuter_migration, old_idx, axis=0)
    mat_commuter_migration = np.delete(mat_commuter_migration, old_idx, axis=1)

    countykey_list = geoger.get_county_ids()
    df_commuter_migration = pd.DataFrame(
        data=mat_commuter_migration, columns=countykey_list)
    df_commuter_migration.index = countykey_list
    commuter_sanity_checks(df_commuter_migration)
    filename = 'migration_bfa_20' + files[0].split(
        '-20')[1][0:2] + '_dim' + str(mat_commuter_migration.shape[0])
    gd.write_dataframe(df_commuter_migration, directory, filename, file_format)
    gd.check_dir(os.path.join(directory.split('pydata')[0], 'mobility'))
    df_commuter_migration.to_csv(
        directory.split('pydata')[0] + 'mobility/commuter_mobility' +
        '_20' + files[0].split('-20')[1][0: 2] + '.txt', sep=' ', index=False,
        header=False)

    # Short end for cases where no center coordinates are given.
    # this is especially important for 'test_epidata_get_vaccination_data'
    if isinstance(center_coordinates, str):
        return df_commuter_migration
        
    commuter_migration_scaled = scale_commuter_mobility(
        mat_commuter_migration, center_coordinates)

    df_commuter_migration_scaled = pd.DataFrame(
        data=commuter_migration_scaled, columns=countykey_list)
    df_commuter_migration_scaled.index = countykey_list
    filename = 'commuter_mobility_scaled_20' + files[0].split(
        '-20')[1][0:2] + '_dim' + str(mat_commuter_migration.shape[0])
    gd.write_dataframe(df_commuter_migration_scaled,
                       directory, filename, file_format)
    gd.check_dir(os.path.join(directory.split('pydata')[0], 'mobility'))
    df_commuter_migration_scaled.to_csv(
        directory.split('pydata')[0] + 'mobility/commuter_mobility_scaled' +
        '_20' + files[0].split('-20')[1][0: 2] + '.txt', sep=' ', index=False,
        header=False)

    return df_commuter_migration_scaled


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
        merge_eisenach=True, out_folder=dd.defaultDict['out_folder'],
        center_coordinates=""):
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
    @return Neighbors of the county with respect to mobility and the number of 
        commuters from and to the neighbors.
    '''
    # This is not very nice either to have the same file with either Eisenach merged or not...
    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)
    try:
        if merge_eisenach:
            commuter = pd.read_json(os.path.join(
                directory, "migration_bfa_2020_dim400.json"))
        else:
            commuter = pd.read_json(os.path.join(
                directory, "migration_bfa_2020_dim401.json"))
    except ValueError:
        print("Commuter data was not found. Download and process it from the internet.")
        commuter = get_commuter_data(
            out_folder=out_folder, center_coordinates=center_coordinates)

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
        merge_eisenach=True, out_folder=dd.defaultDict['out_folder'],
        center_coordinates=""):
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
    @return Neighbors of all counties with respect to mobility.
    '''
    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)
    countyids = geoger.get_county_ids(merge_eisenach=merge_eisenach)
    neighbors_table = []
    for id in countyids:
        neighbors_table.append(
            get_neighbors_mobility(
                id, direction=direction, abs_tol=abs_tol,
                rel_tol=rel_tol, tol_comb=tol_comb,
                merge_eisenach=merge_eisenach,
                out_folder=out_folder, 
                center_coordinates=center_coordinates))

    return dict(zip(countyids, neighbors_table))


def main():
    """! Main program entry."""

    arg_dict = gd.cli("commuter_official")
    ref_year = 2020

    abs_tol = 100  # maximum absolute error allowed per county migration
    rel_tol = 0.01  # maximum relative error allowed per county migration
    path = 'https://statistik.arbeitsagentur.de/Statistikdaten/Detail/' + \
        str(ref_year) + '12/iiia6/beschaeftigung-sozbe-krpend/'

    setup_dict = {'abs_tol': abs_tol,
                  'rel_tol': rel_tol,
                  'path': path}

    memilio_path = up(up(up(up(up(__file__)))))
    path_json = os.path.join(
            memilio_path, 'data', 'mobility' , 'county_centers_dim400.json')
    
    if os.path.isfile(path_json):
        center_coordinates = pd.read_json(path_json)
    else:
        center_coordinates = ""
        


    get_neighbors_mobility(1001, abs_tol=0, rel_tol=0, tol_comb='or',
                           merge_eisenach=True, out_folder=dd.defaultDict['out_folder'],
                           center_coordinates=center_coordinates)

    mat_commuter_migration = get_commuter_data(
        setup_dict, **arg_dict, center_coordinates=center_coordinates)


if __name__ == "__main__":
    main()
