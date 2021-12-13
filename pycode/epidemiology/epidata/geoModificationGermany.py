#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn, Kathrin Rack
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
@file geoModificationGermany.py

@brief Provides methods to return lists of local entities such as federal states
    and counties and geographical merging criteria.
"""

import os
from epidemiology.epidata import defaultDict as dd

import pandas as pd
from epidemiology.epidata import getDataIntoPandasDataFrame as gd

# Merging of Counties that are reported differently, either separatedly or 
# summed, in different data sources
CountyMerging = {
    # Different districts to Berlin; reporting differs according to source
    11000: [11001, 11002, 11003, 11004, 11005, 11006, 11007, 11008, 11009, 
            11010, 11011, 11012],
    # Wartburgkreis and Eisenach to Wartburgkreis (decision from July 1, 2021)
    16063 : [16063, 16056]
}


def get_state_ids(zfill=False):
    """"! Get list of federal state IDs sorted according to state ID.

    @param zfill [Default: False] Defines whether county IDs are zero-filled to
        two digits and returned as a string or returned as an integer.
    @return List of federal IDs sorted according to state ID.
    """         
    unique_geo_entities = sorted(set(dd.State.keys()))
    if zfill:
        unique_geo_entities = [str(id).zfill(2) for id in unique_geo_entities]

    return unique_geo_entities


def get_state_names():
    """"! Get list of federal state names sorted according to state ID.

    @return List of federal names sorted according to state ID.
    """            
    return [dd.State[i] for i in get_state_ids()]

def get_state_names_and_ids(zfill=False):
    """"! Get list of federal state names and IDs sorted according to state ID.

    @param zfill [Default: False] Defines whether county IDs are zero-filled to
        two digits and returned as a string or returned as an integer.
    @return List of federal names and IDs sorted according to state ID.
    """        
    stateids = get_state_ids(zfill=zfill)

    return [[dd.State[int(stateids[i])], stateids[i]] for i in range(len(stateids))]

def get_stateid_to_name():
    """"! Returns a hash map from federal state ID to state name.

    @return hash map from federal state ID to state name.
    """            
    return dd.State

# while reporting for Berlin is just different for different sources, Eisenach
# was merged on political decision with Wartburgkreis on July 1, 2021
def get_county_ids(merge_berlin=True, merge_eisenach=True, zfill=False):
    """"! Get list of county IDs sorted according to county ID.

    @param merge_berlin [Default: True] Defines whether the different districts
        are listed separately or combined as one entity 'Berlin'.
    @param merge_eisenach [Default: True] Defines whether the counties 
        'Wartburgkreis' and 'Eisenach' are listed separately or combined 
        as one entity 'Wartburgkreis'.
    @param zfill [Default: False] Defines whether county IDs are zero-filled to
        five digits and returned as a string or returned as an integer.
    @return List of county IDs sorted according to county ID.
    """        
    unique_geo_entities = sorted(set(dd.County.keys()))
    # disregard different districts of Berlin and only take Berlin as one county
    if merge_berlin:
        for i in CountyMerging[11000]:
            unique_geo_entities.remove(i)
    else:
        unique_geo_entities.remove(11000)
    
    if merge_eisenach:
        unique_geo_entities.remove(CountyMerging[16063][1])

    if zfill:
        return [str(id).zfill(5) for id in unique_geo_entities]
    else:
        return [id for id in unique_geo_entities]


def get_county_names(merge_berlin=True, merge_eisenach=True):
    """"! Get list of county names sorted according to county ID.

    @param merge_berlin [Default: True] Defines whether the different districts
        are listed separately or combined as one entity 'Berlin'.
    @param merge_eisenach [Default: True] Defines whether the counties 
        'Wartburgkreis' and 'Eisenach' are listed separately or combined 
        as one entity 'Wartburgkreis'.
    @return List of county names sorted according to county ID.
    """    
    return [dd.County[i] for i in get_county_ids(merge_berlin, merge_eisenach)]

def get_county_names_and_ids(merge_berlin=True, merge_eisenach=True, zfill=False):
    """"! Get list of county names and IDs sorted according to county ID.

    @param merge_berlin [Default: True] Defines whether the different districts
        are listed separately or combined as one entity 'Berlin'.
    @param merge_eisenach [Default: True] Defines whether the counties 
        'Wartburgkreis' and 'Eisenach' are listed separately or combined 
        as one entity 'Wartburgkreis'.
    @param zfill [Default: False] Defines whether county IDs are zero-filled to
        five digits and returned as a string or returned as an integer.        
    @return List of county names and IDs sorted according to county ID.
    """
    countyids = get_county_ids(merge_berlin, merge_eisenach)

    if zfill:
        return [[dd.County[countyids[i]], str(countyids[i]).zfill(5)] for i in range(len(countyids))]
    else:
        return [[dd.County[countyids[i]], countyids[i]] for i in range(len(countyids))]


def get_countyid_to_name():
    """"! Returns a hash map from county ID to county name.

    @return hash map from county ID to county name.
    """            
    return dd.County


def check_for_all_counties(unique_county_list, merge_berlin=True, merge_eisenach=True):
    """! Checks if all states are mentioned

    This function checks if all counties are available in the list provided.
    If data is incomplete this function returns false and a parent function may 
    try to download from another source.
    Note 1: There is no check if data for every day of every county is available.
    Note 2: If the source data file contains more local entities than the official 
                county list True is returned and the user has to check on its own.

    @param unique_county_list unique county list to check.
    @return Boolean to say if data is complete or not.
    """
    # check for all counties outside Berlin plus districts of Berlin
    missing = len(get_county_ids(merge_berlin, merge_eisenach))-len(unique_county_list)
    if missing != 0:
        print("Downloaded data is not complete. Missing " +
                str(missing) + " counties.")
        if missing < 0:
            # Returning True if source data file contains more counties than list
            print('Source data frame contains more counties than official '
                'county list. This could be OK, please verify yourself.')
            return True
        elif missing < 10:
            print('Missing counties: ' + str(list(set(get_county_ids(merge_berlin, merge_eisenach)
                                                      ).difference(unique_county_list).difference(set({11000})))))
        # Returning False if source data file lacks at least one county
        return False

    # if it seems complete
    return True

def get_countyid_to_stateid_map(merge_eisenach=True, zfill=False):
    """! Creates a hash map from county IDs to state IDs

    @param merge_eisenach [Default: True] Defines whether the counties 
        'Wartburgkreis' and 'Eisenach' are listed separately or combined 
        as one entity 'Wartburgkreis'.
    @param zfill [Default: False]. Defines whether or not all IDs are returned
        as zero-filled strings. By default, integer maps are returned.
    @return County ID to state ID map.
    """
    county_ids = get_county_ids(merge_eisenach=merge_eisenach, zfill=zfill)

    if zfill:
        return {id : id[0:2] for id in county_ids}
    else:
        return {id : int(str(id).zfill(5)[0:2]) for id in county_ids}


def get_stateid_to_countyids_map(merge_eisenach=True, zfill=False):
    """! Creates a hash map from state IDs to lists of county IDs

    @param merge_eisenach [Default: True] Defines whether the counties 
        'Wartburgkreis' and 'Eisenach' are listed separately or combined 
        as one entity 'Wartburgkreis'.
    @param zfill [Default: False]. Defines whether or not all IDs are returned
        as zero-filled strings. By default, integer maps are returned.
    @return State IDs to lists of county IDs map
    """
    county_ids = get_county_ids(merge_eisenach=merge_eisenach, zfill=zfill)
    state_ids = get_state_ids()
    state_to_county_table = [[] for i in range(len(state_ids))]

    for id in county_ids:
        state_to_county_table[int(str(id).zfill(5)[0:2])-1].append(id)

    return dict(zip(state_ids,state_to_county_table))


get_stateid_to_countyids_map()

def get_governing_regions(strict=True):
    """! Creates a sorted list of governing regions which may simply be 
    federal states or intermediate regions which themselves are a real 
    subset of a federal state and to which a certain number of counties 
    is attributed.

    Governing regions are generally denoted by the first three digits of the   
    belonging county IDs. In cases of a trailing zero, only two digits are
    taken and for Rhineland Palatinate and Saxony, the strict definition
    returns the two digit code of the federal state (i.e. 07 and 14).

    Note that this list may include former 'governing regions'. However, this 
    function may only be used to equally distribute information which only exist 
    on the 'governing region' level but not on county level itself or where the 
    county level information seems to be wrong. Then, information is extrapolated 
    with the help of governing regions.

    @param strict [Default: True] Defines whether only regions currently 
        considered as governing regions are returned.
    @return List of governing regions.
    """
    # Take first three digits, apply set() to remove double appearances and
    # sort again.
    if strict == False:
        return sorted(set([id[0:3] for id in get_county_ids(zfill=True)]))
    else:
        # make exceptions for Rhineland Palatinate and Saxony and remove 
        # trailing zeros
        return sorted(
            set(
                [id[0: 3]
                 if
                 not (id[0: 2] in ['07', '14'])
                 and (id[2] != '0') else id[0: 2]
                 for id in get_county_ids(zfill=True)]))


def get_official_county_table():
    """! Downloads a table with information on all German counties from destatis.

    @return County table with essential columns.
    """
    path_counties = 'https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales/Gemeindeverzeichnis/Administrativ/'
    county_table = gd.loadExcel(targetFileName='04-kreise.xlsx?__blob=publicationFile', apiUrl=path_counties, extension='',
                            param_dict={'sheet_name' : 1, 'header' : 5, 'engine' : 'openpyxl'})
    rename_kreise_deu_dict = {
        1: dd.EngEng['idCounty'],
        '2': "type", # name not important, column not used so far
        3: dd.EngEng['county'],
        4: dd.EngEng['nuts3'],
        5: dd.EngEng['area'],
        6: dd.EngEng['population'],
        7: "population_male", # name not important, column not used so far
        8: "population_female", # name not important, column not used so far
        9: "population_per_km2" # name not important, column not used so far
    }
    # rename columns
    county_table.rename(columns = rename_kreise_deu_dict, inplace = True)

    return county_table

def get_nuts3_county_id_map(merge_eisenach=True):
    """! Downloads county list file from destatis and creates hash map from 
    NUTS3 ID to county ID.

    @return Hash map of NUTS3 ID to county ID
    """
    # download county list
    county_table = get_official_county_table()
    # delete rows with nuts3 = NaN
    # take just columns with name dd.EngEng['idCounty'] and dd.EngEng['nuts3']
    key_nuts3 = county_table.dropna(subset=[dd.EngEng['nuts3']])[
        [dd.EngEng['idCounty'], dd.EngEng['nuts3']]]
    # convert ID data types
    key_nuts3 = key_nuts3.astype({dd.EngEng['idCounty'] : int})

    # check for completeness
    if not check_for_all_counties(key_nuts3[dd.EngEng['idCounty']].unique()):
        return dict()

    # make dictionary / hash map from data frame
    nuts3_key_dict = dict(
        zip(key_nuts3[dd.EngEng['nuts3']],
            key_nuts3[dd.EngEng['idCounty']]))

    return nuts3_key_dict



def merge_df_counties(df, merged_id, separated_ids, sorting=[dd.EngEng['date']], columns=dd.EngEng['date'], method='sum'):
    """! Merges the data frame data of different local entities such as the districts of Berlin or
    counties that were merged on political decision in between.

    @param df Original pandas dataframe.
    @param merged_id One new ID or old ID that is part of the list in 
        separated_ids that will be used in the returned data frame.
    @param separated_ids List of old IDs that will be merged.
    @param sorting Column or criterion on how to sort the rearranged frame.
    @param columns columns to be grouped by, Default: 'Date'.
    @param method Method of merging ('sum' [default], 'mean', 'median', 'min', 
        'max')
    @return Reduced data frame with separated_ids information merged to 
        merged_id rows.
    """
    # extract rows of IDs that will be merged
    rows_merged = df[dd.EngEng['idCounty']].isin(separated_ids)
    df_merged = df[rows_merged].copy()
    if not df_merged.empty:
        # set merged ID and county name
        if dd.EngEng['idCounty'] == columns:
            df_merged[dd.EngEng['idCounty']] = merged_id        
        df_merged = df_merged.groupby(columns).agg(method)
        # bring 'columns' which have been transfered to 'index' back as real 
        # columns 
        df_merged.reset_index(inplace=True)

        # reset countyID if not in columns
        if not dd.EngEng['idCounty'] in columns:
            df_merged[dd.EngEng['idCounty']] = merged_id
        # check if column available
        if dd.EngEng['county'] in df.columns:
            df_merged[dd.EngEng['county']] = dd.County[merged_id]

        # reset state ID and state name (only possible if idState and state
        # were already columns in the input data frame)
        if dd.EngEng['idState'] in df.columns:
            df_merged[dd.EngEng['idState']] = df[rows_merged][dd.EngEng['idState']].unique()[0]
        if dd.EngEng['state'] in df.columns:
            df_merged[dd.EngEng['state']] = df[rows_merged][dd.EngEng['state']].unique()[0]

        # remove unmerged rows from old data frame
        df = df[~rows_merged]
        df = pd.concat([df, df_merged], axis=0)

        # resort that final sorting is according to the date
        df.sort_values(sorting, inplace=True)
        df.reset_index(inplace=True, drop=True)

    return df

def merge_df_counties_all(df, sorting=[dd.EngEng['date']], columns=dd.EngEng['date'], method='sum'):
    """! Merges the data frame data of different local entities such as the districts
    of Berlin or counties that were merged on political decision in between according 
    to the lists provided in the dictionary geoModificationGermany.CountyMerging.

    @param df Original pandas dataframe.
    @param columns columns to be grouped by, Default: 'Date'.
    @param sorting Column or criterion on how to sort the rearranged frame. 
        Default: ['Date']   
    @param method Method of merging ('sum' [default], 'mean', 'median', 'min', 
        'max')
    @return Reduced data frame with IDs merged as given in CountyMerging dict.
    """
    for key, val in CountyMerging.items():
        df = merge_df_counties(df, key, val, sorting, columns, method)

    return df
