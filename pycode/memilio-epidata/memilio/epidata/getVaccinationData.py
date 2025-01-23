#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Martin J. Kuehn, Wadim Koslow
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
import itertools
import os
from datetime import datetime, date
from typing import Tuple

import numpy as np
import pandas as pd

from memilio.epidata import progress_indicator
from memilio.epidata import customPlot
from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getCommuterMobility as gcm
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import getPopulationData as gpd
from memilio.epidata import modifyDataframeSeries as mdfs

# activate CoW for more predictable behaviour of pandas DataFrames
pd.options.mode.copy_on_write = True


def download_vaccination_data(read_data, filename, directory, interactive):
    url = "https://raw.githubusercontent.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland/master/Deutschland_Landkreise_COVID-19-Impfungen.csv"
    path = os.path.join(directory + filename + ".json")
    df_data = gd.get_file(path, url, read_data, param_dict={'dtype': {
        'LandkreisId_Impfort': "string", 'Altersgruppe': "string", 'Impfschutz': int, 'Anzahl': int}},
        interactive=interactive)

    return df_data


def sanity_checks(df):
    # test if dataframe is empty
    if df.empty:
        raise gd.DataError(
            "Download of Vaccination Data failed. File is empty.")

    actual_strings_list = df.columns.tolist()
    # check number of data categories
    if len(actual_strings_list) != 5:
        raise gd.DataError("Error: Number of data categories changed.")

    # these strings need to be in the header
    test_strings = {"Impfdatum", "LandkreisId_Impfort",
                    "Altersgruppe", "Impfschutz", "Anzahl"}
    # check if headers are those we want
    for name in test_strings:
        if (name not in actual_strings_list):
            raise gd.DataError("Error: Data categories have changed.")


def compute_vaccination_ratios(
        age_group_list, vaccinations_table, vacc_column, region_column,
        population, merge_2022=True):
    """! Computes vaccination ratios based on the number of vaccinations
    and the corresponding population data

    @param age_group_list List of age groups considered.
    @param vaccinations_table Table of vaccinations (possible multiple columns
        for different number of doses)
    @param vacc_column Column name of vaccinations_table to be considered.
    @param region_column Column of regions in vaccinations table, e.g., ID_County or ID_State.
    @param population Table of population data for the given regions and considered age groups.
    @param merge_2022 [Default: False] Defines whether population data has to be merged to counties as of 2022.
    @return All vaccination ratios per region and age group.
    """
    # create new data frame and reshape it
    df_vacc_ratios = pd.DataFrame(columns=[region_column] + age_group_list)
    df_vacc_ratios[region_column] = vaccinations_table[region_column].unique()
    df_vacc_ratios[age_group_list] = np.array(
        vaccinations_table[vacc_column]).reshape(
        len(vaccinations_table[region_column].unique()),
        len(age_group_list))
    # compute county and age-group-specific vaccination ratios
    if merge_2022:
        population = geoger.merge_df_counties_all(
            population, sorting=[region_column],
            columns=region_column)
    df_vacc_ratios[['r' + age for age in age_group_list]
                   ] = df_vacc_ratios[age_group_list] / population[age_group_list].values

    return df_vacc_ratios


def sanitizing_average_regions(
        df, to_county_map, age_groups, column_names, age_population):
    """! Vaccinations in all regions are split up per population of its counties.
    This is done by summing up all vaccinations in this region and divide this by the population ratios.
    This is done for every age group and number of vaccination seperately.
    A new dataframme is created where the new data is stored.

    @param df DataFrame with Data to compute.
    @param to_county_map dict with regions as keys and countyIDs as values.
    @param age_groups list of all age groups as in df.
    @param column_names list of columns to compute.
    @param age_population Dataframe with number of population per age group and county.
    @return New DataFrame with sanitized data.
    """

    # create list of DataFrames, later to be merged
    df_total = []
    # computation is done for all Age Groups seperately
    for age in age_groups:
        # create subframe with specific age group
        df_age = df[df[dd.EngEng['ageRKI']] == age]
        # loop over all regions
        for counties_list in to_county_map.values():
            vacc_sums = df_age.loc[df_age[dd.EngEng['idCounty']].isin(
                counties_list)].groupby(dd.EngEng['date'])[column_names].sum()
            # get sums of all vaccinations in this region and agegroup
            vacc_sums = pd.concat(
                [vacc_sums] * len(counties_list),
                ignore_index=True)
            # get population ratios per county
            population_ratios = pd.merge(
                df_age.loc
                [df_age[dd.EngEng['idCounty']].isin(counties_list)]
                [dd.EngEng['idCounty']],
                age_population.loc
                [age_population[dd.EngEng['idCounty']].isin(counties_list)]
                [[dd.EngEng['idCounty'],
                  age]])[age] / age_population.loc[
                age_population[dd.EngEng['idCounty']].isin(counties_list)][
                age].sum()
            # for each column: vaccinations = all vaccinations * population_ratios
            for column in column_names:
                df_age.loc[df_age[dd.EngEng['idCounty']].isin(
                    counties_list), column] = vacc_sums[column].values * population_ratios.values

        df_total.append(df_age)

    df_total = pd.concat(df_total, ignore_index=True, sort=False)

    return df_total


def sanitizing_extrapolation_mobility(
        df, age_groups, column_names, age_population, neighbors_mobility):
    """! ATTENTION: DO NOT USE! ONLY FOR BACKWARD STABILITY AND DEVELOPMENT PURPOSES.
    Distributes vaccinations of a county to connected counties if a lot more vaccinations than the federal state average were reported at the newest date.
    Thus for different max dates data for a specific date can be different.
    The average vaccination ratio per age group is only computed for completed vaccinations.
    The average vaccination ratios for partially and refreshed vaccinations are not computed. Those vaccinations are also distributed by the ratios of completed vaccinations.
    Since the distribution is done for one county after another a different order of counties may result in different data.

    @param df DataFrame with Data to compute.
    @param age_groups list of all age groups as in df.
    @param column_names list of columns to compute.
    @param age_population Dataframe with number of population per age group and county.
    @param neighbors_mobility dict with counties as keys and commuter mobility to other counties as values.
    @return New DataFrame with sanitized data.
    """
    max_sanit_threshold_arr = np.zeros(len(age_groups))

    # compute average vaccination ratio per age group for full vaccinations
    aver_ratio = df.groupby(dd.EngEng['ageRKI']).agg({column_names[1]: "sum"})[
        column_names[1]].values / age_population[age_groups].sum().values

    # compute maximum_sanitizing threshold per age group as maxmimum of country-wide ratio + 10%
    # and predefined maximum value; threshold from becoming larger than 1
    for kk in range(len(age_groups)):
        max_sanit_threshold_arr[kk] = min(1, aver_ratio[kk] + 0.1)

    # create copy of dataframe
    df_san = df[:]

    # aggregate total number of vaccinations per county and age group
    vacc_sums_nonsanit = df.groupby(
        [dd.EngEng['idCounty'],
         dd.EngEng['ageRKI']]).agg(
        {column_names[1]: "sum"}).reset_index()
    # create new data frame and reshape it
    df_fullsum = compute_vaccination_ratios(
        age_groups, vacc_sums_nonsanit, column_names[1],
        dd.EngEng['idCounty'],
        age_population)

    # compute average federal state and age-group-specific vaccination ratios
    state_to_county = geoger.get_stateid_to_countyids_map(
        merge_eisenach=True)
    sanitizing_thresholds = []
    stateidx = 0
    for key in state_to_county.keys():
        vaccsum_state = df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
            state_to_county[key]), age_groups].sum()
        popsum_state = age_population.loc[age_population[dd.EngEng['idCounty']].isin(
            state_to_county[key]), age_groups].sum()
        sanitizing_thresholds.append(
            vaccsum_state.values / popsum_state.values)
        # set sanitizing threshold on federal state level as average + 5% or
        # as max_sanit_threshold (e.g., 95%) if average + 5% exceeds this
        # value; computation done separately for each age group
        for ii in range(len(sanitizing_thresholds[stateidx])):
            sanitizing_thresholds[stateidx][ii] = min(
                sanitizing_thresholds[stateidx][ii] + 0.05,
                max_sanit_threshold_arr[ii])

        # compute capacity weight for all counties of the federal state by
        # federal state and age-group-specific sanitizing threshold minus
        # current vaccination ratio
        df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
            state_to_county[key]), ['cw' + age for age in age_groups]] = sanitizing_thresholds[stateidx] - df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                state_to_county[key]), ['r' + age for age in age_groups]].values
        # replace negative numbers by zero, i.e., take maximum of 0 and value
        df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
            state_to_county[key]), ['cw' + age for age in age_groups]] = df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                state_to_county[key]), ['cw' + age for age in age_groups]].mask(
            df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                state_to_county[key]), ['cw' + age for age in age_groups]] < 0, 0)

        # compute equally the vaccination amount that is considered to be
        # distributed (or that can be accepted)
        df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
            state_to_county[key]), ['vd' + age for age in age_groups]] = df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                state_to_county[key]), [age for age in age_groups]].values - sanitizing_thresholds[stateidx] * age_population.loc[age_population[dd.EngEng['idCounty']].isin(
                    state_to_county[key]), age_groups].values

        stateidx += 1

    # go through all counties and age groups and distribute vaccinations
    # to connected counties (i.e., neighbors in mobility definition) if
    # a lot more than the federal state average were reported to be
    # vaccinated here
    for ii in range(len(df_fullsum)):
        # get county id
        id = df_fullsum.loc[ii, dd.EngEng['idCounty']]

        # get number of vaccinations to distribute of considered age group
        vacc_dist = df_fullsum.loc[ii,
                                   ['vd' + age
                                    for age in age_groups]]

        # get capacity weights for neihboring counties to distributed reported
        # vaccinations: access only rows which belong to neighboring
        # counties
        cap_weight = df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
            neighbors_mobility[id][0]), ['cw' + age for age in age_groups]]

        # iterate over age groups
        for ageidx in range(len(age_groups)):
            # check if vaccinations have to be distributed
            if vacc_dist[ageidx] > 1e-10:
                cap_chck = np.zeros(len(neighbors_mobility[id][0])) - 1
                chk_err_idx = 0
                while (chk_err_idx == 0) or (len(np.where(cap_chck > 1e-10)[0]) > 0):
                    neighb_cap_reached = np.where(cap_chck > -1e-10)[0]
                    neighb_open = np.where(cap_chck <= -1e-10)[0]
                    # maximum the neighbor takes before exceeding
                    # the average
                    vacc_nshare_max = df_fullsum.loc[
                        df_fullsum[dd.EngEng['idCounty']].isin(neighbors_mobility[id][0]), [
                            'vd' + age_groups[ageidx]]].reset_index().loc[:, 'vd' + age_groups[ageidx]].values
                    vacc_nshare_max[vacc_nshare_max > 0] = 0
                    # multiply capacity weight with commuter mobility weight and divide
                    # by sum of single products such that the sum of these weights
                    # is equal to 1. In here, adapt for cases where this would
                    # exceed the neighbor's capacity of vaccinations to take
                    # without exceeding the federal states-vaccination average+x%
                    # 1th step: initialize
                    vacc_dist_weight = np.zeros(
                        len(neighbors_mobility[id][0]))
                    # 2th step: adjust for maximum where necessary
                    vacc_dist_weight[neighb_cap_reached] = abs(
                        vacc_nshare_max[neighb_cap_reached] / vacc_dist[ageidx])
                    # 3th step: compute initial weights for all other counties
                    vacc_dist_weight[neighb_open] = neighbors_mobility[id][1][neighb_open] * cap_weight.values[
                        neighb_open,
                        ageidx] / sum(neighbors_mobility[id][1][neighb_open] * cap_weight.values[neighb_open, ageidx])
                    # 4th step: scale according to non-distributed vaccinations
                    vacc_dist_weight[neighb_open] = (
                        1 - sum(vacc_dist_weight[neighb_cap_reached])) * vacc_dist_weight[neighb_open]

                    # potential vaccination distribution after iteration step
                    vacc_nshare_pot = vacc_dist_weight * vacc_dist[ageidx]

                    # check capacity (remove zeros from denominator since there
                    # we would have 0/0 and the value there should be 0)
                    vacc_nshare_pot_denom = vacc_nshare_pot.copy()
                    vacc_nshare_pot_denom[vacc_nshare_pot_denom == 0] = 1
                    cap_chck = (
                        vacc_nshare_max + vacc_nshare_pot) / vacc_nshare_pot_denom

                    chk_err_idx += 1
                    if chk_err_idx > len(neighbors_mobility[id][0]):
                        raise gd.DataError(
                            'Error in functionality of vaccine distribution, exiting.')

                if abs(vacc_dist[ageidx] - sum(vacc_nshare_pot)) > 0.01 * vacc_dist[ageidx]:
                    raise gd.DataError(
                        'Error in functionality of vaccine distribution, exiting.')
                # sort weights from max to min and return indices
                order_dist = np.argsort(vacc_dist_weight)[::-1]
                vacc_nshare_pot = vacc_nshare_pot[order_dist]
                vacc_nshare_pot = vacc_nshare_pot[np.where(
                    vacc_nshare_pot > 0)]

                # get vaccination shares to distribute
                reduc_shares = vacc_nshare_pot / \
                    df_fullsum.loc[ii, age_groups[ageidx]]

                # sum up new additions
                df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(neighbors_mobility[id][0]), [
                    'vd' + age_groups[ageidx]]] = np.array(
                    df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(neighbors_mobility[id][0]), [
                        'vd' + age_groups[ageidx]]]).flatten() + vacc_dist[ageidx] * vacc_dist_weight

                # iterate over neighbors and add potential share of vaccionations
                # to neighbor and decrease local vaccionations at the end
                nidx = 0
                for neighb_id in neighbors_mobility[id][0][order_dist][
                        np.where(vacc_nshare_pot > 0)]:
                    df_san.loc[
                        (df_san[dd.EngEng['idCounty']] == neighb_id)
                        & (df_san[dd.EngEng['ageRKI']] == age_groups[ageidx]),
                        column_names] += (reduc_shares[nidx]
                                          * df_san.loc[(df_san[dd.EngEng['idCounty']] == id)
                                                       & (df_san[dd.EngEng['ageRKI']] == age_groups[ageidx]), column_names]
                                          .values
                                          )
                    nidx += 1
                df_san.loc[(df_san[dd.EngEng['idCounty']] == id)
                           & (df_san[dd.EngEng['ageRKI']] == age_groups[ageidx]),
                           column_names] -= (sum(reduc_shares)
                                             * df_san.loc[(df_san[dd.EngEng['idCounty']] == id)
                                                          & (df_san[dd.EngEng['ageRKI']] == age_groups[ageidx]), column_names]
                                             .values
                                             )

        if len(
                np.where(np.isnan(df_san[column_names]) == True)[0]) > 0:
            raise gd.DataError(
                'Error in functionality of vaccine distribution, NaN found after county '
                + str(id) + '. Exiting program.')

    # create copy only for possible comparison reasons now
    df = df_san[:]

    # create cumulative sum
    groupby_list = [dd.EngEng['date'], dd.EngEng['idState'],
                    dd.EngEng['idCounty'], dd.EngEng['ageRKI']]
    df = df.groupby(
        groupby_list).agg(
        {column_new: "sum" for column_new in column_names})
    df = df.groupby(
        level=[groupby_list.index(dd.EngEng['idCounty']),
               groupby_list.index(dd.EngEng['ageRKI'])]).cumsum().reset_index()

    # check that cumulative sum was done correctly (TODO: to be removed in production)
    for id in geoger.get_county_ids():
        for age in age_groups:
            a = df_san[(df_san.ID_County == id) & (
                df_san.Age_RKI == age)][column_names].sum()
            b = df[(df.ID_County == id) & (
                df.Age_RKI == age)].loc[:, column_names].iloc[-1]
            if sum(a - b) > 1e-8:
                gd.default_print(
                    "Error", "Cumulative sum error in: " + str(id) + " " + str(age))
    ### end of to be removed ###

    return df


def extrapolate_age_groups_vaccinations(
        df_data, population_all_ages, unique_age_groups_old,
        unique_age_groups_new, column_names, age_old_to_all_ages_indices,
        min_all_ages, all_ages_to_age_new_share):
    """! Original age groups (05-11, 12-17, 18-59, 60+) are replaced by infection data age groups
    (0-4, 5-14, 15-34, 35-59, 60-79, 80+). For every county the vacinations of old age groups are split to infection
    data age groups by its population ratio.
    For every age group and county a new dataframe is created. After the extrapolation all subframes are merged together.

    @param age_old_to_all_ages_indices List. List of original ages
    @param df_data DataFrame with Data to compute.
    @param population_all_ages Dataframe with number of population for every age group and county.
    @param unique_age_groups_old List of original age groups.
    @param unique_age_groups_new List of infection data age groups.
    @param column_names List of columns to compute.
    @param age_old_to_age_new_indices Defines in which new age group data from old age group is stored.
    @param min_all_ages List of lower age from all age groups
    @param all_ages_to_age_new_share Age groups indices of all age groups in every new age group
    @return New DataFrame with new age groups.
    """

    # create list of dataframes, later to be merged
    df_data_ageinf_county_cs = []
    groupby_list = list(df_data.columns)
    for i in range(len(column_names)):
        groupby_list.remove(column_names[i])
    # extrapolate age groups for one county at a time
    for countyID in df_data[dd.EngEng['idCounty']].unique():

        # get copy of global dataframe with only entries for current county
        vacc_df = df_data[df_data[dd.EngEng['idCounty']] == countyID]
        # get population data for all age groups in current county
        pop_state = population_all_ages[population_all_ages[dd.EngEng['idCounty']] == countyID]
        # create empty dataframe for each county
        total_county_df = []

        for i in range(0, len(unique_age_groups_old)):

            # get dataframe with only specific agegroup
            county_age_df = vacc_df[vacc_df[dd.EngEng['ageRKI']]
                                    == unique_age_groups_old[i]]
            # create new Dataframe with dates, ids, etc.
            info_df = county_age_df.drop(
                column_names, axis=1).drop(
                dd.EngEng['ageRKI'], axis=1)
            # create new dataframe for vaccination data
            vacc_data_df = []

            # get total population in old agegroup
            total_pop = 0
            for j in age_old_to_all_ages_indices[i]:
                total_pop += float(pop_state[str(min_all_ages[j])].iloc[0])
            # get population ratios in old agegroup
            ratios = [0 for zz in range(0, len(unique_age_groups_new))]
            for j in age_old_to_all_ages_indices[i]:
                ratios[all_ages_to_age_new_share[j][0][1]
                       ] += float(pop_state[str(min_all_ages[j])].iloc[0]) / total_pop
            # split vaccinations in old agegroup to new agegroups
            for j in range(0, len(ratios)):
                new_dataframe = county_age_df[column_names] * ratios[j]
                new_dataframe[dd.EngEng['ageRKI']] = unique_age_groups_new[j]
                vacc_data_df.append(pd.concat(
                    [info_df, new_dataframe],
                    axis=1))
            vacc_data_df = pd.concat(vacc_data_df)

            # merge all dataframes for each age group into one dataframe
            total_county_df.append(vacc_data_df)

        total_county_df = pd.concat(total_county_df).groupby(
            groupby_list).sum().reset_index()

        df_data_ageinf_county_cs.append(total_county_df)

        # test if number of vaccinations in current county are equal in old and new dataframe for random chosen date
        for vacc in column_names:
            if total_county_df[total_county_df[dd.EngEng['date']] == '2022-05-10'][vacc].sum() - vacc_df[vacc_df[dd.EngEng['date']] == '2022-05-10'][vacc].sum() > 1e-5:
                gd.default_print("Error", "Error in transformation...")

    # merge all county specific dataframes
    df_data_ageinf_county_cs = pd.concat(df_data_ageinf_county_cs)

    return df_data_ageinf_county_cs


# gets rki vaccination monitoring data for all states and extrapolates the values for counties according to their
# population Missing ratio values for the two different age groups are also estimated

def fetch_vaccination_data(
        conf_obj,
        filename: str,
        directory: str,
        read_data: str = dd.defaultDict['read_data'],
) -> pd.DataFrame:
    """ Downloads or reads the vaccination data and writes the RKIVaccFull dataset

    @param directory: str
        Path to the output directory
    @param filename: str
        Name of the full dataset filename
    @param conf_obj
        configuration object
    @param read_data bool True or False. Defines if data is read from file or downloaded. Default defined in defaultDict.

    @return pd.DataFrame fetched vaccination data
    """
    out_folder = conf_obj.path_to_use
    no_raw = conf_obj.no_raw

    df_data = download_vaccination_data(
        read_data, filename, directory, conf_obj.interactive)
    if conf_obj.checks:
        sanity_checks(df_data)

    if not no_raw:
        gd.write_dataframe(df_data, directory, filename, "json")
    return df_data


def process_vaccination_data(
        df_data: pd.DataFrame,
        conf_obj,
        directory: str,
        file_format: str = dd.defaultDict['file_format'],
        start_date: date = dd.defaultDict['start_date'],
        end_date: date = dd.defaultDict['end_date'],
        moving_average: int = dd.defaultDict['moving_average'],
        sanitize_data: int = dd.defaultDict['sanitize_data']
) -> dict:
    """! Processes downloaded raw data
    While working with the data
    - the column names are changed to English depending on defaultDict
    - The column "Date" provides information on the date of each data point given in the corresponding columns.

    @param df_data pd.DataFrame a Dataframe containing processed vaccination data
    @param directory: str
        Path to the output directory
    @param conf_obj
        configuration object
    @param file_format str. File format which is used for writing the data. Default defined in defaultDict.
    @param start_date Date of first date in dataframe. Default defined in defaultDict.
    @param end_date Date of last date in dataframe. Default defined in defaultDict.
    @param moving_average int. Integers >=0. Applies an 'moving_average'-days moving average on all time series
        to smooth out effects of irregular reporting. Default defined in defaultDict.
    @param sanitize_data int. Value in {0,1,2,3}; Default: 1. For many counties,
        vaccination data is not correctly attributed to home locations of
        vaccinated persons. If 'sanitize_data' is set to larger 0, this is corrected.
        0: No sanitizing applied.
        1: Averaging ratios over federal states.
        2: Averaging ratios over intermediate regions.
        3: All counties with vaccination quotas of
        more than 'sanitizing_threshold' will be adjusted to the average of its
        federal state and remaining vaccinations will be distributed to
        closely connected neighboring regions using commuter mobility networks.
        The sanitizing threshold will be defined by the age group-specific
        average on the corresponding vaccination ratios on county and federal
        state level.

    @return tuple and DataFrame
    """
    out_folder = conf_obj.path_to_use
    no_raw = conf_obj.no_raw

    with progress_indicator.Spinner(message='Preparing DataFrame'):
        df_data.rename(dd.GerEng, axis=1, inplace=True)

        try:
            df_data[dd.EngEng['date']] = pd.to_datetime(
                df_data[dd.EngEng['date']], format="ISO8601")
        except ValueError:
            try:
                df_data[dd.EngEng['date']] = pd.to_datetime(
                    df_data[dd.EngEng['date']], format="%Y-%m-%d")
            except:
                raise gd.DataError(
                    "Time data can't be transformed to intended format")

        # remove unknown locations if only modest number (i.e. less than 0.1%)
        if df_data[
                df_data[dd.EngEng['idCounty']] == 'u'].agg(
                {'Number': "sum"}).Number / df_data.agg(
                {'Number': "sum"}).Number < 0.001:
            df_data = df_data[df_data[dd.EngEng['idCounty']] != 'u']
        else:
            raise gd.DataError(
                'Too many data items with unknown vaccination location, '
                'please check source data.')

        if df_data[
                df_data[dd.EngEng['ageRKI']] == 'u'].agg(
                {'Number': "sum"}).Number / df_data.agg(
                {'Number': "sum"}).Number < 0.001:
            df_data = df_data[df_data[dd.EngEng['ageRKI']] != 'u']
        else:
            raise gd.DataError(
                'Too many data items with unknown vaccination age, '
                'please check source data.')

        # remove leading zeros for ID_County (if not yet done)
        try:
            df_data[dd.EngEng['idCounty']
                    ] = df_data[dd.EngEng['idCounty']].astype(int)
        except ValueError:
            gd.default_print("Error", 'Data items in ID_County could not be converted to integer. '
                                      'Imputation and/or moving_average computation will FAIL.')

        # NOTE: the RKI vaccination table contains about
        # 180k 'complete' vaccinations in id 17000 Bundesressorts, which
        # can not be attributed to any county, so these are currently ignored!
        # for spatially resolved data, we remove and ignore it.
        df_data = df_data[df_data[dd.EngEng['idCounty']] != 17000]
        # get unique age groups
        unique_age_groups_old = sorted(df_data[dd.EngEng['ageRKI']].unique())

        ##### get population data to sanitize data or to create new age groups #####
        # TODO: The following functionality may be generalized and outsourced to
        # getPopulationData....
        # reasonable max age; defines the extrapolation factor for the two
        # oldest age groups
        max_age_all = 100
        # get age groups separators of original vaccination table
        min_age_old = []
        extrapolate_agegroups = True
        for age in unique_age_groups_old:
            if '-' in age:
                min_age_old.append(int(age.split('-')[0]))
            elif '+' in age:
                min_age_old.append(int(age.split('+')[0]))
            else:
                extrapolate_agegroups = False
                gd.default_print(
                    "Error", "can not extrapolate provided age groups from vaccination data to infection number age groups.")
        min_age_old.append(max_age_all)

    # get population data for all countys (TODO: better to provide a corresponding method for the following lines in getPopulationData itself)

    try:
        population = pd.read_json(
            directory + "county_current_population.json")
    # pandas>1.5 raise FileNotFoundError instead of ValueError
    except (ValueError, FileNotFoundError):
        gd.default_print(
            "Info", "Population data was not found. Download it from the internet.")
        population = gpd.get_population_data(
            read_data=False, file_format=file_format, out_folder=out_folder,
            no_raw=no_raw, merge_eisenach=True)

    with progress_indicator.Spinner(message='Preparing Population data and age groups'):
        min_age_pop = []
        unique_age_groups_pop = list(population.columns)[2:]
        for age in unique_age_groups_pop:
            age = age.split()[0]  # remove " years" from string
            if '-' in age:
                min_age_pop.append(int(age.split('-')[0]))
            elif '>' in age:
                min_age_pop.append(int(age.split('>')[1]) + 1)
            elif '<' in age:
                min_age_pop.append(0)
            else:
                extrapolate_agegroups = False
                gd.default_print(
                    "Error", "can not extrapolate provided age groups from vaccination data to infection number age groups.")
        min_age_pop.append(max_age_all)

        # new age groups, here taken from definition of RKI infection data
        min_age_new = [0, 5, 15, 35, 60, 80, max_age_all]

        # combine all age group breaks
        min_all_ages = sorted(pd.unique(np.array(list(
            itertools.chain(min_age_old, min_age_pop, min_age_new)))))

        # check if the vaccinated age groups in old age groups start at zero
        if min_age_old[0] == 0:
            old_age_not_vacc = 0
            index_shift = []
        else:
            old_age_not_vacc = 1
            index_shift = [0]

        # compute share of all_ages intervals from population intervals
        population_to_all_ages_share = mdfs.create_intervals_mapping(
            min_age_pop, min_all_ages)

        # compute mappings from all ages to old and new intervals
        all_ages_to_age_old_share = mdfs.create_intervals_mapping(
            min_all_ages, index_shift + min_age_old)
        all_ages_to_age_new_share = mdfs.create_intervals_mapping(
            min_all_ages, min_age_new)

        # get interval indices from all age groups that correspond to old age group
        age_old_to_all_ages_indices = [[]
                                       for zz in range(0, len(min_age_old) - 1)]
        for i in range(0, len(unique_age_groups_old)):
            for k in range(0, len(all_ages_to_age_old_share)):
                if all_ages_to_age_old_share[k][0][1] == i + old_age_not_vacc:
                    age_old_to_all_ages_indices[i].append(k)
                elif k == len(all_ages_to_age_old_share) \
                        or all_ages_to_age_old_share[k][0][1] == i + old_age_not_vacc + 1:
                    break

        # get interval indices from all age groups that correspond to new age group
        age_new_to_all_ages_indices = [[]
                                       for zz in range(0, len(min_age_new) - 1)]
        for i in range(0, len(min_age_new)):
            for k in range(0, len(all_ages_to_age_new_share)):
                if all_ages_to_age_new_share[k][0][1] == i:
                    age_new_to_all_ages_indices[i].append(k)
                elif k == len(all_ages_to_age_new_share) \
                        or all_ages_to_age_new_share[k][0][1] == i + 1:
                    break

        # create new data frame and add zero to all new age group columns
        population_all_ages = pd.DataFrame(population[dd.EngEng['idCounty']])
        for i in min_all_ages:
            population_all_ages[str(i)] = 0

        # iterate over all original age groups
        for i in range(0, len(population_to_all_ages_share)):
            # iterate over intervals where population shares are assigned to
            for assign_share in population_to_all_ages_share[i]:
                # assign_share[0]: share / factor, assign_share[1]: column / age group
                population_all_ages[str(min_all_ages[assign_share[1]])
                                    ] += assign_share[0] * population[unique_age_groups_pop[i]]

        # rename last column and save total number per county
        population_all_ages.rename(
            columns={str(min_all_ages[-1]): 'Total'}, inplace=True)
        # remove last entry from all ages to call remaining columns
        min_all_ages = min_all_ages[:-1]
        population_all_ages['Total'] = population_all_ages[[
            str(i) for i in min_all_ages]].sum(axis=1)

        # TODO: a similar functionality has to be implemented as unit test
        if max(
                abs(
                    population[unique_age_groups_pop].sum(axis=1) -
                    population_all_ages[[str(i) for i in min_all_ages]].sum(
                        axis=1))) > 1e-8:
            gd.default_print("Error", "Population does not match expectations")

        population_old_ages = pd.DataFrame(population[dd.EngEng['idCounty']])
        for i in range(len(age_old_to_all_ages_indices)):
            # access columns + start_age_data since county_ID (and maybe other)
            # is in first place
            start_age_data = list(population_all_ages.columns).index('0')
            population_old_ages[unique_age_groups_old[i]] = population_all_ages.iloc[:, np.array(
                age_old_to_all_ages_indices[i]) + start_age_data].sum(axis=1)

        ## only for output meta information purposes ##
        # create hashes to access columns of new age group intervals
        unique_age_groups_new = [
            str(min_age_new[i]) + "-" + str(min_age_new[i + 1] - 1)
            for i in range(len(min_age_new) - 1)]
        population_new_ages = pd.DataFrame(population[dd.EngEng['idCounty']])
        for i in range(len(age_new_to_all_ages_indices)):
            # access columns + start_age_data since county_ID (and maybe other)
            # is in first place
            start_age_data = list(population_all_ages.columns).index('0')
            population_new_ages[unique_age_groups_new[i]] = population_all_ages.iloc[:, np.array(
                age_new_to_all_ages_indices[i]) + start_age_data].sum(axis=1)
        # end of output meta information purposes

        vacc_column_names = [
            dd.EngEng['vaccPartial'],
            dd.EngEng['vaccComplete'],
            dd.EngEng['vaccRefresh']]

        groupby_list = [
            dd.EngEng['date'],
            dd.EngEng['idCounty'],
            dd.EngEng['ageRKI']]

        if sanitize_data < 3:
            moving_average_sanit = moving_average
            impute_sanit = 'forward'
            compute_cumsum = True
        else:
            # do not sum vaccinations for sanitize_data = 3 since these have to
            # be redistributed in sanitizing data approach. Do not impute 'forward' used for
            # cumulative sums but impute 'zeros', i.e., zero vaccinations for
            # missing dates.
            compute_cumsum = False
            moving_average_sanit = 0
            impute_sanit = 'zeros'

        # define new column names
        column_names_dict = {
            1: dd.EngEng['vaccPartial'],
            2: dd.EngEng['vaccComplete'],
            3: dd.EngEng['vaccRefresh'],
            11: dd.EngEng['vaccNotComplete'],
            'additional identifiers': dd.EngEng['vaccRefresh']}

        vacc_column_names, df_data_joined = mdfs.split_column_based_on_values(
            df_data, "Impfschutz", "Number", groupby_list, column_names_dict,
            compute_cumsum)

        ######## data with age resolution as provided in original frame ########
        # transform and write data frame resolved per county and age (with age
        # classes as provided in vaccination tables: 05-11, 12-17, 18-59, 60+)

        df_data_agevacc_county_cs = mdfs.impute_and_reduce_df(
            df_data_joined,
            {dd.EngEng['idCounty']: df_data_joined[dd.EngEng['idCounty']].unique(),
             dd.EngEng['ageRKI']: unique_age_groups_old},
            vacc_column_names,
            impute=impute_sanit, moving_average=moving_average_sanit,
            min_date=start_date, max_date=end_date)

    df_data_agevacc_county_cs = geoger.merge_df_counties_all(
        df_data_agevacc_county_cs,
        sorting=[dd.EngEng["idCounty"],
                 dd.EngEng["ageRKI"],
                 dd.EngEng["date"]],
        columns=[dd.EngEng["ageRKI"],
                 dd.EngEng["date"],
                 ])

    # remove merged information and convert ID from float
    df_data_agevacc_county_cs = df_data_agevacc_county_cs.astype(
        {dd.EngEng['idCounty']: int})

    # insert county names
    df_data_agevacc_county_cs.insert(
        loc=2, column=dd.EngEng["county"],
        value=df_data_agevacc_county_cs[dd.EngEng["idCounty"]].replace(
            dd.County))

    # insert corresponding state IDs
    df_data_agevacc_county_cs.insert(
        loc=1, column=dd.EngEng["idState"],
        value=df_data_agevacc_county_cs[dd.EngEng["idCounty"]])
    county_to_state = geoger.get_countyid_to_stateid_map(merge_eisenach=True)
    for countyid in df_data_agevacc_county_cs[dd.EngEng["idState"]].unique():
        df_data_agevacc_county_cs.loc[df_data_agevacc_county_cs[dd.EngEng["idState"]]
                                      == countyid, dd.EngEng["idState"]] = county_to_state[countyid]

    with progress_indicator.Spinner(message='Sanitizing') as spinner:
        if sanitize_data == 1 or sanitize_data == 2:

            if sanitize_data == 1:
                gd.default_print(
                    'Info', 'Sanitizing activated: Using federal state average values.')
                to_county_map = geoger.get_stateid_to_countyids_map(
                    merge_eisenach=True)
            elif sanitize_data == 2:
                gd.default_print(
                    'Info', 'Sanitizing activated: Using intermediate region average values.')
                to_county_map = geoger.get_intermediateregionid_to_countyids_map(
                    merge_eisenach=True)

            df_data_agevacc_county_cs = sanitizing_average_regions(
                df_data_agevacc_county_cs, to_county_map,
                unique_age_groups_old, vacc_column_names, population_old_ages)

        elif sanitize_data == 3:
            gd.default_print('Info',
                             'Sanitizing activated: Using mobility-based vaccination redistribution approach.')
            # get neighbors based on mobility pattern and store
            # commuter inflow from other counties as first weight to distribute
            # vaccinations from vaccination county to extrapolated home counties
            neighbors_mobility = gcm.get_neighbors_mobility_all(
                direction='in', abs_tol=10, out_folder=out_folder)
            df_data_agevacc_county_cs = sanitizing_extrapolation_mobility(
                df_data_agevacc_county_cs, unique_age_groups_old,
                vacc_column_names, population_old_ages, neighbors_mobility)
            # compute the moving average
            df_data_agevacc_county_cs = mdfs.impute_and_reduce_df(
                df_data_agevacc_county_cs,
                {dd.EngEng['idCounty']: df_data_agevacc_county_cs[dd.EngEng['idCounty']].unique(),
                 dd.EngEng['ageRKI']: unique_age_groups_old},
                vacc_column_names,
                impute='forward', moving_average=moving_average,
                min_date=start_date, max_date=end_date)
        else:
            spinner.stop()
            gd.default_print('Info', 'Sanitizing deactivated.')

    dict_of_data = {
        "df_data_agevacc_county_cs": df_data_agevacc_county_cs,
        "vacc_column_names": vacc_column_names,
        "unique_age_groups_old": unique_age_groups_old,
        "population_old_ages": population_old_ages,
        "extrapolate_agegroups": extrapolate_agegroups,
        "population_all_ages": population_all_ages,
        "unique_age_groups_new": unique_age_groups_new,
        "age_old_to_all_ages_indices": age_old_to_all_ages_indices,
        "population_new_ages": population_new_ages,
        "all_ages_to_age_new_share": all_ages_to_age_new_share,
        "min_all_ages": min_all_ages,
    }
    return dict_of_data


def write_vaccination_data(dict_data: dict,
                           conf_obj,
                           directory: str,
                           file_format: str = dd.defaultDict['file_format'],
                           impute_dates: bool = True,
                           moving_average: int = dd.defaultDict['moving_average'],
                           ) -> None or Tuple:
    """! Writes the vaccination data
    The data is exported in three different ways:
        - all_county_vacc: Resolved per county by grouping all original age groups (05-11, 12-17, 18-59, 60+)
        - all_county_agevacc_vacc: Resolved per county and original age group (05-11, 12-17, 18-59, 60+)
        - all_county_ageinf_vacc: Resolved per county and infection data age group (0-4, 5-14, 15-34, 35-59, 60-79, 80+)
            - To do so getPopulationData is used and age group specific date from the original source
                is extrapolated on the new age groups on county level.

    - Missing dates are imputed for all data frames ('fillDates' is not optional but always executed).
    - A central moving average of N days is optional.

    - Start and end dates can be provided to define the length of the returned data frames.
    Parameters
    ----------
    @param dict_data dict. Contains various datasets or values
        - df_data_agevacc_county_cs: pd.DataFrame a Dataframe containing processed vaccination data
        - vacc_column_names
        - unique_age_groups_old
        - population_old_ages
        - extrapolate_agegroups
        - population_all_ages
        - unique_age_groups_new
        - age_old_to_all_ages_indices
        - min_all_ages
        - all_ages_to_age_new_share
        - population_new_ages

    @param directory: str
        Path to the output directory
    @param conf_obj
        configuration object
    @param file_format: str. File format which is used for writing the data. Default defined in defaultDict.
    @param impute_dates: bool. True or False. Defines if values for dates without new information are imputed. Default defined in defaultDict.
    @param moving_average: int. Integers >=0. Applies an 'moving_average'-days moving average on all time series to smooth out effects of irregular reporting. Default defined in defaultDict.
        sanitize_data: Value in {0,1,2,3}; Default: 1. For many counties, vaccination data is not correctly attributed to home locations of
        vaccinated persons. If 'sanitize_data' is set to larger 0, this is corrected.
        0: No sanitizing applied.
        1: Averaging ratios over federal states.
        2: Averaging ratios over intermediate regions.
        3: All counties with vaccination quotas of more than 'sanitizing_threshold' will be adjusted to the average of its
        federal state and remaining vaccinations will be distributed to closely connected neighboring regions using commuter mobility networks.
        The sanitizing threshold will be defined by the age group-specific average on the corresponding vaccination ratios on county and federal
        state level.
    @return: none
    """

    df_data_agevacc_county_cs = dict_data["df_data_agevacc_county_cs"]
    vacc_column_names = dict_data["vacc_column_names"]
    unique_age_groups_old = dict_data["unique_age_groups_old"]
    population_old_ages = dict_data["population_old_ages"]
    extrapolate_agegroups = dict_data["extrapolate_agegroups"]
    population_all_ages = dict_data["population_all_ages"]
    unique_age_groups_new = dict_data["unique_age_groups_new"]
    age_old_to_all_ages_indices = dict_data["age_old_to_all_ages_indices"]
    min_all_ages = dict_data["min_all_ages"]
    all_ages_to_age_new_share = dict_data["all_ages_to_age_new_share"]
    population_new_ages = dict_data["population_new_ages"]

    # data for all dates is automatically added
    if not impute_dates:
        gd.default_print(
            'Warning', 'Setting impute_dates = True as data for all dates is automatically added.')
        impute_dates = True

    if conf_obj.plot:
        # have a look extrapolated vaccination ratios (TODO: create plotting for production)
        # aggregate total number of vaccinations per county and age group
        latest_date = df_data_agevacc_county_cs[dd.EngEng["date"]][len(
            df_data_agevacc_county_cs.index) - 1].strftime("%Y-%m-%d")
        vacc_sums_nonsanit = df_data_agevacc_county_cs.loc[(
            df_data_agevacc_county_cs.Date == latest_date), ['ID_County', vacc_column_names[1]]]
        # create new data frame and reshape it
        df_fullsum = compute_vaccination_ratios(
            unique_age_groups_old, vacc_sums_nonsanit, vacc_column_names[1],
            dd.EngEng['idCounty'],
            population_old_ages, merge_2022=True)

    df_data_agevacc_state_cs = df_data_agevacc_county_cs.groupby(
        [dd.EngEng['date'],
         dd.EngEng['idState'],
         dd.EngEng['ageRKI']]).agg(
        {column: "sum" for column in vacc_column_names}).reset_index()

    # make plot of absolute numbers original age resolution
    if conf_obj.plot:
        # extract (dummy) date column to plt
        date_vals = df_data_agevacc_county_cs.loc[
            (df_data_agevacc_county_cs[dd.EngEng['ageRKI']] ==
             unique_age_groups_old[0]) &
            (df_data_agevacc_county_cs[dd.EngEng['idCounty']] ==
             geoger.get_county_ids()[0])][dd.EngEng['date']]

        # plot partial vaccination curves for different age groups
        yvals = [
            df_data_agevacc_county_cs.loc
            [df_data_agevacc_county_cs[dd.EngEng['ageRKI']] == age,
             [dd.EngEng['date'],
             dd.EngEng['vaccPartial']]].groupby(dd.EngEng['date']).sum()
            for age in unique_age_groups_old]
        customPlot.plot_multiple_series(
            date_vals, yvals, [age for age in unique_age_groups_old],
            title='Partial vaccination over different age groups',
            xlabel=dd.EngEng['date'],
            ylabel='Number', fig_name="Germany_PartialVacination_Absolute")

        # plot full vaccination curves for different age groups
        yvals = [
            df_data_agevacc_county_cs.loc
            [df_data_agevacc_county_cs[dd.EngEng['ageRKI']] == age,
             [dd.EngEng['date'],
             dd.EngEng['vaccComplete']]].groupby(dd.EngEng['date']).sum()
            for age in unique_age_groups_old]
        customPlot.plot_multiple_series(
            date_vals, yvals, [age for age in unique_age_groups_old],
            title='Full vaccination over different age groups',
            xlabel=dd.EngEng['date'],
            ylabel='Number', fig_name="Germany_FullVacination_Absolute")

    ######## data without age resolution ###########
    # write data frame resolved per county
    df_data_county_cs = df_data_agevacc_county_cs.groupby(
        [dd.EngEng['date'],
         dd.EngEng['idState'],
         dd.EngEng['idCounty']]).agg(
        {col_new: "sum" for col_new in vacc_column_names}).reset_index()

    df_data_state_cs = df_data_county_cs.groupby(
        [dd.EngEng['date'],
         dd.EngEng['idState']]).agg(
        {column: "sum" for column in vacc_column_names}).reset_index()

    ####### age resolved with extrapolation to other age groups #######
    # write data frame resolved per county and age (with age classes as
    # provided in RKI infection tables: 0-4, 5-14, 15-34, 35-59, 60-79, 80+)

    if extrapolate_agegroups:
        with progress_indicator.Spinner(message='Extrapolating age groups'):
            df_data_ageinf_county_cs = extrapolate_age_groups_vaccinations(
                df_data_agevacc_county_cs, population_all_ages,
                unique_age_groups_old, unique_age_groups_new, vacc_column_names,
                age_old_to_all_ages_indices, min_all_ages,
                all_ages_to_age_new_share)

    df_data_ageinf_county_cs.reset_index(drop=True, inplace=True)

    df_data_ageinf_state_cs = df_data_ageinf_county_cs.groupby(
        [dd.EngEng['date'],
         dd.EngEng['idState'],
         dd.EngEng['ageRKI']]).agg(
        {column: "sum" for column in vacc_column_names}).reset_index()

    # make plot of relative numbers of original and extrapolated age resolution
    if conf_obj.plot:
        # merge Eisenach...
        population_new_ages = geoger.merge_df_counties_all(
            population_new_ages, sorting=[dd.EngEng["idCounty"]],
            columns=dd.EngEng["idCounty"])

        # have a look extrapolated vaccination ratios (TODO: create plotting for production)
        # aggregate total number of vaccinations per county and age group
        latest_date = datetime.strptime(
            df_data_ageinf_county_cs[dd.EngEng["date"]]
            [len(df_data_ageinf_county_cs.index) - 1],
            "%Y-%m-%d")
        vacc_sums = df_data_ageinf_county_cs.loc[(df_data_ageinf_county_cs.Date == latest_date.strftime(
            "%Y-%m-%d")), [dd.EngEng['idCounty'], vacc_column_names[1]]]
        df_fullsum_county = compute_vaccination_ratios(
            unique_age_groups_new, vacc_sums, vacc_column_names[1],
            dd.EngEng['idCounty'],
            population_new_ages, merge_2022=True)

        # TODO: Plot information on vaccination ratios of df_fullsum_county

        # have a look extrapolated vaccination ratios (TODO: create plotting for production)
        # aggregate total number of vaccinations per county and age group
        latest_date = datetime.strptime(
            df_data_ageinf_state_cs[dd.EngEng["date"]]
            [len(df_data_ageinf_state_cs.index) - 1],
            "%Y-%m-%d")
        vacc_sums = df_data_ageinf_state_cs.loc[(df_data_ageinf_state_cs.Date == latest_date.strftime(
            "%Y-%m-%d")), [dd.EngEng['idState'], vacc_column_names[1]]]

        # get population per federal state
        population_new_ages.insert(
            loc=1, column=dd.EngEng["idState"],
            value=population_new_ages[dd.EngEng["idCounty"]])
        county_to_state = geoger.get_countyid_to_stateid_map(
            merge_eisenach=True)
        for countyid in population_new_ages[dd.EngEng["idState"]].unique():
            population_new_ages.loc[population_new_ages[dd.EngEng["idState"]]
                                    == countyid, dd.EngEng["idState"]] = county_to_state[countyid]
        population_new_ages_states = population_new_ages.groupby(dd.EngEng['idState']).agg(
            {age_group: "sum" for age_group in unique_age_groups_new}).reset_index()

        df_fullsum_state = compute_vaccination_ratios(
            unique_age_groups_new, vacc_sums, vacc_column_names[1],
            dd.EngEng['idState'],
            population_new_ages_states)

        # TODO: Plot information on vaccination ratios of df_fullsum_state

        # extract (dummy) date column to plt
        date_vals = df_data_agevacc_county_cs.loc[
            (df_data_agevacc_county_cs[dd.EngEng['ageRKI']] ==
             unique_age_groups_old[0]) &
            (df_data_agevacc_county_cs[dd.EngEng['idCounty']] ==
             geoger.get_county_ids()[0])][dd.EngEng['date']]

        # consider first vaccination for new age groups
        yvals = [
            df_data_ageinf_county_cs.loc
            [df_data_ageinf_county_cs[dd.EngEng['ageRKI']] == age,
             [dd.EngEng['date'],
             dd.EngEng['vaccPartial']]].groupby(dd.EngEng['date']).sum()
            for age in unique_age_groups_new]

        customPlot.plot_multiple_series(
            date_vals, yvals, [age for age in unique_age_groups_new],
            title='Partial vaccination over different age groups',
            xlabel=dd.EngEng['date'],
            ylabel='Number',
            fig_name="Germany_PartialVacination_AgeExtr_Absolute")

        # consider full vaccination for new age groups
        yvals = [
            df_data_ageinf_county_cs.loc
            [df_data_ageinf_county_cs[dd.EngEng['ageRKI']] == age,
             [dd.EngEng['date'],
             dd.EngEng['vaccComplete']]].groupby(dd.EngEng['date']).sum()
            for age in unique_age_groups_new]

        customPlot.plot_multiple_series(
            date_vals, yvals, [age for age in unique_age_groups_new],
            title='Full vaccination over different age groups',
            xlabel=dd.EngEng['date'],
            ylabel='Number',
            fig_name="Germany_FullVacination_AgeExtr_Absolute")

    if not conf_obj.to_dataset:
        # store data for all counties
        filename = 'vacc_county_agevacc'
        filename = gd.append_filename(filename, impute_dates, moving_average)
        gd.write_dataframe(df_data_agevacc_county_cs,
                           directory, filename, file_format)

        # store data for all federal states: group information on date, state and age level
        # (i.e., aggregate information of all counties per federal state)
        filename = 'vacc_states_agevacc'
        filename = gd.append_filename(filename, impute_dates, moving_average)
        gd.write_dataframe(df_data_agevacc_state_cs,
                           directory, filename, file_format)

        # store data for all counties
        filename = 'vacc_county'
        filename = gd.append_filename(filename, impute_dates, moving_average)
        gd.write_dataframe(df_data_county_cs, directory, filename, file_format)

        # store data for all federal states: group information on date, state and age level
        # (i.e., aggregate information of all counties per federal state)
        filename = 'vacc_states'
        filename = gd.append_filename(filename, impute_dates, moving_average)
        gd.write_dataframe(df_data_state_cs, directory, filename, file_format)

        ####### age resolved with extrapolation to other age groups #######
        # write data frame resolved per county and age (with age classes as
        # provided in RKI infection tables: 0-4, 5-14, 15-34, 35-59, 60-79, 80+)

        # store data for all counties
        filename = 'vacc_county_ageinf'
        filename = gd.append_filename(filename, impute_dates, moving_average)
        gd.write_dataframe(df_data_ageinf_county_cs,
                           directory, filename, file_format)

        # store data for all federal states: group information on date, state and age level
        # (i.e., aggregate information of all counties per federal state)
        filename = 'vacc_states_ageinf'
        filename = gd.append_filename(filename, impute_dates, moving_average)
        gd.write_dataframe(df_data_ageinf_state_cs,
                           directory, filename, file_format)
        return None
    else:
        return (df_data_agevacc_county_cs, df_data_agevacc_state_cs,
                df_data_county_cs, df_data_state_cs,
                df_data_ageinf_county_cs, df_data_ageinf_state_cs)


def get_vaccination_data(
        read_data: str = dd.defaultDict['read_data'],
        file_format: str = dd.defaultDict['file_format'],
        out_folder: str = dd.defaultDict['out_folder'],
        start_date: date = dd.defaultDict['start_date'],
        end_date: date = dd.defaultDict['end_date'],
        moving_average: int = dd.defaultDict['moving_average'],
        sanitize_data: int = dd.defaultDict['sanitize_data'],
        impute_dates: bool = True,
        **kwargs
):
    """! Downloads the RKI vaccination data and provides different kind of structured data.

    The data is read from the internet.
    The file is read in or stored at the folder "out_folder"/Germany/.
    To store and change the data we use pandas.

    While working with the data
    - the column names are changed to English depending on defaultDict
    - The column "Date" provides information on the date of each data point given in the corresponding columns.

    - The data is exported in three different ways:
        - all_county_vacc: Resolved per county by grouping all original age groups (05-11, 12-17, 18-59, 60+)
        - all_county_agevacc_vacc: Resolved per county and original age group (05-11, 12-17, 18-59, 60+)
        - all_county_ageinf_vacc: Resolved per county and infection data age group (0-4, 5-14, 15-34, 35-59, 60-79, 80+)
            - To do so getPopulationData is used and age group specific date from the original source
                is extrapolated on the new age groups on county level.

    - Missing dates are imputed for all data frames ('fillDates' is not optional but always executed).
    - A central moving average of N days is optional.

    - Start and end dates can be provided to define the length of the returned data frames.

    @param read_data [Currently not used] True or False. Defines if data is read from file or downloaded.
        Here Data is always downloaded from the internet.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    @param out_folder Folder where data is written to. Default defined in defaultDict.
    @param start_date Date of first date in dataframe. Default defined in defaultDict.
    @param end_date Date of last date in dataframe. Default defined in defaultDict.
    @param moving_average Integers >=0. Applies an 'moving_average'-days moving average on all time series
        to smooth out effects of irregular reporting. Default defined in defaultDict.
    @param sanitize_data: Value in {0,1,2,3}; Default: 1. For many counties,
        vaccination data is not correctly attributed to home locations of
        vaccinated persons. If 'sanitize_data' is set to larger 0, this is
        corrected.
        0: No sanitizing applied.
        1: Averaging ratios over federal states.
        2: Averaging ratios over intermediate regions.
        3: All counties with vaccination quotas of more than
        'sanitizing_threshold' will be adjusted to the average of its
        federal state and remaining vaccinations will be distributed to
        closely connected neighboring regions using commuter mobility networks.
        The sanitizing threshold will be defined by the age group-specific
        average on the corresponding vaccination ratios on county and federal
        state level.
    @param impute_dates bool True or False. Defines if values for dates without new information are imputed. Default defined in defaultDict.
    @param to_dataset bool True or False. Whether to return the dataframe as an object instead of json file.
        If True - returns objects with dataframes
        If False - write dataframes into files
        Default defined in defaultDict.

    @return None
    """
    conf = gd.Conf(out_folder, **kwargs)
    out_folder = conf.path_to_use

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    filename = "RKIVaccFull"
    raw_df = fetch_vaccination_data(
        conf_obj=conf,
        filename=filename,
        directory=directory,
        read_data=read_data,
    )
    process_dict_df = process_vaccination_data(
        df_data=raw_df,
        conf_obj=conf,
        directory=directory,
        start_date=start_date,
        end_date=end_date,
        file_format=file_format,
        moving_average=moving_average,
        sanitize_data=sanitize_data
    )
    silver_datasets = write_vaccination_data(dict_data=process_dict_df,
                                             conf_obj=conf,
                                             directory=directory,
                                             file_format=file_format,
                                             impute_dates=impute_dates,
                                             moving_average=moving_average,
                                             )
    if conf.to_dataset:
        return silver_datasets


def main():
    """! Main program entry."""

    arg_dict = gd.cli("vaccination")
    get_vaccination_data(**arg_dict)


if __name__ == "__main__":
    main()
