#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
from datetime import date
import itertools
import pandas as pd
import time
import numpy as np
import os

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from epidemiology.epidata import getPopulationData

# Downloads vaccine data from RKI


def download_vaccine_data():
    # RKI content from github
    url = 'https://raw.githubusercontent.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland/master/Aktuell_Deutschland_Landkreise_COVID-19-Impfungen.csv'
    # empty data frame to return if not read correctly
    df = pd.DataFrame()
    # try to read csv
    try:
        df = pd.read_csv(url, dtype={'LandkreisId_Impfort': "string",
                         'Altersgruppe': "string", 'Impfschutz': int, 'Anzahl': int})
    except:
        print("Error in reading csv. Returning empty data frame.")

    return df


def fill_df(df_old, group_by_items, avg_by, moving_average, min_date='', max_date=''):
    """! Fills missing dates of df and optionally calculates the the 7 day moving average of the data

    @param df_old old pandas dataframe
    @param group_by_items Column names for grouping by and items of particular group specification (e.g., for region: list of county IDs)
    @param avg_by List of columns for which the 7-day moving average should be computed (ICU, ICU_invasive)
    @param moving_average Flag which indicates whether to compute the moving average
    @param min_date [Default: ''] If set, minimum date to be set in new data frame for all items in group_by
    @param max_date [Default: ''] If set, maximum date to be set in new data frame for all items in group_by
    @return dataframe with imputed dates (and moving average if requested)
    """
    df_new = []

    # derive time from date
    try:
        df_old.Date = df_old.Date.dt.date
    except:
        df_old['Date'] = pd.to_datetime(df_old['Date'])
        df_old.Date = df_old.Date.dt.date

    # create empty copy of the df
    df_new = pd.DataFrame(columns=df_old.columns)

    # remove 'index' column if available
    try:
        df_new = df_new.drop(columns='index')
    except:
        pass

    # range of dates which should be filled
    if min_date == '':
        min_date = min(df_old.Date)
    if max_date == '':
        max_date = max(df_old.Date)
    idx = pd.date_range(min_date, max_date)

    # create list of all possible groupby columns combinations
    unique_ids = []
    for group_key in group_by_items.keys():
        unique_ids.append(group_by_items[group_key])
    unique_ids_comb = list(itertools.product(*unique_ids))
    # create list of keys/group_by column names
    group_by = list(group_by_items.keys())

    # loop over all regions/ages/gender
    for ids in unique_ids_comb:
        df_local = df_old.copy()
        counter = 0
        # filter df
        while counter < len(ids):
            df_local = df_local[df_local[group_by[counter]] == ids[counter]]
            counter += 1

        # create missing dates
        df_local.index = df_local.Date
        df_local_new = df_local.reindex(idx)
        df_local_new.Date = idx
        df_local_new.index = (range(len(idx)))

        if len(df_local) > 0:
            # create values for first date
            values = {}
            for column in df_local.columns:
                values[column] = df_local[column][0]
            for avg in avg_by:
                values[avg] = 0

            # fill values of missing dates based on last entry
            df_local_new.fillna(method='ffill', inplace=True)
            # fill value of the first date, if it doesn't exist yet
            # has to be conducted in second step to not impute 'value' at first missing value if start is present
            df_local_new.fillna(values, limit=1, inplace=True)
            # fill remaining values (between first date and first reported date of the df_local)
            df_local_new.fillna(method='ffill', inplace=True)

            # compute 7 day moving average
            if moving_average:
                for avg in avg_by:
                    # compute moving average in new column
                    df_local_new['MA' + avg] = df_local_new[avg].rolling(
                        window=7, min_periods=1, center=True).mean()
                    df_local_new['MA' + avg] = df_local_new['MA' +
                                                            avg].fillna(df_local_new[avg])
                    # overwrite daily values by moving averages
                    df_local_new[avg] = df_local_new['MA' + avg]
                    # remove helper column 'MA'+column_name
                    df_local_new.drop('MA' + avg, axis=1, inplace=True)

        else:
            # Decide whether to activate the warning or not. Currently, it happens that certain counties
            # have not had any kind of refreshing vaccinations. Then, the warning is misleading.
            # print('Warning: Tuple ' + str(ids) + ' not found in local data frame. Imputing zeros.')
            # create zero values for numbers of non-mentioned regions, age groups etc.
            values = {}
            counter = 0
            while counter < len(ids):
                values[group_by[counter]] = ids[counter]
                counter += 1
            for avg in avg_by:
                values[avg] = 0

            df_local_new.fillna(values, inplace=True)

        # append current local entity (i.e., county or state)
        df_new = df_new.append(df_local_new)
        # rearrange indices from 0 to N
        df_new.index = (range(len(df_new)))

    return df_new

# creates a mapping from given intervals to new desired intervals


def create_intervals_mapping(from_lower_bounds, to_lower_bounds):
    """! Creates a mapping from given intervals to new desired intervals

    @param from_lower_bounds lower bounds of original intervals
    @param to_lower_bounds desired lower bounds of new intervals
    @return mapping from intervals to intervals
    """
    # compute share of all_ages intervals from population intervals
    from_to_mapping = [[] for i in range(0, len(from_lower_bounds)-1)]
    j = 0  # iterator over all age breaks
    for i in range(0, len(from_lower_bounds)-1):
        # check if lower bound is larger in lower resolved data
        if from_lower_bounds[i] >= to_lower_bounds[j]:
            # Example: min_age_pop[i]=3 and min_age_pop[i+1]=6 shall be mapped on all_ages[j]=3 and all_ages[j+1]=5
            #   Then the all ages interval from j to j+1 will obtain the share
            #       x = (all_ages[j+1] - all_ages[j]) / (min_age_pop[i+1] - min_age_pop[i])
            #   of population age group 3-6
            #   in the next step, check if 6 is larger than all_ages[j+2]=Y
            #       if no: add the remaining part 1-x to j+1
            #       if yes: compute the corresponding share and go through it iteratively
            share = 0
            # if not, the remaining share will be assigned all ages j
            while from_lower_bounds[i+1] > to_lower_bounds[j+1]:
                share += (to_lower_bounds[j+1] - to_lower_bounds[j]) / \
                    (from_lower_bounds[i+1] - from_lower_bounds[i])
                from_to_mapping[i].append([share, j])
                j += 1
            from_to_mapping[i].append([1-share, j])
            # if both upper bounds are equal, then all ages j will not get any more share from any old group
            if from_lower_bounds[i+1] == to_lower_bounds[j+1]:
                j += 1

    return from_to_mapping

# splits a column based on its values to multiple columns


def split_column_based_on_values(df_global, column_ident, column_vals, new_column_labels):
    """! Fills missing dates of df and optionally calculates the the 7 day moving average of the data

    @param df_global global pandas dataframe
    @param column_ident identifier to split accordingly
    @param column_vals values to be put in new columns
    @param new_column_labels new labels for resulting columns
    @return dataframe with imputed dates (and moving average if requested)
    """
    # check number of given names and correct if necessary
    column_identifiers = df_global[column_ident].unique()
    if len(column_identifiers) != len(new_column_labels):
        print("Wrong number of new labels handed over. Naming new columns 'new_column_i' for i=1,2,...")
        new_column_labels = []
        for i in column_identifiers:
            new_column_labels.append('new_column_'+str(i))

    df_subset = []
    for i in range(0, len(column_identifiers)):
        df_subset.append(
            df_global[df_global[column_ident] == column_identifiers[i]])
        df_subset[i] = df_subset[i].drop(columns=column_ident)
        df_subset[i] = df_subset[i].rename(
            columns={column_vals: new_column_labels[i]})

    return df_subset


# gets rki vaccine monitoring data for all states and extrapolates the values for counties according to their population
# Missing ratio values for the two different age groups are also estimated
def get_vaccine_data(read_data=dd.defaultDict['read_data'],
                     file_format=dd.defaultDict['file_format'],
                     out_folder=dd.defaultDict['out_folder'],
                     no_raw=dd.defaultDict['no_raw'],
                     start_date=dd.defaultDict['start_date'],
                     end_date=dd.defaultDict['end_date'],
                     moving_average=dd.defaultDict['moving_average']
                     ):

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    df_data = download_vaccine_data()

    df_data.rename(dd.GerEng, axis=1, inplace=True)

    # remove unknown locations if only modest number
    if df_data[df_data['ID_County'] == 'u'].agg({'Number': sum}).Number < 10000:
        df_data = df_data[df_data['ID_County'] != 'u']
    else:
        print('Too many data items with unknown vaccination location, please check source data.')

    if df_data[df_data['Age_RKI'] == 'u'].agg({'Number': sum}).Number < 10000:
        df_data = df_data[df_data['Age_RKI'] != 'u']
    else:
        print('Too many data items with unknown vaccination age, please check source data.')

    # remove leading zeros for ID_County (if not yet done)
    try:
        df_data['ID_County'] = df_data['ID_County'].astype(int)
    except:
        print('Data items in ID_County could not be converted to integer. fill_data or moving_average computation will FAIL.')

    # get county ids # NOTE: the RKI vaccination table contains about 180k 'complete' vaccinations in id 17000 Bundesressorts,
    # which could not be attributed to any county, so these are currently ignored!
    unique_geo_entities = sorted(set(dd.County.keys()))
    # disregard different districts of Berlin and only take Berlin as one county
    for i in range(11001, 11013):
        unique_geo_entities.remove(i)

    unique_age_groups_old = sorted(df_data['Age_RKI'].unique())

    # df_data now becomes an array
    column_names_new = ["Vacc_partially", "Vacc_completed", "Vacc_refreshed"]
    df_data = split_column_based_on_values(
        df_data, "Impfschutz", "Number", column_names_new)

    # extract min, max dates and all possible values for groups from all data frames
    min_date = date.today()
    max_date = date(1, 1, 1)
    for i in range(0, len(df_data)):
        # pandas internally creates timestamps with to_datetime; to avoid warnings, convert to datetime again
        if pd.to_datetime(min(df_data[i].Date)).date() < min_date:
            min_date = pd.to_datetime(min(df_data[i].Date)).date()
        if pd.to_datetime(min(df_data[i].Date)).date() > max_date:
            max_date = pd.to_datetime(max(df_data[i].Date)).date()

    # write data frame resolved per county
    start_time = time.perf_counter()
    df_data_county_cs = []
    for i in range(0, len(df_data)):
        df_data[i] = df_data[i].reset_index()
        # group by date and county
        df_data_reduced = df_data[i].groupby(
            ['Date', 'ID_County']).agg({column_names_new[i]: sum})
        # compute cummulative sum over group index 1 (ID_County)
        df_data_reduced = df_data_reduced.groupby(
            level=1).cumsum().reset_index()
        df_data_county_cs.append(fill_df(df_data_reduced, {"ID_County": unique_geo_entities}, [
            column_names_new[i]], moving_average, min_date, max_date))
        # merge all data of vaccinations
        if i > 0:
            df_data_county_cs[0] = df_data_county_cs[0].merge(
                df_data_county_cs[i], on=['Date', 'ID_County'])

    # remove merged information and convert ID from float
    df_data_county_cs = df_data_county_cs[0].astype({'ID_County': int})

    filename = 'all_county_vacc'
    if moving_average:
        filename += '_ma'

    gd.write_dataframe(df_data_county_cs, directory, filename, file_format)
    end_time = time.perf_counter()
    print("Time needed: " + str(int(end_time - start_time)) + " sec")

    # write data frame resolved per county and age (with age classes as provided in vaccination tables: 12-17, 18-59, 60+)
    start_time = time.perf_counter()
    df_data_agevacc_county_cs = []
    for i in range(0, len(df_data)):
        df_data[i] = df_data[i].reset_index()
        # group by date and county
        df_data_reduced = df_data[i].groupby(
            ['Date', 'ID_County', 'Age_RKI']).agg({column_names_new[i]: sum})
        # compute cummulative sum over group index 1 (ID_County)
        df_data_reduced = df_data_reduced.groupby(
            level=[1, 2]).cumsum().reset_index()
        df_data_agevacc_county_cs.append(fill_df(df_data_reduced, {'ID_County': unique_geo_entities, 'Age_RKI': unique_age_groups_old}, [
            column_names_new[i]], moving_average, min_date, max_date))
        # merge all data of vaccinations
        if i > 0:
            df_data_agevacc_county_cs[0] = df_data_agevacc_county_cs[0].merge(
                df_data_agevacc_county_cs[i], on=['Date', 'ID_County', 'Age_RKI'])

    df_data_agevacc_county_cs = df_data_agevacc_county_cs[0].astype(
        {'ID_County': int})  # remove merged information and convert ID from float

    filename = 'all_county_agevacc_vacc'
    if moving_average:
        filename += '_ma'

    gd.write_dataframe(df_data_agevacc_county_cs,
                       directory, filename, file_format)
    end_time = time.perf_counter()
    print("Time needed: " + str(int(end_time - start_time)) + " sec")

    # write data frame resolved per county and age (with age classes as provided in RKI infection tables: 0-4, 5-14, 15-34, 35-59, 60-79, 80+)
    start_time = time.perf_counter()
    # reasonable max age; defines the extrapolation factor for the two oldest age groups
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
            print("Error in provided age groups from vaccination data; can not extrapolate to infection number age groups.")
    min_age_old.append(max_age_all)

    # get population data for all countys
    population = getPopulationData.get_age_population_data(write_df=False)
    min_age_pop = []
    extrapolate_agegroups = True
    unique_age_groups_pop = list(population.columns)[2:]
    for age in unique_age_groups_pop:
        age = age.split()[0]  # remove " years" from string
        if '-' in age:
            min_age_pop.append(int(age.split('-')[0]))
        elif '>' in age:
            min_age_pop.append(int(age.split('>')[1]))
        elif '<' in age:
            min_age_pop.append(0)
        else:
            extrapolate_agegroups = False
            print("Error in provided age groups from population data; can not extrapolate to infection number age groups.")
    min_age_pop.append(max_age_all)

    # new age groups, here taken from definition of RKI infection data
    min_age_new = [0, 5, 15, 35, 60, 80, max_age_all]

    # combine all age group breaks
    min_all_ages = sorted(pd.unique(list(itertools.chain(
        min_age_old, min_age_pop, min_age_new))))

    # get number of new age groups that are not vaccinated at all
    j = 0
    new_age_not_vacc = 0
    while min_age_new[j+1] < min_age_old[0]:
        new_age_not_vacc += 1
        j += 1

    # compute share of all_ages intervals from population intervals
    population_to_all_ages_share = create_intervals_mapping(
        min_age_pop, min_all_ages)

    # compute mappings from all ages to old and new intervals
    all_ages_to_age_old_share = create_intervals_mapping(
        min_all_ages, [0] + min_age_old)
    all_ages_to_age_new_share = create_intervals_mapping(
        min_all_ages, min_age_new)

    # compute indices of (partially) shared intervals from old to new
    age_old_to_age_new_share = create_intervals_mapping(
        [0] + min_age_old, min_age_new)
    age_old_to_age_new_share[0] = []
    for i in range(1, len(age_old_to_age_new_share)):
        age_old_to_age_new_share[i] = [x[1]
                                       for x in age_old_to_age_new_share[i]]

    # get interval indices from all age groups that correspond to old age group
    age_old_to_all_ages_indices = [[]
                                   for zz in range(0, len(min_age_old)-1)]
    for i in range(0, len(unique_age_groups_old)):
        for k in range(0, len(all_ages_to_age_old_share)):
            if all_ages_to_age_old_share[k][0][1] == i + new_age_not_vacc:
                age_old_to_all_ages_indices[i].append(k)
            elif k == len(all_ages_to_age_old_share) or all_ages_to_age_old_share[k][0][1] == i + new_age_not_vacc + 1:
                break

    # get interval indices from all age groups that correspond to new age group
    age_new_to_all_ages_indices = [[]
                                   for zz in range(0, len(min_age_new)-1)]
    for i in range(0, len(min_age_new)):
        for k in range(0, len(all_ages_to_age_new_share)):
            if all_ages_to_age_new_share[k][0][1] == i:
                age_new_to_all_ages_indices[i].append(k)
            elif k == len(all_ages_to_age_new_share) or all_ages_to_age_new_share[k][0][1] == i + 1:
                break

    # create new data frame and add zero to all new age group columns
    all_ages_populations = pd.DataFrame(population["ID_County"])
    for i in min_all_ages:
        all_ages_populations[str(i)] = 0

    # iterate over all original age groups
    for i in range(0, len(population_to_all_ages_share)):
        # iterate over intervals where population shares are assigned to
        for assign_share in population_to_all_ages_share[i]:
            # assign_share[0]: share / factor, assign_share[1]: column / age group
            all_ages_populations[str(
                min_all_ages[assign_share[1]])] += assign_share[0] * population[unique_age_groups_pop[i]]

    # rename last column and save total number per county
    all_ages_populations.rename(
        columns={str(min_all_ages[-1]): 'Total'}, inplace=True)
    # remove last entry from all ages to call remaining columns
    min_all_ages = min_all_ages[:-1]
    all_ages_populations['Total'] = all_ages_populations[[
        str(i) for i in min_all_ages]].sum(axis=1)

    # TODO: a similar function has to be implemented as unit test
    if max(abs(population[unique_age_groups_pop].sum(axis=1)-all_ages_populations[[str(i) for i in min_all_ages]].sum(axis=1))) > 1e-8:
        print("ERROR")

    err_chk_time = 0
    if extrapolate_agegroups:

        # create new, empty data frame
        df_data_ageinf_county_cs = pd.DataFrame(
            columns=df_data_agevacc_county_cs.columns)

        # loop over all counties
        for county in unique_geo_entities:
            # get population data of current state
            pop_state = all_ages_populations[all_ages_populations.ID_County == county]

            # access entries for the current county
            df_local_old = df_data_agevacc_county_cs[(
                df_data_agevacc_county_cs.ID_County == county)]
            # new data frame for the county
            df_local_new = pd.DataFrame(columns=['Date', 'ID_County'])

            # create dataframes for new age groups
            df_local_age_new = []
            for i in range(0, len(min_age_new)-1):
                # copy days and county ID (by using a dummy age group accessor)
                df_local_age_new.append(df_local_old.loc[df_data_agevacc_county_cs.Age_RKI ==
                                                         unique_age_groups_old[0], ["Date", "ID_County"]].copy())
                # add age group column
                df_local_age_new[i]["Age_RKI"] = str(
                    min_age_new[i])+"-"+str(min_age_new[i+1]-1)
                # add specific zero columns
                df_local_age_new[i][column_names_new] = 0
                # reset indices to 0,...,length to be able to sum up columns correctly
                df_local_age_new[i] = df_local_age_new[i].reset_index().drop(
                    columns='index')

            # iterate over old data frame and assign with correct shares to new age ranges
            j = 0
            for i in range(0, len(unique_age_groups_old)):

                # get data specific for one county and one age group of the original vaccination age group table
                df_local_age_old = df_local_old[df_local_old.Age_RKI ==
                                                unique_age_groups_old[i]]

                # get all age groups which are younger than the minimum age of vaccinated persons
                while min_age_new[j+1] < min_age_old[i]:
                    # append zero vaccination group to local, county-specific data frame
                    df_local_new = df_local_new.append(df_local_age_new[j])
                    j += 1
                df_local_new = df_local_new.reset_index().drop(columns='index')

                # copy shares of vaccinations from old age group columns to new age group columns
                for k in age_old_to_age_new_share[i + new_age_not_vacc]:

                    common_ages_all_ages = list(set(age_new_to_all_ages_indices[k]).intersection(
                        age_old_to_all_ages_indices[i]))
                    # splitting values of old column by comparing subset population sizes
                    # first, get population share that moves to new column
                    share_new = 0
                    total_old = 0
                    for age_group in common_ages_all_ages:
                        share_new += float(
                            pop_state[str(min_all_ages[age_group])])
                    for age_group in age_old_to_all_ages_indices[i]:
                        total_old += float(
                            pop_state[str(min_all_ages[age_group])])
                    # second add share to new column
                    dummy_frame = df_local_age_old[column_names_new].copy()
                    dummy_frame = dummy_frame.reset_index()
                    df_local_age_new[k][column_names_new] += share_new / \
                        total_old * dummy_frame[column_names_new]

            for k in range(j, len(min_age_new)-1):
                # append summed vaccination age group to local, county-specific data frame
                df_local_new = df_local_new.append(df_local_age_new[k])

            # TODO: unit test for some randomly selected counties and dates (here only stupid testing)
            start_time_err = time.perf_counter()
            if(abs(df_local_new.groupby("Date").agg({"Vacc_completed": sum}) - df_data_agevacc_county_cs[(
                    df_data_agevacc_county_cs.ID_County == county)].groupby("Date").agg({"Vacc_completed": sum})).sum()["Vacc_completed"]) > 1e-5:
                print("Error in transformation...")
            err_chk_time += time.perf_counter() - start_time_err

            df_data_ageinf_county_cs = df_data_ageinf_county_cs.append(
                df_local_new)

    filename = 'all_county_ageinf_vacc'
    if moving_average:
        filename += '_ma'

    gd.write_dataframe(df_data_ageinf_county_cs,
                       directory, filename, file_format)
    end_time = time.perf_counter()
    print("Time needed: " + str(int(end_time - start_time)) + " sec")
    print("Error check time " + str(err_chk_time))


def main():
    """! Main program entry."""

    arg_dict = gd.cli("vaccine")
    get_vaccine_data(**arg_dict)


if __name__ == "__main__":

    main()
