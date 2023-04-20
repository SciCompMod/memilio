#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn
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
@file modifyDataframeSeries.py
@brief Tools for modifying data frame series like imputing zeros for unknown dates,
    copying previous values, and/or computing moving averages
"""
import itertools
from datetime import datetime, timedelta

import numpy as np
import pandas as pd

from memilio.epidata import defaultDict as dd


def impute_and_reduce_df(
        df_old, group_by_cols, mod_cols, impute='forward', moving_average=0,
        min_date='', max_date='', start_w_firstval=False):
    """! Impute missing dates of dataframe time series and optionally calculates a moving average of the data.
    Extracts Dates between min and max date.

    @param df_old old pandas dataframe
    @param group_by_cols Column names for grouping by and items of particular group specification (e.g., for region: list of county oder federal state IDs)
    @param mod_cols List of columns for which the imputation and/or moving average is conducted (e.g., Confirmed or ICU)
    @param impute [Default: 'forward'] imputes either based on older values ('forward') or zeros ('zeros')
    @param moving_average [Default: 0, no averaging] Number of days over which to compute the moving average
    @param min_date [Default: '', taken from df_old] If set, minimum date to be set in new data frame for all items in group_by
    @param max_date [Default: '', taken from df_old] If set, maximum date to be set in new data frame for all items in group_by
    @param start_w_firstval: [Default: False] If True and min_date < first date in dataframe, then between min_date and first date, the value
        of the first date will be repeated backwards. If False, then zero is set there.
    @return dataframe with imputed dates (and moving average if requested)
    """
    # derive date from time
    try:
        df_old.Date = df_old.Date.dt.date.astype(df_old.dtypes.Date)
    except AttributeError:
        df_old[dd.EngEng['date']] = pd.to_datetime(df_old[dd.EngEng['date']])
        df_old.Date = df_old.Date.dt.date.astype(df_old.dtypes.Date)

    # create empty copy of the df
    df_new = pd.DataFrame(columns=df_old.columns)
    # make pandas use the same data types....
    df_new = df_new.astype(dtype=dict(zip(df_old.columns, df_old.dtypes)))

    # remove 'index' column if available
    try:
        df_new = df_new.drop(columns='index')
    except KeyError:
        pass

    # range of dates which should be filled
    if min_date == '':
        min_date = min(df_old[dd.EngEng['date']])
    if max_date == '':
        max_date = max(df_old[dd.EngEng['date']])

    # Transform dates to datetime
    if isinstance(min_date, str) == True:
        min_date = datetime.strptime(min_date, "%Y-%m-%d")
    if isinstance(max_date, str) == True:
        max_date = datetime.strptime(max_date, "%Y-%m-%d")

    start_date = min_date
    end_date = max_date
    # shift start and end date for relevant dates to compute moving average.
    # if moving average is odd, both dates are shifted equaly.
    # if moving average is even, start date is shifted one day more than end date.
    if moving_average > 0:
        end_date = end_date + timedelta(int(np.floor((moving_average-1)/2)))
        start_date = start_date - timedelta(int(np.ceil((moving_average-1)/2)))

    idx = pd.date_range(start_date, end_date)

    # create list of all possible groupby columns combinations
    unique_ids = []
    for group_key in group_by_cols.keys():
        unique_ids.append(group_by_cols[group_key])
    unique_ids_comb = list(itertools.product(*unique_ids))
    # create list of keys/group_by column names
    group_by = list(group_by_cols.keys())

    # loop over all items in columns that are given to group by (i.e. regions/ages/gender)
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
            # depending on 'start_w_firstval', missing values at the beginning
            # of the frame will either be set to zero or to the first available
            # value in the data frame
            if not start_w_firstval:
                for avg in mod_cols:
                    values[avg] = 0

            if impute == 'zeros':
                # impute zeros at missing dates
                for keys, vals in values.items():
                    df_local_new[keys].fillna(vals, inplace=True)
            else:
                # fill values of missing dates based on last entry
                df_local_new.fillna(method='ffill', inplace=True)
                # fill value of the first date, if it doesn't exist yet
                # has to be conducted in second step to not impute 'value'
                # at first missing value if start is present
                for keys, vals in values.items():
                    df_local_new[keys].fillna(vals, limit=1, inplace=True)
                # fill remaining values (between first date and first
                # reported date of the df_local)
                df_local_new.fillna(method='ffill', inplace=True)

            # compute 'moving average'-days moving average
            if moving_average > 0:
                for avg in mod_cols:
                    # compute moving average in new column
                    df_local_new['MA' + avg] = df_local_new[avg].rolling(
                        window=moving_average, min_periods=moving_average,
                        center=True).mean().fillna(df_local_new[avg])
                    df_local_new['MA' + avg] = df_local_new['MA' +
                                                            avg].fillna(df_local_new[avg])
                    # overwrite daily values by moving averages
                    df_local_new[avg] = df_local_new['MA' + avg]
                    # remove helper column 'MA'+column_name
                    df_local_new.drop('MA' + avg, axis=1, inplace=True)

        else:
            # Decide whether to activate the following warning or not.
            # It can happen that certain local entities do not have certain data
            # at all (e.g., many counties do not have had any kind of
            # refreshing vaccinations so far.) Then, the following warning
            # is misleading.
            # print('Warning: Tuple ' + str(ids) + ' not found in local data frame. Imputing zeros.')
            # create zero values for non-existent time series
            values = {}
            counter = 0
            while counter < len(ids):
                values[group_by[counter]] = ids[counter]
                counter += 1
            for avg in mod_cols:
                values[avg] = 0
            # TODO: by this the corresponding columns will be zero-filled
            #       other entries such as names etc will get lost here
            #       any idea to introduce these names has to be found.
            for keys, vals in values.items():
                df_local_new[keys].fillna(vals, inplace=True)

        # append current local entity (i.e., county or state)
        df_new = pd.concat([df_new, df_local_new])
        # rearrange indices from 0 to N
        df_new.index = (range(len(df_new)))

    # extract min and max date
    df_new = extract_subframe_based_on_dates(df_new, min_date, max_date)

    return df_new


def split_column_based_on_values(
        df_to_split, column_to_split, column_vals_name, groupby_list,
        column_identifiers_to_names_dict, compute_cumsum):
    """! Splits a column in a dataframe into separate columns. For each unique value that appears in a selected column,
    all corresponding values in another column are transfered to a new column. If required, cumulative sum is calculated in new generated columns.

    @param df_to_split global pandas dataframe
    @param column_to_split identifier of the column for which separate values will define separate dataframes
    @param column_vals_name The name of the original column which will be split into separate columns named according to new_column_labels.
    @param groupby_list The name of the original columns with which data of new_column_labels can be joined.
    @param column_identifiers_to_names_dict Dict for new labels of resulting columns.
    @param compute_cumsum Computes cumulative sum in new generated columns
    @return a dataframe with the new splitted columns
    """
    column_identifiers = sorted(df_to_split[column_to_split].unique())
    new_column_labels = []
    # for column identifiers not in column_identifiers_to_names_dict.keys()
    # use the dict value with key 'additional identifiers' and add _2,_3,_4...
    key_to_start_count = [
        key for key, value in column_identifiers_to_names_dict.items()
        if value == column_identifiers_to_names_dict['additional identifiers']][0]

    for i in column_identifiers:
        if i in column_identifiers_to_names_dict.keys():
            new_column_labels.append(column_identifiers_to_names_dict[i])
        else:
            new_column_labels.append(
                column_identifiers_to_names_dict['additional identifiers'] + '_' +
                str(i - key_to_start_count + 1))

    # create empty copy of the df_to_split
    df_joined = pd.DataFrame(
        columns=df_to_split.columns).drop(
        columns=[column_to_split, column_vals_name])

    for i in range(0, len(column_identifiers)):
        df_reduced = df_to_split[df_to_split[column_to_split] == column_identifiers[i]].rename(
            columns={column_vals_name: new_column_labels[i]}).drop(columns=column_to_split)
        df_reduced = df_reduced.groupby(
            groupby_list).agg({new_column_labels[i]: sum})
        if compute_cumsum:
            # compute cummulative sum over level index of ID_County and level
            # index of Age_RKI
            df_reduced = df_reduced.groupby(
                level=[groupby_list.index(dd.EngEng['idCounty']),
                       groupby_list.index(dd.EngEng['ageRKI'])]).cumsum()
        # joins new generated column to dataframe
        df_joined = df_reduced.reset_index().join(
            df_joined.set_index(groupby_list),
            on=groupby_list, how='outer')

    return new_column_labels, df_joined


def extract_subframe_based_on_dates(df, start_date, end_date):
    """! Removes all data with date lower than start date or higher than end date.

    Returns the Dataframe with only dates between start date and end date.
    Resets the Index of the Dataframe.

    @param df The dataframe which has to be edited
    @param start_date Date of first date in dataframe
    @param end_date Date of last date in dataframe

    @return a dataframe with the extracted dates
    """

    upperdate = datetime.strftime(end_date, '%Y-%m-%d')
    lowerdate = datetime.strftime(start_date, '%Y-%m-%d')

    # Removes dates higher than end_date
    df_new = df[df[dd.EngEng['date']] <= upperdate]
    # Removes dates lower than start_date
    df_new = df_new[df_new[dd.EngEng['date']] >= lowerdate]

    df_new.reset_index(drop=True, inplace=True)

    return df_new


def insert_column_by_map(df, col_to_map, new_col_name, map):
    """! Adds a column to a given dataframe based on a mapping of values of a given column

    The mapping is defined by a list containing tupels of the form (new_value, old_value)
    where old_value is a value in the col_to_map and new_value the value
    that is added in the new column if col_to_map contains the old_value.
    @param df dataframe to modify
    @param col_to_map column containing values to be mapped
    @param new_col_name name of the new column containing the mapped values
    @param map List of tuples of values in the column to be added and values in the given column
    @return dataframe df with column of state names correspomding to state ids
    """
    df_new = df.copy()
    loc_new_col = df_new.columns.get_loc(col_to_map)+1
    df_new.insert(loc=loc_new_col, column=new_col_name,
                  value=df_new[col_to_map])
    for item in map:
        df_new.loc[df_new[col_to_map] == item[1], [new_col_name]] = item[0]
    return df_new


def create_intervals_mapping(from_lower_bounds, to_lower_bounds):
    """! Creates a mapping from given intervals to new desired intervals

    @param from_lower_bounds lower bounds of original intervals
    @param to_lower_bounds desired lower bounds of new intervals
    @return mapping from intervals to intervals
        The mapping is given as a list of tupels for every original interval.
        The list contains a tuple for every new interval intersecting the
        original interval. Each tuple defines the share of the original interval
        that is mapped to the new interval and the index of the new interval. We
        assume that the range of the intervals mapped from is contained in the
        range of the intervals mapped to.
        For example for from_lower_bounds = [5,20,30,80,85,90] and
        to_lower_bounds=[0,15,20,60,100] given the mapping would be
        [[[2/3,0], [1/3,1]],
         [[1,2]],
         [[3/5,2], [2/5,3]],
         [[1,3]],
         [[1,3]]]
    """
    if (from_lower_bounds[0] != to_lower_bounds[0] or
            from_lower_bounds[-1] != to_lower_bounds[-1]):
        print("Range of intervals mapped from is different than range of " +
              "intervals mapped to. Therefore, empty entries are possible.")

    extended_begin = False
    extended_end = False
    if (from_lower_bounds[0] < to_lower_bounds[0]):
        to_lower_bounds.insert(0, from_lower_bounds[0])
        extended_begin = True

    if (from_lower_bounds[-1] > to_lower_bounds[-1]):
        to_lower_bounds.append(from_lower_bounds[-1])
        extended_end = True

    # compute the shares of the original intervals mapped to the new intervals
    from_to_mapping = [[] for i in range(0, len(from_lower_bounds)-1)]
    to_idx = 0  # iterator over new intervals
    # iterate over original intervals
    for from_idx in range(0, len(from_lower_bounds) - 1):
        remaining_share = 1  # share of original interval to be distributed
        # position up to which the distribution is already computed
        pos = from_lower_bounds[from_idx]
        len_orig_interval = from_lower_bounds[from_idx+1] - pos
        # find first new interval intersecting the original interval
        while pos >= to_lower_bounds[to_idx+1]:
            to_idx += 1
        while from_lower_bounds[from_idx+1] > to_lower_bounds[to_idx+1]:
            # compute share of original interval that is send to new interval
            share = (to_lower_bounds[to_idx+1] - pos) / len_orig_interval
            if extended_begin:
                from_to_mapping[from_idx].append([share, to_idx-1])
            else:
                from_to_mapping[from_idx].append([share, to_idx])
            remaining_share -= share
            pos = to_lower_bounds[to_idx+1]
            to_idx += 1
        # if upper bound of the new interval is not smaller than upper bound of
        # the original interval assign remaining share of the original interval
        # to the new interval
        if (extended_begin and to_idx == 0) or (extended_end and to_idx == len(from_lower_bounds)):
            continue
        else:
            from_to_mapping[from_idx].append([remaining_share, to_idx])

    return from_to_mapping


def fit_age_group_intervals(
        df_age_in, age_out, df_population=None, max_age=100):
    """! Creates a mapping from given intervals to new desired intervals. Provide all intervals as "x-y".
         Boundary age groups can be provided with "<x" or ">y". Minimum and maximum are then taken as 0 and 99, respectively.
        Example:
        If df_population is set, we can use this data set to best interpolate
        @df_age_in to the desired age stratification of @age_out.
        Where this data is not finely enough resolved or if this data set is not
        provided, we assume the population to be equally distributed.
        Ex.
        df_age_in = ["1-10": 4, "11-60": 10, "61-99": 8]
        age_out = ["1-5", "6-10", "11-50", "51-99"]
        returns
        ["1-5": 2, "6-10": 2, "11-50": 8, "51-99": 10]
        if no population data is provided.

        If we also provide the population data
        population = ["1-5": 40, "6-7": 5, "8-10": 5, "11-60": 25, "61-99": 25],
        The output is:
        ["1-5": 3.2, "6-10": 0.8, "11-50": 8., "51:99": 10.]

    @param df_age_in Dataframe with columns of different age intervals and one row for subpopulation sizes for an arbitrary feature.
    @param age_out Desired age group distribution in list of strings. 
    @param df_population Total population data of the same structure as df_age_in used to inter- or extrapolate date of @df_age_in.
    @return Subpopulations of @df_age_in inter- or extrapolated to age stratification as required by @age_out.

    """
    # First check if the input is valid.

    # get minimum ages of each group
    age_in_min = []
    max_entry_in = 0
    min_entry_in = max_age
    for age in df_age_in.columns:
        if "year" in age:
            age = age.split()[0]  # remove " years" from string
        if '-' in age:
            age_in_min.append(int(age.split('-')[0]))
            min_entry_in = np.minimum(min_entry_in, int(age.split('-')[0]))
            max_entry_in = np.maximum(max_entry_in, int(age.split('-')[1]))
        elif '>' in age:
            age_in_min.append(int(age.split('>')[1])+1)
            max_entry_in = np.maximum(max_entry_in, max_age)
        elif '<' in age:
            age_in_min.append(0)
            min_entry_in = 0
            max_entry_in = np.maximum(max_entry_in, int(age.split('<')[1]))
        else:
            raise ValueError("Undefined entry for one age group in df_age_in")
    # if the max_entry is already an element of the vector, then we dont need to add it.
    # Otherwise we add the element to the end of the list and also set interset_max to true.
    inserted_max = False
    if max_entry_in not in age_in_min:
        age_in_min.append(max_entry_in)
        inserted_max = True

    age_out_min = []
    max_entry_out = 0
    min_entry_out = max_age
    for age in age_out:
        if "year" in age:
            age = age.split()[0]  # remove " years" from string
        if '-' in age:
            age_out_min.append(int(age.split('-')[0]))
            min_entry_out = np.minimum(min_entry_out, int(age.split('-')[0]))
            max_entry_out = np.maximum(max_entry_out, int(age.split('-')[1]))
        elif '>' in age:
            age_out_min.append(int(age.split('>')[1])+1)
            max_entry_out = np.maximum(max_entry_out, int(age.split('>')[1])+1)
        elif '<' in age:
            age_out_min.append(0)
            min_entry_out = 0
        else:
            raise ValueError("Undefined entry for one age group in age_out")
    if min_entry_out < min_entry_in or max_entry_out > max_entry_in:
        print(
            "Data from input data frame does not fit to desired output. Required data that is missing is interpreted as zero.")
    if max_entry_in not in age_out_min:
        age_out_min.append(max_entry_in)

    # if the first min age from age_out is greater than the first element in age_in,
    # we need to add the element and delete it later in the return.
    inserted_min = False
    if age_in_min[0] != age_out_min[0]:
        age_out_min.insert(0, age_in_min[0])
        inserted_min = True

    # when no df_population is given, we assume the data is equally distributed
    if df_population is None:
        population = df_age_in.iloc[0].to_numpy()
        age_shares = create_intervals_mapping(age_in_min, age_out_min)
        ans = np.zeros(len(age_out_min))
        population_indx = 0
        for age_entry_from in age_shares:
            for share_new in age_entry_from:
                ans[share_new[1]] += share_new[0] * population[population_indx]
            population_indx += 1
    else:
        age_pop_min = []
        max_entry_pop = 0
        min_entry_pop = max_age
        for age in df_population.columns:
            if "year" in age:
                age = age.split()[0]  # remove " years" from string
            if '-' in age:
                age_pop_min.append(int(age.split('-')[0]))
                min_entry_pop = np.minimum(
                    min_entry_pop, int(age.split('-')[0]))
                max_entry_pop = np.maximum(
                    max_entry_pop, int(age.split('-')[1]))
            elif '>' in age:
                age_pop_min.append(int(age.split('>')[1])+1)
                max_entry_pop = np.maximum(
                    max_entry_pop, int(age.split('>')[1])+1)
            elif '<' in age:
                age_pop_min.append(0)
                min_entry_pop = 0
            else:
                raise ValueError(
                    "Undefined entry for one age group in population data")
        if min_entry_out < min_entry_in or max_entry_out > max_entry_in:
            print(
                "Data from input data frame does not fit to population data. Required data that is missing is interpreted as zero.")
        if max_entry_in not in age_pop_min:
            age_pop_min.append(max_entry_in)

        # As already done with age_out, we need to check the min values for the population data again.
        if age_in_min[0] != age_pop_min[0]:
            age_pop_min.insert(0, age_in_min[0])

        # get weights from population file
        pop_data = df_population.iloc[0].to_numpy()
        new_pop = np.zeros(df_population.iloc[0].shape[0])
        age_shares_pop = create_intervals_mapping(age_in_min, age_pop_min)
        population_indx = 0
        for age_share in age_shares_pop:
            sum_pop = sum(df_population.iloc[0].to_numpy()[
                          age_share[0][1]:age_share[-1][1]+1])
            for age in age_share:
                new_pop[age[1]] += pop_data[age[1]
                                            ] / sum_pop * df_age_in.iloc[0][population_indx]

            population_indx += 1

        # population is now stored in new_pop for the population age groups
        # now, we need to distribute them to the aim age group
        age_shares = create_intervals_mapping(age_pop_min, age_out_min)
        ans = np.zeros(len(age_pop_min))
        population_indx = 0
        for age_entry_from in age_shares:
            for share_new in age_entry_from:
                ans[share_new[1]] += share_new[0] * new_pop[population_indx]
            population_indx += 1

    if inserted_min:
        ans = ans[1:]
    if inserted_max:
        ans = ans[:-1]

    return ans
