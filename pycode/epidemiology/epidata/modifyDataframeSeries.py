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
@brief Tools modify data frame series like imputing zeros for unknown dates, 
    copying previous values, and/or computing moving averages
"""
import pandas as pd
import itertools

from epidemiology.epidata import defaultDict as dd


def impute_and_reduce_df(df_old, group_by_cols, mod_cols, impute='forward', moving_average=0, min_date='', max_date='', start_w_firstval=False):
    """! Impute missing dates of dataframe time series and optionally calculates a moving average of the data

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
    # derive time from date
    try:
        df_old.Date = df_old.Date.dt.date
    except:
        df_old[dd.EngEng['date']] = pd.to_datetime(df_old[dd.EngEng['date']])
        df_old.Date = df_old.Date.dt.date

    # create empty copy of the df
    df_new = pd.DataFrame(columns=df_old.columns)
    # make pandas use the same data types....
    df_new = df_new.astype(dtype=dict(zip(df_old.columns, df_old.dtypes)))

    # remove 'index' column if available
    try:
        df_new = df_new.drop(columns='index')
    except:
        pass

    # range of dates which should be filled
    if min_date == '':
        min_date = min(df_old[dd.EngEng['date']])
    if max_date == '':
        max_date = max(df_old[dd.EngEng['date']])
    idx = pd.date_range(min_date, max_date)

    # create list of all possible groupby columns combinations
    unique_ids = []
    for group_key in group_by_cols.keys():
        unique_ids.append(group_by_cols[group_key])
    unique_ids_comb = list(itertools.product(*unique_ids))
    # create list of keys/group_by column names
    group_by = list(group_by_cols.keys())

    #loop over all regions/ages/gender
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
                df_local_new.fillna(values, inplace=True)
            else:
                # fill values of missing dates based on last entry
                df_local_new.fillna(method='ffill', inplace=True)
                # fill value of the first date, if it doesn't exist yet
                # has to be conducted in second step to not impute 'value' 
                # at first missing value if start is present
                df_local_new.fillna(values, limit=1, inplace=True)
                # fill remaining values (between first date and first 
                # reported date of the df_local)
                df_local_new.fillna(method='ffill', inplace=True)

            # compute 'moving average'-days moving average
            if moving_average > 0:
                for avg in mod_cols:
                    # compute moving average in new column
                    df_local_new['MA' + avg] = df_local_new[avg].rolling(
                        window=moving_average, min_periods=moving_average, center=True).mean().fillna(df_local_new[avg])
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
            df_local_new.fillna(values, inplace=True)

        # append current local entity (i.e., county or state)
        df_new = df_new.append(df_local_new)
        # rearrange indices from 0 to N
        df_new.index = (range(len(df_new)))

    return df_new
