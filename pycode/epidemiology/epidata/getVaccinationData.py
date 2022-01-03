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
import sys
import time
import os
import itertools
import pandas as pd
import numpy as np

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from epidemiology.epidata import getPopulationData
from epidemiology.epidata import modifyDataframeSeries
from epidemiology.epidata import customPlot
from epidemiology.epidata import geoModificationGermany as geoger
from epidemiology.epidata import getCommuterMobility as gcm

# Downloads vaccination data from RKI


def download_vaccination_data():
    # RKI content from github
    url = 'https://raw.githubusercontent.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland/master/Aktuell_Deutschland_Landkreise_COVID-19-Impfungen.csv'
    # empty data frame to return if not read correctly
    df = pd.DataFrame()
    # try to read csv
    try:
        df = pd.read_csv(
            url,
            dtype={'LandkreisId_Impfort': "string", 'Altersgruppe': "string",
                   'Impfschutz': int, 'Anzahl': int})
    except:
        print("Error in reading csv. Returning empty data frame.")

    return df

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
                share += (to_lower_bounds[j+1] - to_lower_bounds[j]
                         ) / (from_lower_bounds[i+1] - from_lower_bounds[i])
                from_to_mapping[i].append([share, j])
                j += 1
            from_to_mapping[i].append([1-share, j])
            # if both upper bounds are equal, then all ages j will not get any more share from any old group
            if from_lower_bounds[i+1] == to_lower_bounds[j+1]:
                j += 1

    return from_to_mapping

# splits a column based on its values to multiple columns


def split_column_based_on_values(
        df_to_split, column_to_split, column_vals_name, new_column_labels):
    """! Splits a dataframe into seperate dataframes for each unique value in a column. 

    @param df_to_split global pandas dataframe
    @param column_to_split identifier of the column in which seperate values are split into seperate dataframes
    @param column_vals_name name of the column where the values of column_to_split are stored. This column will be renamed to new_column_labels
    @param new_column_labels new labels for resulting columns. There have to be the same amount of names and unique values.
    @return an array with the new splitted dataframes
    """
    # TODO: Maybe we should input a map e.g. 1->Vacc_partially, 2->Vacc_completed etc. This is more future proof.
    # check number of given names and correct if necessary
    column_identifiers = df_to_split[column_to_split].unique()
    if len(column_identifiers) != len(new_column_labels):
        print("Wrong number of new labels handed over. Naming new columns 'new_column_i' for i=1,2,...")
        new_column_labels = []
        for i in column_identifiers:
            new_column_labels.append('new_column_'+str(i))

    df_subset = []
    for i in range(0, len(column_identifiers)):
        df_subset.append(
            df_to_split[df_to_split[column_to_split] == column_identifiers[i]].copy())
        df_subset[i] = df_subset[i].drop(columns=column_to_split)
        df_subset[i] = df_subset[i].rename(
            columns={column_vals_name: new_column_labels[i]})

    return df_subset


# gets rki vaccination monitoring data for all states and extrapolates the values for counties according to their population
# Missing ratio values for the two different age groups are also estimated
def get_vaccination_data(read_data=dd.defaultDict['read_data'],
                         file_format=dd.defaultDict['file_format'],
                         out_folder=dd.defaultDict['out_folder'],
                         no_raw=dd.defaultDict['no_raw'],
                         start_date=dd.defaultDict['start_date'],
                         end_date=dd.defaultDict['end_date'],
                         make_plot=dd.defaultDict['make_plot'],
                         moving_average=dd.defaultDict['moving_average'],
                         sanitize_data=True,
                         max_sanit_threshold=0.95):
    """! Downloads the RKI vaccination data and provides different kind of structured data.

    The data is read from the internet.
    The file is read in or stored at the folder "out_folder"/Germany/.
    To store and change the data we use pandas.

    While working with the data
    - the column names are changed to English depending on defaultDict
    - The column "Date" provides information on the date of each data point given in the corresponding columns.

    - The data is exported in three different ways:
        - all_county_vacc: Resolved per county by grouping all original age groups (12-17, 18-59, 60+)
        - all_county_agevacc_vacc: Resolved per county and original age group (12-17, 18-59, 60+)
        - all_county_ageinf_vacc: Resolved per county and infection data age group (0-4, 5-14, 15-34, 35-59, 60-79, 80+)
            - To do so getPopulationData is used and age group specific date from the original source
                is extrapolated on the new age groups on county level.

    - Missing dates are imputed for all data frames ('fillDates' is not optional but always executed). 
    - A central moving average of N days is optional.

    - Start and end dates can be provided to define the length of the returned data frames.

    @param read_data False [Default]. Data is always downloaded from the internet.
    @param file_format File format which is used for writing the data. Default defined in defaultDict.
    @param out_folder Path to folder where data is written in folder out_folder/Germany.
    @param no_raw True or False [Default]. Defines if raw data is saved or not.
    @param start_date [Default = '', taken from read data] Start date of stored data frames.
    @param end_date [Default = '', taken from read data] End date of stored data frames.
    @param make_plot False [Default] or True. Defines if plots are generated with matplotlib.
    @param sanitize_data True [Default] or False. For certain federal states and counties, vaccination data
        is not correctly attributed to home locations of vaccinated persons. If 'sanitize_data' is set, 
        all counties with vaccination quotas of more than 'sanitizing_threshold' will be adjusted to the average of its
        federal state and remaining vaccinations will be distributed to closely connected neighboring regions 
        using commuter mobility networks. The sanitizing threshold will be defined by the age group-specific average
        on the corresponding federal state level + 5 % or max_sanit_threshold if this number exceeds max_sanit_threshold.
    @param max_sanit_threshold [Default: 0.95] Maximum relative number of vaccinations per county and age group.
    @param moving_average 0 [Default] or Number>0. Defines the number of days for which a centered moving average is computed.
    """
    start_time = time.perf_counter()
    # data for all dates is automatically added
    impute_dates = True

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    df_data = download_vaccination_data()

    df_old = df_data.copy()

    if not no_raw:
        gd.write_dataframe(df_data, directory, "RKIVaccFull", "json")

    df_data.rename(dd.GerEng, axis=1, inplace=True)

    # remove unknown locations if only modest number (i.e. less than 0.1%)
    if df_data[df_data[dd.EngEng['idCounty']] == 'u'].agg({'Number': sum}).Number \
        / df_data[df_data[dd.EngEng['idCounty']] != 'u'].agg({'Number': sum}).Number < 0.001:
        df_data = df_data[df_data[dd.EngEng['idCounty']] != 'u']
    else:
        sys.exit('Too many data items with unknown vaccination location, '
        'please check source data.')

    if df_data[df_data[dd.EngEng['ageRKI']] == 'u'].agg({'Number': sum}).Number \
         / df_data[df_data[dd.EngEng['ageRKI']] != 'u'].agg({'Number': sum}).Number < 0.001:
        df_data = df_data[df_data[dd.EngEng['ageRKI']] != 'u']
    else:
        sys.exit('Too many data items with unknown vaccination age, '
        'please check source data.')

    # remove leading zeros for ID_County (if not yet done)
    try:
        df_data[dd.EngEng['idCounty']
                ] = df_data[dd.EngEng['idCounty']].astype(int)
    except:
        print('Data items in ID_County could not be converted to integer. '
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
            print("Error in provided age groups from vaccination data; "
                  "can not extrapolate to infection number age groups.")
    min_age_old.append(max_age_all)

    # get population data for all countys (TODO: better to provide a corresponding method for the following lines in getPopulationData itself)
    try:
        population = pd.read_json(
            directory + "county_current_population_dim401.json")
    except:
        print("Population data was not found. Download it from the internet.")
        population = getPopulationData.get_age_population_data(read_data=False,
                            file_format = file_format,
                            out_folder = out_folder,
                            no_raw=no_raw,
                            write_df=True,
                            merge_eisenach=False)

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
            print("Error in provided age groups from population data;"
                  " can not extrapolate to infection number age groups.")
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
    age_old_to_all_ages_indices = [[] for zz in range(0, len(min_age_old)-1)]
    for i in range(0, len(unique_age_groups_old)):
        for k in range(0, len(all_ages_to_age_old_share)):
            if all_ages_to_age_old_share[k][0][1] == i + new_age_not_vacc:
                age_old_to_all_ages_indices[i].append(k)
            elif k == len(all_ages_to_age_old_share) \
                    or all_ages_to_age_old_share[k][0][1] == i + new_age_not_vacc + 1:
                break

    # get interval indices from all age groups that correspond to new age group
    age_new_to_all_ages_indices = [[] for zz in range(0, len(min_age_new)-1)]
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
        print("ERROR")

    population_old_ages = pd.DataFrame(population[dd.EngEng['idCounty']])
    for i in range(len(age_old_to_all_ages_indices)):
        # access columns + start_age_data since county_ID (and maybe other)
        # is in first place
        start_age_data = list(population_all_ages.columns).index('0')
        population_old_ages[unique_age_groups_old[i]] = population_all_ages.iloc[:, np.array(
            age_old_to_all_ages_indices[i])+start_age_data].sum(axis=1)
    ############## end of potential outsourcing #####################

    # df_data now becomes an array
    column_names_new = [
        dd.EngEng['vaccPartial'],
        dd.EngEng['vaccComplete'],
        dd.EngEng['vaccRefresh']]
    df_data_array = split_column_based_on_values(
        df_data, "Impfschutz", "Number", column_names_new)

    # extract min, max dates and all possible values for groups
    # from all data frames
    min_date = date.today()
    max_date = date(1, 1, 1)
    for i in range(0, len(df_data_array)):
        # pandas internally creates timestamps with to_datetime; to
        # avoid warnings, convert to datetime again
        if pd.to_datetime(min(df_data_array[i].Date)).date() < min_date:
            min_date = pd.to_datetime(min(df_data_array[i].Date)).date()
        if pd.to_datetime(min(df_data_array[i].Date)).date() > max_date:
            max_date = pd.to_datetime(max(df_data_array[i].Date)).date()

    ######## data with age resolution as provided in original frame ########
    # transform and write data frame resolved per county and age (with age
    # classes as provided in vaccination tables: 12-17, 18-59, 60+)
    df_data_agevacc_county_cs = []
    for i in range(0, len(df_data_array)):
        df_data_array[i] = df_data_array[i].reset_index()
        # group by date and county
        df_data_reduced = df_data_array[i].groupby(
            [dd.EngEng['date'],
             dd.EngEng['idCounty'],
             dd.EngEng['ageRKI']]).agg(
            {column_names_new[i]: sum})
        if not sanitize_data:
            # compute cummulative sum over level index 1 (ID_County) and level
            # index 2 (Age_RKI)
            df_data_reduced = df_data_reduced.groupby(
                level=[1, 2]).cumsum().reset_index()
            df_data_agevacc_county_cs.append(
                modifyDataframeSeries.impute_and_reduce_df(
                    df_data_reduced,
                    {dd.EngEng['idCounty']: df_data_reduced[dd.EngEng['idCounty']].unique(),
                     dd.EngEng['ageRKI']: unique_age_groups_old},
                    [column_names_new[i]],
                    impute='forward', moving_average=moving_average,
                    min_date=min_date, max_date=max_date))
        else:
            # do not sum vaccinations since these have to be redistributed in
            # sanitizing data approach. Also do not imput 'forward' used for
            # cumulative sums but impute 'zeros', i.e., zero vaccinations for
            # missing dates.
            df_data_reduced = df_data_reduced.reset_index()
            df_data_agevacc_county_cs.append(
                modifyDataframeSeries.impute_and_reduce_df(
                    df_data_reduced,
                    {dd.EngEng['idCounty']: df_data_reduced[dd.EngEng['idCounty']].unique(),
                     dd.EngEng['ageRKI']: unique_age_groups_old},
                    [column_names_new[i]],
                    impute='zeros', moving_average=0,
                    min_date=min_date, max_date=max_date))
        # merge all data of vaccinations
        if i > 0:
            df_data_agevacc_county_cs[0] = df_data_agevacc_county_cs[0].merge(
                df_data_agevacc_county_cs[i],
                on=[dd.EngEng['date'],
                    dd.EngEng['idCounty'],
                    dd.EngEng['ageRKI']])

    # remove merged information and convert ID from float
    df_data_agevacc_county_cs = df_data_agevacc_county_cs[0].astype(
        {dd.EngEng['idCounty']: int})

    # insert county names
    df_data_agevacc_county_cs.insert(
        loc=2, column=dd.EngEng["county"],
        value=df_data_agevacc_county_cs[dd.EngEng["idCounty"]].replace(
            dd.County))

    end_time = time.perf_counter()
    print("Time needed initial phase: " + str(int(end_time - start_time)) + " sec")

    if sanitize_data:
        start_time = time.perf_counter()
        # create copy of dataframe
        df_san = df_data_agevacc_county_cs.copy()

        # aggregate total number of vaccinations per county and age group
        vacc_sums_nonsanit = df_data_agevacc_county_cs.groupby(
            [dd.EngEng['idCounty'],
             dd.EngEng['ageRKI']]).agg(
            {column_names_new[1]: sum}).reset_index()
        # create new data frame and reshape it
        df_fullsum = pd.DataFrame(
            columns=[dd.EngEng['idCounty']] + unique_age_groups_old)
        df_fullsum[dd.EngEng['idCounty']
                   ] = vacc_sums_nonsanit[dd.EngEng['idCounty']].unique()
        df_fullsum[unique_age_groups_old] = np.array(
            vacc_sums_nonsanit[column_names_new[1]]).reshape(
            len(vacc_sums_nonsanit[dd.EngEng['idCounty']].unique()),
            len(unique_age_groups_old))
        # compute county and age-group-specific vaccination ratios
        df_fullsum[['r'+age for age in unique_age_groups_old]
                   ] = df_fullsum[unique_age_groups_old] / population_old_ages[unique_age_groups_old].values

        # compute average federal state and age-group-specific vaccination ratios
        state_to_county = geoger.get_stateid_to_countyids_map(
            merge_eisenach=False)
        sanitizing_thresholds = []
        stateidx = 0
        for key in state_to_county.keys():
            vaccsum_state = df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                state_to_county[key]), unique_age_groups_old].sum()
            popsum_state = population_old_ages.loc[population_old_ages[dd.EngEng['idCounty']].isin(
                state_to_county[key]), unique_age_groups_old].sum()
            sanitizing_thresholds.append(
                vaccsum_state.values / popsum_state.values)
            # set sanitizing threshold on federal state level as average + 5% or
            # as max_sanit_threshold (e.g., 95%) if average + 5% exceeds this
            # value
            for ii in range(len(sanitizing_thresholds[stateidx])):
                sanitizing_thresholds[stateidx][ii] = min(
                    sanitizing_thresholds[stateidx][ii] + 0.05,
                    max_sanit_threshold)

            # compute capacity weight for all counties of the federal state by
            # federal state and age-group-specific sanitizing threshold minus
            # current vaccination ratio
            df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                state_to_county[key]), ['cw'+age for age in unique_age_groups_old]] = sanitizing_thresholds[stateidx] - df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                    state_to_county[key]), ['r'+age for age in unique_age_groups_old]].values
            # replace negative numbers by zero, i.e., take maximum of 0 and value
            df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                state_to_county[key]), ['cw'+age for age in unique_age_groups_old]] = df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                    state_to_county[key]), ['cw'+age for age in unique_age_groups_old]].mask(df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                        state_to_county[key]), ['cw'+age for age in unique_age_groups_old]] < 0, 0)

            # compute equally the vaccination amount that is considered to be
            # distributed (or that can be accepted)
            df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                state_to_county[key]), ['vd'+age for age in unique_age_groups_old]] = df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                    state_to_county[key]), [age for age in unique_age_groups_old]].values - sanitizing_thresholds[stateidx] * population_old_ages.loc[population_old_ages[dd.EngEng['idCounty']].isin(
                        state_to_county[key]), unique_age_groups_old].values

            stateidx += 1

        # get neighbors based on mobility pattern and store
        # commuter inflow from other counties as first weight to distribute
        # vaccinations from vaccination county to extrapolated home counties 
        neighbors_mobility = gcm.get_neighbors_mobility_all(
            direction='in', abs_tol=10, merge_eisenach=False, 
            out_folder = out_folder)

        end_time = time.perf_counter()
        print("Time needed for preparing sanitizing: " +
              str(int(end_time - start_time)) + " sec")

        # go through all counties and age groups and distribute vaccinations
        # to connected counties (i.e., neighbors in mobility definition) if
        # a lot more than the federal state average were reported to be
        # vaccinated here
        for ii in range(len(df_fullsum)):
            start_time_county = time.perf_counter()
            # get county id
            id = df_fullsum.loc[ii, dd.EngEng['idCounty']]

            # get number of vaccinations to distribute of considered age group
            vacc_dist = df_fullsum.loc[ii,
                                       ['vd' + age
                                        for age in unique_age_groups_old]]

            # get capacity weights for neihboring counties to distributed reported
            # vaccinations: access only rows which belong to neighboring
            # counties
            cap_weight = df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(
                neighbors_mobility[id][0]), ['cw'+age for age in unique_age_groups_old]]

            # iterate over age groups
            for ageidx in range(len(unique_age_groups_old)):
                # check if vaccinations have to be distributed
                if vacc_dist[ageidx] > 1e-10:
                    cap_chck = np.zeros(len(neighbors_mobility[id][0]))-1
                    chk_err_idx = 0
                    while (chk_err_idx == 0) or (len(np.where(cap_chck > 1e-10)[0]) > 0):
                        neighb_cap_reached = np.where(cap_chck > -1e-10)[0]
                        neighb_open = np.where(cap_chck < -1e-10)[0]
                        # maximum the neighbor takes before exceeding
                        # the average
                        vacc_nshare_max = df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(neighbors_mobility[id][0]), [
                            'vd'+unique_age_groups_old[ageidx]]].reset_index().loc[:, 'vd'+unique_age_groups_old[ageidx]].values
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
                        vacc_dist_weight[neighb_open] = neighbors_mobility[id][1][neighb_open] * cap_weight.values[neighb_open,
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
                            vacc_nshare_max + vacc_nshare_pot)/vacc_nshare_pot_denom

                        chk_err_idx += 1
                        if chk_err_idx > len(neighbors_mobility[id][0]):
                            sys.exit(
                                'Error in functionality of vaccine distribution, exiting.')

                    if abs(vacc_dist[ageidx] - sum(vacc_nshare_pot)) > 0.01 * vacc_dist[ageidx]:
                        sys.exit(
                            'Error in functionality of vaccine distribution, exiting.')
                    # sort weights from max to min and return indices
                    order_dist = np.argsort(vacc_dist_weight)[::-1]
                    vacc_nshare_pot = vacc_nshare_pot[order_dist]
                    vacc_nshare_pot = vacc_nshare_pot[np.where(
                        vacc_nshare_pot > 0)]

                    # get vaccination shares to distribute
                    reduc_shares = vacc_nshare_pot / \
                        df_fullsum.loc[ii, unique_age_groups_old[ageidx]]

                    # sum up new additions
                    df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(neighbors_mobility[id][0]), [
                        'vd'+unique_age_groups_old[ageidx]]] = np.array(df_fullsum.loc[df_fullsum[dd.EngEng['idCounty']].isin(neighbors_mobility[id][0]), [
                            'vd'+unique_age_groups_old[ageidx]]]).flatten() + vacc_dist[ageidx]*vacc_dist_weight

                    # iterate over neighbors and add potential share of vaccionations
                    # to neighbor and decrease local vaccionations at the end
                    nidx = 0
                    for neighb_id in neighbors_mobility[id][0][order_dist][
                            np.where(vacc_nshare_pot > 0)]:
                        df_san.loc[(df_san[dd.EngEng['idCounty']] == neighb_id) & (df_san[dd.EngEng['ageRKI']] == unique_age_groups_old[ageidx]), column_names_new] += reduc_shares[nidx]*df_san.loc[(
                            df_san[dd.EngEng['idCounty']] == id) & (df_san[dd.EngEng['ageRKI']] == unique_age_groups_old[ageidx]), column_names_new].values
                        nidx += 1
                    df_san.loc[(df_san[dd.EngEng['idCounty']] == id) & (df_san[dd.EngEng['ageRKI']] == unique_age_groups_old[ageidx]), column_names_new] -= sum(
                        reduc_shares)*df_san.loc[(df_san[dd.EngEng['idCounty']] == id) & (df_san[dd.EngEng['ageRKI']] == unique_age_groups_old[ageidx]), column_names_new].values

            if ii < 401:
                end_time = time.perf_counter()
                print("Time needed for sanitizing county " + str(id) + ": " +
                      str(int(end_time - start_time_county)) + " sec")

            if len(
                    np.where(np.isnan(df_san[column_names_new]) == True)[0]) > 0:
                sys.exit(
                    'Error in functionality of vaccine distribution, NaN found after county '
                    + str(id) + '. Exiting program.')

        # create copy only for possible comparison reasons now
        dummy = df_data_agevacc_county_cs.copy()
        df_data_agevacc_county_cs = df_san.copy()

        # create cumulative sum
        df_data_agevacc_county_cs = df_data_agevacc_county_cs.groupby(
            [dd.EngEng['date'],
             dd.EngEng['idCounty'],
             dd.EngEng['ageRKI']]).agg(
            {column_new: sum for column_new in column_names_new})
        df_data_agevacc_county_cs = df_data_agevacc_county_cs.groupby(
            level=[1, 2]).cumsum().reset_index()

        # check that cumulative sum was done correctly (TODO: to be removed in production)
        for id in geoger.get_county_ids():
            for age in unique_age_groups_old:
                a = df_san[(df_san.ID_County == id) & (
                    df_san.Age_RKI == age)].sum()[column_names_new]
                b = df_data_agevacc_county_cs[(df_data_agevacc_county_cs.ID_County == id) & (
                    df_data_agevacc_county_cs.Age_RKI == age)].loc[:, column_names_new].iloc[-1]
                if sum(a-b) > 1e-8:
                    print("Error in: " + str(id) + " " + str(age))
        ### end of to be removed ###

        # compute the moving average
        df_data_agevacc_county_cs = modifyDataframeSeries.impute_and_reduce_df(
            df_data_agevacc_county_cs,
            {dd.EngEng['idCounty']: df_data_agevacc_county_cs[dd.EngEng['idCounty']].unique(),
             dd.EngEng['ageRKI']: unique_age_groups_old},
            [column_names_new[i]],
            impute='forward', moving_average=moving_average,
            min_date=min_date, max_date=max_date)

        end_time = time.perf_counter()
        print(
            "Time needed for sanitizing: " + str(int(end_time - start_time)) +
            " sec")

    # merge Eisenach and Wartburgkreis from vaccination data
    df_data_agevacc_county_cs = geoger.merge_df_counties_all(
        df_data_agevacc_county_cs,
        sorting=[dd.EngEng["idCounty"],
                 dd.EngEng["ageRKI"],
                 dd.EngEng["date"]],
        columns=[dd.EngEng["ageRKI"],
                 dd.EngEng["date"]])

   
    if make_plot:
        # have a look extrapolated vaccination ratios (TODO: to be removed or outsourced with plotting for production)
        # aggregate total number of vaccinations per county and age group
        latest_date=df_data_agevacc_county_cs[dd.EngEng["date"]][len(df_data_agevacc_county_cs.index)-1].strftime("%Y-%m-%d")
        vacc_sums_nonsanit = df_data_agevacc_county_cs.loc[(
            df_data_agevacc_county_cs.Date == latest_date), ['ID_County', column_names_new[1]]]
        # create new data frame and reshape it
        df_fullsum = pd.DataFrame(columns=[dd.EngEng['idCounty']] + unique_age_groups_old)
        df_fullsum[dd.EngEng['idCounty']
                        ] = vacc_sums_nonsanit[dd.EngEng['idCounty']].unique()
        df_fullsum[unique_age_groups_old] = np.array(
            vacc_sums_nonsanit[column_names_new[1]]).reshape(
            len(vacc_sums_nonsanit[dd.EngEng['idCounty']].unique()),
            len(unique_age_groups_old))
        # compute county and age-group-specific vaccination ratios
        df_fullsum[['r'+age for age in unique_age_groups_old] ] = df_fullsum[unique_age_groups_old] / geoger.merge_df_counties_all(
            population_old_ages, sorting=[dd.EngEng["idCounty"]], columns=dd.EngEng["idCounty"])[unique_age_groups_old].values
        ## end of the be removed ###

    # store data
    filename = 'all_county_agevacc_vacc'
    filename = gd.append_filename(filename, impute_dates, moving_average)
    gd.write_dataframe(df_data_agevacc_county_cs,
                       directory, filename, file_format)
    # make plot of absolute numbers original age resolution
    if make_plot:
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
        customPlot.plotList(
            date_vals, yvals, [age for age in unique_age_groups_old],
            'Partial vaccination over different age groups', dd.EngEng['date'],
            'Number', "Germany_PartialVacination_Absolute")

        # plot full vaccination curves for different age groups
        yvals = [
            df_data_agevacc_county_cs.loc
            [df_data_agevacc_county_cs[dd.EngEng['ageRKI']] == age,
             [dd.EngEng['date'],
              dd.EngEng['vaccComplete']]].groupby(dd.EngEng['date']).sum()
            for age in unique_age_groups_old]
        customPlot.plotList(
            date_vals, yvals, [age for age in unique_age_groups_old],
            'Full vaccination over different age groups', dd.EngEng['date'],
            'Number', "Germany_FullVacination_Absolute")

    ######## data without age resolution ###########
    # write data frame resolved per county
    start_time = time.perf_counter()
    df_data_county_cs = df_data_agevacc_county_cs.groupby(
        [dd.EngEng['date'],
         dd.EngEng['idCounty']]).agg(
        {col_new: sum for col_new in column_names_new}).reset_index()

    # store data
    filename = 'all_county_vacc'
    filename = gd.append_filename(filename, impute_dates, moving_average)
    gd.write_dataframe(df_data_county_cs, directory, filename, file_format)
    end_time = time.perf_counter()
    print("Time needed: " + str(int(end_time - start_time)) + " sec")

    ####### age resolved with extrapolation to other age groups #######
    # write data frame resolved per county and age (with age classes as
    # provided in RKI infection tables: 0-4, 5-14, 15-34, 35-59, 60-79, 80+)
    start_time = time.perf_counter()

    err_chk_time = 0
    if extrapolate_agegroups:
        # merge Eisenach (merge_eisenach=True)....
        population_all_ages = geoger.merge_df_counties_all(
            population_all_ages, sorting=[dd.EngEng["idCounty"]],
            columns=dd.EngEng["idCounty"])

        # create new, empty data frame
        df_data_ageinf_county_cs = pd.DataFrame(
            columns=df_data_agevacc_county_cs.columns)

        # loop over all counties
        for county in df_data_agevacc_county_cs[dd.EngEng['idCounty']].unique():
            # get population data of current state
            pop_state = population_all_ages[population_all_ages[dd.EngEng['idCounty']] == county]

            # access entries for the current county
            df_local_old = df_data_agevacc_county_cs[(
                df_data_agevacc_county_cs[dd.EngEng['idCounty']] == county)]
            # new data frame for the county
            df_local_new = pd.DataFrame(
                columns=[dd.EngEng['date'], dd.EngEng['idCounty']])

            # create dataframes for new age groups
            df_local_age_new = []
            for i in range(0, len(min_age_new)-1):
                # copy days and county ID (by using a dummy age group accessor)
                df_local_age_new.append(
                    df_local_old.loc
                    [df_data_agevacc_county_cs[dd.EngEng['ageRKI']] ==
                     unique_age_groups_old[0],
                     [dd.EngEng['date'],
                      dd.EngEng['idCounty']]].copy())
                # add age group column
                df_local_age_new[i][dd.EngEng['ageRKI']] = str(
                    min_age_new[i])+"-"+str(min_age_new[i+1]-1)
                # add specific zero columns
                df_local_age_new[i][column_names_new] = 0
                # reset indices to 0,...,length to be able to sum up columns correctly
                df_local_age_new[i] = df_local_age_new[i].reset_index(
                    drop=True)

            # iterate over old data frame and assign with correct shares
            # to new age ranges
            j = 0
            for i in range(0, len(unique_age_groups_old)):

                # get data specific for one county and one age group
                # of the original vaccination age group table
                df_local_age_old = df_local_old[df_local_old.Age_RKI ==
                                                unique_age_groups_old[i]]

                # get all age groups which are younger than the
                # minimum age of vaccinated persons
                while min_age_new[j+1] < min_age_old[i]:
                    # append zero vaccination group to local,
                    # county-specific data frame
                    df_local_new = df_local_new.append(
                        df_local_age_new[j].copy())
                    j += 1
                df_local_new = df_local_new.reset_index(drop=True)

                # copy shares of vaccinations from old age group columns
                # to new age group columns
                for k in age_old_to_age_new_share[i + new_age_not_vacc]:

                    common_ages_all_ages = list(
                        set(age_new_to_all_ages_indices[k]).intersection(
                            age_old_to_all_ages_indices[i]))
                    # splitting values of old column by comparing subset
                    # population sizes first, get population share
                    # that moves to new column
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
                # append summed vaccination age group to local,
                # county-specific data frame
                df_local_new = df_local_new.append(df_local_age_new[k].copy())

            # insert county name
            df_local_new.insert(
                loc=2, column=dd.EngEng["county"],
                value=df_local_new[dd.EngEng["idCounty"]].replace(
                    dd.County))

            # TODO: unit test for some randomly selected countiesand dates (here only stupid testing)
            start_time_err = time.perf_counter()
            if(abs(df_local_new.groupby(dd.EngEng['date']).agg({dd.EngEng['vaccComplete']: sum}) - df_data_agevacc_county_cs[(
                    df_data_agevacc_county_cs[dd.EngEng['idCounty']] == county)].groupby(dd.EngEng['date']).agg({dd.EngEng['vaccComplete']: sum})).sum()[dd.EngEng['vaccComplete']]) > 1e-5:
                print("Error in transformation...")
            err_chk_time += time.perf_counter() - start_time_err

            df_data_ageinf_county_cs = df_data_ageinf_county_cs.append(
                df_local_new.copy())

    # check for now, to be removed for production...
    for id in geoger.get_county_ids():
         if abs(df_data_ageinf_county_cs[df_data_ageinf_county_cs.ID_County==1001].groupby(['Date']).sum()-df_data_agevacc_county_cs[df_data_agevacc_county_cs.ID_County==1001].groupby(['Date']).sum()[column_names_new]).max().max() > 1e-8:
             print('Error in ' + str(id))
    ### end of to be removed ###

    filename = 'all_county_ageinf_vacc'
    filename = gd.append_filename(filename, impute_dates, moving_average)

    gd.write_dataframe(df_data_ageinf_county_cs,
                       directory, filename, file_format)
    end_time = time.perf_counter()
    print("Time needed for age extrapolation: " +
          str(int(end_time - start_time)) + " sec")
    print("Error check time " + str(err_chk_time))

    # make plot of relative numbers of original and extrapolated age resolution
    if make_plot:
        # extract (dummy) date column to plt
        date_vals = df_data_agevacc_county_cs.loc[
            (df_data_agevacc_county_cs[dd.EngEng['ageRKI']] ==
             unique_age_groups_old[0]) &
            (df_data_agevacc_county_cs[dd.EngEng['idCounty']] ==
             geoger.get_county_ids()[0])][dd.EngEng['date']]

        # create hashes to access columns of new age group intervals
        unique_age_groups_new = [
            str(min_age_new[i]) + "-" + str(min_age_new[i + 1] - 1)
            for i in range(len(min_age_new) - 1)]

        # consider first vaccination for new age groups
        yvals = [
            df_data_ageinf_county_cs.loc
            [df_data_ageinf_county_cs[dd.EngEng['ageRKI']] == age,
             [dd.EngEng['date'],
             dd.EngEng['vaccPartial']]].groupby(dd.EngEng['date']).sum()
            for age in unique_age_groups_new]

        customPlot.plotList(
            date_vals, yvals, [age for age in unique_age_groups_new],
            'Partial vaccination over different age groups', dd.EngEng['date'],
            'Number', "Germany_PartialVacination_AgeExtr_Absolute")

        # consider full vaccination for new age groups
        yvals = [
            df_data_ageinf_county_cs.loc
            [df_data_ageinf_county_cs[dd.EngEng['ageRKI']] == age,
             [dd.EngEng['date'],
             dd.EngEng['vaccComplete']]].groupby(dd.EngEng['date']).sum()
            for age in unique_age_groups_new]

        customPlot.plotList(
            date_vals, yvals, [age for age in unique_age_groups_new],
            'Full vaccination over different age groups', dd.EngEng['date'],
            'Number', "Germany_FullVacination_AgeExtr_Absolute")


def main():
    """! Main program entry."""

    arg_dict = gd.cli("vaccination")
    get_vaccination_data(**arg_dict)


if __name__ == "__main__":

    main()
