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
from datetime import datetime, timedelta
import sys
import time
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from scipy.cluster import hierarchy

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import defaultDict as dd
from memilio.epidata import customPlot


def validate(df_npis_old, df_npis, df_infec_rki, countyID, npiCode,
             start_npi_cols, npi_incid_start, start_date_validation,
             end_date_validation, fine_resolution):
    """! Validates the transformed NPI data based on read in NPI data list.
    Also works for incidence-dependent NPIs as long as no activation or lifting
    delay is used.
    @param df_npis_old Original data frame.
    @param df_npis New NPI data frame with potential activation of incidence-
        dependent NPIs.
    @param df_infec_rki Case data for incidence comparison.
    @param countyID CountyID of county to be validated.
    @@param npiCode NPI Code of code to be validated.
    @param start_npi_cols Index of column where NPI information start in
        original data frame.
    @param npi_incid_start Minimum incidence for activation of NPI codes.
    @param start_date_validation Start date for validation.
    @param end_date_validation End date for validation.
    """

    if fine_resolution == 1:
        npiCodes = [npiCode + code
                    for code in [''] + ['_' + str(i) for i in range(1, 6)]]
    else:
        npiCodes = [npiCode]
    for npiCode in npiCodes:
        dummy_old_rows = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == countyID) & (df_npis_old[dd.EngEng['npiCode']].isin(
            npiCodes))].iloc[:, :list(df_npis_old.columns).index(end_date_validation.strftime('d%Y%m%d'))+1]
    dummy_old = np.zeros(
        (dummy_old_rows.shape[1]-start_npi_cols, dummy_old_rows.shape[0]))
    for i in range(dummy_old.shape[1]):
        dummy_old[:, i] = dummy_old_rows.values[i][start_npi_cols:]

    dummy_new = df_npis.loc[(df_npis[dd.EngEng['idCounty']] == countyID) & (
        df_npis[dd.EngEng['date']] <= end_date_validation), npiCodes[0]].values
    incid = df_infec_rki[(df_infec_rki[dd.EngEng['idCounty']] == countyID) &
                         (df_infec_rki[dd.EngEng['date']] <=
                             end_date_validation) &
                         (df_infec_rki[dd.EngEng['date']] >=
                             start_date_validation)]['Incidence'].values

    if fine_resolution == 1:
        for col in range(dummy_old.shape[1]):
            npi_index = np.where(dummy_old[:, col] >= 1)[0]
            incid_index = np.where(incid >= npi_incid_start[npiCodes[col]])[0]
            active_index = np.sort(
                list(set(npi_index).intersection(incid_index)))
            nonactive_index = np.sort(
                list(set(npi_index).difference(active_index)))

            # deactivate values based on incidence before taking
            # the maximum over NPI group
            dummy_old[:, col][list(nonactive_index)] = 0

        # the intermediate step is necessary since the array has a 2nd dimension
        # of size zero if we directly set dummy_old = dummy_old.max(axis=1)
        dummy_old[:, 0] = dummy_old.max(axis=1)
        dummy_old = dummy_old[:, 0:1]

    npi_index = np.where(dummy_old[:, 0] >= 1)[0]
    incid_index = np.where(incid >= npi_incid_start[npiCodes[0]])[0]
    active_index = np.sort(
        list(set(npi_index).intersection(incid_index)))
    nonactive_index = np.sort(
        list(set(npi_index).difference(active_index)))

    valdiff = 0
    for i in range(1, 6):
        # these values >=1 are set to 1
        valdiff += (i-1)*len(
            np.where(dummy_old[list(active_index)] == i)[0])
        # these values >=1 are set to 0
        valdiff += i * len(
            np.where(dummy_old[list(nonactive_index)] == i)[0])
    # -99 always set to 0
    valdiff += 99*len(np.where(dummy_old == -99)[0])
    return [abs(dummy_old[:, 0]-dummy_new).sum(), valdiff, dummy_old, dummy_new]


def print_manual_download(filename, url):
    """! Print download message to ask the user manually download a file. """

    print(
        'This script needs manual downloading of files. Please register'
        ' at corona-datenplatform.com and download ' + filename + ' from ' + url +
        '. Then move it to a folder named raw_data in this directory.')


def read_files(directory, fine_resolution):
    """! Reads files from local directory and returns data in dataframes
    
    @param directory Directory where data is loaded from.
    @param fine_resolution 2 [Default] or 0 or 1. Defines which categories
        are considered.
        If '2' is set, all the subcategories (~1200) are considered.
        If '1' is set, all incidence levels of subcategories are merged and
            ~200 NPIs are considered.
        If '0' is chosen only the main, summarizing categories (~20) are used.
    @return Data frames df_npis_old (Decreed, encoded NPIs for all German
        counties) and df_npis_desc (Description of NPIs)
    """
    if fine_resolution > 0:
        try:
            df_npis_old = pd.read_csv(
                os.path.join(directory, 'kr_massnahmen_unterkategorien.csv'),
                sep=',')  # , nrows=1248)  # 1248 for debugging, only reading Flensburg
        except FileNotFoundError:
            print_manual_download(
                'kr_massnahmen_unterkategorien.csv',
                'https://www.corona-datenplattform.de/dataset/massnahmen_unterkategorien_kreise')
            raise FileNotFoundError
        df_npis_old.rename(dd.GerEng, axis=1, inplace=True)

        # check if rows hospitals and geriatric care are still empty;
        # these fields have been empty so far and are thus not used
        test_codes = ['M23_010', 'M23_020', 'M23_030', 'M23_040',
                      'M23_050', 'M23_060', 'M24_010', 'M24_020',
                      'M24_030', 'M24_040', 'M24_050', 'M24_060']
        for tcode in test_codes:
            for i in [''] + ["_" + str(i) for i in range(1, 6)]:
                if (df_npis_old[df_npis_old[dd.EngEng['npiCode']] == tcode+i].iloc[:, 6:].max().max() > 0):
                    print(tcode+i + " used.")
        # end check

    else:  # read aggregated NPIs

        try:
            df_npis_old = pd.read_csv(os.path.join(
                directory, 'kr_massnahmen_oberkategorien.csv'))
        except FileNotFoundError:
            print_manual_download(
                'datensatzbeschreibung_massnahmen.xlsx',
                'https://www.corona-datenplattform.de/dataset/massnahmen_oberkategorien_kreise')
            raise FileNotFoundError
        df_npis_old.rename(dd.GerEng, axis=1, inplace=True)

    # read data frame of variable names and descriptions
    try:
        if fine_resolution > 0:
            df_npis_desc = pd.read_excel(
                os.path.join(
                    directory, 'datensatzbeschreibung_massnahmen.xlsx'),
                sheet_name=2, engine = 'openpyxl')
        else:
            df_npis_desc = pd.read_excel(
                os.path.join(
                    directory, 'datensatzbeschreibung_massnahmen.xlsx'),
                sheet_name=3, engine = 'openpyxl')
    except FileNotFoundError:
        print_manual_download(
            'datensatzbeschreibung_massnahmen.xlsx',
            'https://www.corona-datenplattform.de/dataset/massnahmen_unterkategorien_kreise')
        raise FileNotFoundError

    return df_npis_old, df_npis_desc


def activate_npis_based_on_incidence(local_incid, npi_lifting_delay, npi_activation_delay, incid_threshold):
    """! 
    Warning: delay = 0 does NOT work!!!
    Computes an activation vector according to a given incidence threshold,
    observed incidence and activation or lifting delays.

    NPI can only be activated or lifted the day AFTER 
    incidence is below/over threshold for N days. The 
    incidence on day N only effects the NPI on day N+1 and
    NOT ON day N. Therefore we shift the incidence one day forward
    to match the indices of our dataframe df_local_new so that
    the NPIs can be calculated on the respective day.
    
    Please also note that the first column will always returned as false 
    so the dataframe should not start with dates where NPIs are implemented.
    For the Corona Datenplattform frame which starts from 2020-03-01
    this is no problem for the first days as there were no NPIs.
    @param local_incid Incidence observed in local region.
    @param npi_lifting_delay Delay in number of days when to lift a NPI if
        threshold is no longer exceeded.
    @param npi_activation_delay Delay in number of days when to implement a
        NPI if threshold is exceeded.
    @param incid_threshold Threshold to be considered.
    """
    # 
    # Example (threshold=3.5): 
    # local_incid=pd.Series([2,4,2,4,2,2,4,4,2,4,2,2,2,2])
    # Yesterdays incidence is over the threshold on following days:
    # [?,0,1,0,1,0,0,1,1,0,1,0,0,0]
    # Yesterday for first day is not known. Thus, first day's boolean
    # is always set to the same boolean as second day's boolean.
    # [0,0,1,0,1,0,0,1,1,0,1,0,0,0]
    
    # First get a Series with 0 for yesterdays incidence
    # is below threshold and 1 for incidence over threshold
    yesterdays_incid_over_threshold = (local_incid.shift(
        1).fillna(local_incid[0]) > incid_threshold).astype(int)            # maybe >=

    # If incidence is above threshold for 
    # 1+npi_activation_delay days, the NPI gets activated.
    # Similarly, if incidence is below threshold for 
    # 1+npi_lifting_delay days, the NPI is lifted.
    # 
    # Example (continued):
    # With yesterdays incidence over threshold on days:
    # [0,0,1,0,1,0,0,1,1,0,1,0,0,0]
    # npi_lifting_delay=2, npi_activation_delay=1
    # NPI should be activated on day 9 and lifted on day 14
    # int_active should then be:
    # [0,0,0,0,0,0,0,0,1,1,1,1,1,0]
    #
    # Another example: 
    # With yesterday's incidence over threshold on days:
    # [1,1,0,1,0,0,1,1,0,1,0,0,0,0]
    # npi_lifting_delay=3, npi_activation_delay=1
    # NPI should be activated on day 2 and lifted on day 14
    # int_active should then be:
    # [0,1,1,1,1,1,1,1,1,1,1,1,1,0]

    # get a zero filled Series with same length to be
    # filled with ones where NPI is active
    int_active = pd.Series(np.zeros(len(local_incid), dtype=int))
    # loop over every day
    for i in range(len(yesterdays_incid_over_threshold)):
        # Set int_active=0 where last npi_lifting_delay+1 days are 0
        if yesterdays_incid_over_threshold[max(0,i-npi_lifting_delay):i+1].values.sum() == 0:
            int_active[i] = 0
        # Set int_active=1 where last npi_activation_delay+1 days are all 1
        elif yesterdays_incid_over_threshold[max(0,i-npi_activation_delay):i+1].values.sum() == npi_activation_delay+1:
            int_active[i] = 1
        # If no condition applies, set int_active to the value of the previous day
        elif i>0: # if 
            int_active[i] = int_active[i-1]

    return int_active


def drop_codes_and_categories(npi_codes_prior, npi_codes_prior_desc, df_npis_old, fine_resolution):
    """! Drops codes and categories from original data frame if they are not 
    used.

    @param npi_codes_prior NPI codes read from description sheet.
    @param npi_codes_prior_desc NPI code descriptions read from description sheet.
    @param df_npis_old Original data frame with encoding (decreed, yes or no)
        for different NPIs
    @param fine_resolution 2 [Default] or 0 or 1. Defines which categories
        are considered.
        If '2' is set, all the subcategories (~1200) are considered.
        If '1' is set, all incidence levels of subcategories are merged and
            ~200 NPIs are considered.
        If '0' is chosen only the main, summarizing categories (~20) are used.
    @return Returns dropped codes, prior codes and reduced original data frame.
    """
    if fine_resolution > 0:
        # correct M04_N codes to M_04_M_N, N in {1,...,5}, M in {120,110,100,130,140}
        # (M04_1, i.e. i=1, has been corrected in original file but not for i>1)
        for i in range(2, 6):
            npi_codes_prior[npi_codes_prior == 'M04_'+str(i)] = ['M04_120_'+str(
                i), 'M04_110_'+str(i), 'M04_100_'+str(i), 'M04_130_'+str(i), 'M04_140_'+str(i)]

        # correct M05_N codes to M_05_M_N, N in {1,...,5}, M in {130,150,120,140,110,100,160}
        for i in range(1, 6):
            npi_codes_prior[npi_codes_prior == 'M05_'+str(i)] = ['M05_130_'+str(i), 'M05_150_'+str(
                i), 'M05_120_'+str(i), 'M05_140_'+str(i), 'M05_110_'+str(i), 'M05_100_'+str(i), 'M05_160_'+str(i)]

        # correct 'M16_200_2' to missing 'M16_100_2'
        npi_codes_prior[npi_codes_prior == 'M16_200_2'] = 'M16_100_2'

        # check for missing codes
        npi_codes_prior_data = df_npis_old[dd.EngEng['npiCode']].unique()

        missing_codes = list(set(npi_codes_prior).difference(
            npi_codes_prior_data))
        if len(missing_codes) > 0:
            # if incidence is grouped, only search for grouping codes without
            # having a detailed "_DETAIL" naming as of MCODE_NUMBER_DETAIL
            if fine_resolution == 1:
                missing_grouped_codes = []
                for mcode in missing_codes:
                    # only consider incidence independent npis
                    # only exit if one of these (i.e., MCODE_NUMBER) is missing
                    if len(mcode.split('_')) != 3:
                        missing_grouped_codes.append(mcode)
                if len(missing_grouped_codes) > 0:  # only MCODE_NUMBER codes
                    raise gd.DataError('Missing NPI codes: ' +
                                str(missing_grouped_codes))
            else:
                raise gd.DataError('Missing NPI codes: ' + str(missing_codes))

        # we dont have any explanations from "datensatzbeschreibung_massnahmen"
        # on these codes, so drop the rows.
        codes_dropped = list(set(npi_codes_prior_data).difference(
            npi_codes_prior))
        # also remove dummy 'Platzhalter' categories
        dummy_categories = []
        for i in range(len(npi_codes_prior)):
            if 'Platzhalter' in npi_codes_prior_desc[i]:
                dummy_categories.append(npi_codes_prior[i])
        # codes without explanation and dummy categories
        # sorting done for consistenty, maybe not necessary
        codes_dropped = list(np.sort(codes_dropped + dummy_categories))
        if len(codes_dropped) > 0:
            df_npis_old = df_npis_old[~df_npis_old[dd.EngEng['npiCode']].isin(
                codes_dropped)].reset_index(drop=True)
            # for every main code removed, all 5 subcodes have to be removed;
            # if this is not the case, the naming of them is wrong/not consistent
            if (len(codes_dropped) % 6) != 0:
                raise gd.DataError('Error in NPI names, please check.')
    else:
        # no dropping for fine_resolution == 0
        codes_dropped=[]
    
    return codes_dropped, npi_codes_prior, df_npis_old



def get_npi_data(fine_resolution=2,
                 file_format=dd.defaultDict['file_format'],
                 out_folder=dd.defaultDict['out_folder'],
                 start_date=dd.defaultDict['start_date'],
                 end_date=dd.defaultDict['end_date'],
                 counties_considered=geoger.get_county_ids()
                 ):
    """! Loads a certain resolution of recorded NPI data from
    the Corona Datenplattform and extracts the counties asked for and
    activates the NPIs if they are incidence dependent.

    Results' data frames will be stored in the directory as:
        -fine_resolution=2: germany_counties_npi_subcat
        -fine_resolution=1: germany_counties_npi_subcat_incgrouped
        -fine_resolution=0: germany_counties_npi_maincat

    Needs the files 'cases_all_county_all_dates_repdate.json' and 
    'county_current_population.json' which can be created by the functions
    getCasesData.py (with argument --rep-date) and getPopulationData.py.

    Please manually download
    - kr_massnahmen_unterkategorien.csv
    - datensatzbeschreibung_massnahmen.xlsx
    from https://www.corona-datenplattform.de/dataset/massnahmen_unterkategorien_kreise
    and
    - kr_massnahmen_oberkategorien.csv
    from https://www.corona-datenplattform.de/dataset/massnahmen_oberkategorien_kreise
    and move it to the *directory*-path mentioned in the beginning of the function.

    @param fine_resolution 2 [Default] or 0 or 1. Defines which categories
        are considered.
        If '2' is set, all the subcategories (~1200) are considered.
        If '1' is set, all incidence levels of subcategories are merged and
            ~200 NPIs are considered.
        If '0' is chosen only the main, summarizing categories (~20) are used.
    @param file_format File format which is used for writing the data.
        Default defined in defaultDict.
    @param out_folder Path to folder where data is written in folder
        out_folder/Germany.
    @param start_date [Default = '', taken from read data] Start date
        of stored data frames.
    @param end_date [Default = '', taken from read data] End date of
        stored data frames.
    @param counties_considered [Default: 'All']. Either 'All' or a list of
        county IDs from 1001 to 16xxx.
    """

    if counties_considered == 'All':
        counties_considered = geoger.get_county_ids()

    directory = out_folder
    directory = os.path.join(directory, 'Germany/')
    gd.check_dir(directory)

    if fine_resolution > 0:
        # defines delay in number of days between exceeding
        # incidence threshold and NPI getting active
        # delay = 0 means only one day is considered (=no delay)
        npi_activation_delay = 3
        npi_lifting_delay = 5
        # depending on the federal state and time period, there are 
        # huge deviations of the lifting and activation delay which was usually
        # between 1 and 14 days
        # we use npi_lifting_delay = 5 and npi_activation_delay = 3 
        # as this is the most common and has at some point been used in almost every county
        print('Using a delay of NPI activation of ' +
              str(npi_activation_delay) + ' days.')
        print('Using a delay of NPI lifting of ' +
              str(npi_lifting_delay) + ' days.')

    # read manual downloaded files from directory
    df_npis_old, df_npis_desc = read_files(directory, fine_resolution)

    # get existing codes that are used (in df_npis_old M22-M24 are empty)
    npi_codes_prior = df_npis_desc['Variablenname']
    npi_codes_prior_desc = df_npis_desc['Variable']

    # for fine_resolution > 0 deactivation of non-combinable
    # incidence-dependent NPIs has to be conducted; therefore we defined a
    # matrix of possible combinations of NPIs (marked with an X if combinable)
    # NPIs of different main category (e.g., M01a and M04) can always be
    # combined; only those of, e.g., M01a_010_3 and M01a_080_4 can exclude each
    # other
    if fine_resolution > 0:
        df_npis_combinations_pre = pd.read_excel(
            os.path.join(
                directory, 'combination_npis.xlsx'), engine = 'openpyxl')

        # rename essential columns and throw away others
        column_names = ['Unnamed: ' + str(i) for i in range(3, 19)]
        rename_columns = {column_names[i]: i for i in range(len(column_names))}
        df_npis_combinations_pre.rename(columns=rename_columns, inplace=True)
        df_npis_combinations_pre = df_npis_combinations_pre[[
            'Variablenname'] + [i for i in range(0, 16)]]
        # replace empty cells by zeros and x-marked cells by ones
        df_npis_combinations_pre = df_npis_combinations_pre.replace(np.nan, 0)
        df_npis_combinations_pre = df_npis_combinations_pre.replace('x', 1)

        # extract different NPI groups and store indices of NPIs belonging
        # to the different groups
        npi_groups_combinations = pd.Series(
            code.split('_')[0]
            for code in df_npis_combinations_pre['Variablenname'])
        npi_groups_combinations_unique = npi_groups_combinations.unique()
        npi_groups_idx = []
        for code in npi_groups_combinations_unique:
            npi_groups_idx.append(
                list(
                    npi_groups_combinations
                    [npi_groups_combinations == code].index))
        # create hash table of main code to contained codes and combination matrix
        df_npis_combinations = {
            npi_groups_combinations_unique[i]:
            [
                list(
                    df_npis_combinations_pre['Variablenname']
                    [npi_groups_idx[0]]),
                np.eye(len(npi_groups_idx[i]))]
            for i in range(len(npi_groups_combinations_unique))}

        # run through all groups and set possible combinations according to
        # read combination matrix
        start_comb_matrix = list(
            df_npis_combinations_pre.columns).index('Variablenname')+1
        for i in range(len(npi_groups_idx)):
            codes_local = df_npis_combinations_pre.loc[npi_groups_idx[i],
                                                       'Variablenname'].values
            df_npis_combinations[npi_groups_combinations_unique[i]][1] = df_npis_combinations_pre.iloc[npi_groups_idx[i],
                                                                                                       start_comb_matrix:start_comb_matrix+len(npi_groups_idx[i])].values
            if (df_npis_combinations[npi_groups_combinations_unique[i]][1]-np.transpose(df_npis_combinations[npi_groups_combinations_unique[i]][1])).max() > 0:
                print('Error in input file: Please correct combination matrix input.')
            # make it a dataframe to allow easy removal of code lines and rows
            # if they are not used later on
            df_npis_combinations[npi_groups_combinations_unique[i]][1] = pd.DataFrame(
                df_npis_combinations[npi_groups_combinations_unique[i]][1],
                columns=codes_local)
            df_npis_combinations[npi_groups_combinations_unique[i]][1].insert(
                0, 'Code', codes_local)

        del df_npis_combinations_pre

        # use to_excel function and specify the sheet_name and index
        # to store the dataframe in specified sheet if file not yet existent
        # otherwise just valid results against stored sheets
        write_file = False
        if not os.path.exists(os.path.join(
                directory,
                'combinations_npis_cleanoutput.xlsx')):
            writer = pd.ExcelWriter(os.path.join(
                directory, 'combinations_npis_cleanoutput.xlsx'))
            write_file = True
        for i in range(len(npi_groups_combinations_unique)):
            codes_local = df_npis_combinations[npi_groups_combinations_unique[i]
                                               ][1].columns[1:]
            df_out = df_npis_combinations[npi_groups_combinations_unique[i]][
                1].copy()
            df_out.insert(
                0, 'Description (German)',
                [desc
                 for desc in npi_codes_prior_desc
                 [npi_codes_prior.isin(codes_local)].values])
            try:
                # store verified output
                df_in_valid = pd.read_excel(
                    os.path.join(
                        directory, 'combinations_npis_cleanoutput.xlsx'),
                    sheet_name=i, engine = 'openpyxl')
                if not df_in_valid.drop(columns='Unnamed: 0').equals(df_out):
                    print('Error in combination matrix.')
                del df_in_valid
            except:
                pass

            if write_file:
                df_out.to_excel(
                    writer, sheet_name=npi_groups_combinations_unique[i])
            del df_out
        if write_file:
            writer.save()

    # correct differences in codes between data sheet and explanation sheet
    codes_dropped, npi_codes_prior, df_npis_old = drop_codes_and_categories(
        npi_codes_prior, npi_codes_prior_desc, df_npis_old, fine_resolution)

    # sort NPI codes according to numeric values (argsort gives indices
    # in input list to be used for sorted array)
    npi_codes_sorting = np.argsort(npi_codes_prior.values)
    npi_codes = list(npi_codes_prior[npi_codes_sorting])
    if fine_resolution > 0:
        # for subcategories, description is in "Beschreibung" column; The
        # "Variable" column is repeated after the ";" sign
        # (except for 6 first rows where there is probably some error)
        npi_desc = list(df_npis_desc["Beschreibung"][npi_codes_sorting])

        # Check for consistent naming in descriptions
        # Errors are known for the first 6 rows
        dummy_a = list(df_npis_desc["Variable"][npi_codes_sorting])
        dummy_b = df_npis_desc["Beschreibung"][npi_codes_sorting]
        dummy_c = [str(x).split("; ")[1] for x in dummy_b]
        errors = []
        for i in range(len(dummy_a)):
            if not dummy_a[i] == dummy_c[i]:
                errors.append(i)
        if not errors == [0, 1, 2, 3, 4, 5]:
            print("Additional errors in consistent naming.")
        # End of check

        # correct for consistent naming (mainly done for plotting reasons,
        # otherwise naming column not strictly necessary)
        for i in range(errors[0], errors[-1]+1):
            npi_desc[i] = npi_desc[i].split("; ")[0] + "; " + dummy_a[i]

    else:
        # extract variable names for main categories
        npi_desc = list(df_npis_desc["Variable"][npi_codes_sorting])

    del df_npis_desc

    # combine NPI codes and descriptions to ensure that both are ordered
    # the same way; delete npi_codes or npi_desc for not using hereafter
    idx_codes_retained = ~pd.Series(npi_codes).isin(codes_dropped)
    npis = pd.DataFrame({
        dd.EngEng['npiCode']: list(pd.Series(npi_codes)[idx_codes_retained]),
        dd.EngEng['desc']: list(pd.Series(npi_desc)[idx_codes_retained])})
    del npi_codes
    del npi_desc
    # remove rows and columns of unused codes
    if fine_resolution > 0:
        for code in df_npis_combinations.keys():
            local_codes_used_rows = df_npis_combinations[code][1].Code.isin(
                npis.NPI_code)
            local_codes_used_cols = df_npis_combinations[code][1].columns.isin(
                npis.NPI_code)

            # overwrite item 0 since codes are stored in *.columns
            df_npis_combinations[code] = df_npis_combinations[code][1].loc[local_codes_used_rows,
                                                                        local_codes_used_cols].reset_index(drop=True).copy()

    # prepare grouping of NPIs to reduce product space of
    # NPI x active_from_inc (with values "incidence does not matter", and
    # incidence 0, 10, 35, 50, 100) to NPI
    if fine_resolution == 1:
        # create hash table from parental or main code/main category
        # to list of subcodes/subcategories
        maincode_to_npicodes_map = dict()
        major_code = npis.iloc[:, 0][0]
        maincode_to_npicodes_map[major_code] = []
        for code in npis.iloc[:, 0]:
            if major_code in code:
                maincode_to_npicodes_map[major_code].append(code)
            else:
                major_code = code
                maincode_to_npicodes_map[major_code] = [code]

        npi_codes_aggregated = []
        for main_code in maincode_to_npicodes_map.keys():
            if main_code.count('_') > 1:
                raise gd.DataError('Error. Subcode assigned as main code.')
            npi_codes_aggregated.append(main_code)

        npis_final = npis[npis[dd.EngEng['npiCode']].isin(
            npi_codes_aggregated)].reset_index()
    else:
        npis_final = npis

    # extract incidence-threshold for NPIs
    if fine_resolution > 0:
        npi_incid_start = dict()
        for i in range(len(npis)):
            incid_threshold = 1e5
            if npis.loc[i, dd.EngEng['desc']].split(' ')[0] == 'Unabhängig':
                # set -1 for incidence-independent NPIs
                incid_threshold = -1
            elif npis.loc[i, dd.EngEng['desc']].split(' ')[0] == 'Ab':
                incid_threshold = int(
                    npis.loc[i, dd.EngEng['desc']].split(' ')[1])
            else:
                raise gd.DataError(
                    'Error in description file. NPI activation can not '
                    'be computed. Exiting.')
            npi_incid_start[npis.loc[i, dd.EngEng['npiCode']]
                            ] = incid_threshold

        # get all incidence thresholds (This list has to be sorted)
        incidence_thresholds = []
        for code, threshold in npi_incid_start.items():
            if len(code.split('_')) < 3:
                if not (threshold, '') in incidence_thresholds:
                    incidence_thresholds.append((threshold, ''))
            else:
                if not (threshold, '_' + code.split('_')[2]) in incidence_thresholds:
                    incidence_thresholds.append(
                        (threshold, '_' + code.split('_')[2]))
        for i in range(len(incidence_thresholds)-1):
            if incidence_thresholds[i][0] > incidence_thresholds[i+1][0]:
                raise gd.DataError('List needs to be sorted.')

        # create hash map from thresholds to NPI indices
        incidence_thresholds_to_npis = dict(
            zip(incidence_thresholds, [[] for i in range(len(incidence_thresholds))]))
        for i in range(len(npis)):
            code_considered = npis.loc[i, dd.EngEng['npiCode']]
            incval = npi_incid_start[code_considered]
            if len(code_considered.split('_')) < 3:
                incidence_thresholds_to_npis[(incval, '')].append(i)
            else:
                incidence_thresholds_to_npis[(
                    incval, '_' + code_considered.split('_')[2])].append(i)

    # check if more than the county of Eisenach would be removed with
    # current county list
    counties_removed = df_npis_old[
        ~df_npis_old[dd.EngEng['idCounty']].isin(counties_considered)][
        dd.EngEng['idCounty']].unique()
    if list(counties_removed) != [16056]:
        raise gd.DataError('Error. Other counties than that of Eisenach were removed.')
    # remove rows for Eisenach
    df_npis_old = df_npis_old[df_npis_old[dd.EngEng['idCounty']].isin(
        counties_considered)].reset_index(drop=True)

    start_npi_cols = list(
        df_npis_old.columns).index(
        dd.EngEng['npiCode']) + 1

    # store string dates 'dYYYYMMDD' in list before parsing
    str_dates = list(df_npis_old.iloc[:, start_npi_cols:].columns)
    # convert string dates into other format
    dates_new = [datetime.strptime(old_date, "d%Y%m%d")
                 for old_date in str_dates]

    # check for missing dates
    date_diff = [
        (dates_new[i + 1] - dates_new[i]).days
        for i in range(len(dates_new) - 1)]
    date_diff_idx = np.where(np.array(date_diff) > 1)[0]
    if max(date_diff) > 1:
        print("Error. Dates missing in data frame:")
        for i in date_diff_idx:
            print(
                "\t - From " + str(dates_new[i] + timedelta(1)) + " until " +
                str(dates_new[i] + timedelta(date_diff[i] - 1)))
        raise gd.DataError('Exiting. Dates missing in data frame.')

    min_date = []
    max_date = []

    # get RKI infectious numbers to find dates where incidence-dependent
    # NPIs were active
    if fine_resolution > 0:
        df_infec_rki = pd.read_json(os.path.join(
            directory, 'cases_all_county_all_dates_repdate.json'))
        df_infec_rki[dd.EngEng['date']] = pd.to_datetime(
            df_infec_rki[dd.EngEng['date']])
        df_population = pd.read_json(
            directory + "county_current_population.json")
        min_date.append(
            df_infec_rki[dd.EngEng['date']].min().to_pydatetime())
        max_date.append(
            df_infec_rki[dd.EngEng['date']].max().to_pydatetime())

    # adapt time series according to available dates and start_date,
    # end_date input parameter
    start_date_new = max(
        min_date + [min(dates_new), pd.to_datetime(start_date)])
    end_date_new = min(
        max_date + [max(dates_new),
                    pd.to_datetime(end_date)])

    # create new data frame for all NPIs given in the columns,
    # resolved by county and day
    df_npis = pd.DataFrame(
        columns=[dd.EngEng['date']] + [dd.EngEng['idCounty']] +
        list(npis_final[dd.EngEng['npiCode']]))
    # convert NPI data from object to int such that correlations can be
    # computed
    df_npis = df_npis.astype(dict(
        zip(
            [dd.EngEng['date']] + [dd.EngEng['idCounty']] +
            list(npis_final[dd.EngEng['npiCode']]), ['str', 'int'] +
            ['int' for i in npis_final[dd.EngEng['npiCode']]])))

    # iterate over countyIDs
    counters = np.zeros(4)  # time counter for output only
    countyidx = 0
    # replace -99 ("not used anymore") by 0 ("not used")
    # replace 2,3,4,5 ("mentioned in ...") by 1 ("mentioned")
    df_npis_old.replace([-99, 2, 3, 4, 5], [0, 1, 1, 1, 1], inplace=True)

    for countyID in counties_considered:
        cid = 0
        countyidx += 1

        if fine_resolution > 0:
            # compute incidence based on previous data frames
            df_infec_local = df_infec_rki[df_infec_rki[dd.EngEng['idCounty']] == countyID].copy(
            )
            pop_local = df_population.loc[df_population[dd.EngEng['idCounty']]
                                          == countyID, dd.EngEng['population']].values[0]
            # consider difference between current day and day-7 to compute incidence
            df_infec_local['Incidence'] = df_infec_local[dd.EngEng['confirmed']].diff(
                periods=7).fillna(df_infec_local[dd.EngEng['confirmed']]) / pop_local * 100000

            # set to main data frame
            df_infec_rki.loc[df_infec_rki[dd.EngEng['idCounty']] ==
                             countyID, 'Incidence'] = df_infec_local['Incidence'].values

            # cut infection information at start_date_new and end_date_new
            df_infec_local = df_infec_local[(df_infec_local[dd.EngEng['date']] >= start_date_new) & (
                df_infec_local[dd.EngEng['date']] <= end_date_new)].reset_index()

        # get county-local data frame
        start_time = time.perf_counter()
        df_local_old = df_npis_old[df_npis_old[dd.EngEng['idCounty']]
                                   == countyID].copy()

        # potentially remove rows if they are not in npis dict
        npi_rows = [i in npis[dd.EngEng['npiCode']].values
                    for i in df_local_old[dd.EngEng['npiCode']]]


        # create columns for date, county ID 
        df_local_new = pd.DataFrame(
            columns=[dd.EngEng['date']] + [dd.EngEng['idCounty']])

        counters[cid] += time.perf_counter()-start_time
        cid += 1

        start_time = time.perf_counter()

        # old dataframe has npi codes as columns and date values as rows
        # new dataframe should be transposed
        df_local_new = df_local_old.iloc[npi_rows, start_npi_cols-1:].set_index(
            dd.EngEng['npiCode']).transpose().copy()
        # get datetime as a column (previously index after transposing)
        df_local_new = df_local_new.reset_index(
            drop=False).rename(
            columns={'index': dd.EngEng['date']})
        # reset index name (which is dd.EngEng['npiCode'] after transposing)
        df_local_new.rename_axis('', axis=1, inplace=True)
        # change time format from 'dYYYYMMDD' to datetime timestamps
        df_local_new[dd.EngEng['date']] = pd.to_datetime(
            df_local_new[dd.EngEng['date']], format='d%Y%m%d')
        # fill in column for county ID
        df_local_new[dd.EngEng['idCounty']] = countyID
        # sort columns as to {Date, ID_County, npi_codes...}
        # for now this can be done alphabetically
        df_local_new.sort_index(axis=1, inplace=True)

        counters[cid] += time.perf_counter()-start_time
        cid += 1

        ### evaluate NPIs mentioned with respect to confirmed cases ###
        # values > 0
        #   - for NPIs independent of new infections mean "mentioned" = "active"
        #   - for NPIs dependent on incidence "mentioned" does not mean
        #       active and evaluation has to be conducted against confirmed
        #       infections to determine whether the NPI was active
        start_time = time.perf_counter()
        if fine_resolution > 0:
            # cut NPI information at start_date_new and end_date_new
            df_local_new = df_local_new.loc[(df_local_new[dd.EngEng['date']] >= start_date_new) & (
                df_local_new[dd.EngEng['date']] <= end_date_new), :].reset_index()
            try:
                df_local_new = df_local_new.drop(columns='index')
            except:
                pass
            # get index of first NPI column in local data frame
            npis_idx_start = list(
                df_local_new.columns).index(
                npis[dd.EngEng['npiCode']][0])

            # iterate through all NPIs and activate if incidence threshold
            # is exceeded
            for level, npi_indices in incidence_thresholds_to_npis.items():
                if level[0] >= 0:  # level[0] = incidvalthrsh
                    local_incid = df_infec_local['Incidence'].copy()

                    # get days where npis are active as int (1/0)
                    int_active = activate_npis_based_on_incidence(local_incid, npi_lifting_delay, npi_activation_delay, level[0])

                    # multiply rows of data frame by either 1 if threshold
                    # passed (i.e., mentioned NPI is active) or zero
                    # (i.e., mentioned NPI is not active)
                    # 'mul' multiplies the original data frame row by row
                    # with the respective value in int_active
                    df_local_new.iloc[:, npis_idx_start + np.array(npi_indices)] \
                        = df_local_new.iloc[:, npis_idx_start + np.array(npi_indices)].mul(int_active, axis=0)

            # if new, dynamic NPIs for higher incidence (more restrictions,
            # i.e., stricter) cannot be combined with previous, dynamic
            # NPIs for lower indices (less restrictions, less strict),
            # the latter have to be deactivated
            # (incidence_thresholds_to_npis.keys() has to be sorted !)
            levels_exclusion = list(reversed(incidence_thresholds_to_npis.keys()))[
                0:-1]  # level<0 means non-incidence dependent and always active
            for level in levels_exclusion:
                level_lower = [lev for lev in levels_exclusion
                               if lev[0] < level[0]]
                for code in df_npis_combinations.keys():
                    code_cols = df_npis_combinations[code].columns
                    # iterate over subcode indices
                    for scidx in range(len(code_cols)-1):
                        # check if code was used, otherwise nothing to
                        # exclude, i.e. no combination possible anyway.
                        indicator_code_active = df_local_new.loc[:,
                                                                 code_cols
                                                                 [scidx] +
                                                                 level
                                                                 [1]]
                        indicator_code_active_idx = np.where(
                            indicator_code_active > 0)[0]
                        if len(indicator_code_active_idx) > 0:
                            # extract codes
                            subcodes_nocombi = df_npis_combinations[code].loc[scidx, :]
                            # only consider those codes which cannot be
                            # combined; for these values of 1 have to be
                            # set to 0
                            subcodes_nocombi = list(
                                subcodes_nocombi
                                [subcodes_nocombi == 0].index)
                            # iterate over exclusive subcodes
                            for subcode_excl in subcodes_nocombi:
                                # iterate over less strict dynamic NPIs
                                # i.e., where threshold is higher
                                for level_other in level_lower:
                                    # deactivate potential NPIs (with code:
                                    # subcode_excl + level_other[1]) on days
                                    # where NPI code_cols[scidx] + level[1]
                                    # is active
                                    if df_local_new.loc[indicator_code_active_idx, subcode_excl + level_other[1]].any():
                                        # print('Deactivating for ' + 'County ' + str(countyID) + ' and Incidence ' + str(level_other[0]))
                                        # print('\t' + npi_codes_prior_desc[npi_codes_prior[npi_codes_prior==subcode_excl + level_other[1]].index].values[0])
                                        # print('Due to Incidence > ' + str(level[0]) + ' and NPI ')
                                        # print('\t' + npi_codes_prior_desc[npi_codes_prior[npi_codes_prior==code_cols[scidx] + level[1]].index].values[0])
                                        # print(list(df_local_new.loc[indicator_code_active_idx,'Date']))
                                        # print('\n')
                                        df_local_new.loc[indicator_code_active_idx,
                                                         subcode_excl + level_other[1]] = 0

            # reduction of factor space NPI x incidence threshold to NPI
            # by max aggregation of all incidence threshold columns per NPI
            if fine_resolution == 1:
                for main_code, codes_group in maincode_to_npicodes_map.items():
                    # group by incidence (former codes X1_Y, X1_Z were transformed
                    # to X1, X2) and write max value to main code column
                    df_local_new.loc[:, main_code] = df_local_new.loc[:, codes_group].max(
                        axis=1)
                # remove subcategory columns
                df_local_new = df_local_new.loc[:, [
                    dd.EngEng['date'], dd.EngEng['idCounty']] + npi_codes_aggregated].copy()

        counters[cid] += time.perf_counter()-start_time
        cid += 1
        ### ###

        start_time = time.perf_counter()

        df_npis = pd.concat([df_npis.copy(), df_local_new.copy()], ignore_index=True)
        counters[cid] += time.perf_counter()-start_time
        cid += 1

        # divide working time by completed number of counties and multiply
        # by remaining number of counties to estimate time remaining
        time_remain = sum(
            counters) / countyidx * (len(counties_considered) - countyidx)
        # print progress
        if countyidx == 1 or countyidx % int(
                len(counties_considered) / 10) == 0:
            print('Progress ' + str(countyidx) + ' / ' +
                  str(len(counties_considered)) +
                  '. Estimated time remaining: ' +
                  str(int(time_remain / 60)) + ' min.')

    # print sub counters
    print('Sub task counters are: ')
    print(counters)

    # reset index and drop old index column
    df_npis.reset_index(inplace=True)
    try:
        df_npis = df_npis.drop(columns='index')
    except:
        pass
    try:
        df_npis = df_npis.drop(columns='level_0')
    except:
        pass

    #### start validation ####
    if fine_resolution == 2 and (npi_activation_delay + npi_lifting_delay == 0):
        start_date_validation = datetime(2020, 3, 1)
        end_date_validation = datetime(2022, 2, 15)

        for countyID in counties_considered:
            for npiCode in [
                'M01a_010', 'M01a_150', 'M05_120', 'M01a_010',
                    'M18_030', 'M01b_020', 'M02b_035', 'M16_050']:
                for subcode in [''] + ['_'+str(i) for i in range(1, 6)]:
                    [
                        a, b, oldf, newf] = validate(
                        df_npis_old, df_npis, df_infec_rki, countyID,
                        npiCode + subcode, start_npi_cols, npi_incid_start,
                        start_date_validation, end_date_validation,
                        fine_resolution)
                    if (a != b):
                        print('Error in NPI activation computation')
                        print(a, b, a - b)

    elif fine_resolution == 1:
        start_date_validation = datetime(2020, 3, 1)
        end_date_validation = datetime(2022, 2, 15)

        for countyID in counties_considered:
            for npiCode in [
                'M01a_010', 'M01a_150', 'M05_120', 'M01a_010',
                    'M18_030', 'M01b_020', 'M02b_035', 'M16_050']:
                [a, b, oldf, newf] = validate(df_npis_old, df_npis,
                                              df_infec_rki, countyID, npiCode, start_npi_cols,
                                              npi_incid_start, start_date_validation,
                                              end_date_validation, fine_resolution)
                if (a != b):
                    print('Error in NPI activation computation')
                    print(a, b, a == b)
    #### end validation ####

    if fine_resolution > 0:
        if fine_resolution == 1:
            filename = 'germany_counties_npi_subcat_incgrouped'
        else:
            filename = 'germany_counties_npi_subcat'
    else:
        filename = 'germany_counties_npi_maincat'
    gd.write_dataframe(df_npis, directory, filename, file_format)


def main():
    """! Main program entry."""

    # arg_dict = gd.cli("testing")
    get_npi_data(fine_resolution=2)


if __name__ == "__main__":

    main()
