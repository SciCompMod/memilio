#############################################################################
# Copyright (C) 2020-2026 MEmilio
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
from memilio.epidata import getPopulationData as gpd
from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd
from datetime import date, datetime, timedelta
import time
import os
import copy
import pandas as pd
import numpy as np
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

# activate CoW for more predictable behaviour of pandas DataFrames
pd.options.mode.copy_on_write = True


def validate(df_npis_old, df_npis, df_infec_rki, countyID, npiCode,
             start_npi_cols, npi_incid_start, start_date_validation,
             end_date_validation, fine_resolution):
    """ Validates the transformed NPI data based on read in NPI data list.
    Also works for incidence-dependent NPIs as long as no activation or lifting
    delay is used.

    :param df_npis_old: Original data frame.
    :param df_npis: New NPI data frame with potential activation of incidence-
        dependent NPIs.
    :param df_infec_rki: Case data for incidence comparison.
    :param countyID: CountyID of county to be validated.
    :param npiCode: NPI Code of code to be validated.
    :param start_npi_cols: Index of column where NPI information start in
        original data frame.
    :param npi_incid_start: Minimum incidence for activation of NPI codes.
    :param start_date_validation: Start date for validation.
    :param end_date_validation: End date for validation.
    :param npiCode: 
    :param fine_resolution: 

    """

    if fine_resolution == 2:
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
    """ Print download message to ask the user manually download a file.

    :param filename: 
    :param url: 

    """

    gd.default_print("Error",
                     'This script needs manual downloading of files. Please register'
                     ' at corona-datenplatform.com and download ' + filename + ' from ' + url +
                     '. Then move it to a folder named raw_data in this directory.')


def read_files(directory, fine_resolution, run_checks):
    """Reads files from local directory and returns data in dataframes.

    :param directory: Directory where data is loaded from.
    :param fine_resolution: 2 [Default] or 0 or 1. Defines which categories are considered.
        If '2' is set, all the subcategories (~1200) are considered.
        If '1' is set, all incidence levels of subcategories are merged and ~200 NPIs are considered.
        If '0' is chosen only the main, summarizing categories (~20) are used.
    :param run_checks: 
    :returns: Data frames df_npis_old (Decreed, encoded NPIs for all German
        counties) and df_npis_desc (Description of NPIs).
    """
    if fine_resolution > 0:
        try:
            try:
                codelist = [
                    'm01a', 'm01b', 'm02a', 'm02b', 'm03', 'm04', 'm05', 'm06',
                    'm07', 'm08', 'm09', 'm10', 'm11', 'm12', 'm13', 'm14',
                    'm15', 'm16', 'm17', 'm18', 'm19', 'm20', 'm21']
                counter_codes = 0
                for code in codelist:
                    df_npis_per_code = pd.read_csv(
                        os.path.join(directory,
                                     f'kr_massn_unterkat_{code}.csv'),
                        sep=',')

                    # set some parameters for dataframe
                    if counter_codes == 0:
                        counties = np.sort(df_npis_per_code.ags5.unique())
                        num_counties = len(df_npis_per_code.ags5.unique())

                        # extract dates from data
                        dates = df_npis_per_code.iloc[:int(
                            df_npis_per_code.shape[0]/num_counties), 5]
                        # rename dates so that they match dates from other npi dataframe
                        dates_new = [
                            'd' + date.replace('-', '') for date in dates]

                        df_local = [pd.DataFrame()
                                    for i in range(num_counties)]

                    #  set df for all counties
                    for i in range(0, num_counties):
                        if counter_codes == 0:
                            df_local[i] = pd.DataFrame(
                                columns=list(df_npis_per_code.columns[0: 5]) +
                                ['code'] + dates_new)

                        dummy_to_append = pd.DataFrame(
                            columns=['code'] + dates_new,
                            data=copy.deepcopy(df_npis_per_code
                                               [df_npis_per_code.ags5 == counties[i]].
                                               iloc[:, 6:].T.reset_index().values))

                        df_local[i] = pd.concat([df_local[i], dummy_to_append])

                        if df_npis_per_code.iloc[i * len(dates): (i + 1) *
                                                 len(dates),
                                                 3].nunique() > 1:
                            raise gd.DataError(
                                'Dates are not sorted as expected.')

                        # Set first five columns so that they match old format of data frame (from kr_massnahmen_unterkategorien.csv)
                        if counter_codes == len(codelist)-1:
                            df_local[i][df_local[i].columns[0:5]
                                        ] = df_npis_per_code.iloc[i*len(dates), 0:5].values

                    counter_codes += 1
                df_npis_old = pd.concat([df_local[i]
                                        for i in range(num_counties)])
                df_npis_old.rename(dd.GerEng, axis=1, inplace=True)
                df_npis_old['NPI_code'] = df_npis_old['NPI_code'].str.replace(
                    'code_m', 'M')
            except FileNotFoundError:
                df_npis_old = pd.read_csv(
                    os.path.join(
                        directory, 'kr_massnahmen_unterkategorien.csv'),
                    sep=',')
                df_npis_old.rename(dd.GerEng, axis=1, inplace=True)
        except FileNotFoundError:
            print_manual_download(
                'kr_massnahmen_unterkategorien.csv',
                'https://www.corona-datenplattform.de/dataset/massnahmen_unterkategorien_kreise')
            raise FileNotFoundError
        # check if rows hospitals and geriatric care are still empty;
        # these fields have been empty so far and are thus not used
        test_codes = ['M23_010', 'M23_020', 'M23_030', 'M23_040',
                      'M23_050', 'M23_060', 'M24_010', 'M24_020',
                      'M24_030', 'M24_040', 'M24_050', 'M24_060']
        for tcode in test_codes:
            for i in [''] + ["_" + str(i) for i in range(1, 6)]:
                if (df_npis_old[df_npis_old[dd.EngEng['npiCode']] == tcode+i].iloc[:, 6:].max().max() > 0):
                    gd.default_print("Debug", tcode+i + " used.")
        # end check

    else:  # read aggregated NPIs

        try:
            df_npis_old = pd.read_csv(os.path.join(
                directory, 'kr_massnahmen_oberkategorien.csv'))
        except FileNotFoundError:
            print_manual_download(
                'kr_massnahmen_oberkategorien.csv',
                'https://www.corona-datenplattform.de/dataset/massnahmen_oberkategorien_kreise')
            raise FileNotFoundError
        df_npis_old.rename(dd.GerEng, axis=1, inplace=True)

    # read data frame of variable names and descriptions
    try:
        if fine_resolution > 0:
            df_npis_desc = pd.read_excel(
                os.path.join(
                    directory, 'datensatzbeschreibung_massnahmen.xlsx'),
                sheet_name=2, engine=gd.Conf.excel_engine)
        else:
            df_npis_desc = pd.read_excel(
                os.path.join(
                    directory, 'datensatzbeschreibung_massnahmen.xlsx'),
                sheet_name=3, engine=gd.Conf.excel_engine)
    except FileNotFoundError:
        print_manual_download(
            'datensatzbeschreibung_massnahmen.xlsx',
            'https://www.corona-datenplattform.de/dataset/massnahmen_unterkategorien_kreise')
        raise FileNotFoundError

    # download combinations of npis
    try:
        fname = 'combination_npis_incl_ranking.xlsx'
        if fine_resolution > 0:
            df_npis_combinations_pre = pd.read_excel(
                os.path.join(
                    directory, fname), engine=gd.Conf.excel_engine)
    except FileNotFoundError:
        raise FileNotFoundError('File ' + fname + ' not found.')

    if run_checks:
        npi_sanity_check(df_npis_old, df_npis_desc, df_npis_combinations_pre)
    else:
        gd.default_print(
            'Warning', "Sanity checks for NPI data have not been executed.")

    return df_npis_old, df_npis_desc, df_npis_combinations_pre


def activate_npis_based_on_incidence(
        local_incid, npi_lifting_days_threshold, npi_activation_days_threshold,
        incid_threshold):
    """
    Computes an activation vector according to a given incidence threshold,
    observed incidence and activation or lifting delays.

    In order for incidence-dependent NPIs to become active, the incidence
    has to exceed the threshold for npi_activation_days_threshold many days.
    For a formerly active NPI to be lifted, the incidence has to be below
    the threshold for npi_lifting_days_threshold many days.

    If one of the former cases holds true, then the activation or lifting happens
    two days after the satisfaction of the criterion. This is in accordance with
    case reporting that can only happen after the day has finished (a 'delay' of one
    day is introduced here) and as these reports generally appeared in the morning
    for the previous day, the NPI was not directly be activated or lifted that same
    day but only on the next day (another delay of one day). Hence,
    the incidence-dependent NPI is activated or lifted two days after the threshold
    is/is not anymore exceeded, additionally considering the number of consecutive
    days to implement or lift (see second paragraph above).

    Please see the examples for a better understanding.

    Example 1 (Threshold=3.5):
    local_incid=pd.Series([2, 4, 2, 4, 2, 2, 4, 4, 2, 4, 2, 2, 2, 2])
    Yesterdays incidence is over the threshold on following days:
    [?, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0]

    Yesterday for first day is not known. Thus, first day's boolean
    is always set to the same boolean as second day's boolean.
    [0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0]

    With yesterdays incidence over threshold on days:
    [0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0]
    Example 1a) ... and npi_lifting_days_threshold=2, npi_activation_days_threshold=1,
    the NPI will be activated on days 4 and 9 and lifted on days 8 and 14, i.e.,
    int_active then is:
    [0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0]
    Example 1b) ... and npi_lifting_days_threshold=3, npi_activation_days_threshold=2,
    the NPI will be activated on day 10 (and would be lifted on day 15;
    which is not in the vector anymore), i.e., int_active then is:
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]

    Example 2:
    With yesterday's incidence over threshold on days:
    [1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0]
    and npi_lifting_days_threshold=3, npi_activation_days_threshold=1,
    the NPI will be activated on day 2 and lifted on day 14, i.e.,
    int_active then is:
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]

    Please also note that the first column will always returned as false
    so the dataframe should not start with dates where NPIs are implemented.
    For the Corona Datenplattform frame which starts from 2020-03-01
    this is no problem for the first days as there were no NPIs.

    :param local_incid: Incidence observed in local region.
    :param npi_lifting_days_threshold: Number of days for the incidence to be
        lower than threshold to lift the NPI.
    :param npi_activation_days_threshold: Number of days the incidence threshold
        is exceeded before activation of the NPI.
    :param incid_threshold: Threshold to be considered.

    """

    if npi_lifting_days_threshold < 1 or npi_activation_days_threshold < 1:
        raise ValueError(
            'Activation or lifting day variables need to be 1 or larger')

    # First get a Series with 0 for yesterdays incidence
    # is below threshold and 1 for incidence over threshold
    yesterdays_incid_over_threshold = (local_incid.shift(
        1).fillna(local_incid[0]) > incid_threshold).astype(int)

    # get a zero filled Series with same length to be
    # filled with ones where NPI is active
    int_active = pd.Series(np.zeros(len(local_incid), dtype=int))
    # loop over every day
    for i in range(len(yesterdays_incid_over_threshold)):
        # Set int_active=0 where last npi_lifting_days_threshold+1 days did not exceed
        # the threshold. Look only until the day before yesterday (max(...):i is exclusive sum for i)
        # as we assume a necessary delay of 24h to implement an intervention (see explanation above)
        if yesterdays_incid_over_threshold[max(0, i-npi_lifting_days_threshold):i].values.sum() == 0:
            int_active[i] = 0
        # Set int_active=1 where last npi_activation_days_threshold+1 days did
        # all exceed the threshold
        elif yesterdays_incid_over_threshold[max(0, i-npi_activation_days_threshold):i].values.sum() == npi_activation_days_threshold:
            int_active[i] = 1
        # If no condition applies, set int_active to the value of the previous day
        # for i=0, int_active always will be zero (see comment above)
        elif i > 0:
            int_active[i] = int_active[i-1]
        # elif i==0 int active is 0

    return int_active


def drop_codes_and_categories(
        npi_codes_prior, npi_codes_prior_desc, df_npis_old, fine_resolution):
    """Drops codes and categories from original data frame if they are not used.

    :param npi_codes_prior: NPI codes read from description sheet.
    :param npi_codes_prior_desc: NPI code descriptions read from description sheet.
    :param df_npis_old: Original data frame with encoding (decreed, yes or no)
        for different NPIs.
    :param fine_resolution: 2 [Default] or 0 or 1. Defines which categories are considered.
        If '2' is set, all the subcategories (~1200) are considered.
        If '1' is set, all incidence levels of subcategories are merged and ~200 NPIs are considered.
        If '0' is chosen only the main, summarizing categories (~20) are used.
    :returns: Returns dropped codes, prior codes and reduced original data frame.
    """
    if fine_resolution > 0:

        for i in range(1, 6):
            # correct M04_N codes to M_04_M_N, N in {1,...,5}, M in {120,110,100,130,140}
            # (M04_1, i.e. i=1, has been corrected in original file but not for i>1)
            if i != 1:
                npi_codes_prior[npi_codes_prior == 'M04_'+str(i)] = ['M04_120_'+str(
                    i), 'M04_110_'+str(i), 'M04_100_'+str(i), 'M04_130_'+str(i), 'M04_140_'+str(i)]
        # correct M05_N codes to M_05_M_N, N in {1,...,5}, M in {130,150,120,140,110,100,160}
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
        codes_dropped = []

    return codes_dropped, npi_codes_prior, df_npis_old


def npi_sanity_check(df_npis_old, df_npis_desc, df_npis_combinations_pre):
    """

    :param df_npis_old: 
    :param df_npis_desc: 
    :param df_npis_combinations_pre: 

    """
    # Check if all counties are in df_npis_old
    if not np.array_equal(df_npis_old.ID_County.unique().astype(int), np.array(geoger.get_county_ids(merge_eisenach=False)).astype(int)):
        raise gd.DataError('Not all counties found in DataFrame.')
    # Check if all NPIs are in df_npis_old
    if len(df_npis_old.NPI_code.unique()) != 1152:
        raise gd.DataError('Missing NPIs in DataFrame.')
    # check columns of df_npis_old (6 columns with info, 883 dates until 2022/07/31)
    if len(df_npis_old.columns) != 889:
        raise gd.DataError('Unexpected length of DataFrame.')
    # check if Variablenname, Variable and Beschreibung are in df_npis_desc columns
    if not ('Variablenname' in df_npis_desc.columns):
        raise gd.DataError('Column Variablenname not found.')
    if not ('Variable' in df_npis_desc.columns):
        raise gd.DataError('Column Variable not found.')
    if not ('Beschreibung' in df_npis_desc.columns):
        raise gd.DataError('Column Beschreibung not found.')
    # df_npis_desc should have 1224 rows
    if len(df_npis_desc) != 1224:
        raise gd.DataError('Unexpected length of description DataFrame.')
    # df_npis_combinations_pre should habe 204 rows (1224/6)
    if len(df_npis_combinations_pre) != 204:
        raise gd.DataError('Unexpected length of combination DataFrame.')
    # combination part should have values NaN and x
    for column in df_npis_combinations_pre.columns[5:]:
        if not set(df_npis_combinations_pre[column].unique().astype(str)).issubset({'nan', 'x'}):
            raise gd.DataError('Unexpected values in combination matrix.')


def get_npi_data(fine_resolution=2,
                 file_format=dd.defaultDict['file_format'],
                 out_folder=dd.defaultDict['out_folder'],
                 start_date=dd.defaultDict['start_date'],
                 end_date=dd.defaultDict['end_date'],
                 counties_considered=geoger.get_county_ids(),
                 npi_activation_days_threshold=3,
                 npi_lifting_days_threshold=5,
                 **kwargs
                 ):
    """Loads a certain resolution of recorded NPI data from the Corona Datenplattform.

    Extracts the counties asked for and activates the NPIs if they are incidence dependent.

    Results' data frames will be stored in the directory as:

    - fine_resolution=2: germany_counties_npi_subcat
    - fine_resolution=1: germany_counties_npi_subcat_incgrouped
    - fine_resolution=0: germany_counties_npi_maincat

    Needs the files 'cases_all_county_all_dates_repdate.json' and
    'county_current_population.json' which can be created by the functions
    getCasesData.py (with argument --rep-date) and getPopulationData.py.

    Please manually download kr_massnahmen_unterkategorien.csv and
    datensatzbeschreibung_massnahmen.xlsx from
    https://www.corona-datenplattform.de/dataset/massnahmen_unterkategorien_kreise
    and kr_massnahmen_oberkategorien.csv from
    https://www.corona-datenplattform.de/dataset/massnahmen_oberkategorien_kreise
    and move it to the directory-path mentioned in the beginning of the function.

    :param fine_resolution: 2 [Default] or 0 or 1. Defines which categories are considered.
        If '2' is set, all the subcategories (~1200) are considered.
        If '1' is set, all incidence levels of subcategories are merged and ~200 NPIs are considered.
        If '0' is chosen only the main, summarizing categories (~20) are used.
    :param file_format: File format which is used for writing the data. Default defined in defaultDict.
    :param out_folder: Path to folder where data is written in folder
        out_folder/Germany. (Default value = dd.defaultDict['out_folder'])
    :param start_date: Start date of stored data frames. Default = '', taken from read data.
    :param end_date: End date of stored data frames. Default = '', taken from read data.
    :param counties_considered: Either 'All' or a list of county IDs from 1001 to 16xxx. Default: 'All'.
    :param npi_activation_days_threshold: Defines necessary number of days exceeding case
        incidence threshold to activate NPIs. Default: 3.
    :param npi_lifting_days_threshold: Defines necessary number of days below case incidence
        threshold to lift NPIs. (Default value = 5)
    :param kwargs: Additional keyword arguments.
    """
    conf = gd.Conf(out_folder, **kwargs)
    out_folder = conf.path_to_use

    # Depending on the federal state and time period, there are
    # huge differences for number of days before the lifting and activation.
    # It was usually between 1 and 14 days. We use npi_lifting_days_threshold = 5
    # and npi_activation_days_threshold = 3 as default averaged value.

    if counties_considered == 'All':
        counties_considered = geoger.get_county_ids()
    try:
        counties_considered.remove(16056)
    except ValueError:
        pass

    # if only one county is considered, it should be a list with one element
    if not isinstance(counties_considered, list):
        counties_considered = [counties_considered]

    directory = out_folder
    directory = os.path.join(directory, 'Germany', 'pydata')
    gd.check_dir(directory)

    # read manual downloaded files from directory
    df_npis_old, df_npis_desc, df_npis_combinations_pre = read_files(
        directory, fine_resolution, conf.checks)

    gd.default_print('Debug', 'Download completed.')

    # Compute column index of NPI start (columns with NPIs start with days
    # which are provided in format dYYYYMMDD).
    npi_start_col = np.where(
        df_npis_old.columns.str.contains('d2') == True)[0][0]

    # get existing codes that are used; for fine resolution we don't
    # have codes M22 - M24 but these are still listed in description
    if fine_resolution > 0:
        # count how many codes contain M22, M23 or M24
        num_nonexistent_codes = df_npis_desc['Variablenname'].str.count(
            "M22|M23|M24").sum()
        # do not include these nonexistent codes
        if num_nonexistent_codes != 0:
            npi_codes_prior = df_npis_desc['Variablenname'].iloc[: -
                                                                 num_nonexistent_codes]
            npi_codes_prior_desc = df_npis_desc['Variable'].iloc[: -
                                                                 num_nonexistent_codes]
        else:
            npi_codes_prior = df_npis_desc['Variablenname']
            npi_codes_prior_desc = df_npis_desc['Variable']
    # for fine_resolution = 0 df_npis_old M22-M24 are empty)
    else:
        npi_codes_prior = df_npis_desc['Variablenname']
        npi_codes_prior_desc = df_npis_desc['Variable']

    # For fine_resolution > 0 deactivation of non-combinable / conflicting
    # NPIs has to be conducted.
    #
    # NPIs of different main categories (e.g., M01a and M04) can always be
    # prescribed together as they target different locations and settings.
    #
    # NPIs with the same main code (i.e., targeting the same location, e.g.,
    # schools, or the same set of NPIs, e.g., masks) can exclude each other.
    # Exclusion happens based on table provided in xlsx or csv format.
    #
    # In first place, NPIs of higher stringency index as defined by the Corona-
    # Datenplattform deactivate NPIs with lower stringency index.
    # NPIs of the same main code and with the same stringency index may or
    # may not exclude each other according to the threshold they were
    # prescribed with. Prescribed and active NPIs with high incidence thresholds
    # deactivate conflicting NPIs with lower thresholds (as the latter are
    # considered to be less strict).
    if fine_resolution > 0:
        num_nonexistent_codes_pre = df_npis_combinations_pre['Variablenname'].str.count(
            "M22|M23|M24").sum()
        if num_nonexistent_codes_pre != 0:
            df_npis_combinations_pre = df_npis_combinations_pre.iloc[: -
                                                                     num_nonexistent_codes_pre, :]

        # drop 0 column if existent
        try:
            df_npis_combinations_pre.drop(columns=0, inplace=True)
        except KeyError:
            pass
        # rename essential columns and throw away others
        columns_used = np.where(
            (df_npis_combinations_pre == 'x').any() == True)[0]
        column_names = list(
            df_npis_combinations_pre.columns[[i for i in columns_used]])
        rename_columns = {column_names[i]: i for i in range(len(column_names))}
        df_npis_combinations_pre.rename(columns=rename_columns, inplace=True)
        df_npis_combinations_pre = df_npis_combinations_pre[[
            'Variablenname', 'Massnahmenindex'] + [i for i in range(0, len(columns_used))]]
        # replace empty cells by zeros and x-marked cells by ones
        # This has to be done by replacing the values with the same dtype and then changing the dtype
        # Pandas 3.0 will not allow downcasting with replace operations
        df_npis_combinations_pre = df_npis_combinations_pre.replace(
            np.nan, '0')
        df_npis_combinations_pre = df_npis_combinations_pre.replace('x', '1')
        df_npis_combinations_pre[df_npis_combinations_pre.columns[2:]
                                 ] = df_npis_combinations_pre[df_npis_combinations_pre.columns[2:]].astype(int)

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

        # create hash table of main code to strictness rankings inside main
        # code and combination matrix inside the same strictness rank
        df_npis_combinations = {
            npi_groups_combinations_unique[i]:
            [
                {df_npis_combinations_pre['Variablenname'][npi_groups_idx[i]].values[j]:
                 df_npis_combinations_pre['Massnahmenindex'][npi_groups_idx[i]].values[j] for j in range(
                    len(df_npis_combinations_pre['Variablenname'][npi_groups_idx[i]]))},
                np.eye(len(npi_groups_idx[i]))]
            for i in range(len(npi_groups_combinations_unique))}

        # run through all groups and set possible combinations according to
        # read combination matrix
        # find begin of combination matrix (find first '0' column)
        start_comb_matrix = list(df_npis_combinations_pre.columns).index(0)
        for i in range(len(npi_groups_idx)):
            codes_local = df_npis_combinations_pre.loc[npi_groups_idx[i],
                                                       'Variablenname'].values
            npic_uniq = npi_groups_combinations_unique[i]  # reduce code length
            df_npis_combinations[npic_uniq][1] = df_npis_combinations_pre.iloc[np.array(npi_groups_idx[i]),
                                                                               start_comb_matrix:start_comb_matrix+len(npi_groups_idx[i])].values
            if (df_npis_combinations[npic_uniq][1]-np.transpose(df_npis_combinations[npic_uniq][1])).max() > 0:
                gd.default_print(
                    'Error', 'Input file does not match with data. Please correct combination matrix input.')
            # make it a dataframe to allow easy removal of code lines and rows
            # if they are not used later on
            df_npis_combinations[npic_uniq][1] = pd.DataFrame(
                df_npis_combinations[npic_uniq][1],
                columns=codes_local)
            df_npis_combinations[npic_uniq][1].insert(
                0, 'Code', codes_local)

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
            df_out = copy.deepcopy(df_npis_combinations[npi_groups_combinations_unique[i]][
                1])
            df_out.insert(
                0, 'Description (German)',
                [desc
                 for desc in npi_codes_prior_desc
                 [npi_codes_prior.isin(codes_local)].values])
            # validate if combinations cleanout file exists, else write this file
            if write_file == False:
                # store verified output
                df_in_valid = pd.read_excel(
                    os.path.join(
                        directory, 'combinations_npis_cleanoutput.xlsx'),
                    sheet_name=i, engine=gd.Conf.excel_engine)
                if not df_in_valid.drop(columns='Unnamed: 0').equals(df_out):
                    gd.default_print('Error', 'Error in combination matrix.')
                del df_in_valid
            else:
                df_out.to_excel(
                    writer, sheet_name=npi_groups_combinations_unique[i])
            del df_out
        if write_file:
            writer.close()

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
            gd.default_print(
                "Error", "Additional errors in consistent naming.")
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
                npis['NPI_code'])
            local_codes_used_cols = df_npis_combinations[code][1].columns.isin(
                npis['NPI_code'])

            # remove strictness indices of unused codes
            df_npis_combinations[code][0] = {
                key: val for key, val in df_npis_combinations[code][0].items()
                if key in npis['NPI_code'].values}
            # remove columns of combinations
            df_npis_combinations[code][1] = df_npis_combinations[code][1].loc[local_codes_used_rows,
                                                                              local_codes_used_cols].reset_index(drop=True)

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

    # Remove counties which are not considered. Check if all considered counties are in the dataframe
    counties_removed = df_npis_old[
        ~df_npis_old[dd.EngEng['idCounty']].isin(counties_considered)][
        dd.EngEng['idCounty']].unique()

    if set(counties_considered).difference(counties_removed) == set(
            counties_considered) and np.array_equal(
            sorted(np.append(counties_considered, counties_removed)),
            sorted(df_npis_old[dd.EngEng['idCounty']].unique())):
        pass
    else:
        raise gd.DataError('Error. Considered counties hae been removed.')
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
        gd.default_print("Error", "Dates missing in data frame:")
        for i in date_diff_idx:
            gd.default_print("Debug",
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
        try:
            df_population = pd.read_json(
                directory + "county_current_population.json")
        except FileNotFoundError:
            df_population = gpd.get_population_data()
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

    # iterate over countyIDs
    counters = np.zeros(4)  # time counter for output only
    countyidx = 0

    # Infer type of columns to be able to use replace with ints without downcasting.
    df_npis_old = df_npis_old.infer_objects()
    df_npis_old.replace([-99, 2, 3, 4, 5],
                        [0, 1, 1, 1, 1], inplace=True)

    counter_cases_start = 0

    # setup dataframe for each maingroup, same format as df_npi_combinations
    # used to count codes that occur simultaneously now (before any (de-)activation)
    df_count_joint_codes = copy.deepcopy(df_npis_combinations)
    for maincode in df_count_joint_codes.keys():
        df_count_joint_codes[maincode][1] *= 0
    df_counted_joint_codes = count_code_multiplicities(df_npis_old, df_count_joint_codes,
                                                       counties_considered=counties_considered)
    save_interaction_matrix(df_counted_joint_codes, 'joint_codes', directory)
    plot_interaction_matrix('joint_codes', directory)

    # create dataframe to count multiple codes after incidence dependent (de-)activation
    df_incid_depend = pd.DataFrame()
    df_count_incid_depend = copy.deepcopy(df_npis_combinations)
    for maincode in df_count_incid_depend.keys():
        df_count_incid_depend[maincode][1] *= 0

    # create dataframe to count multiple codes after strictness deactivation
    df_count_active = copy.deepcopy(df_npis_combinations)
    for maincode in df_count_active.keys():
        df_count_active[maincode][1] *= 0

    # setup dataframe for each maingroup, same format as df_npi_combinations
    # used to count number of codes that are deactivated
    df_count_deactivation = copy.deepcopy(df_npis_combinations)
    for maincode in df_count_deactivation.keys():
        df_count_deactivation[maincode][1] *= 0

    all_subcodes = []
    for maincode in df_npis_combinations.keys():
        all_subcodes += df_npis_combinations[maincode][1].columns.to_list()
        # check (and validate) that element 0 and 1 in df_npis_combination match.
        if df_npis_combinations[maincode][1].columns.to_list() != list(
                df_npis_combinations[maincode][0].keys()):
            raise gd.DataError('Error. Description and table do not match.')

    # create new data frame for all NPIs
    df_npis = pd.DataFrame()

    for countyID in counties_considered:
        cid = 0
        countyidx += 1

        # compute incidence for given county and store in other data frame
        if fine_resolution > 0:
            # compute incidence based on previous data frames
            df_infec_local = copy.deepcopy(
                df_infec_rki[df_infec_rki[dd.EngEng['idCounty']] == countyID])
            pop_local = df_population.loc[df_population[dd.EngEng['idCounty']]
                                          == countyID, dd.EngEng['population']].values[0]

            # consider difference between current day and day-7 to compute incidence
            # As a helper, repeat first entry seven times, incidence then always starts with 0.
            cases_first_value = df_infec_local[dd.EngEng['confirmed']].values[0]
            df_infec_local_repeat_first_entry = [
                cases_first_value for i in range(7)] + list(
                df_infec_local[dd.EngEng['confirmed']].values.transpose())

            df_infec_local['Incidence'] = (pd.Series(
                df_infec_local_repeat_first_entry).diff(periods=7) /
                pop_local * 100000)[7:].values

            # set to main data frame
            df_infec_rki.loc[df_infec_rki[dd.EngEng['idCounty']] ==
                             countyID, 'Incidence'] = df_infec_local['Incidence'].values

            # cut infection information at start_date_new and end_date_new
            df_infec_local = df_infec_local[(df_infec_local[dd.EngEng['date']] >= start_date_new) & (
                df_infec_local[dd.EngEng['date']] <= end_date_new)].reset_index()

            # Count counties with start cases >= 1:
            # In this case NPI activation cannot be ensured to work as expected
            if cases_first_value >= 1:
                counter_cases_start += 1

        # get county-local data frame
        start_time = time.perf_counter()
        df_local_old = copy.deepcopy(df_npis_old[df_npis_old[dd.EngEng['idCounty']]
                                                 == countyID])

        # get number of codes of one NPI (incidence indep. + dep.)
        # for fine_resolution=1, inc_codes=1, for fine_res=2, inc_codes=6
        inc_codes = len(np.where(npis_final.NPI_code.str.contains(
            npis[dd.EngEng['npiCode']][0]))[0])

        # Consistency of incidence independent and dependent NPIs:
        # The same NPI should not be prescribed multiple times at the same day
        # for different incidence-dependent thresholds or incidence-independently.
        # In order to avoid contradictions, only retain the strictest mentioned
        # implementation. Incidence-independent is always stricter than any
        # incidence-dependent implementation.
        for i in range(int(len(df_local_old)/inc_codes)):

            # check if access is correct
            if not all(
                [npis[dd.EngEng['npiCode']][i * inc_codes] in npi_code_test
                 for npi_code_test in df_local_old.iloc
                 [inc_codes * i: inc_codes * (i + 1),
                  npi_start_col - 1].to_list()]):
                raise gd.DataError('Wrong NPI rows aggregated.')

            sum_npi_inc = np.where(
                df_local_old.iloc[inc_codes*i:inc_codes*(i+1), npi_start_col:].sum() > 1)
            if (len(sum_npi_inc[0]) > 0):
                gd.default_print("Trace",
                                 'Reduce multiple prescription in county ' + str(countyID) +
                                 ' for NPI ' + str(npis.loc[inc_codes*i, 'Description']))
                for j in sum_npi_inc[0]:
                    # get lowest index (i.e., strictest implementation of NPI).
                    idx_start = np.where(
                        df_local_old.iloc[inc_codes*i:inc_codes*(i+1), npi_start_col+j])[0].min()
                    # Remove less strict and thus contradictory
                    # implementations of the same NPI the same day.
                    df_local_old.iloc[inc_codes*i+idx_start +
                                      1:inc_codes*(i+1), npi_start_col+j] = 0

                if not all(
                    df_local_old.iloc
                    [inc_codes * i: inc_codes * (i + 1),
                     npi_start_col + sum_npi_inc[0]].sum() == 1):
                    raise gd.DataError('Consistency correction failed.')

        ## end of consistency correction ##

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
        df_local_new = copy.deepcopy(df_local_old.iloc[npi_rows, start_npi_cols-1:].set_index(
            dd.EngEng['npiCode']).transpose())
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
            except KeyError:
                pass
            # get index of first NPI column in local data frame
            npis_idx_start = list(
                df_local_new.columns).index(
                npis[dd.EngEng['npiCode']][0])

            # extract local incidence from local frame
            local_incid = copy.deepcopy(df_infec_local['Incidence'])

            # iterate through all NPIs and activate if incidence threshold
            # is exceeded
            for level, npi_indices in incidence_thresholds_to_npis.items():
                if level[0] >= 0:  # level[0] = incidvalthrsh

                    # get days where npis are active as int (1/0)
                    int_active = activate_npis_based_on_incidence(
                        local_incid, npi_lifting_days_threshold,
                        npi_activation_days_threshold, level[0])

                    # multiply rows of data frame by either 1 if threshold
                    # passed (i.e., mentioned NPI is active) or zero
                    # (i.e., mentioned NPI is not active)
                    # 'mul' multiplies the original data frame row by row
                    # with the respective value in int_active
                    df_local_new.iloc[:, npis_idx_start + np.array(npi_indices)] \
                        = df_local_new.iloc[:, npis_idx_start + np.array(npi_indices)].mul(int_active, axis=0)

            # merge incidence dependent NPIs to have only one column for each subcode
            df_local_new_merged = copy.deepcopy(df_local_new.iloc[:, :2])
            for subcode in all_subcodes:
                # extract columns which have the subcode as part of the column
                # name and sum over all these subcodes
                df_local_new_merged[subcode] = df_local_new.filter(
                    regex=subcode).sum(axis=1)
            # strictness deactivation is done with this merged dataframe

            df_incid_depend = pd.concat(
                [df_incid_depend, copy.deepcopy(df_local_new_merged)])

            if df_local_new_merged.iloc[:, 2:].max().max() > 1:
                raise gd.DataError('Error in merging...')

            # Remove conflicting NPIs according to strictness index of Corona-
            # Datenplattform and exclusion criteria defined in df_npis_combinations
            for maincode in df_npis_combinations.keys():
                # get all subcodes
                subcodes = list(df_npis_combinations[maincode][0].keys())
                subcodes_strictness_values = list(
                    df_npis_combinations[maincode][0].values())
                if len(subcodes) != len(subcodes_strictness_values):
                    raise gd.DataError(
                        'Subcode and strictness array inconsistent.')
                # sort index reversely with the strictest (highest) index first
                idx_strictness_sorted_rev = np.argsort(
                    subcodes_strictness_values)[::-1]
                for i in range(len(idx_strictness_sorted_rev)-1):
                    # get index of NPI of a certain strictness
                    idx_strictness = idx_strictness_sorted_rev[i]
                    # get code of corresponding NPI
                    subcode = subcodes[idx_strictness]

                    # get indices of days where subcode is active
                    subcode_active = np.where(
                        df_local_new_merged.loc[:, subcode] > 0)[0]

                    if len(subcode_active) > 0:
                        # get indices of less strict NPIs
                        codes_less_strict = df_npis_combinations[maincode][1].columns[np.sort(
                            idx_strictness_sorted_rev[i+1:])]

                        subcodes_nocombi = df_npis_combinations[maincode][1].loc[idx_strictness, :]
                        # only consider those codes which cannot be combined;
                        # for these, values of 1 have to be set to 0
                        subcodes_nocombi = list(
                            subcodes_nocombi[subcodes_nocombi == 0].index)

                        # intersect non-combinable subcodes with less strict subcodes
                        subcodes_deactivation = np.sort(
                            list(set(codes_less_strict).intersection(subcodes_nocombi)))

                        for nocombi_code in subcodes_deactivation:
                            # check where the less strict NPI is mentioned, only
                            # considering rows where the stricter NPI is mentioned.
                            days_deact = np.where(
                                df_local_new_merged.loc[subcode_active, nocombi_code] > 0)[0]
                            if len(days_deact) > 0:
                                gd.default_print("Trace", 'Deactivating for ' +
                                                 'County ' + str(countyID)+'\t' + str(nocombi_code) + ' due to ' +
                                                 str(subcode) + ' on ' + str(len(days_deact)) + ' days.\n')
                                # take subcode_active rows as days_deact is
                                # numbering inside subcode_active rows only,
                                # not numbering on the whole df_local_new_merged
                                # data frame
                                df_local_new_merged.loc[subcode_active,
                                                        nocombi_code] = 0
                                df_count_deactivation[maincode][1].loc[idx_strictness,
                                                                       nocombi_code] += len(days_deact)

            # count joint codes from after strictness based deactivation
            df_count_active = count_code_multiplicities(
                df_local_new_merged, df_count_active, [countyID], False)

            # count joint codes from after incidence based activation
            df_count_incid_depend = count_code_multiplicities(
                df_incid_depend, df_count_incid_depend, [countyID], False)

            # for fine resolution = 1 only consider merged dataframe
            if fine_resolution == 1:
                df_local_new = copy.deepcopy(df_local_new_merged)
            else:
                # multiply subcode columns with incidence dependent subcode columns in df_local_new
                for subcode in all_subcodes:
                    for incidcode in ['', '_1', '_2', '_3', '_4', '_5']:
                        df_local_new[subcode +
                                     incidcode] *= df_local_new_merged[subcode]

        counters[cid] += time.perf_counter()-start_time
        cid += 1
        ### ###

        start_time = time.perf_counter()

        df_npis = copy.deepcopy(pd.concat(
            [df_npis, df_local_new],
            ignore_index=True))
        counters[cid] += time.perf_counter()-start_time
        cid += 1

        # divide working time by completed number of counties and multiply
        # by remaining number of counties to estimate time remaining
        time_remain = sum(
            counters) / countyidx * (len(counties_considered) - countyidx)
        # print progress
        if countyidx == 1 or countyidx % int(
                len(counties_considered) / 10) == 0:
            gd.default_print('Debug', 'Progress ' + str(countyidx) + ' / ' +
                             str(len(counties_considered)) +
                             '. Estimated time remaining: ' +
                             str(int(time_remain / 60)) + ' min.')

    save_interaction_matrix(df_count_deactivation,
                            'count_deactivation', directory)
    plot_interaction_matrix('count_deactivation', directory)

    if counter_cases_start >= len(counties_considered)*0.05:
        gd.default_print('Warning', 'DataFrame starts with reported cases > 0 '
                         'for more than 5 percent of the counties to be considered. '
                         'In this case, incidence computation and activation of '
                         'incidence-dependent NPIs cannot be ensured to work correctly. '
                         'Please consider a start date of some weeks ahead of the '
                         'time window to be analyzed for NPI\'s effects.')

    save_interaction_matrix(df_count_incid_depend,
                            'joint_codes_incid_depend', directory)
    plot_interaction_matrix('joint_codes_incid_depend', directory)

    save_interaction_matrix(df_count_active, 'joint_codes_active', directory)
    plot_interaction_matrix('joint_codes_active', directory)

    # print sub counters
    gd.default_print('Debug', 'Sub task counters are: '+str(counters))

    # reset index and drop old index column
    df_npis.reset_index(inplace=True)
    try:
        df_npis = df_npis.drop(columns='index')
    except KeyError:
        pass
    try:
        df_npis = df_npis.drop(columns='level_0')
    except KeyError:
        pass

    #### start validation ####
    if fine_resolution > 0 and npi_activation_days_threshold == 1 and npi_lifting_days_threshold == 1:

        for countyID in counties_considered:
            for npiCode in [
                'M01a_010', 'M01a_150', 'M05_120', 'M01a_010',
                    'M18_030', 'M01b_020', 'M02b_035', 'M16_050']:
                [
                    a, b, oldf, newf] = validate(
                    df_npis_old, df_npis, df_infec_rki, countyID,
                    npiCode, start_npi_cols, npi_incid_start,
                    start_date_new, end_date_new,
                    fine_resolution)
                if (a != b):
                    gd.default_print('Error', 'Error in NPI activation computation' +
                                     str(a) + str(b) + str(a - b))

    #### end validation ####

    if fine_resolution > 0:
        if fine_resolution == 1:
            filename = 'germany_counties_npi_subcat_incgrouped'
        else:
            filename = 'germany_counties_npi_subcat'
    else:
        filename = 'germany_counties_npi_maincat'
    gd.write_dataframe(df_npis, directory, filename, file_format)

    return df_npis


def count_code_multiplicities(df_npis_input, df_count, counties_considered, initial_data_frame=True):
    """ Count for all pairs of NPI codes how many times they were
    mentioned at the same day in the initial data frame.

    :param df_npis_input: 
    :param df_count: 
    :param counties_considered: 
    :param initial_data_frame:  (Default value = True)

    """
    for county in counties_considered:
        df_local = df_npis_input[df_npis_input[dd.EngEng['idCounty']] == county]
        if initial_data_frame:
            # get column where dates start
            npi_start_col = np.where(
                df_local.columns.str.startswith('d2') == True)[0][0]
        # prepare dictionnary for dates when code was mentioned
        code_dates = {}
        # run through all maincodes (i.e., first 3-4 characters like M01a or M11)
        # diagonal entries
        for maincode in df_count.keys():
            code_list = df_count[maincode][1].columns
            # iterate over code/row indices 0 to n
            for code_idx in range(len(code_list)):

                # initial data frame (df_npis_old) and reworked new data frames
                # are transposed (NPIs and dates in rows and columns switched)
                if initial_data_frame:
                    # get dates where NPI is mentioned as existing in potential intervention set
                    npi_rows = df_local.NPI_code.str.contains(
                        code_list[code_idx])
                    npi_dates_in_df = np.where(
                        df_local[npi_rows].iloc[:, npi_start_col:].max() > 0)[0]
                    # store non-transformed dates in code_dict
                    code_dates[code_list[code_idx]] = df_local.iloc[:,
                                                                    npi_start_col + npi_dates_in_df].columns
                    # count number of multiply mentionned NPIs with different incidence thresholds for the same day
                    df_count[maincode][1].iloc[code_idx, code_idx] += df_local[npi_rows].iloc[:,
                                                                                              npi_start_col + npi_dates_in_df].sum().sum() - len(npi_dates_in_df)
                else:
                    # get dates where NPI is mentioned as existing in potential intervention set
                    npi_cols = df_local.columns.str.contains(
                        code_list[code_idx])
                    npi_dates_in_df = np.where(
                        df_local.loc[:, npi_cols].max(axis=1) > 0)[0]
                    # store transformed dates in code_dict
                    code_dates[code_list[code_idx]
                               ] = df_local.iloc[npi_dates_in_df, 0].to_list()

        # offdiagonal entries (as before, use that code_dates has been filled for all diagonal entries, i.e., all codes)
        for maincode in df_count.keys():
            code_list = df_count[maincode][1].columns
            # iterate over rows in matrix df_count with code/row indices 0 to n
            for code_idx in range(len(code_list)):
                # iterate over code/column indices 0 to code_idx-1 (not filling diagonal)
                # Note that the upper diagonal part of the matrix does not
                # need to be considered as matrix is symmetric.
                for code_idx_other in range(code_idx):
                    df_count[maincode][1].iloc[code_idx, code_idx_other] += len(set(
                        code_dates[code_list[code_idx]]).intersection(set(code_dates[code_list[code_idx_other]])))

    return df_count


def save_interaction_matrix(df_interactions, filename, directory):
    """ Saves interaction matrices for all subcodes in provided main codes.

    :param df_interactions: 
    :param filename: 
    :param directory: 

    """

    writer = pd.ExcelWriter(
        os.path.join(directory, filename + '.xlsx'),
        engine='xlsxwriter')
    for code in df_interactions.keys():
        df_interactions[code][1].to_excel(writer, sheet_name=code)
    writer.close()

# saves plot in folder directory/heatmaps_filename


def plot_interaction_matrix(filename, directory):
    """ Reads interaction matrices from hard drive and writes heatmap plots
         to hard drive. Separates diagonal and offdiagonal entries as
         interactions inside one NPI are counted for all incidence dependent
         sublevels while between NPIs only one interaction is counted if more
         than one sublevel is mentioned on each of the sides.

    :param filename: 
    :param directory: 

    """
    target_directory = os.path.join(directory, 'heatmaps_' + filename)
    if not os.path.exists(target_directory):
        os.makedirs(target_directory)

    try:
        codelist = pd.ExcelFile(os.path.join(
            directory, filename + '.xlsx'), engine=gd.Conf.excel_engine).sheet_names
    except FileNotFoundError:
        raise FileNotFoundError('File ' + filename + ' not found.')

    # invert color map elements for tab20c such that subcolors are shown
    # from light to dark
    cmap = copy.copy(plt.get_cmap('tab20b'))
    colors = [
        cmap(i)
        for i in np.array(
            [list(range(4 * (i + 1) - 1, 4 * i - 1, -1)) for i in
             range(5)]).flatten()]
    colors = colors + [(0.6, 0.6, 0.6), (0.4, 0.4, 0.4),
                       (0.2, 0.2, 0.2), (0, 0, 0)]
    cmap = mpl.colors.ListedColormap(colors)

    for code in codelist:
        df = pd.read_excel(
            os.path.join(directory, filename + '.xlsx'),
            sheet_name=code, engine=gd.Conf.excel_engine)

        # remove first column and convert to numpy array
        array_exclusion = df.iloc[:, 1:].to_numpy()

        # separate diag and offdiag
        array_exclusion_diag = copy.deepcopy(array_exclusion.diagonal())
        # set diag = 0
        for i in range(array_exclusion.shape[0]):
            array_exclusion[i, i] = 0

        if filename != 'count_deactivation':
            # for count deactivation xlabel != ylabel
            # else matrix is of squared form and symmetric
            array_exclusion += np.transpose(array_exclusion)
        positions = [i for i in range(len(df.columns)-1)]
        plt.xticks(positions, df.columns.to_list()[1:], rotation='vertical')
        plt.yticks(positions, df.columns.to_list()[1:])

        # use different labels and title for each filename
        if filename == 'count_deactivation':
            plt.xlabel('Second NPI')
            plt.ylabel('First NPI')
            plt.title('NPI deactivations')
        elif filename == 'joint_codes_incid_depend':
            plt.xlabel('NPI')
            plt.ylabel('NPI')
            plt.title('Joint NPI prescriptions (incidence dependent)')
        elif filename == 'joint_codes_active':
            plt.xlabel('NPI')
            plt.ylabel('NPI')
            plt.title('Joint NPI implementations')
        elif filename == 'joint_codes':
            plt.xlabel('NPI')
            plt.ylabel('NPI')
            plt.title('Joint NPI prescriptions')
        else:
            raise gd.DataError('Unknown filename: ' + filename)

        # plot offdiagonal (interactions between NPIs)
        # Set vmin = 1 so that only combinations that are simultaneously active
        # at least on one day are in color, else use white.
        # Set vmax = 1e6 to be adjusted with colormap, this value is larger
        # than the maximum in all dataframes, this way colors of heatmaps are
        # comparable across different visualizations
        # (e.g. between codes or between joint_codes and exclusions)
        plt.imshow(array_exclusion, cmap=cmap,
                   norm=mpl.colors.LogNorm(vmin=1, vmax=1e6))
        plt.colorbar()
        plt.tight_layout()
        plt.savefig(
            os.path.join(target_directory, 'InterNPIs_' + filename + '_{}'.format(
                code)), dpi=300)
        plt.close()

        # plot diagonal (interactions between incidence levels of one NPIs)
        plt.figure()
        positions = [i for i in range(len(df.columns)-1)]
        plt.yticks(positions, df.columns.to_list()[1:])
        plt.xticks([])
        plt.imshow(np.array([array_exclusion_diag]).T,
                   cmap=cmap, norm=mpl.colors.LogNorm(vmin=1, vmax=1e6))
        plt.colorbar()
        plt.title('Intra-NPI duplicates')
        plt.tight_layout()
        plt.savefig(
            os.path.join(target_directory, 'IntraNPIs_' + filename + '_{}'.format(
                code)), dpi=300)
        plt.close()


def main():
    """ Main program entry."""

    # arg_dict = gd.cli("testing")
    df = get_npi_data(start_date=date(2020, 1, 1),
                      fine_resolution=2, file_format='csv')


if __name__ == "__main__":

    main()
