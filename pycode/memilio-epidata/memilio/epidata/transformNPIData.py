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
from datetime import datetime, timedelta, date
import sys
import time
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.spatial.distance import pdist
from scipy.cluster import hierarchy
from sklearn.cluster import KMeans

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import defaultDict as dd
from memilio.epidata import customPlot


def evaluate_clustering(corr_mat, idx_to_cluster_idx, indices_all):
    """! Computes a score for a particular clustering based on the 
    correlation matrix. The score is computed as the percentage of 'higher' 
    or 'high' values  (e.g., between 0.5 and 0.75 or 0.75 and 1) of the 
    correlation matrix that are to be found in the diagonal blocks of the 
    clustered correlation matrix vs these values in the offdiagonal blocks.

    @param corr_mat correlation matrix between the features / data set items  
        that were clustered.
    @param idx_to_cluster_idx Mapping of data item to cluster index.
    @param indices_all List of indices of all data items.

    @return Scores for the provided clustering.
    """

    if idx_to_cluster_idx.min() == 1:
        idx_to_cluster_idx -= 1

    # store indices of clusters
    clusters = [[] for i in range(idx_to_cluster_idx.max()+1)]
    for ii in range(len(idx_to_cluster_idx)):
        clusters[idx_to_cluster_idx[ii]].append(ii)
    # store remaining/perpendicular indices for all clusters
    clusters_perp = [[] for i in range(idx_to_cluster_idx.max()+1)]
    for ii in range(len(clusters)):
        clusters_perp[ii] = list(indices_all.difference(set(clusters[ii])))
    # extract correlation values of block diagonals and offdiagonals separ.
    corr_diag = []
    corr_offdiag = []
    for ii in range(len(clusters)):
        corr_diag = np.append(corr_diag, abs(
            corr_mat[np.ix_(clusters[ii], clusters[ii])].flatten()))
        corr_offdiag = np.append(corr_offdiag, abs(
            corr_mat[np.ix_(clusters[ii], clusters_perp[ii])].flatten()))

    corr_thresholds = [0.25, 0.5, 0.75]
    cluster_quantification = np.zeros(6)
    for ii in range(len(corr_thresholds)):
        num_diag = len(np.where(corr_diag > corr_thresholds[ii])[0])
        num_offdiag = len(np.where(corr_offdiag > corr_thresholds[ii])[0])
        if ii < len(corr_thresholds)-1:
            num_diag -= len(np.where(corr_diag > corr_thresholds[ii+1])[0])
            num_offdiag -= len(np.where(corr_offdiag >
                               corr_thresholds[ii+1])[0])
        cluster_quantification[2*ii] = num_diag / (num_diag+num_offdiag)
        cluster_quantification[2*ii+1] = (
            num_diag+num_offdiag) / (len(indices_all)**2)

    # print scores on clustering
    print("Number of clusters: " + str(len(clusters)) +
          ", shares diag/all between [0.25, 0.5, 0.75]: %.4f" %
          cluster_quantification[0] + " (%.4f" % cluster_quantification[1] +
          "), " + " %.4f " % cluster_quantification[2] + " (%.4f" %
          cluster_quantification[3] + "), " + " %.4f " %
          cluster_quantification[4] + " (%.4f" % cluster_quantification[5] + ")")

    return cluster_quantification


def compute_hierarch_clustering(corr_mat, corr_pairwdist,
                                metrics=['single', 'complete', 'average',
                                         'weighted', 'centroid', 'median',
                                         'ward']):
    """! Computes a hierarchical clustering for a (list of) metric(s) and 
    provides the maximum cophenetic distance(s) as well as a score for the 
    clustering (see @method evaluate_clustering(...)).

    @param corr_mat correlation matrix between the features / data set items  
        to be clustered hierarchically.
    @param corr_pairwdist Computed pairwise distance between the features / data
        set items.
    @param metric Metric or list of metrics to compute the hierarchical 
        clustering.

    @return (List of) hierarchical clustering(s), maximum cophenetic distance(s)
        and scores of the hierarchical clustering.
    """
    # NOTE: if changing metric, pay attention to linkage methods;
    #       'centroid', 'median', and 'ward' are correctly defined only if
    #       Euclidean pairwise metric is used.
    # Based on the distances, we compute an hierarchical clustering for
    # different metrics
    max_coph_corr = 0
    scores = dict()
    # allow single entry
    if not isinstance(metrics, list):
        metrics = [metrics]
    # iterate over list
    for metric in metrics:
        cluster_hierarch = hierarchy.linkage(corr_pairwdist, method=metric)
        # compute cophentic correlation distance
        coph_corr, coph_dists = hierarchy.cophenet(
            cluster_hierarch, pdist(corr_mat))
        scores[metric] = coph_corr
        if coph_corr > max_coph_corr:
            max_coph_corr = coph_corr
            max_metric = metric
            max_coph_dist = coph_dists

    cluster_hierarch = hierarchy.linkage(corr_pairwdist, method=max_metric)

    print(
        "Cophentic correlation distance for metric " + max_metric + ": " +
        str(max_coph_corr))

    return cluster_hierarch, max_coph_dist, scores


def flatten_hierarch_clustering(corr_mat, cluster_hierarch, weights):
    """! Flattens a hierarchical clustering for a (list of) maximum cophenetic 
    distance(s) in the flat clusters and evaluates the resulting clustering with
    respect to the corresponding correlation matrix.

    @param corr_mat correlation matrix between the features / data set items  
        clustered hierarchically.
    @param cluster_hierarch hierarchical clustering of given features  / data 
        set items.
    @param weigths Maximum cophenetic distance or list of maximum cophenetic 
        distances to compute the flat clustering(s).

    @return flat clustering(s) according to the (list of) maximum distance(s).
    """

    # all indices in npis_corr from 0 to n-1
    npi_indices_all = set(range(corr_mat.shape[0]))
    npi_idx_to_cluster_idx_list = []
    # allow single entries
    if not isinstance(weights, list):
        weights = [weights]
    # iterate over weights
    for weight in weights:
        # use the given weight to flatten the dendrogram
        npi_idx_to_cluster_idx = hierarchy.fcluster(
            cluster_hierarch, weight, criterion='distance')

        # evaluate clustering
        evaluate_clustering(corr_mat, npi_idx_to_cluster_idx, npi_indices_all)

        # append new npi_idx to cluster_idx assignment to list of assignments
        npi_idx_to_cluster_idx_list.append(npi_idx_to_cluster_idx)

    return npi_idx_to_cluster_idx_list


def print_manual_download(filename, url):

    print(
        'This script needs manual downloading of files. Please register'
        ' at corona-datenplatform.com and download ' + filename + ' from ' + url +
        '. Then move it to a folder named raw_data in this directory.')


def transform_npi_data(fine_resolution=2,
                       read_data=dd.defaultDict['read_data'],
                       file_format=dd.defaultDict['file_format'],
                       out_folder=dd.defaultDict['out_folder'],
                       start_date=dd.defaultDict['start_date'],
                       end_date=dd.defaultDict['end_date'],
                       make_plot=dd.defaultDict['make_plot'],
                       ):
    """! Loads a certain resolution of recorded NPI data from
    the Corona Datenplattform and transforms it according to the 
    arguments given. 

    For full functionality, please manually download
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
    @param make_plot False [Default] or True. Defines if plots are
        generated with matplotlib.
    @param moving_average 0 [Default] or Number>0. Defines the number of
        days for which a centered moving average is computed.
    """

    directory = out_folder
    directory = os.path.join(directory, 'Germany/')
    gd.check_dir(directory)

    if not read_data:

        if fine_resolution > 0:
            # defines delay in number of days between exceeding
            # incidence threshold and NPI getting active
            npi_activation_delay = 0
            print('Using a delay of NPI activation of ' +
                  str(npi_activation_delay) + ' days.')

            try:
                df_npis_old = pd.read_csv(
                    os.path.join(
                        directory, 'kr_massnahmen_unterkategorien.csv'),
                    sep=',')  # , nrows=1248)  # 1248 for debugging, only reading Flensburg
            except FileNotFoundError:
                print_manual_download(
                    'kr_massnahmen_unterkategorien.csv',
                    'https://www.corona-datenplattform.de/dataset/massnahmen_unterkategorien_kreise')
                raise FileNotFoundError
            df_npis_old.rename(dd.GerEng, axis=1, inplace=True)

            # check if rows hospitals and geriatric care are still empty
            # these fields have been empty so far and are thus not used
            test_codes = ['M23_010', 'M23_020', 'M23_030', 'M23_040',
                          'M23_050', 'M23_060', 'M24_010', 'M24_020',
                          'M24_030', 'M24_040', 'M24_050', 'M24_060']
            for tcode in test_codes:
                for i in [''] + ["_" + str(i) for i in range(1, 6)]:
                    if(df_npis_old[df_npis_old[dd.EngEng['npiCode']] == tcode+i].iloc[:, 6:].max().max() > 0):
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

    else:  # read formatted file

        if fine_resolution > 0:
            if fine_resolution == 1:
                filename = 'germany_counties_npi_subcat_incgrouped'
            else:
                filename = 'germany_counties_npi_subcat'
        else:
            filename = 'germany_counties_npi_maincat'
        df_npis = pd.read_json(directory + filename + ".json")

    # read data frame of variable names and descriptions
    try:
        if fine_resolution > 0:
            df_npis_desc = pd.read_excel(
                os.path.join(
                    directory, 'datensatzbeschreibung_massnahmen.xlsx'),
                sheet_name=2)
        else:
            df_npis_desc = pd.read_excel(
                os.path.join(
                    directory, 'datensatzbeschreibung_massnahmen.xlsx'),
                sheet_name=3)
    except FileNotFoundError:
        print_manual_download(
            'datensatzbeschreibung_massnahmen.xlsx',
            'https://www.corona-datenplattform.de/dataset/massnahmen_unterkategorien_kreise')
        raise FileNotFoundError

    # get existing codes that are used (in df_npis_old M22-M24 are empty)
    npi_codes_prior = df_npis_desc['Variablenname']

    # correct differences in codes between data sheet and explanation sheet
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
        if not read_data:
            npi_codes_prior_data = df_npis_old[dd.EngEng['npiCode']].unique()
        else:
            npi_codes_prior_data = list(df_npis.columns[2:])

        missing_codes = list(set(npi_codes_prior).difference(
            npi_codes_prior_data))
        if len(missing_codes) > 0:
            # if incidence is grouped, only search for grouping codes without
            # having a detailed "_DETAIL" naming as of MCODE_NUMBER_DETAIL
            if fine_resolution == 1:
                missing_grouped_codes = []
                for mcode in missing_codes:
                    if len(mcode.split('_')) < 2:
                        missing_grouped_codes.append(mcode)
                if len(missing_grouped_codes) > 0:  # only MCODE_NUMBER codes
                    sys.exit('Missing NPI codes: ' +
                             str(missing_grouped_codes))
            else:
                sys.exit('Missing NPI codes: ' + str(missing_codes))

        # we dont have any explanations from "datensatzbeschreibung_massnahmen"
        # on these codes, so drop the rows.
        codes_dropped = list(set(npi_codes_prior_data).difference(
            npi_codes_prior))
        if len(codes_dropped) > 0:
            df_npis_old = df_npis_old[~df_npis_old[dd.EngEng['npiCode']].isin(
                codes_dropped)].reset_index()

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

    # NPIs groups codes and description to ensure that both are ordered
    # the same way
    npis_dummy = {dd.EngEng['npiCode']: npi_codes, dd.EngEng['desc']: npi_desc}
    npis = pd.DataFrame(npis_dummy)

    # transform data from original format to desired format
    if not read_data:
        # prepare grouping of NPIs to reduce product space of
        # NPI x active_from_inc (with values "incidence does not matter", and
        # incidence 0, 10, 35, 50, 100) to NPI
        if fine_resolution == 1:
            # create hash table from parental or main code/main category
            # to list of subcodes/subcategories
            maincode_to_npicodes_map = dict()
            npicodes_to_maincode_map = dict()
            major_code = npi_codes[0]
            maincode_to_npicodes_map[major_code] = []
            for code in npi_codes:
                npicodes_to_maincode_map[code] = major_code
                if major_code in code:
                    maincode_to_npicodes_map[major_code].append(code)
                else:
                    major_code = code
                    maincode_to_npicodes_map[major_code] = [code]

            npi_codes_aggregated = []
            for main_code in maincode_to_npicodes_map.keys():
                if main_code.count('_') > 1:
                    sys.exit('Error. Subcode assigned as main code.')
                npi_codes_aggregated.append(main_code)

            npis_final = npis[npis[dd.EngEng['npiCode']].isin(
                npi_codes_aggregated)].reset_index()
        else:
            npis_final = npis

        # extract incidence-threshold for NPIs
        if fine_resolution > 0:
            npi_incid_start = dict()
            for i in range(len(npis)):
                incid_threshold = 1e10
                if npis.loc[i, dd.EngEng['desc']].split(' ')[0] == 'Unabhängig':
                    # set -1 for incidence-independent NPIs
                    incid_threshold = -1
                elif npis.loc[i, dd.EngEng['desc']].split(' ')[0] == 'Ab':
                    incid_threshold = int(
                        npis.loc[i, dd.EngEng['desc']].split(' ')[1])
                else:
                    sys.exit(
                        'Error in description file. NPI activation can not '
                        'be computed. Exiting.')
                npi_incid_start[npis.loc[i, dd.EngEng['npiCode']]
                                ] = incid_threshold

            # get all incidence thresholds
            incidence_thresholds = sorted(set(npi_incid_start.values()))

        # create hash map from thresholds to NPI indices
        incidence_thresholds_to_npis = dict(
            zip(incidence_thresholds, [[] for i in range(len(incidence_thresholds))]))
        for i in range(len(npis)):
            incval = npi_incid_start[npis.loc
                                     [i, dd.EngEng['npiCode']]]
            incidence_thresholds_to_npis[incval].append(i)

        # get county ids
        unique_geo_entities = geoger.get_county_ids()
        # check if more than the county of Eisenach would be removed with
        # current county list
        counties_removed = df_npis_old[
            ~df_npis_old[dd.EngEng['idCounty']].isin(unique_geo_entities)][
            dd.EngEng['idCounty']].unique()
        if len(counties_removed) == 1 and counties_removed[0] != 16056:
            sys.exit('Error. Other counties than that of Eisenach were removed.')
        # remove rows for Eisenach
        df_npis_old = df_npis_old[df_npis_old[dd.EngEng['idCounty']].isin(
            unique_geo_entities)]

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
            sys.exit('Exiting. Dates missing in data frame.')

        min_date = []
        max_date = []

        # get RKI infectious numbers to find dates where incidence-dependent
        # NPIs were active
        if fine_resolution > 0:
            df_infec_rki = pd.read_json(os.path.join(
                directory, 'all_county_all_dates_repdate_rki.json'))
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
        counters = np.zeros(5)  # time counter for output only
        countyidx = 0
        # unique_geo_entities:
        for countyID in [5315, 1001, 9162, 16071, 11000, 1060, 5566]:
            cid = 0
            countyidx += 1

            if fine_resolution > 0:
                # compute incidence based on previous data frames
                df_infec_local = df_infec_rki[df_infec_rki[dd.EngEng['idCounty']] == countyID].copy(
                )
                pop_local = df_population.loc[df_population[dd.EngEng['idCounty']]
                                              == countyID, dd.EngEng['population']].values[0]
                incidence_local = df_infec_local[dd.EngEng['confirmed']].diff(
                    periods=7).fillna(df_infec_local[dd.EngEng['confirmed']])
                df_infec_local['Incidence'] = incidence_local / pop_local * 100000

                # set to main data frame
                df_infec_rki.loc[df_infec_rki[dd.EngEng['idCounty']] ==
                                 countyID, 'Incidence'] = df_infec_local['Incidence'].values

                # cut infection information at start_date_new and end_date_new
                df_infec_local = df_infec_local.loc[(df_infec_local[dd.EngEng['date']] >= start_date_new) & (
                    df_infec_local[dd.EngEng['date']] <= end_date_new), :].reset_index()

            # get county-local data frame
            start_time = time.perf_counter()
            df_local_old = df_npis_old[df_npis_old[dd.EngEng['idCounty']]
                                       == countyID].copy()

            # remove potential rows of which codes are not in npi_codes_considered
            npi_rows = [i in npis[dd.EngEng['npiCode']].values
                        for i in df_local_old[dd.EngEng['npiCode']]]

            # get list of NPI codes, ordered as the rows in the current data frame
            npi_codes_ordered_as_rows = df_local_old[dd.EngEng['npiCode']][
                npi_rows].to_list()

            # get indices of rows for the NPI codes as in the sorted npi_codes list
            # may be superfluous if NPI code rows are sorted correctly
            npi_code_rows_to_sorted = [
                npi_codes_ordered_as_rows.index(i) for i in
                npis[dd.EngEng['npiCode']].values]

            # access NPI values matrix and store it as integers
            npi_vals = df_local_old.iloc[npi_rows, start_npi_cols:].astype(int)

            # create columns for date, county ID and NPI code
            df_local_new = pd.DataFrame(
                columns=[dd.EngEng['date']] + [dd.EngEng['idCounty']] +
                list(npis[dd.EngEng['npiCode']]))

            counters[cid] += time.perf_counter()-start_time
            cid += 1

            start_time = time.perf_counter()

            # fill in NPI values by transposing from columns to rows
            df_local_new[dd.EngEng['date']] = dates_new
            df_local_new[dd.EngEng['idCounty']] = countyID
            # possible resorting of rows such that they are sorted according to
            # a literal sorting of the code strings
            df_local_new[npis[dd.EngEng['npiCode']]] = np.transpose(
                npi_vals.iloc[npi_code_rows_to_sorted, :].values)

            counters[cid] += time.perf_counter()-start_time
            cid += 1

            start_time = time.perf_counter()

            # replace -99 ("not used anymore") by 0 ("not used")
            df_local_new[npis[dd.EngEng['npiCode']]
                         ] = df_local_new[npis[dd.EngEng['npiCode']]].replace(-99, 0)
            # replace 2,3,4,5 ("mentioned in ...") by 1 ("mentioned")
            df_local_new[npis[dd.EngEng['npiCode']]
                         ] = df_local_new[npis[dd.EngEng['npiCode']]].replace([2, 3, 4, 5], 1)

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
                for incidvalthrsh, npi_indices in incidence_thresholds_to_npis.items():
                    if incidvalthrsh >= 0:
                        local_incid = df_infec_local['Incidence'].copy()
                        if npi_activation_delay > 0:
                            # shift values to npi_activation_delay days later
                            local_incid.iloc[npi_activation_delay:
                                             ] = local_incid.iloc[0:-npi_activation_delay].values
                            # take constant value of day 0 for first delay days
                            local_incid.iloc[:npi_activation_delay] = local_incid.iloc[0]
                        # compare incidence against threshold
                        int_active = (
                            local_incid >= incidvalthrsh).astype(int)
                        # multiply rows of data frame by either 1 if threshold
                        # passed (i.e., mentioned NPI is active) or zero
                        # (i.e., mentioned NPI is not active)
                        # 'mul' multiplies the original data frame row by row
                        # with the respective value in int_active
                        df_local_new.iloc[:, npis_idx_start + np.array(npi_indices)] \
                            = df_local_new.iloc[:, npis_idx_start + np.array(npi_indices)].mul(int_active, axis=0)

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
            df_npis = df_npis.append(df_local_new,
                                     ignore_index=True).copy()
            counters[cid] += time.perf_counter()-start_time
            cid += 1

            # divide working time by completed number of counties and multiply
            # by remaining number of counties to estimate time remaining
            time_remain = sum(
                counters) / countyidx * (len(unique_geo_entities) - countyidx)
            # print progress
            if countyidx == 1 or countyidx % int(
                    len(unique_geo_entities) / 10) == 0:
                print('Progress ' + str(countyidx) + ' / ' +
                      str(len(unique_geo_entities)) +
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
        if fine_resolution == 2:
            # the following validation is only valid for end_date May 14, 2021
            # and data frames of mentioned NPIs, not active NPIs
            start_date_validation = datetime(2020, 3, 1)
            end_date_validation = datetime(2022, 2, 15)
            end_date_validation_str = 'd20220215'

            # Cologne for M01a_010 and all dates (no changes, 0 to 0)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 5315) & (df_npis_old[dd.EngEng['npiCode']] == 'M01a_010')].iloc[:, :list(
                df_npis_old.columns).index(end_date_validation_str)+1].values[0][start_npi_cols:]
            dummy_new = df_npis.loc[(df_npis[dd.EngEng['idCounty']] == 5315) & (
                df_npis[dd.EngEng['date']] <= end_date_validation), 'M01a_010'].values
            valdiff = 0
            for i in range(2, 6):
                valdiff += (i-1)*len(np.where(dummy_old == i)[0])
            valdiff += 99*len(np.where(dummy_old == -99)[0])
            print(abs(dummy_old-dummy_new).sum() == valdiff)

            # Cologne for M01a_150_2 and all dates (no changes, 0 to 0, 1 to 1)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 5315) & (df_npis_old[dd.EngEng['npiCode']] == 'M01a_150_2')].iloc[:, :list(
                df_npis_old.columns).index(end_date_validation_str)+1].values[0][start_npi_cols:]
            dummy_new = df_npis.loc[(df_npis[dd.EngEng['idCounty']] == 5315) & (
                df_npis[dd.EngEng['date']] <= end_date_validation), 'M01a_150_2'].values
            incid = df_infec_rki[(df_infec_rki[dd.EngEng['idCounty']] == 5315) &
                                 (df_infec_rki[dd.EngEng['date']] <=
                                  end_date_validation) &
                                 (df_infec_rki[dd.EngEng['date']] >=
                                  start_date_validation)]['Incidence'].values
            npi_index = np.where(dummy_old >= 1)[0]
            incid_index = np.where(incid >= 10)[0]
            active_index = np.sort(
                list(set(npi_index).intersection(incid_index)))
            nonactive_index = np.sort(
                list(set(npi_index).difference(active_index)))
            valdiff = 0
            for i in range(2, 6):
                # these values >=1 are set to 1
                valdiff += (i-1)*len(
                    np.where(dummy_old[list(active_index)] == i)[0])
                # these values >=1 are set to 0
                valdiff += i * len(
                    np.where(dummy_old[list(nonactive_index)] == i)[0])
            # -99 always set to 0
            valdiff += 99*len(np.where(dummy_old == -99)[0])
            print(abs(dummy_old-dummy_new).sum() == valdiff)

            # Flensburg for M05_120 and all dates ('2's become '1's)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 1001) & (df_npis_old[dd.EngEng['npiCode']] == 'M05_120')].iloc[:, :list(
                df_npis_old.columns).index(end_date_validation_str)+1].values[0][start_npi_cols:]
            dummy_new = df_npis.loc[(df_npis[dd.EngEng['idCounty']] == 1001) & (
                df_npis[dd.EngEng['date']] <= end_date_validation), 'M05_120'].values
            incid = df_infec_rki[(df_infec_rki[dd.EngEng['idCounty']] == 1001) &
                                 (df_infec_rki[dd.EngEng['date']] <=
                                  end_date_validation) &
                                 (df_infec_rki[dd.EngEng['date']] >=
                                  start_date_validation)]['Incidence'].values
            npi_index = np.where(dummy_old >= 1)[0]
            incid_index = np.where(incid >= 10)[0]
            active_index = np.sort(
                list(set(npi_index).intersection(incid_index)))
            nonactive_index = np.sort(
                list(set(npi_index).difference(active_index)))
            valdiff = 0
            for i in range(2, 6):
                # these values >=1 are set to 1
                valdiff += (i-1)*len(
                    np.where(dummy_old[list(active_index)] == i)[0])
                # these values >=1 are set to 0
                valdiff += i * len(
                    np.where(dummy_old[list(nonactive_index)] == i)[0])
            # -99 always set to 0
            valdiff += 99*len(np.where(dummy_old == -99)[0])
            print(abs(dummy_old-dummy_new).sum() == valdiff)

            # Munich for M01a_010_4 and all dates (-99 becomes 0, 0 stays 0)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 9162) & (df_npis_old[dd.EngEng['npiCode']] == 'M01a_010_4')].iloc[:, :list(
                df_npis_old.columns).index(end_date_validation_str)+1].values[0][start_npi_cols:]
            dummy_new = df_npis.loc[(df_npis[dd.EngEng['idCounty']] == 9162) & (
                df_npis[dd.EngEng['date']] <= end_date_validation), 'M01a_010_4'].values
            incid = df_infec_rki[(df_infec_rki[dd.EngEng['idCounty']] == 9162) &
                                 (df_infec_rki[dd.EngEng['date']] <=
                                  end_date_validation) &
                                 (df_infec_rki[dd.EngEng['date']] >=
                                  start_date_validation)]['Incidence'].values
            npi_index = np.where(dummy_old >= 1)[0]
            incid_index = np.where(incid >= 10)[0]
            active_index = np.sort(
                list(set(npi_index).intersection(incid_index)))
            nonactive_index = np.sort(
                list(set(npi_index).difference(active_index)))
            valdiff = 0
            for i in range(2, 6):
                # these values >=1 are set to 1
                valdiff += (i-1)*len(
                    np.where(dummy_old[list(active_index)] == i)[0])
                # these values >=1 are set to 0
                valdiff += i * len(
                    np.where(dummy_old[list(nonactive_index)] == i)[0])
            # -99 always set to 0
            valdiff += 99*len(np.where(dummy_old == -99)[0])
            print(abs(dummy_old-dummy_new).sum() == valdiff)

            # Weimar for M12_030_3 and all dates (-99 becomes 0, 1 stays 1, 0 stays 0)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 16071) & (df_npis_old[dd.EngEng['npiCode']] == 'M18_030_4')].iloc[:, :list(
                df_npis_old.columns).index(end_date_validation_str)+1].values[0][start_npi_cols:]
            dummy_new = df_npis.loc[(df_npis[dd.EngEng['idCounty']] == 16071) & (
                df_npis[dd.EngEng['date']] <= end_date_validation), 'M18_030_4'].values
            incid = df_infec_rki[(df_infec_rki[dd.EngEng['idCounty']] == 16071) &
                                 (df_infec_rki[dd.EngEng['date']] <=
                                  end_date_validation) &
                                 (df_infec_rki[dd.EngEng['date']] >=
                                  start_date_validation)]['Incidence'].values
            npi_index = np.where(dummy_old >= 1)[0]
            incid_index = np.where(incid >= 50)[0]
            active_index = np.sort(
                list(set(npi_index).intersection(incid_index)))
            nonactive_index = np.sort(
                list(set(npi_index).difference(active_index)))
            valdiff = 0
            for i in range(2, 6):
                # these values >=1 are set to 1
                valdiff += (i-1)*len(
                    np.where(dummy_old[list(active_index)] == i)[0])
                # these values >=1 are set to 0
                valdiff += i * len(
                    np.where(dummy_old[list(nonactive_index)] == i)[0])
            # -99 always set to 0
            valdiff += 99*len(np.where(dummy_old == -99)[0])
            print(abs(dummy_old-dummy_new).sum() == valdiff)  # TODO

            # Berlin for M01b_020 and all dates (2 becomes 1, 1 stays 1, 0 stays 0)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 11000) & (df_npis_old[dd.EngEng['npiCode']] == 'M01b_020')].iloc[:, :list(
                df_npis_old.columns).index(end_date_validation_str)+1].values[0][start_npi_cols:]
            dummy_new = df_npis.loc[(df_npis[dd.EngEng['idCounty']] == 11000) & (
                df_npis[dd.EngEng['date']] <= end_date_validation), 'M01b_020'].values
            incid = df_infec_rki[(df_infec_rki[dd.EngEng['idCounty']] == 11000) &
                                 (df_infec_rki[dd.EngEng['date']] <=
                                  end_date_validation) &
                                 (df_infec_rki[dd.EngEng['date']] >=
                                  start_date_validation)]['Incidence'].values
            npi_index = np.where(dummy_old >= 1)[0]
            incid_index = np.where(incid >= 10)[0]
            active_index = np.sort(
                list(set(npi_index).intersection(incid_index)))
            nonactive_index = np.sort(
                list(set(npi_index).difference(active_index)))
            valdiff = 0
            for i in range(2, 6):
                # these values >=1 are set to 1
                valdiff += (i-1)*len(
                    np.where(dummy_old[list(active_index)] == i)[0])
                # these values >=1 are set to 0
                valdiff += i * len(
                    np.where(dummy_old[list(nonactive_index)] == i)[0])
            # -99 always set to 0
            valdiff += 99*len(np.where(dummy_old == -99)[0])
            print(abs(dummy_old-dummy_new).sum() == valdiff)  # TODO

            # Segeberg for M02b_035 and all dates (2 -> 1, 3 -> 1, 5 -> 1)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 1060) & (df_npis_old[dd.EngEng['npiCode']] == 'M02b_035')].iloc[:, :list(
                df_npis_old.columns).index(end_date_validation_str)+1].values[0][start_npi_cols:]
            dummy_new = df_npis.loc[(df_npis[dd.EngEng['idCounty']] == 1060) & (
                df_npis[dd.EngEng['date']] <= end_date_validation), 'M02b_035'].values
            incid = df_infec_rki[(df_infec_rki[dd.EngEng['idCounty']] == 1060) &
                                 (df_infec_rki[dd.EngEng['date']] <=
                                  end_date_validation) &
                                 (df_infec_rki[dd.EngEng['date']] >=
                                  start_date_validation)]['Incidence'].values
            npi_index = np.where(dummy_old >= 1)[0]
            incid_index = np.where(incid >= 10)[0]
            active_index = np.sort(
                list(set(npi_index).intersection(incid_index)))
            nonactive_index = np.sort(
                list(set(npi_index).difference(active_index)))
            valdiff = 0
            for i in range(2, 6):
                # these values >=1 are set to 1
                valdiff += (i-1)*len(
                    np.where(dummy_old[list(active_index)] == i)[0])
                # these values >=1 are set to 0
                valdiff += i * len(
                    np.where(dummy_old[list(nonactive_index)] == i)[0])
            # -99 always set to 0
            valdiff += 99*len(np.where(dummy_old == -99)[0])
            print(abs(dummy_old-dummy_new).sum() == valdiff)  # TODO

            # Steinfurt for M16_050 and all dates (4 -> 1, ...)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 5566) & (df_npis_old[dd.EngEng['npiCode']] == 'M16_050')].iloc[:, :list(
                df_npis_old.columns).index(end_date_validation_str)+1].values[0][start_npi_cols:]
            dummy_new = df_npis.loc[(df_npis[dd.EngEng['idCounty']] == 5566) & (
                df_npis[dd.EngEng['date']] <= end_date_validation), 'M16_050'].values
            incid = df_infec_rki[(df_infec_rki[dd.EngEng['idCounty']] == 5566) &
                                 (df_infec_rki[dd.EngEng['date']] <=
                                  end_date_validation) &
                                 (df_infec_rki[dd.EngEng['date']] >=
                                  start_date_validation)]['Incidence'].values
            npi_index = np.where(dummy_old >= 1)[0]
            incid_index = np.where(incid >= 10)[0]
            active_index = np.sort(
                list(set(npi_index).intersection(incid_index)))
            nonactive_index = np.sort(
                list(set(npi_index).difference(active_index)))
            valdiff = 0
            for i in range(2, 6):
                # these values >=1 are set to 1
                valdiff += (i-1)*len(
                    np.where(dummy_old[list(active_index)] == i)[0])
                # these values >=1 are set to 0
                valdiff += i * len(
                    np.where(dummy_old[list(nonactive_index)] == i)[0])
            # -99 always set to 0
            valdiff += 99*len(np.where(dummy_old == -99)[0])
            print(abs(dummy_old-dummy_new).sum() == valdiff)

            x = 15
        elif fine_resolution == 1:
            # replace more detailed code names X_Y with major code X
            df_npis_old[dd.EngEng['npiCode']] = df_npis_old[dd.EngEng[
                'npiCode']].replace(npicodes_to_maincode_map)

            # Cologne for M01a_010 and all dates (no changes)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 5315) & (
                df_npis_old[dd.EngEng['npiCode']] == 'M01a_010')].values[0][start_npi_cols:]
            for subcode in range(1, 6):  # add subcode values
                dummy_old = np.maximum(dummy_old, df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 5315) & (
                    df_npis_old[dd.EngEng['npiCode']] == 'M01a_010')].values[subcode][start_npi_cols:])
            dummy_new = df_npis.loc[df_npis[dd.EngEng['idCounty']] ==
                                    5315, 'M01a_010'].values
            print(abs(dummy_old-dummy_new).sum() == 0)

            # Flensburg for M05_120 and all dates ('2's become '1's)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 1001) & (
                df_npis_old[dd.EngEng['npiCode']] == 'M05_120')].values[0][start_npi_cols:]
            for subcode in range(1, 6):  # add subcode values
                dummy_old = np.maximum(dummy_old,
                                       df_npis_old
                                       [(df_npis_old
                                         [dd.EngEng['idCounty']] == 1001) &
                                        (df_npis_old[dd.EngEng['npiCode']] ==
                                         'M05_120')].values[subcode]
                                       [start_npi_cols:])
            dummy_new = df_npis.loc[df_npis[dd.EngEng['idCounty']] ==
                                    1001, 'M05_120'].values
            print(abs(dummy_old-dummy_new).sum() == 5)

            # Munich for M01a_010 and all dates (-99 becomes 0, 0 stays 0)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 9162) & (
                df_npis_old[dd.EngEng['npiCode']] == 'M01a_010')].values[0][start_npi_cols:]
            for subcode in range(1, 6):  # add subcode values
                dummy_old = np.maximum(dummy_old, df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 9162) & (
                    df_npis_old[dd.EngEng['npiCode']] == 'M01a_010')].values[subcode][start_npi_cols:])
            dummy_new = df_npis.loc[df_npis[dd.EngEng['idCounty']] ==
                                    9162, 'M01a_010'].values
            print(abs(dummy_old-dummy_new).sum() == 0)

            # Weimar for M12_030_3 and all dates (-99 becomes 0, 1 stays 1, 0 stays 0)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 16071) & (
                df_npis_old[dd.EngEng['npiCode']] == 'M12_030')].values[0][start_npi_cols:]
            for subcode in range(1, 6):  # add subcode values
                dummy_old = np.maximum(dummy_old,
                                       df_npis_old
                                       [(df_npis_old
                                         [dd.EngEng['idCounty']] == 16071) &
                                        (df_npis_old[dd.EngEng['npiCode']] ==
                                         'M12_030')].values[subcode]
                                       [start_npi_cols:])
            dummy_new = df_npis.loc[df_npis[dd.EngEng['idCounty']] ==
                                    16071, 'M12_030'].values
            print(abs(dummy_old-dummy_new).sum() == 19)

            # Berlin for M01b_020 and all dates (2 becomes 1, 1 stays 1, 0 stays 0)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 11000) & (
                df_npis_old[dd.EngEng['npiCode']] == 'M01b_020')].values[0][start_npi_cols:]
            for subcode in range(1, 6):  # add subcode values
                dummy_old = np.maximum(dummy_old, df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 11000) & (
                    df_npis_old[dd.EngEng['npiCode']] == 'M01b_020')].values[subcode][start_npi_cols:])
            dummy_new = df_npis.loc[df_npis[dd.EngEng['idCounty']] ==
                                    11000, 'M01b_020'].values
            print(abs(dummy_old-dummy_new).sum() == 82)

            # Segeberg for M02b_035 and all dates (2 -> 1, 3 -> 1, 5 -> 1)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 1060) & (
                df_npis_old[dd.EngEng['npiCode']] == 'M02b_035')].values[0][start_npi_cols:]
            for subcode in range(1, 6):  # add subcode values
                dummy_old = np.maximum(dummy_old, df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 1060) & (
                    df_npis_old[dd.EngEng['npiCode']] == 'M02b_035')].values[subcode][start_npi_cols:])
            dummy_new = df_npis.loc[df_npis[dd.EngEng['idCounty']] ==
                                    1060, 'M02b_035'].values
            print(abs(dummy_old-dummy_new).sum() == 151+2*53+4*22)

            # Steinfurt for M16_050 and all dates (4 -> 1, ...)
            dummy_old = df_npis_old[(df_npis_old[dd.EngEng['idCounty']] == 5566) & (
                df_npis_old[dd.EngEng['npiCode']] == 'M16_050')].values[0][start_npi_cols:]
            for subcode in range(1, 6):  # add subcode values
                dummy_old = np.maximum(dummy_old,
                                       df_npis_old
                                       [(df_npis_old
                                         [dd.EngEng['idCounty']] == 5566) &
                                        (df_npis_old[dd.EngEng['npiCode']] ==
                                         'M16_050')].values[subcode]
                                       [start_npi_cols:])
            dummy_new = df_npis.loc[df_npis[dd.EngEng['idCounty']] ==
                                    5566, 'M16_050'].values
            print(abs(dummy_old-dummy_new).sum() == 32+2*20+3*22)
        #### end validation ####

        if fine_resolution > 0:
            if fine_resolution == 1:
                filename = 'germany_counties_npi_subcat_incgrouped'
            else:
                filename = 'germany_counties_npi_subcat'
        else:
            filename = 'germany_counties_npi_maincat'
        gd.write_dataframe(df_npis, directory, filename, file_format)

        # stupid validation
        # df_validation = pd.read_json(directory + filename + ".json")
        # if len(
        #     np.where(
        #         df_validation.iloc[:, start_npi_cols - 1:] != df_npis.iloc
        #         [:, start_npi_cols - 1:])[0]) > 0:
        #     print('Error in file writing/reading')

    # get code levels (main/subcodes) and position of main codes
    # code_level = [i.count('_') for i in npi_codes]
    # main_code_pos = [i for i in range(len(code_level)) if code_level[i] == 1]

    # check if any other integer than 0: not implemented or 1: implemented is
    # used (maybe to specify the kind of implementation)
    if len(np.where(df_npis[npi_codes_considered] > 1)[0]) > 0:

        print("Info: Please ensure that NPI information is only boolean.")

    else:
        # sum over different NPIs and plot share of countires implementing
        # these NPIs versus counties without corresponding actions
        df_npis_aggregated = df_npis.groupby(
            dd.EngEng['date']).agg(
            {i: sum for i in npi_codes_considered}).copy()
        npis_total_sum = df_npis_aggregated.sum()

        npi_codes_empty = list(np.array(npi_codes_considered)[
                               np.where(npis_total_sum == 0)[0]])

        npi_unused_indices_all = []
        npi_used_indices_all = []
        npi_unused_indices = []
        npi_used_indices = []
        for i in range(len(npi_codes_considered)):
            if npi_codes_considered[i] in npi_codes_empty:
                npi_unused_indices.append(i)
                npi_unused_indices_all.append(
                    npis[dd.EngEng['npiCode']].index(npi_codes_considered[i]))
            else:
                npi_used_indices.append(i)
                npi_used_indices_all.append(
                    npis[dd.EngEng['npiCode']].index(npi_codes_considered[i]))

        npis_unused = np.array(npis[dd.EngEng['desc']])[npi_unused_indices_all]
        npis_used = np.array(npis[dd.EngEng['desc']])[npi_used_indices_all]
        npi_codes_used = list(np.array(npi_codes_considered)[npi_used_indices])
        npi_codes_unused = list(
            np.array(npi_codes_considered)[npi_unused_indices])

        # open file to write unused categories
        if fine_resolution > 0:
            if fine_resolution == 1:
                filename = 'unused_subcats_incgrouped.txt'
            else:
                filename = 'unused_subcats.txt'
        else:
            filename = 'unused_maincats.txt'
        file_npi = open(directory + filename, 'w')
        # Writing unused NPIs
        for i in range(len(npis_unused)):
            file_npi.write(npi_codes_unused[i] + ": " + npis_unused[i])
            file_npi.write("\n")
        # Closing file
        file_npi.close()

        # open file to write unused categories
        if fine_resolution > 0:
            if fine_resolution == 1:
                filename = 'used_subcats_incgrouped.txt'
            else:
                filename = 'used_subcats.txt'
        else:
            filename = 'used_maincats.txt'
        file_npi = open(directory + filename, 'w')
        # Writing unused NPIs
        for i in range(len(npis_used)):
            file_npi.write(npi_codes_used[i] + ": " + npis_used[i])
            file_npi.write("\n")
        # Closing file
        file_npi.close()

        df_npis_used = df_npis[[dd.EngEng['date'],
                                dd.EngEng['idCounty']] + npi_codes_used].copy()
        if fine_resolution > 0:
            if fine_resolution == 1:
                filename = 'germany_counties_npi_subcat_used_incgrouped'
            else:
                filename = 'germany_counties_npi_subcat_used'
        else:
            filename = 'germany_counties_npi_maincat_used'
        gd.write_dataframe(df_npis_used, directory, filename, file_format)

        # compute correlations
        npis_corr = df_npis_used.iloc[:, 2:].corr().values
        # plot log-colored correlations
        plt.imshow(abs(npis_corr), cmap='gray_r')
        # plot histogram
        plt.figure()
        plt.hist(npis_corr.flatten(), bins=50)
        plt.title("Correlation histogram", fontsize=18)
        plt.xlabel("Correlation", fontsize=12)
        plt.ylabel("Number of values", fontsize=12)

        # We understand the rows of npis_corr, the correlations of one NPI
        # to the others as one node in the #NPIs-used-dimensional space.
        # We compute the pairwise distances of these nodes. Then, nodes with
        # similar correlations towards all other nodes exhibit small distances
        corr_pairwdist = hierarchy.distance.pdist(
            npis_corr, metric='euclidean')

        # compute hierarchical clustering (via best-suited metric)
        compare_metrics = True
        if compare_metrics:
            # centroid
            metric = 'centroid'
            cluster_hierarch, coph_dist, scores = compute_hierarch_clustering(
                abs(npis_corr),
                corr_pairwdist,
                metric)
            # # plot dendrogram
            plt.figure()
            plt.title(metric)
            hierarchy.dendrogram(cluster_hierarch)
            plt.show()
            max_coph_dist = coph_dist.max()
            flatten_hierarch_clustering(
                abs(npis_corr), cluster_hierarch,
                [wg * max_coph_dist
                 for wg in [0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75]])
            # ward
            metric = 'ward'
            cluster_hierarch, coph_dist, scores = compute_hierarch_clustering(
                npis_corr,
                corr_pairwdist,
                metric)
            # # plot dendrogram
            # plt.figure()
            # plt.title(metric)
            # hierarchy.dendrogram(cluster_hierarch)
            # plt.show()
            max_coph_dist = coph_dist.max()
            flatten_hierarch_clustering(
                abs(npis_corr), cluster_hierarch,
                [wg * max_coph_dist for wg in [0.1, 0.125, 0.15, 0.175, 0.2]])
            # average
            metric = 'average'
            cluster_hierarch, coph_dist, scores = compute_hierarch_clustering(
                abs(npis_corr),
                corr_pairwdist,
                metric)
            # # plot dendrogram
            # plt.figure()
            # plt.title(metric)
            # hierarchy.dendrogram(cluster_hierarch)
            # plt.show()
            max_coph_dist = coph_dist.max()
            flatten_hierarch_clustering(
                npis_corr, cluster_hierarch,
                [wg * max_coph_dist
                 for wg in [0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65]])

        metric = 'centroid'
        cluster_hierarch, coph_dist, scores = compute_hierarch_clustering(
            npis_corr,
            corr_pairwdist,
            metric)
        # # plot dendrogram
        # plt.figure()
        # plt.title(metric)
        # hierarchy.dendrogram(cluster_hierarch)
        # plt.show()
        max_coph_dist = coph_dist.max()
        npi_idx_to_cluster_idx = flatten_hierarch_clustering(
            npis_corr, cluster_hierarch,
            [wg * max_coph_dist
             for wg in [0.65]])

        cluster_dict = dict()
        cluster_codes = [[] for i in range(npi_idx_to_cluster_idx[0].max()+1)]
        cluster_desc = [[] for i in range(npi_idx_to_cluster_idx[0].max()+1)]
        for i in range(len(npi_idx_to_cluster_idx[0])):
            cluster_dict[npi_codes_used[i]
                         ] = "CM_" + str(npi_idx_to_cluster_idx[0][i]).zfill(3)
            cluster_codes[npi_idx_to_cluster_idx[0]
                          [i]].append(npi_codes_used[i])
            cluster_desc[npi_idx_to_cluster_idx[0]
                         [i]].append(str(npis_used[i]))

        # create clustered dataframe
        df_npis_clustered = df_npis[[
            dd.EngEng['date'], dd.EngEng['idCounty']]].copy()

        for i in range(len(cluster_codes)):
            df_npis_clustered["CM_" + str(i).zfill(3)
                              ] = df_npis[cluster_codes[i]].max(axis=1).copy()

        npis_corr_cluster = df_npis_clustered.corr()
        # npis_corr_cluster[abs(npis_corr_cluster)<0.25] = 0
        plt.imshow(abs(npis_corr_cluster), cmap='gray_r')
        plt.title('Absolute correlation>0.25 of clustered NPIs')
        plt.xlabel('NPI cluster')
        plt.ylabel('NPI cluster')
        plt.colorbar()

        # open file to write unused categories
        if fine_resolution > 0:
            if fine_resolution == 1:
                filename = 'clusters_subcats_incgrouped.txt'
            else:
                filename = 'clusters_subcats.txt'
        else:
            filename = 'clusters_maincats.txt'
        file_npi = open(directory + filename, 'w')
        # Writing unused NPIs
        for i in range(len(cluster_codes)):
            file_npi.write("Cluster " + str(i) + "\n")
            for j in range(len(cluster_codes[i])):
                file_npi.write(cluster_codes[i][j] + ": " + cluster_desc[i][j])
                file_npi.write("\n")
            file_npi.write("\n")
        # Closing file
        file_npi.close()

        npi_idx_new = np.argsort(npi_idx_to_cluster_idx[0])
        npis_corr_reorder = npis_corr[npi_idx_new, :][:, npi_idx_new]

        plt.imshow(abs(npis_corr_reorder), cmap='gray_r')
        plt.colorbar()

        # npi_indices_all = set(range(npis_corr.shape[0]))
        # for i in [40]:#[10, 20, 40, 80, 160]:
        #     kmeans_npis = KMeans(n_clusters=i).fit(df_npis_used.iloc[:,2:].T)
        #     evaluate_clustering(npis_corr, kmeans_npis.labels_, npi_indices_all)

        # for i in [40]:#[10, 20, 40, 80, 160]:
        #     kmeans_corr = KMeans(n_clusters=i).fit(npis_corr)
        #     evaluate_clustering(npis_corr, kmeans_corr.labels_, npi_indices_all)

        # corr_threshold = 0.5
        # corr_indices_threshold = np.where(npis_corr > corr_threshold)
        # npis_corr_threshold = np.zeros(npis_corr.shape)
        # npis_corr_threshold[corr_indices_threshold] = npis_corr[corr_indices_threshold]
        # plt.imshow(npis_corr_threshold, cmap='gray_r')

        # plot share of counties that implement the main categories
        if make_plot:
            # plot four different subsets of curves for better distinction
            j = 0
            if fine_resolution > 0:
                num_images = 15
            else:
                num_images = 1
            for i in [
                slice(
                    int(len(npi_codes_used) / num_images) * i,
                    min(
                        int(len(npi_codes_used) / num_images) *
                        (i + 1),
                        len(npis_used))) for i in range(
                    num_images + 1)]:
                customPlot.plotList(df_npis_aggregated.index,
                                    [df_npis_aggregated[code]
                                     for code in npi_codes_used[i]],
                                    npis_used[i],
                                    'Counties implementing NPI main categories',
                                    'Date', 'Number', "Counties_NPI_main_" +
                                    str(j) + "_of_"+str(num_images))
                j += 1


def main():
    """! Main program entry."""

    # arg_dict = gd.cli("testing")
    transform_npi_data(fine_resolution=2, read_data=False, make_plot=True)


if __name__ == "__main__":

    main()
