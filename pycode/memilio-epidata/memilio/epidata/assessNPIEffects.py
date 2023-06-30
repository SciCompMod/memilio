#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from memilio.epidata import customPlot
from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist


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
        # compute cophenetic correlation distance
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


def analyze_npi_data(
        read_data, make_plot, fine_resolution, npis, directory, file_format,
        npi_codes_considered):

    if not read_data:
        raise FileNotFoundError('')
        # transform_npi_data(fine_resolution=2,
        #                file_format=dd.defaultDict['file_format'],
        #                out_folder=dd.defaultDict['out_folder'],
        #                start_date=dd.defaultDict['start_date'],
        #                end_date=dd.defaultDict['end_date'],
        #                make_plot=dd.defaultDict['make_plot'],
        #                )

    else:  # read formatted file
        npis = pd.read_excel(
            os.path.join(
                directory, 'datensatzbeschreibung_massnahmen.xlsx'),
            sheet_name=3, engine='openpyxl')
        if fine_resolution > 0:
            if fine_resolution == 1:
                filename = 'germany_counties_npi_subcat_incgrouped'
            else:
                filename = 'germany_counties_npi_subcat'
        else:
            filename = 'germany_counties_npi_maincat'
        df_npis = pd.read_csv(directory + filename + ".csv")
    try:
        df_npis = df_npis.drop('Unnamed: 0', axis=1)
    except KeyError:
        pass
    npis = pd.read_json(os.path.join(directory, 'npis.json'))
    # merge subcodes of considered npis (Do we want this? To discuss...)
        # get code levels (main/subcodes) and position of main codes
        # code_level = [i.count('_') for i in npi_codes]
        # main_code_pos = [i for i in range(len(code_level)) if code_level[i] == 1]

    # check if any other integer than 0: not implemented or 1: implemented is
    # used (maybe to specify the kind of implementation)

    npi_codes_considered = npis.NPI_code.values.tolist()#[x for x in npis.NPI_code if len(x) <=8]

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
                    np.where(npis[dd.EngEng['npiCode']] == 'M01a_010')[0][0])
            else:
                npi_used_indices.append(i)
                npi_used_indices_all.append(
                    np.where(npis[dd.EngEng['npiCode']] == 'M01a_010')[0][0])

        npis_unused = np.array(npis['Description'])[npi_unused_indices_all] # for all fr1 codes: len=153
        npis_used = np.array(npis['Description'])[npi_used_indices_all] # for all fr1 codes: len=14
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
            npis_corr, # same result as abs(npis_corr)
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
             for wg in [0.50]]) # 10 cluster

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
    fine_resolution = 2
    npis_final = []
    directory = os.path.join(dd.defaultDict['out_folder'], 'Germany/')
    file_format = 'json'
    npi_codes_considered = ['M01a_010', 'M01a_020']
    analyze_npi_data(True, True, fine_resolution, npis_final,
                     directory, file_format, npi_codes_considered)


if __name__ == "__main__":

    main()
