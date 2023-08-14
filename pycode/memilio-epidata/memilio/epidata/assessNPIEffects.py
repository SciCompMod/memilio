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
import matplotlib.cm as cm
import numpy as np
import pandas as pd
from memilio.epidata import customPlot
from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from scipy.cluster import hierarchy
from sklearn.metrics import silhouette_samples, silhouette_score

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

    ### old method
    '''
    # extract correlation values of block diagonals and offdiagonals separ.
    corr_diag = []
    corr_offdiag = []
    for ii in range(len(clusters)):
        corr_diag = np.append(corr_diag, abs(
            corr_mat[np.ix_(clusters[ii], clusters[ii])].flatten()))
        corr_offdiag = np.append(corr_offdiag, abs(
            corr_mat[np.ix_(clusters[ii], clusters_perp[ii])].flatten()))

    step_width=0.02
    corr_thresholds = np.linspace(0+step_width, 1, int(1/step_width))
    p_diag = np.array([len(np.where(corr_offdiag <= corr_thresholds[i])[0])
                      for i in range(len(corr_thresholds))])/len(corr_offdiag)*step_width
    p_offdiag = np.array([len(np.where(corr_diag < corr_thresholds[i])[0])
                         for i in range(len(corr_thresholds))])/len(corr_diag)*step_width


    return clusters, p_diag.sum()+p_offdiag.sum()
    '''
    ### new method
    if idx_to_cluster_idx.max() > 0:
        sample_silhouette_values = silhouette_samples(corr_mat, idx_to_cluster_idx)

    if (idx_to_cluster_idx.max() < 15):
        return clusters, -1
    else:
        return clusters, sample_silhouette_values.min()

# TODO: Used name 'methods' instead of 'metrics'. This is conform with documentation of hierarchy.linkage
    # and does not lead to confusion with metric used in pdist. To discuss!


def compute_hierarch_clustering(corr_pairwdist,
                                methods=['single', 'complete', 'average',
                                         'weighted', 'centroid', 'median',
                                         'ward']):
    """! Computes a hierarchical clustering for a (list of) method(s) and
    provides the maximum cophenetic distance(s) as well as a score for the
    clustering (see @method evaluate_clustering(...)).

    @param corr_mat correlation matrix between the features / data set items
        to be clustered hierarchically.
    @param corr_pairwdist Computed pairwise distance between the features / data
        set items.
    @param method method or list of methods to compute the hierarchical
        clustering.

    @return (List of) hierarchical clustering(s), maximum cophenetic distance(s)
        and scores of the hierarchical clustering.
    """

    # NOTE: if changing method, pay attention to linkage methods;
    #       'centroid', 'median', and 'ward' are correctly defined only if
    #       Euclidean pairwise metric is used in distance matrix that we used as input.

    # Based on the distances, we compute an hierarchical clustering for
    # different methods
    max_coph_corr_dist = 0
    scores = dict()
    # allow single entry
    if not isinstance(methods, list):
        methods = [methods]
    # iterate over list
    for method in methods:
        cluster_hierarch = hierarchy.linkage(corr_pairwdist, method=method)
        # compute cophenetic correlation distance and cophenetic distance matrix
        # TODO: Why was pdist(corr_mat) used as input for hierarchy.cophenet?
        # Shouldn't we use the same distance matrix as input both for linkage and cophenet?
        # To discuss!
        # TODO: < or > max_coph_corr_dist?
        coph_corr_dist, coph_dist_mat = hierarchy.cophenet(
            cluster_hierarch, corr_pairwdist)
        scores[method] = coph_corr_dist
        if coph_corr_dist > max_coph_corr_dist:
            max_coph_corr_dist = coph_corr_dist
            max_method = method
            max_coph_dist_mat = coph_dist_mat

    cluster_hierarch = hierarchy.linkage(corr_pairwdist, method=max_method)

    print(
        "Cophenetic correlation distance for method " + max_method + ": " +
        str(max_coph_corr_dist))

    return cluster_hierarch, max_coph_dist_mat


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

    total_eval_number = [[] for weight in weights]
    n=0
    # iterate over weights
    for weight in weights:
        # use the given weight to flatten the dendrogram
        npi_idx_to_cluster_idx = hierarchy.fcluster(
            cluster_hierarch, weight, criterion='distance')

        # evaluate clustering
        clusters, total_eval_number[n]=evaluate_clustering(corr_mat, npi_idx_to_cluster_idx, npi_indices_all)

        # append new npi_idx to cluster_idx assignment to list of assignments
        npi_idx_to_cluster_idx_list.append(npi_idx_to_cluster_idx)
        n+=1
    
    # print scores on clustering
    print("Number of clusters: " + str(len(clusters)) + "; evaluation number: " + str(round(np.nanmax(np.array(total_eval_number)), 4)))

    return npi_idx_to_cluster_idx_list[total_eval_number.index(np.nanmax(np.array(total_eval_number)))]

def silhouette(X, cluster_sizes, cluster_labels, label):

    if not isinstance(cluster_sizes, list):
        cluster_sizes = [cluster_sizes]

    plt.figure()

    for n_clusters in cluster_sizes:

        sample_silhouette_values = silhouette_samples(X, cluster_labels)

        y_lower = 10
        for i in range(n_clusters):
            ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

            ith_cluster_silhouette_values.sort()

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.hsv(float(i) / n_clusters)
            plt.fill_betweenx(
                np.arange(y_lower, y_upper),
                0,
                ith_cluster_silhouette_values,
                facecolor=color,
                edgecolor=color,
                alpha=0.7,
            )

            plt.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

            y_lower = y_upper + 10

        plt.xlabel("The silhouette coefficient values")
        plt.ylabel("Cluster label")

        plt.axvline(x=0, c='black', alpha = 0.3)
        plt.axvline(x=0.25, c='black', alpha = 0.3)
        plt.axvline(x=0.5, c='black', alpha = 0.3)
        plt.axvline(x=0.75, c='black', alpha = 0.3)

        plt.yticks([])  # Clear the yaxis labels / ticks
        plt.xticks([-0.4,-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])

        plt.title(
            "Silhouette analysis for " + label + " clustering with n_clusters = " + str(n_clusters),
        )

        plt.show(block=False)

        return sample_silhouette_values

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
    #df_npis = mdfs.extract_subframe_based_on_dates(df_npis, date(2021,1,1), date(2021,6,1))
    npis = pd.read_json(os.path.join(directory, 'npis.json'))
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
        # TODO: difference between scipy.cluster.hierarchy.distance.pdist and scipy.spatial.distance.pdist?
        corr_pairwdist = hierarchy.distance.pdist(
            npis_corr, metric='euclidean')
        
        npi_codes_used = np.asarray(df_npis_used.iloc[:, 2:].columns)

        # compute hierarchical clustering (via best-suited method)
        compare_methods = True
        if compare_methods:

            # centroid
            method = 'centroid'
            cluster_hierarch, coph_dist_mat = compute_hierarch_clustering(
                corr_pairwdist,
                method)
            # plot dendrogram
            plt.figure()
            plt.title(method)
            hierarchy.dendrogram(cluster_hierarch)
            # plt.show()
            max_coph_dist = coph_dist_mat.max()
            # TODO: Discuss why abs(npis_corr) is used as input and not corr_pairwdist
            npi_idx_to_cluster_idx = flatten_hierarch_clustering(
                abs(npis_corr), cluster_hierarch,
                [wg * max_coph_dist
                 for wg in np.linspace(0.01, 1, 500)])
            
            samples = silhouette(npis_corr, npi_idx_to_cluster_idx.max()+1, npi_idx_to_cluster_idx, label = method)

            # ward
            method = 'ward'
            cluster_hierarch, coph_dist_mat = compute_hierarch_clustering(
                corr_pairwdist,
                method)
            # plot dendrogram
            plt.figure()
            plt.title(method)
            hierarchy.dendrogram(cluster_hierarch)
            max_coph_dist = coph_dist_mat.max()
            npi_idx_to_cluster_idx = flatten_hierarch_clustering(
                abs(npis_corr), cluster_hierarch,
                [wg * max_coph_dist for wg in np.linspace(0.01, 1, 500)])
            
            samples = silhouette(npis_corr, npi_idx_to_cluster_idx.max()+1, npi_idx_to_cluster_idx, label = method)

            # average
            method = 'average'
            cluster_hierarch, coph_dist_mat = compute_hierarch_clustering(
                corr_pairwdist,
                method)
            # plot dendrogram
            plt.figure()
            plt.title(method)
            hierarchy.dendrogram(cluster_hierarch)
            max_coph_dist = coph_dist_mat.max()
            npi_idx_to_cluster_idx = flatten_hierarch_clustering(
                npis_corr, cluster_hierarch,
                [wg * max_coph_dist
                 for wg in np.linspace(0.01, 1, 500)])       
            
            samples = silhouette(npis_corr, npi_idx_to_cluster_idx.max()+1, npi_idx_to_cluster_idx, label = method)

        # centroid
        method = 'centroid'
        cluster_hierarch, coph_dist_mat = compute_hierarch_clustering(
            corr_pairwdist,
            method)
        # plot dendrogram
        plt.figure()
        plt.title(method)
        hierarchy.dendrogram(cluster_hierarch)
        plt.show()
        max_coph_dist = coph_dist_mat.max()
        npi_idx_to_cluster_idx = flatten_hierarch_clustering(
            npis_corr, cluster_hierarch,
            [wg * max_coph_dist
             for wg in [0.50]]) # 10 cluster

        cluster_dict = dict()
        cluster_codes = [[] for i in range(npi_idx_to_cluster_idx.max()+1)]
        cluster_desc = [[] for i in range(npi_idx_to_cluster_idx.max()+1)]
        for i in range(len(npi_idx_to_cluster_idx)):
            cluster_dict[npi_codes_used[i]
                         ] = "CM_" + str(npi_idx_to_cluster_idx[i]).zfill(3)
            cluster_codes[npi_idx_to_cluster_idx
                          [i]].append(npi_codes_used[i])
            cluster_desc[npi_idx_to_cluster_idx
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
                plt.tight_layout()
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
