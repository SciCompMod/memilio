import numpy as np
import pandas as pd
import os
from sklearn.preprocessing import FunctionTransformer

from memilio.epidata import transformMobilityData as tmd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.simulation.osecir import (ModelGraph, set_edges)
from memilio.epidata import modifyDataframeSeries as mdfs


def remove_confirmed_compartments(dataset_entries, num_groups):
    """! The compartments which contain confirmed cases are not needed and are 
        therefore omitted by summarizing the confirmed compartment with the 
        original compartment. 
    @param dataset_entries Array that contains the compartmental data with 
            confirmed compartments. 
    @param num_groups Number of age groups.
    @return Array that contains the compartmental data without confirmed compartments. 
   """

    new_dataset_entries = []
    for i in dataset_entries:
        dataset_entries_reshaped = i.reshape(
            [num_groups, int(np.asarray(dataset_entries).shape[1]/num_groups)]
        )
        sum_inf_no_symp = np.sum(dataset_entries_reshaped[:, [2, 3]], axis=1)
        sum_inf_symp = np.sum(dataset_entries_reshaped[:, [4, 5]], axis=1)
        dataset_entries_reshaped[:, 2] = sum_inf_no_symp
        dataset_entries_reshaped[:, 4] = sum_inf_symp
        new_dataset_entries.append(
            np.delete(dataset_entries_reshaped, [3, 5], axis=1).flatten()
        )
    return new_dataset_entries


def getBaselineMatrix():
    """! loads the baselinematrix
    """

    baseline_contact_matrix0 = os.path.join(
        "./data/contacts/baseline_home.txt")
    baseline_contact_matrix1 = os.path.join(
        "./data/contacts/baseline_school_pf_eig.txt")
    baseline_contact_matrix2 = os.path.join(
        "./data/contacts/baseline_work.txt")
    baseline_contact_matrix3 = os.path.join(
        "./data/contacts/baseline_other.txt")

    baseline = np.loadtxt(baseline_contact_matrix0) \
        + np.loadtxt(baseline_contact_matrix1) + \
        np.loadtxt(baseline_contact_matrix2) + \
        np.loadtxt(baseline_contact_matrix3)

    return baseline


def getMinimumMatrix():
    """! loads the minimum matrix
    """

    minimum_contact_matrix0 = os.path.join(
        "./data/contacts/minimum_home.txt")
    minimum_contact_matrix1 = os.path.join(
        "./data/contacts/minimum_school_pf_eig.txt")
    minimum_contact_matrix2 = os.path.join(
        "./data/contacts/minimum_work.txt")
    minimum_contact_matrix3 = os.path.join(
        "./data/contacts/minimum_other.txt")

    minimum = np.loadtxt(minimum_contact_matrix0) \
        + np.loadtxt(minimum_contact_matrix1) + \
        np.loadtxt(minimum_contact_matrix2) + \
        np.loadtxt(minimum_contact_matrix3)

    return minimum


def make_graph(directory, num_regions, countykey_list, models):
    """! 
    @param directory Directory with mobility data. 
    @param num_regions Number (int) of counties that should be added to the 
            grap-ODE model. Equals 400 for whole Germany. 
    @param countykey_list List of keys/IDs for each county. 
    @models models List of osecir Model with one model per population. 
    @return graph Graph-ODE model. 
   """
    graph = ModelGraph()
    for i in range(num_regions):
        graph.add_node(int(countykey_list[i]), models[i])

    num_locations = 4

    set_edges(os.path.abspath(os.path.join(directory, os.pardir)),
              graph, num_locations)
    return graph


def transform_mobility_directory():
    """! Transforms the mobility data by merging Eisenach and Wartburgkreis
    """
    # get mobility data directory
    arg_dict = gd.cli("commuter_official")

    directory = arg_dict['out_folder'].split('/pydata')[0]
    directory = os.path.join(directory, 'mobility/')

    # Merge Eisenach and Wartbugkreis in Input Data
    tmd.updateMobility2022(directory, mobility_file='twitter_scaled_1252')
    tmd.updateMobility2022(
        directory, mobility_file='commuter_migration_scaled')
    return directory


def get_population():
    df_population = pd.read_json(
        'data/pydata/Germany/county_population.json')
    age_groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80-130']

    df_population_agegroups = pd.DataFrame(
        columns=[df_population.columns[0]] + age_groups)
    for region_id in df_population.iloc[:, 0]:
        df_population_agegroups.loc[len(df_population_agegroups.index), :] = [int(region_id)] + list(
            mdfs.fit_age_group_intervals(df_population[df_population.iloc[:, 0] == int(region_id)].iloc[:, 2:], age_groups))

    population = df_population_agegroups.values.tolist()
    return population


def scale_data(data):
    num_groups = int(np.asarray(data['inputs']).shape[2] / 8)
    transformer = FunctionTransformer(np.log1p, validate=True)

    # Scale inputs
    inputs = np.asarray(
        data['inputs']).transpose(2, 0, 1, 3).reshape(num_groups * 8, -1)
    scaled_inputs = transformer.transform(inputs)
    original_shape_input = np.asarray(data['inputs']).shape

    # Reverse the reshape
    reshaped_back = scaled_inputs.reshape(original_shape_input[2],
                                          original_shape_input[0],
                                          original_shape_input[1],
                                          original_shape_input[3])

    # Reverse the transpose
    original_inputs = reshaped_back.transpose(1, 2, 0, 3)
    scaled_inputs = original_inputs.transpose(0, 3, 1, 2)

    # Scale labels
    labels = np.asarray(
        data['labels']).transpose(2, 0, 1, 3).reshape(num_groups * 8, -1)
    scaled_labels = transformer.transform(labels)
    original_shape_labels = np.asarray(data['labels']).shape

    # Reverse the reshape
    reshaped_back = scaled_labels.reshape(original_shape_labels[2],
                                          original_shape_labels[0],
                                          original_shape_labels[1],
                                          original_shape_labels[3])

    # Reverse the transpose
    original_labels = reshaped_back.transpose(1, 2, 0, 3)
    scaled_labels = original_labels.transpose(0, 3, 1, 2)

    return scaled_inputs, scaled_labels
