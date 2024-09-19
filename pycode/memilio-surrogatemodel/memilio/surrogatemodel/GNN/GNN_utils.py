import numpy as np
import pandas as pd
import os

from sklearn.preprocessing import FunctionTransformer

from memilio.epidata import transformMobilityData as tmd
from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.simulation.osecir import (ModelGraph, set_edges)


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
