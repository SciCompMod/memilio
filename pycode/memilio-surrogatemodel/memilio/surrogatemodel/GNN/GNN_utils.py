#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Agatha Schmidt, Henrik Zunker
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
Utility functions for GNN-based surrogate models.

This module provides helper functions for data preprocessing, transformation,
and graph construction used in Graph Neural Network surrogate models for
epidemiological simulations.
"""

import os

import numpy as np
from sklearn.preprocessing import FunctionTransformer

from memilio.simulation.osecir import ModelGraph, set_edges


# Default number of compartments in the ODE-SECIR model (without confirmed compartments)
DEFAULT_NUM_COMPARTMENTS = 8


def remove_confirmed_compartments(dataset_entries, num_groups):
    """Removes confirmed compartments from simulation data by merging them with base compartments.

    The ODE-SECIR model includes separate "confirmed" compartments that track
    diagnosed cases. For GNN training, these are merged back into their base
    compartments (InfectedNoSymptoms + InfectedNoSymptomsConfirmed -> InfectedNoSymptoms,
    InfectedSymptoms + InfectedSymptomsConfirmed -> InfectedSymptoms).

    :param dataset_entries: Compartment data array containing confirmed compartments.
                           Shape: [num_timesteps, num_groups * num_compartments_with_confirmed]
    :param num_groups: Number of age groups in the model.
    :returns: Array with confirmed compartments merged into base compartments.
             Shape: [num_timesteps, num_groups * num_compartments_without_confirmed]

    """

    new_dataset_entries = []
    for timestep_data in dataset_entries:
        # Reshape to separate age groups and compartments
        data_reshaped = timestep_data.reshape(
            [num_groups, int(np.asarray(
                dataset_entries).shape[1] / num_groups)]
        )

        # Merge InfectedNoSymptoms (index 2) with InfectedNoSymptomsConfirmed (index 3)
        sum_infected_no_symptoms = np.sum(data_reshaped[:, [2, 3]], axis=1)

        # Merge InfectedSymptoms (index 4) with InfectedSymptomsConfirmed (index 5)
        sum_infected_symptoms = np.sum(data_reshaped[:, [4, 5]], axis=1)

        # Replace original compartments with merged values
        data_reshaped[:, 2] = sum_infected_no_symptoms
        data_reshaped[:, 4] = sum_infected_symptoms

        # Remove confirmed compartments (indices 3 and 5) and flatten
        new_dataset_entries.append(
            np.delete(data_reshaped, [3, 5], axis=1).flatten()
        )

    return new_dataset_entries


def get_baseline_contact_matrix(data_dir):
    """Loads and sums baseline contact matrices for all location types.

    Loads contact matrices for home, school, work, and other locations, then
    returns their sum as the total baseline contact matrix.

    :param data_dir: Root directory containing contact matrix data.
    :returns: Combined baseline contact matrix as numpy array.
    :raises FileNotFoundError: If any contact matrix file is not found.

    """
    contact_dir = os.path.join(data_dir, "Germany", "contacts")

    contact_files = [
        "baseline_home.txt",
        "baseline_school_pf_eig.txt",
        "baseline_work.txt",
        "baseline_other.txt"
    ]

    baseline_matrix = None

    for filename in contact_files:
        filepath = os.path.join(contact_dir, filename)

        if not os.path.exists(filepath):
            raise FileNotFoundError(
                f"Contact matrix file not found: {filepath}")

        matrix = np.loadtxt(filepath)

        if baseline_matrix is None:
            baseline_matrix = matrix
        else:
            baseline_matrix += matrix

    return baseline_matrix


def create_mobility_graph(mobility_dir, num_regions, county_ids, models):
    """Creates a graph-ODE model with mobility connections between regions.

    Constructs a graph, where each node represents a region (county) with
    its own ODE-SECIR model, and edges represent mobility flows between regions.

    :param mobility_dir: Directory containing mobility data files.
    :param num_regions: Number of regions/counties to include in the graph.
    :param county_ids: List of county IDs/keys for each region.
    :param models: List of ODE-SECIR Model instances, one per county.
    :returns: Configured graph with nodes and mobility edges.

    """
    # Initialize empty graph
    graph = ModelGraph()

    # Add one node per region with its model
    for i in range(num_regions):
        graph.add_node(int(county_ids[i]), models[i])

    # Number of contact locations (home, school, work, other)
    num_contact_locations = 4

    # Add mobility edges between nodes
    parent_dir = os.path.abspath(os.path.join(mobility_dir, os.pardir))
    set_edges(parent_dir, graph, num_contact_locations)

    return graph


def scale_data(
        data, transform=True, num_compartments=DEFAULT_NUM_COMPARTMENTS):
    """Applies logarithmic transformation to simulation data for training.

    Transforms compartment data using log1p (log(1+x)) to stabilize training.
    This transformation helps handle the wide range of population values and
    improves gradient flow during neural network training.

    :param data: Dictionary containing 'inputs' and 'labels' keys with simulation data.
                Inputs shape: [num_samples, time_steps, num_nodes, num_features]
                Labels shape: [num_samples, time_steps, num_nodes, num_features]
    :param transform: Whether to apply log transformation (True) or just reshape (False).
    :param num_compartments: Number of compartments per age group in the model (default: 8).
    :returns: Tuple of (scaled_inputs, scaled_labels), both with shape 
             [num_samples, num_features, time_steps, num_nodes]
    :raises ValueError: If input data is not numeric or has unexpected structure.

    """
    # Validate input data types
    inputs_array = np.asarray(data['inputs'])
    labels_array = np.asarray(data['labels'])

    if not np.issubdtype(inputs_array.dtype, np.number):
        raise ValueError("Input data must be numeric.")
    if not np.issubdtype(labels_array.dtype, np.number):
        raise ValueError("Label data must be numeric.")

    # Calculate number of age groups from data shape
    num_groups = int(inputs_array.shape[-1] / num_compartments)

    # Initialize transformer (log1p for numerical stability)
    transformer = FunctionTransformer(np.log1p, validate=True)

    # Process inputs
    # Reshape: [samples, timesteps, nodes, features] -> [nodes, samples, timesteps, features]
    #       -> [nodes * compartments, samples * timesteps * age_groups]
    inputs_reshaped = inputs_array.transpose(2, 0, 1, 3).reshape(
        num_groups * num_compartments, -1
    )

    if transform:
        inputs_transformed = transformer.transform(inputs_reshaped)
    else:
        inputs_transformed = inputs_reshaped

    original_shape_input = inputs_array.shape

    # Reverse reshape to separate dimensions
    inputs_back = inputs_transformed.reshape(
        original_shape_input[2],
        original_shape_input[0],
        original_shape_input[1],
        original_shape_input[3]
    )

    # Reverse transpose and reorder to [samples, features, timesteps, nodes]
    scaled_inputs = inputs_back.transpose(1, 2, 0, 3).transpose(0, 3, 1, 2)

    # Process labels with same procedure
    labels_reshaped = labels_array.transpose(2, 0, 1, 3).reshape(
        num_groups * num_compartments, -1
    )

    if transform:
        labels_transformed = transformer.transform(labels_reshaped)
    else:
        labels_transformed = labels_reshaped

    original_shape_labels = labels_array.shape

    labels_back = labels_transformed.reshape(
        original_shape_labels[2],
        original_shape_labels[0],
        original_shape_labels[1],
        original_shape_labels[3]
    )

    scaled_labels = labels_back.transpose(1, 2, 0, 3).transpose(0, 3, 1, 2)

    return scaled_inputs, scaled_labels
