#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Henrik Zunker
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
import numpy as np
import h5py
import pandas as pd
import memilio.simulation as mio
from memilio.simulation import (AgeGroup)

from memilio.simulation.osecir import (
    InfectionState, Model, ModelGraph, set_nodes)
from memilio.surrogatemodel.GNN.GNN_utils import remove_confirmed_compartments
from memilio.surrogatemodel.GNN.data_generation import (
    set_contact_matrices, set_covid_parameters, Location)


def extrapolate_data(start_date, num_days, data_dir):
    '''
    Extrapolate data using the C++ backend.
    This function sets up an empty model and calls the C++ backend to proceed the real world
    data. 
    :param start_date: The date to start the simulation from.
    :param num_days: The number of days to simulate.
    :param data_dir: The directory where the data is stored.
    :return: None
    '''
    # Set up the model
    model = Model(6)
    model.parameters.StartDay = start_date.day_in_year
    set_covid_parameters(model)
    set_contact_matrices(model, data_dir)

    graph = ModelGraph()

    # Set the parameters for the nodes
    scaling_factor_infected = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    scaling_factor_icu = 1.0
    tnt_capacity_factor = 7.5 / 100000.

    # Path containing the population data
    data_dir_Germany = os.path.join(data_dir, "Germany")
    pydata_dir = os.path.join(data_dir_Germany, "pydata")

    path_population_data = os.path.join(pydata_dir,
                                        "county_current_population.json")
    # Set the nodes in the model
    set_nodes(
        model.parameters,
        start_date, start_date + num_days,
        pydata_dir,
        path_population_data, True, graph, scaling_factor_infected,
        scaling_factor_icu, tnt_capacity_factor, num_days, True)


def read_results_h5(path, group_key='Total'):
    ''' Reads results from an HDF5 file and processes the data.
    :param path: Path to the HDF5 file.
    :param group_key: The key for the group to read from the HDF5 file.
    :return: A dictionary with processed results.
    '''
    with h5py.File(path, 'r') as f:
        keys = list(f.keys())
        res = {}
        for i, key in enumerate(keys):
            group = f[key]
            total = group[group_key][()]
            # remove confirmed compartments
            sum_inf_no_symp = np.sum(total[:, [2, 3]], axis=1)
            sum_inf_symp = np.sum(total[:, [4, 5]], axis=1)
            total[:, 2] = sum_inf_no_symp
            total[:, 4] = sum_inf_symp
            total = np.delete(total, [3, 5], axis=1)
            res[key] = total
    return res


def get_ground_truth_data(start_date, num_days, data_dir, create_new=False):
    '''
    Generates ground truth data for the specified start date and number of days.
    :param start_date: The date to start the simulation from.
    :param num_days: The number of days to simulate.
    :param data_dir: The directory where the data is stored.
    :param create_new: If True, generates new data; otherwise, reads existing data.
    :return: A numpy array containing the ground truth data.
    '''
    # Check if the path exists and create new data if needed
    path_rki_h5 = os.path.join(data_dir, "Germany", "pydata", "Results_rki.h5")
    if not os.path.isfile(path_rki_h5) or create_new:
        print("Generating real data from C++ backend...")
        extrapolate_data(start_date, num_days, data_dir)
        print("Data generation complete.")

    # age-stratified input.
    num_age_groups = 6
    all_age_data_list = []
    group_keys = ['Group1', 'Group2', 'Group3', 'Group4', 'Group5', 'Group6']
    for age in range(num_age_groups):
        age_data_dict = read_results_h5(path_rki_h5, group_keys[age])
        age_data_np = np.array(list(age_data_dict.values()))
        all_age_data_list.append(age_data_np)

    all_age_data = np.stack(all_age_data_list, axis=-1)
    num_nodes, timesteps, _, _ = all_age_data.shape
    ground_truth_all_nodes_np = np.reshape(
        all_age_data.transpose(0, 1, 3, 2), (num_nodes, timesteps, -1))

    return ground_truth_all_nodes_np


def extrapolate_ground_truth_data(data_dir, num_days=360):
    '''
    Extrapolates ground truth data for the specified number of days.
    :param data_dir: The directory where the data is stored.
    :param num_days: The number of days to simulate.
    :return: None
    '''
    cwd = os.getcwd()
    data_dir = os.path.join(cwd, data_dir)
    start_dates = [
        mio.Date(2020, 6, 1)
    ]

    for start_date in start_dates:
        get_ground_truth_data(
            start_date, num_days, data_dir, create_new=True)
        print(f"Ground truth data for {start_date} generated.")

    path_rki_h5 = data_dir + "/Germany/pydata/Results_rki.h5"
    num_age_groups = 6
    all_age_data_list = []
    group_keys = ['Group1', 'Group2', 'Group3', 'Group4', 'Group5', 'Group6']
    for age in range(num_age_groups):
        age_data_dict = read_results_h5(path_rki_h5, group_keys[age])
        age_data_np = np.array(list(age_data_dict.values()))
        all_age_data_list.append(age_data_np)
    all_age_data = np.stack(all_age_data_list, axis=-1)

    # Save the ground truth data
    output_path = data_dir + "/Germany/pydata/ground_truth_all_nodes.pickle"
    with open(output_path, 'wb') as f:
        pd.to_pickle(all_age_data, f)
    print(f"Ground truth data saved to {output_path}")


def generate_bounds(data_dir):
    path = data_dir + "/Germany/pydata/ground_truth_all_nodes.pickle"
    with open(path, 'rb') as f:
        ground_truth_all_nodes_np = pd.read_pickle(f)

    # Calculate the bounds
    lower_bound = np.min(ground_truth_all_nodes_np, axis=1)
    upper_bound = np.max(ground_truth_all_nodes_np, axis=1)
    path_upper_bound = data_dir + "/Germany/pydata/ground_truth_upper_bound.pickle"
    path_lower_bound = data_dir + "/Germany/pydata/ground_truth_lower_bound.pickle"
    with open(path_upper_bound, 'wb') as f:
        pd.to_pickle(upper_bound, f)
    with open(path_lower_bound, 'wb') as f:
        pd.to_pickle(lower_bound, f)

    print(f"Upper bound saved to {path_upper_bound}")
    print(f"Lower bound saved to {path_lower_bound}")


def main():
    cwd = os.getcwd()
    num_days = 180
    data_dir = os.path.join(cwd, "data")

    extrapolate_ground_truth_data(
        data_dir,
        num_days
    )
    generate_bounds(data_dir)


if __name__ == "__main__":
    main()
