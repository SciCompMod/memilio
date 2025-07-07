#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Manuel Heger, Henrik Zunker
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

import numpy as np
import tensorflow as tf
import os
from memilio.epidata.modifyDataframeSeries import (fit_age_group_intervals)
import pandas as pd


def interpolate_age_groups(data_entry):
    """ Interpolates the age groups from the population data into the age groups used in the simulation.
    We assume that the people in the age groups are uniformly distributed.

    :param data_entry: Data entry containing the population data.
    :returns: List containing the population in each age group used in the simulation.

    """
    df_age_in = pd.DataFrame({
        '<3': 1, '3-5': 1, '6-14': 1, '15-17': 1, '18-24': 1, '25-29': 1,
        '30-39': 1, '40-49': 1, '50-64': 1, '64-74': 1, '>74': 1
    })
    age_out = ["<5", "5-14", "15-34", "35-59", "60-79", ">79"]

    return (fit_age_group_intervals(df_age_in, age_out, data_entry))


def remove_confirmed_compartments(result_array, delete_indices):
    """ Removes the confirmed compartments which are not used in the data generation.

    :param result_array: Array containing the simulation results.
    :param delete_indices: flat indices indicating position containing data from confirmed compartments. 
    :returns: Array containing the simulation results without the confirmed compartments.

    """
    return np.delete(result_array, delete_indices, axis=1)


def normalize_simulation_data(data, transformer, num_runs, num_groups=6, num_compartments=8):
    """ Transforms the data by a logarithmic normalization.
    Reshaping is necessary, because the transformer needs an array with dimension <= 2.

    :param data: Data to be transformed.
    :param transformer: Transformer used for the transformation.
    :param num_runs: Number of samples represented by the data.
    :param num_groups: Number of age groups represented by data.
    :param num_compartments: Number of compartments. 
    :returns: Transformed data.

    """
    data = np.asarray(data).reshape(num_runs, -1)
    scaled_data = transformer.transform(data)
    return tf.convert_to_tensor(scaled_data.reshape(num_runs, -1, num_groups*num_compartments))


def save_model(model, path, modelname):
    """
    Saving a trained model. 

    :param model: trained tensorflow keras model 
    :param path: path where the model should be stored 
    :param modelname: the name of the model 
    """
    if not os.path.isdir(path):
        os.mkdir(path)
    path_to_file = os.path.join(path, modelname + ".keras")
    model.save(path_to_file)
    print("Model successfully saved")


def load_model(path):
    """
    Loading a trained model. 

    :param path: path to the .keras file containing the desired model
    :returns: trained tf.keras model 
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(
            "There is no .keras model stored at the given directory.")
    return tf.keras.models.load_model(path)


def calc_split_index(n, split_train=0.7,
                     split_valid=0.2, split_test=0.1):
    """
    Calculating the indixes for a split_train:split_valid:split_test decomposition of a set with size n 

    It must hold split_train + split_valid + split_test = 1

    :param n: integer value 
    :param split_train: value between 0 and 1
    :param split_valid: value between 0 and 1
    :param split_test: value between 0 and 1
    :returns: a list of the form [i_train, i_valid, i_test]
    """
    if split_train + split_valid + split_test > 1 + 1e-10:
        raise ValueError(
            "Summed data set shares are greater than 1. Please adjust the values.")
    n_train = int(n * split_train)
    n_valid = int(n * split_valid)
    n_test = n - n_train - n_valid

    return [n_train, n_valid, n_test]


def flat_input(input):
    """ Flatten input dimension

    :param input: input array of size (n,k,l)
    :returns: reshaped array of size (n, k*l)

    """
    dim = tf.reduce_prod(tf.shape(input)[1:])
    return tf.reshape(input, [-1, dim])
