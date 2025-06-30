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


def interpolate_age_groups(data_entry):
    """ Interpolates the age groups from the population data into the age groups used in the simulation.
    We assume that the people in the age groups are uniformly distributed.

    :param data_entry: Data entry containing the population data.
    :returns: List containing the population in each age group used in the simulation.

    """
    age_groups = {
        "A00-A04": data_entry['<3 years'] + data_entry['3-5 years'] * 2 / 3,
        "A05-A14": data_entry['3-5 years'] * 1 / 3 + data_entry['6-14 years'],
        "A15-A34": data_entry['15-17 years'] + data_entry['18-24 years'] + data_entry['25-29 years'] + data_entry['30-39 years'] * 1 / 2,
        "A35-A59": data_entry['30-39 years'] * 1 / 2 + data_entry['40-49 years'] + data_entry['50-64 years'] * 2 / 3,
        "A60-A79": data_entry['50-64 years'] * 1 / 3 + data_entry['65-74 years'] + data_entry['>74 years'] * 1 / 5,
        "A80+": data_entry['>74 years'] * 4 / 5
    }
    return [age_groups[key] for key in age_groups]


def remove_confirmed_compartments(result_array):
    """ Removes the confirmed compartments which are not used in the data generation.

    :param result_array: Array containing the simulation results.
    :returns: Array containing the simulation results without the confirmed compartments.

    """
    num_groups = int(result_array.shape[1] / 10)
    delete_indices = [index for i in range(
        num_groups) for index in (3+10*i, 5+10*i)]
    return np.delete(result_array, delete_indices, axis=1)


def transform_data(data, transformer, num_runs, num_groups=6, num_compartments=8):
    """ Transforms the data by a logarithmic normalization.
    Reshaping is necessary, because the transformer needs an array with dimension <= 2.

    :param data: Data to be transformed.
    :param transformer: Transformer used for the transformation.
    :param num_runs: Number of samples represented by the data. 
    :param num_groups: Number of age groups represented by data.
    :param num_compartments: Number of compartments. 
    :returns: Transformed data.

    """

    data = np.asarray(data).transpose(2, 0, 1).reshape(
        num_groups*num_compartments, -1)
    scaled_data = transformer.transform(data)
    return tf.convert_to_tensor(scaled_data.transpose().reshape(num_runs, -1, num_groups*num_compartments))


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
    # Wird noch gelöscht
    print(path_to_file)
    #####
    print("Model successfully saved")


def load_model(path):
    """
    Loading a trained model. 

    :param path: path to the .keras file containing the desired model
    :returns: trained tf.keras model 
    """
    if not os.path.isfile(path):
        raise FileExistsError(
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
