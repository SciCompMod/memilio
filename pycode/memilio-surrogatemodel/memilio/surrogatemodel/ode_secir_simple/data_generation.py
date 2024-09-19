#############################################################################
# Copyright (C) 2020-2024 MEmilio
#
# Authors: Agatha Schmidt, Henrik Zunker, Khoa Nguyen
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
import copy
import os
import pickle
import random
from datetime import date

import numpy as np
import pandas as pd
import tensorflow as tf
from progress.bar import Bar
from sklearn.preprocessing import FunctionTransformer

from memilio.simulation import (AgeGroup, ContactMatrix, Damping, LogLevel,
                                UncertainContactMatrix, set_log_level)
from memilio.simulation.osecir import (Index_InfectionState,
                                       InfectionState, Model, Simulation,
                                       interpolate_simulation_result, simulate)

from memilio.surrogatemodel.utils_surrogatemodel import remove_confirmed_compartments


def run_secir_simple_simulation(days):
    """! Uses an ODE SECIR model allowing for asymptomatic infection. The model is not stratified by region or demographic properties such as age.
    Virus-specific parameters are fixed and initial number of persons in the particular infection states are chosen randomly from defined ranges.

    @param Days Describes how many days we simulate within a single run.
    @return List containing the populations in each compartment for each day of the simulation.
   """
    set_log_level(LogLevel.Off)

    populations = [50_000]
    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1
    num_groups = 1

    # Initialize Parameters
    model = Model(1)

    A0 = AgeGroup(0)

    # Set parameters
    # Compartment transition duration
    model.parameters.TimeExposed[A0] = 3.2
    model.parameters.TimeInfectedNoSymptoms[A0] = 2.
    model.parameters.TimeInfectedSymptoms[A0] = 6.
    model.parameters.TimeInfectedSevere[A0] = 12.
    model.parameters.TimeInfectedCritical[A0] = 8.

    # Initial number of people in each compartment with random numbers
    model.populations[A0, Index_InfectionState(
        InfectionState.Exposed)] = 60 * random.uniform(0.2, 1)
    model.populations[A0, Index_InfectionState(
        InfectionState.InfectedNoSymptoms)] = 55 * random.uniform(0.2, 1)
    model.populations[A0, Index_InfectionState(
        InfectionState.InfectedSymptoms)] = 50 * random.uniform(0.2, 1)
    model.populations[A0, Index_InfectionState(
        InfectionState.InfectedSevere)] = 12 * random.uniform(0.2, 1)
    model.populations[A0, Index_InfectionState(
        InfectionState.InfectedCritical)] = 3 * random.uniform(0.2, 1)
    model.populations[A0, Index_InfectionState(
        InfectionState.Recovered)] = 50 * random.random()
    model.populations[A0, Index_InfectionState(InfectionState.Dead)] = 0
    model.populations.set_difference_from_total(
        (A0, Index_InfectionState(InfectionState.Susceptible)), populations[0])

    # Compartment transition propabilities
    model.parameters.RelativeTransmissionNoSymptoms[A0] = 0.5
    model.parameters.TransmissionProbabilityOnContact[A0] = 0.1
    model.parameters.RecoveredPerInfectedNoSymptoms[A0] = 0.09
    model.parameters.RiskOfInfectionFromSymptomatic[A0] = 0.25
    model.parameters.SeverePerInfectedSymptoms[A0] = 0.2
    model.parameters.CriticalPerSevere[A0] = 0.25
    model.parameters.DeathsPerCritical[A0] = 0.3
    # twice the value of RiskOfInfectionFromSymptomatic
    model.parameters.MaxRiskOfInfectionFromSymptomatic[A0] = 0.5

    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (num_groups, num_groups)) * 10
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
        (num_groups, num_groups)) * 0

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)
    # Interpolate simulation result on days time scale
    result = interpolate_simulation_result(result)

    # Using an array instead of a list to avoid problems with references
    result_array = result.as_ndarray()

    result_array = remove_confirmed_compartments(
        result_array[1:, :].transpose(), 1)

    dataset_entries = copy.deepcopy(result_array)

    return dataset_entries


def generate_data(
        num_runs, path, input_width, label_width, normalize=True,
        save_data=True):
    """! Generate data sets of num_runs many equation-based model simulations and transforms the computed results by a log(1+x) transformation.
    Divides the results in input and label data sets and returns them as a dictionary of two TensorFlow Stacks.

    In general, we have 10 different compartments. However, we aggregate the InfectedNoSymptoms and InfectedSymptomsNoConfirmed compartments. The same
    holds for the InfectedSymptoms and InfectedSymptomsConfirmed compartments. So, we end up with only 8 different compartments. If we choose, 
    input_width = 5 and label_width = 20, the dataset has 
    - input with dimension 5 x 8
    - labels with dimension 20 x 8

   @param num_runs Number of times, the function run_secir_simple_simulation is called.
   @param path Path, where the dataset is saved to.
   @param input_width Int value that defines the number of time series used for the input.
   @param label_width Int value that defines the size of the labels.
   @param normalize [Default: true] Option to transform dataset by logarithmic normalization.
   @param save_data [Default: true] Option to save the dataset.
   @return Data dictionary of input and label data sets.
   """
    data = {
        "inputs": [],
        "labels": []
    }

    # The number of days is the same as the sum of input and label width.
    # Since the first day of the input is day 0, we still need to subtract 1.
    days = input_width + label_width - 1

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)
    for _ in range(0, num_runs):
        data_run = run_secir_simple_simulation(days)
        data['inputs'].append(data_run[:input_width])
        data['labels'].append(data_run[input_width:])
        bar.next()
    bar.finish()

    if normalize:
        # logarithmic normalization
        transformer = FunctionTransformer(np.log1p, validate=True)
        inputs = np.asarray(data['inputs']).transpose(2, 0, 1).reshape(8, -1)
        scaled_inputs = transformer.transform(inputs)
        scaled_inputs = scaled_inputs.transpose().reshape(num_runs, input_width, 8)
        scaled_inputs_list = scaled_inputs.tolist()

        labels = np.asarray(data['labels']).transpose(2, 0, 1).reshape(8, -1)
        scaled_labels = transformer.transform(labels)
        scaled_labels = scaled_labels.transpose().reshape(num_runs, label_width, 8)
        scaled_labels_list = scaled_labels.tolist()

        # cast dfs to tensors
        data['inputs'] = tf.stack(scaled_inputs_list)
        data['labels'] = tf.stack(scaled_labels_list)

    if save_data:
        # check if data directory exists. If necessary, create it.
        if not os.path.isdir(path):
            os.mkdir(path)

        # save dict to json file
        with open(os.path.join(path, 'data_secir_simple.pickle'), 'wb') as f:
            pickle.dump(data, f)
    return data


if __name__ == "__main__":
    # Store data relative to current file two levels higher.
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')

    input_width = 5
    label_width = 30
    num_runs = 10000
    data = generate_data(num_runs, path_data, input_width,
                         label_width, save_data=True)
