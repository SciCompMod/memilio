#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
import copy
import os
import pickle
import random
import json
from datetime import date

import numpy as np
import tensorflow as tf
from progress.bar import Bar
from sklearn.preprocessing import FunctionTransformer

from memilio.simulation import (AgeGroup, Damping, LogLevel, set_log_level)
from memilio.simulation.osecir import (Index_InfectionState,
                                       InfectionState, Model,
                                       interpolate_simulation_result, simulate)

from memilio.surrogatemodel.utils_surrogatemodel import (
    getBaselineMatrix, getMinimumMatrix, get_population)


def interpolate_age_groups(data_entry):
    """! Interpolates the age groups from the population data into the age groups used in the simulation. 
    We assume that the people in the age groups are uniformly distributed.
    @param data_entry Data entry containing the population data.
    @return List containing the population in each age group used in the simulation.
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
    """! Removes the confirmed compartments which are not used in the data generation.
    @param result_array Array containing the simulation results.
    @return Array containing the simulation results without the confirmed compartments.
    """
    num_groups = int(result_array.shape[1] / 10)
    delete_indices = [index for i in range(
        num_groups) for index in (3+10*i, 5+10*i)]
    return np.delete(result_array, delete_indices, axis=1)


def transform_data(data, transformer, num_runs):
    """! Transforms the data by a logarithmic normalization. 
    Reshaping is necessary, because the transformer needs an array with dimension <= 2.
    @param data Data to be transformed.
    @param transformer Transformer used for the transformation.
    @return Transformed data.
    """
    data = np.asarray(data).transpose(2, 0, 1).reshape(48, -1)
    scaled_data = transformer.transform(data)
    return tf.convert_to_tensor(scaled_data.transpose().reshape(num_runs, -1, 48))


def run_secir_groups_simulation(days, damping_day, populations):
    """! Uses an ODE SECIR model allowing for asymptomatic infection with 6 different age groups. The model is not stratified by region. 
    Virus-specific parameters are fixed and initial number of persons in the particular infection states are chosen randomly from defined ranges.
    @param Days Describes how many days we simulate within a single run.
    @param damping_day The day when damping is applied.
    @param populations List containing the population in each age group.
    @return List containing the populations in each compartment used to initialize the run.
   """
    set_log_level(LogLevel.Off)

    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1

    # Define age Groups
    groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']
    num_groups = len(groups)

    # Initialize Parameters
    model = Model(num_groups)

    # Set parameters
    for i in range(num_groups):
        # Compartment transition duration
        model.parameters.TimeExposed[AgeGroup(i)] = 3.2
        model.parameters.TimeInfectedNoSymptoms[AgeGroup(i)] = 2.
        model.parameters.TimeInfectedSymptoms[AgeGroup(i)] = 6.
        model.parameters.TimeInfectedSevere[AgeGroup(i)] = 12.
        model.parameters.TimeInfectedCritical[AgeGroup(i)] = 8.

        # Initial number of people in each compartment with random numbers
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Exposed)] = random.uniform(
            0.00025, 0.0005) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedNoSymptoms)] = random.uniform(
            0.0001, 0.00035) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedNoSymptomsConfirmed)] = 0
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSymptoms)] = random.uniform(
            0.00007, 0.0001) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSymptomsConfirmed)] = 0
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSevere)] = random.uniform(
            0.00003, 0.00006) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedCritical)] = random.uniform(
            0.00001, 0.00002) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Recovered)] = random.uniform(
            0.002, 0.008) * populations[i]
        model.populations[AgeGroup(i),
                          Index_InfectionState(InfectionState.Dead)] = 0
        model.populations.set_difference_from_group_total_AgeGroup(
            (AgeGroup(i), Index_InfectionState(InfectionState.Susceptible)), populations[i])

        # Compartment transition propabilities
        model.parameters.RelativeTransmissionNoSymptoms[AgeGroup(i)] = 0.5
        model.parameters.TransmissionProbabilityOnContact[AgeGroup(i)] = 0.1
        model.parameters.RecoveredPerInfectedNoSymptoms[AgeGroup(i)] = 0.09
        model.parameters.RiskOfInfectionFromSymptomatic[AgeGroup(i)] = 0.25
        model.parameters.SeverePerInfectedSymptoms[AgeGroup(i)] = 0.2
        model.parameters.CriticalPerSevere[AgeGroup(i)] = 0.25
        model.parameters.DeathsPerCritical[AgeGroup(i)] = 0.3
        # twice the value of RiskOfInfectionFromSymptomatic
        model.parameters.MaxRiskOfInfectionFromSymptomatic[AgeGroup(i)] = 0.5

    # StartDay is the n-th day of the year
    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

    # Load baseline and minimum contact matrix and assign them to the model
    baseline = getBaselineMatrix()
    minimum = getMinimumMatrix()

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = minimum

    # Generate a damping matrix and assign it to the model
    damping = np.ones((num_groups, num_groups)
                      ) * np.float16(random.uniform(0, 0.5))

    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=(damping), t=damping_day, level=0, type=0))

    damped_contact_matrix = model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
        damping_day+1)

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)

    # Interpolate simulation result on days time scale
    result = interpolate_simulation_result(result)

    result_array = remove_confirmed_compartments(
        np.transpose(result.as_ndarray()[1:, :]))

    # Omit first column, as the time points are not of interest here.
    dataset_entries = copy.deepcopy(result_array)

    return dataset_entries.tolist(), damped_contact_matrix


def generate_data(
        num_runs, path_out, path_population, input_width, label_width,
        normalize=True, save_data=True):
    """! Generate data sets of num_runs many equation-based model simulations and transforms the computed results by a log(1+x) transformation.
    Divides the results in input and label data sets and returns them as a dictionary of two TensorFlow Stacks.
    In general, we have 8 different compartments and 6 age groups.  If we choose, 
    input_width = 5 and label_width = 20, the dataset has 
    - input with dimension 5 x 8 x 6
    - labels with dimension 20 x 8 x 6
   @param num_runs Number of times, the function run_secir_groups_simulation is called.
   @param path_out Path, where the dataset is saved to.
   @param path_population Path, where we try to read the population data.
   @param input_width Int value that defines the number of time series used for the input.
   @param label_width Int value that defines the size of the labels.
   @param normalize [Default: true] Option to transform dataset by logarithmic normalization.
   @param save_data [Default: true] Option to save the dataset.
   @return Data dictionary of input and label data sets.
   """
    data = {
        "inputs": [],
        "labels": [],
        "contact_matrix": [],
        "damping_day": []
    }

    # The number of days is the same as the sum of input and label width.
    # Since the first day of the input is day 0, we still need to subtract 1.
    days = input_width + label_width - 1

    # Load population data
    # population = get_population(path_population)
    population = get_population()

    # show progess in terminal for longer runs
    # Due to the random structure, there's currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)
    for _ in range(0, num_runs):

        # Generate a random damping day
        damping_day = random.randrange(
            input_width, input_width+label_width)

        data_run, damped_contact_matrix = run_secir_groups_simulation(
            days, damping_day, population[random.randint(0, len(population) - 1)])
        data['inputs'].append(data_run[:input_width])
        data['labels'].append(data_run[input_width:])
        data['contact_matrix'].append(np.array(damped_contact_matrix))
        data['damping_day'].append(damping_day)
        bar.next()
    bar.finish()

    if normalize:
        # logarithmic normalization
        transformer = FunctionTransformer(np.log1p, validate=True)

        # transform inputs and labels
        data['inputs'] = transform_data(data['inputs'], transformer, num_runs)
        data['labels'] = transform_data(data['labels'], transformer, num_runs)
    else:
        data['inputs'] = tf.convert_to_tensor(data['inputs'])
        data['labels'] = tf.convert_to_tensor(data['labels'])

    if save_data:
        # check if data directory exists. If necessary, create it.
        if not os.path.isdir(path_out):
            os.mkdir(path_out)

        # save dict to json file
        with open(os.path.join(path_out, 'data_secir_groups.pickle'), 'wb') as f:
            pickle.dump(data, f)
    return data


if __name__ == "__main__":
    # Store data relative to current file two levels higher.
    path = os.path.dirname(os.path.realpath(__file__))
    path_output = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')

    path_population = os.path.abspath(
        r"data//pydata//Germany//county_population.json")

    input_width = 5
    label_width = 30
    num_runs = 10
    data = generate_data(num_runs, path_output, path_population, input_width,
                         label_width, save_data=True)
