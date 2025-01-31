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
    remove_confirmed_compartments, get_population, getBaselineMatrix)


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
    """! Uses an ODE SECIR model allowing for asymptomatic infection with 6 different age groups. 
    The model is not stratified by region. 
    Virus-specific parameters are fixed and initial number of persons in the particular 
    infection states are chosen randomly from defined ranges.
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
    groups_id = np.arange(num_groups)

    # age specific parameters
    TransmissionProbabilityOnContact = [0.03, 0.06, 0.06, 0.06, 0.09, 0.175]
    RecoveredPerInfectedNoSymptoms = [0.25, 0.25, 0.2, 0.2, 0.2, 0.2]
    SeverePerInfectedSymptoms = [0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225]
    CriticalPerSevere = [0.075, 0.075, 0.075, 0.15, 0.3, 0.4]
    DeathsPerCritical = [0.05, 0.05, 0.14, 0.14, 0.4, 0.6]

    TimeInfectedNoSymptoms = [2.74, 2.74, 2.565, 2.565, 2.565, 2.565]
    TimeInfectedSymptoms = [7.02625, 7.02625,
                            7.0665, 6.9385, 6.835, 6.775]
    TimeInfectedSevere = [5, 5, 5.925, 7.55, 8.5, 11]
    TimeInfectedCritical = [6.95, 6.95, 6.86, 17.36, 17.1, 11.6]

    # Initialize Parameters
    model = Model(num_groups)

    # Set parameters
    for i, rho, muCR, muHI, muUH, muDU, tc, ti, th, tu in zip(groups_id, TransmissionProbabilityOnContact, RecoveredPerInfectedNoSymptoms, SeverePerInfectedSymptoms, CriticalPerSevere, DeathsPerCritical, TimeInfectedNoSymptoms, TimeInfectedSymptoms, TimeInfectedSevere, TimeInfectedCritical):
        # Compartment transition duration
        model.parameters.TimeExposed[AgeGroup(i)] = 3.335
        model.parameters.TimeInfectedNoSymptoms[AgeGroup(i)] = tc
        model.parameters.TimeInfectedSymptoms[AgeGroup(i)] = ti
        model.parameters.TimeInfectedSevere[AgeGroup(i)] = th
        model.parameters.TimeInfectedCritical[AgeGroup(i)] = tu

        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSymptoms)] = populations[i] * random.uniform(0.0001, 0.05)

        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Exposed)] = model.populations[AgeGroup(i),
                                                         InfectionState.InfectedSymptoms].value * random.uniform(0.1, 5)

        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedNoSymptoms)] = model.populations[AgeGroup(i),
                                                                    InfectionState.InfectedSymptoms].value * random.uniform(0.1, 5)

        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSevere)] = model.populations[AgeGroup(i),
                                                                InfectionState.InfectedSymptoms].value * random.uniform(0.001, 1)

        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedCritical)] = model.populations[AgeGroup(i),
                                                                  InfectionState.InfectedSevere].value * random.uniform(0.001, 1)

        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Dead)] = model.populations[AgeGroup(i),
                                                      InfectionState.InfectedCritical].value * random.uniform(0.001, 1)

        subtotal = (model.populations[AgeGroup(i), InfectionState.InfectedSymptoms].value
                    + model.populations[AgeGroup(i),
                                        InfectionState.Exposed].value
                    + model.populations[AgeGroup(i),
                                        InfectionState.InfectedNoSymptoms].value
                    + model.populations[AgeGroup(i), InfectionState.InfectedSevere].value
                    + model.populations[AgeGroup(i), InfectionState.InfectedCritical].value
                    + model.populations[AgeGroup(i), InfectionState.Dead].value
                    )

        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Recovered)] = random.uniform(0, (populations[i]-subtotal))

        # print('get_population_ageGroups: ' +
        #      str(model.populations.get_group_total_AgeGroup(AgeGroup(i))))
        model.populations.set_difference_from_group_total_AgeGroup(
            (AgeGroup(i), Index_InfectionState(InfectionState.Susceptible)), populations[i])

        # Compartment transition propabilities
        model.parameters.RelativeTransmissionNoSymptoms[AgeGroup(i)] = 1
        model.parameters.TransmissionProbabilityOnContact[AgeGroup(i)] = rho
        model.parameters.RecoveredPerInfectedNoSymptoms[AgeGroup(i)] = muCR
        model.parameters.RiskOfInfectionFromSymptomatic[AgeGroup(i)] = 0.25
        model.parameters.SeverePerInfectedSymptoms[AgeGroup(i)] = muHI
        model.parameters.CriticalPerSevere[AgeGroup(i)] = muUH
        model.parameters.DeathsPerCritical[AgeGroup(i)] = muDU
        # twice the value of RiskOfInfectionFromSymptomatic
        model.parameters.MaxRiskOfInfectionFromSymptomatic[AgeGroup(i)] = 0.5

    # StartDay is the n-th day of the year
    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

    # Load baseline and minimum contact matrix and assign them to the model
    baseline = getBaselineMatrix()

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline

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

    return dataset_entries.tolist(), damped_contact_matrix, damping


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
        "damping_day": [],
        "damping_coeff": []
    }

    # The number of days is the same as the sum of input and label width.
    # Since the first day of the input is day 0, we still need to subtract 1.
    days = input_width + label_width - 1

    # Load population data
    population = get_population(path_population)
    population = np.asarray(population).sum(axis=0)

    # show progess in terminal for longer runs
    # Due to the random structure, there's currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)
    for _ in range(0, num_runs):

        # Generate a random damping day
        damping_day = random.randrange(
            input_width, input_width+label_width)

        data_run, damped_contact_matrix, damping_coeff = run_secir_groups_simulation(
            days, damping_day, population)
        data['inputs'].append(data_run[:input_width])
        data['labels'].append(data_run[input_width:])
        data['contact_matrix'].append(np.array(damped_contact_matrix))
        data['damping_day'].append(damping_day)
        data['damping_coeff'].append(damping_coeff)
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
        with open(os.path.join(path_out, f'data_secir_groups_{label_width}days_I_based_Germany_100k.pickle'), 'wb') as f:
            pickle.dump(data, f)
    return data


if __name__ == "__main__":
    # Store data relative to current file two levels higher.
    cwd = os.getcwd()
    path_output = os.path.join(cwd, "saves")

    path_population = cwd + "//data//pydata//Germany//county_current_population.json"

    input_width = 5
    label_width = 90
    num_runs = 100000
    data = generate_data(num_runs, path_output, path_population, input_width,
                         label_width)
