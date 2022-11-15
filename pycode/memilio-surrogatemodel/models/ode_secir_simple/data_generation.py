from cProfile import label
from memilio.simulation import UncertainContactMatrix, ContactMatrix, Damping
from memilio.simulation.secir import Model, simulate, AgeGroup, Index_InfectionState, Simulation, interpolate_simulation_result
from memilio.simulation.secir import InfectionState as State
import numpy as np
import pandas as pd
import pickle
from datetime import date
from math import ceil
import random
import os
import copy
from progress.bar import Bar  # pip install progess
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf
import matplotlib.pyplot as plt
import memilio.simulation as mio
from sklearn.preprocessing import FunctionTransformer


def run_secir_simulation(days):
    """! Slightly modified example of Secir Simple.

   Here, the initial values are choosen randomly, so the model is no longer deterministic.
   Necessary to create the dataset for later training process.
   This method is called within a loop in the function 'generate_data', which also sets the number of runs.

   """
    mio.set_log_level(mio.LogLevel.Off)

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
    model.parameters.IncubationTime[A0] = 5.2
    model.parameters.TimeInfectedSymptoms[A0] = 6.
    model.parameters.SerialInterval[A0] = 4.2
    model.parameters.TimeInfectedSevere[A0] = 12.
    model.parameters.TimeInfectedCritical[A0] = 8.

    # Initial number of people in each compartment with random numbers
    model.populations[A0, Index_InfectionState(
        State.Exposed)] = 60 * random.uniform(0.2, 1)
    model.populations[A0, Index_InfectionState(
        State.InfectedNoSymptoms)] = 55 * random.uniform(0.2, 1)
    model.populations[A0, Index_InfectionState(
        State.InfectedSymptoms)] = 50 * random.uniform(0.2, 1)
    model.populations[A0, Index_InfectionState(
        State.InfectedSevere)] = 12 * random.uniform(0.2, 1)
    model.populations[A0, Index_InfectionState(
        State.InfectedCritical)] = 3 * random.uniform(0.2, 1)
    model.populations[A0, Index_InfectionState(
        State.Recovered)] = 50 * random.random()
    model.populations[A0, Index_InfectionState(State.Dead)] = 0
    model.populations.set_difference_from_total(
        (A0, Index_InfectionState(State.Susceptible)), populations[0])

    # Compartment transition propabilities
    model.parameters.RelativeTransmissionNoSymptoms[A0] = 0.5
    model.parameters.TransmissionProbabilityOnContact[A0] = 1.0
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
        (num_groups, num_groups)) * 1
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
        (num_groups, num_groups)) * 0

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)
    result = interpolate_simulation_result(result)

    # when directly using a list instead of an array, we get problems with references
    result_array = result.as_ndarray()
    dataset = []

    # Omit first column, as the time points are not of interest here.
    for row in np.transpose(result_array[1:, :]):
        dataset.append(copy.deepcopy(row))

    return dataset


def generate_data(num_runs, path, input_width, label_width, save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often

    In general, we have 8 different compartments. If we choose, 
    input_width = 5 and label_width = 20, the dataset has 
    - input with dimension 5 x 8
    - labels with dimension 20 x 8

   @param num_runs Number of times, the function run_secir_simulation is called.
   @param path Path, where the dataset is saved to.
   @param input_width Int value that defines the number of time series used for the input.
   @param label_width Int value that defines the size of the labels.
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """
    data = {
        "inputs": [],
        "labels": []
    }

    days = input_width + label_width - 1

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)
    for _ in range(0, num_runs):
        data_run = run_secir_simulation(days)
        data["inputs"].append(data_run[:input_width])
        data['labels'].append(data_run[input_width:])
        bar.next()
    bar.finish()

    if save_data:

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

        # check if data directory exists. If necessary, create it.
        if not os.path.isdir(path):
            os.mkdir(path)

        # save dict to json file
        with open(os.path.join(path, 'data_secir_simple.pickle'), 'wb') as f:
            pickle.dump(data, f)


def splitdata(inputs, labels, split_train=0.7,
              split_valid=0.2, split_test=0.1):
    """! Split data in train, valid and test

   @param inputs input dataset
   @param labels label dataset
   @param split_train ratio of train datasets
   @param split_valid ratio of validation datasets
   @param split_test ratio of test datasets
   """

    if split_train + split_valid + split_test != 1:
        ValueError("summed Split ratios not equal 1! Please adjust the values")
    elif inputs.shape[0] != labels.shape[0] or inputs.shape[2] != labels.shape[2]:
        ValueError("Number of batches or features different for input and labels")

    n = inputs.shape[0]
    n_train = int(n * split_train)
    n_valid = int(n * split_valid)
    n_test = int(n * split_test)

    if n_train + n_valid + n_test != n:
        n_test = n - n_train - n_valid

    inputs_train, inputs_valid, inputs_test = tf.split(
        inputs, [n_train, n_valid, n_test], 0)
    labels_train, labels_valid, labels_test = tf.split(
        labels, [n_train, n_valid, n_test], 0)

    data = {
        "train_inputs": inputs_train,
        "train_labels": labels_train,
        "valid_inputs": inputs_valid,
        "valid_labels": labels_valid,
        "test_inputs": inputs_test,
        "test_labels": labels_test
    }

    return data


if __name__ == "__main__":
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')

    input_width = 5
    label_width = 30
    num_runs = 1000
    generate_data(num_runs, path_data, input_width,
                  label_width)
