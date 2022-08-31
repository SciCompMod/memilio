from cProfile import label
import json
import pickle
from memilio.simulation import UncertainContactMatrix, ContactMatrix, Damping
from memilio.simulation.secir import SecirModel, simulate, AgeGroup, Index_InfectionState, SecirSimulation
from memilio.simulation.secir import InfectionState as State
import numpy as np
import pandas as pd
from datetime import date
from math import ceil
import random
import os
from progress.bar import Bar  # pip install progess
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf
import matplotlib.pyplot as plt
import memilio.simulation as mio
from datetime import datetime, date

import seaborn as sns  # plot after normalization


def run_secir_groups_simulation(days, damping_data, populations):
    """
    Runs the c++ secir model using mulitple age groups
    and plots the results
    """

    # Define Comartment names
    compartments = ['Susceptible', 'Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead']
    # Define age Groups
    groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']

    days = days  # number of days to simulate
    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1
    num_groups = len(groups)
    num_compartments = len(compartments)

    # set contact frequency matrix
    data_dir = "data"
    baseline_contact_matrix0 = os.path.join(
        data_dir, "contacts/baseline_home.txt")
    baseline_contact_matrix1 = os.path.join(
        data_dir, "contacts/baseline_school_pf_eig.txt")
    baseline_contact_matrix2 = os.path.join(
        data_dir, "contacts/baseline_work.txt")
    baseline_contact_matrix3 = os.path.join(
        data_dir, "contacts/baseline_other.txt")

    # Initialize Parameters
    model = SecirModel(len(populations))

    # set parameters
    for i in range(num_groups):
        # Compartment transition duration

        model.parameters.IncubationTime[AgeGroup(i)] = 5.2  # R_2^(-1)+R_3^(-1)
        model.parameters.InfectiousTimeMild[AgeGroup(
            i)] = 6.  # 4-14  (=R4^(-1))
        # 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        model.parameters.SerialInterval[AgeGroup(i)] = 4.2
        model.parameters.HospitalizedToHomeTime[AgeGroup(
            i)] = 12.  # 7-16 (=R5^(-1))
        model.parameters.HomeToHospitalizedTime[AgeGroup(
            i)] = 5.  # 2.5-7 (=R6^(-1))
        model.parameters.HospitalizedToICUTime[AgeGroup(
            i)] = 2.  # 1-3.5 (=R7^(-1))
        model.parameters.ICUToHomeTime[AgeGroup(i)] = 8.  # 5-16 (=R8^(-1))
        model.parameters.ICUToDeathTime[AgeGroup(i)] = 5.  # 3.5-7 (=R5^(-1))

        # Initial number of peaople in each compartment
        model.populations[AgeGroup(i), State.Exposed] = random.uniform(
            0.00025, 0.0005) * populations[i]
        model.populations[AgeGroup(i), State.Carrier] = random.uniform(
            0.0001, 0.00035) * populations[i]
        model.populations[AgeGroup(i), State.Infected] = random.uniform(
            0.00007, 0.0001) * populations[i]
        model.populations[AgeGroup(
            i), State.Hospitalized] = random.uniform(
            0.00003, 0.00006) * populations[i]
        model.populations[AgeGroup(i), State.ICU] = random.uniform(
            0.00001, 0.00002) * populations[i]
        model.populations[AgeGroup(i), State.Recovered] = random.uniform(
            0.002, 0.008) * populations[i]
        model.populations[AgeGroup(i), State.Dead] = 0
        model.populations.set_difference_from_group_total_AgeGroup(
            (AgeGroup(i), State.Susceptible), populations[i])

        # Compartment transition propabilities

        model.parameters.RelativeCarrierInfectability[AgeGroup(i)] = 0.5
        model.parameters.InfectionProbabilityFromContact[AgeGroup(i)] = 0.7
        model.parameters.AsymptomaticCasesPerInfectious[AgeGroup(
            i)] = 0.1  # 0.01-0.16
        model.parameters.RiskOfInfectionFromSymptomatic[AgeGroup(
            i)] = 0.1  # 0.05-0.5
        model.parameters.HospitalizedCasesPerInfectious[AgeGroup(
            i)] = 0.05  # 0.1-0.35
        model.parameters.ICUCasesPerHospitalized[AgeGroup(
            i)] = 0.1  # 0.15-0.4
        model.parameters.DeathsPerICU[AgeGroup(i)] = 0.3  # 0.15-0.77
        # twice the value of RiskOfInfectionFromSymptomatic
        model.parameters.MaxRiskOfInfectionFromSymptomatic[AgeGroup(i)] = 0.5

    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

    # set contact rates and emulate some mitigations
    # set contact frequency matrix
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.loadtxt(baseline_contact_matrix0) \
        + np.loadtxt(baseline_contact_matrix1) + \
        np.loadtxt(baseline_contact_matrix2) + \
        np.loadtxt(baseline_contact_matrix3)
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
        (num_groups, num_groups)) * 0

    # add one random damping matrix and date
    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=np.array(damping_data[0]), t=damping_data[1], level=0, type=0))

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)
    # print(result.get_last_value())

    num_time_points = result.get_num_time_points()
    result_array = result.as_ndarray()
    t = result_array[0, :]
    group_data = np.transpose(result_array[1:, :])

    # sum over all groups
    data = np.zeros((num_time_points, num_compartments))
    for i in range(num_groups):
        data += group_data[:, i * num_compartments: (i + 1) * num_compartments]

    # Run Simulation
    data = np.zeros((days, len(compartments) * num_groups))
    dataset = []
    data[0] = model.populations.get_compartments()
    # dataset.append(model.populations.get_compartments())
    for day in range(1, days):
        result = simulate(0, day, dt, model)
        data[day] = result.get_last_value()[:]
        # dataset.append(val)

    for row in data:
        dataset.append(row)

    return dataset


def generate_data(num_runs, path, input_width, label_width, normalize_labels=False, save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often

   @param num_runs Number of times, the function run_secir_simulation is called.
   @param path Path, where the datasets are stored.
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """
    data = {
        "inputs": [],
        "labels": [],
        "dampings": []
    }

    days = input_width + label_width

    # create 100 different dampings
    # for _ in range(100):
    #     # create random damping matrix and date in a own list
    #     # and save in list dampings
    #     damping_single = []

    #     damping_single.append(damping_matrix.tolist())
    #     damping_single.append()
    #     data["dampings"].append(damping_single)

    # get population sizes
    population = get_population()

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)
    for _ in range(num_runs):
        damping_matrix = np.zeros((6, 6))  # no of groups = 6
        for idx in range(6):
            for i in range(idx+1):
                damp = round(random.random(), 4)
                damping_matrix[idx][i] = damp
                damping_matrix[i][idx] = damp
        data["dampings"].append(damping_matrix)
        damping_date = float(random.randint(10, days))
        data_run = run_secir_groups_simulation(
            days, [damping_matrix, damping_date],
            population[random.randint(0, len(population) - 1)])

        # # drop columns susceptible und recovered
        # compartments = ['Susceptible', 'Exposed', 'Carrier',
        #                 'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead']
        # df_data = pd.DataFrame(data=data_run, columns=compartments)
        # df_data.drop(['Susceptible', 'Recovered'], axis=1, inplace=True)

        # data_run = df_data.values  # update data_run

        inputs = data_run[:input_width]
        # add damping date to input
        for i in range(input_width):
            inputs[i] = np.append(inputs[i], damping_date)

        data["inputs"].append(inputs)
        data["labels"].append(data_run[input_width:])
        bar.next()

    bar.finish()

    if save_data:

        # cast dfs to tensors
        data["inputs"] = tf.stack(data["inputs"])
        data["labels"] = tf.stack(data["labels"])
        # data["dampings"] = tf.stack(data["dampings"])

        # normalize data
        # data["inputs"] = normalize(data["inputs"])
        # if normalize_labels:
        #     data["labels"] = normalize(data["labels"])

        # data = splitdata(data["inputs"], data["labels"])

        # check if data directory exists. If necessary create it.
        if not os.path.isdir(path):
            os.mkdir(path)

        # save dict to json file
        with open(os.path.join(path, 'data_secir_age_groups.pickle'), 'wb') as f:
            pickle.dump(data, f)


def get_population(path="data\pydata\Germany\county_population_dim401.json"):
    with open(path) as f:
        data = json.load(f)
    population = []
    for data_entry in data:
        population_county = []
        population_county.append(
            data_entry['<3 years'] + data_entry['3-5 years'] / 2)
        population_county.append(data_entry['6-14 years'])
        population_county.append(
            data_entry['15-17 years'] + data_entry['18-24 years'] + data_entry['25-29 years'] + data_entry['30-39 years']/2)
        population_county.append(
            data_entry['30-39 years']/2 + data_entry['40-49 years'] + data_entry['50-64 years'] * 2/3)
        population_county.append(
            data_entry['65-74 years'] + data_entry['>74 years'] * 0.2 + data_entry['50-64 years'] * 1/3)
        population_county.append(
            data_entry['>74 years'] * 0.8)

        population.append(population_county)
    return population

    return 0


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


def splitdampings(damping, split_train=0.7,
                  split_valid=0.2, split_test=0.1):
    """! Split dampings in train, valid and test

   @param damping damping matrices
   @param labels label dataset
   @param split_train ratio of train datasets
   @param split_valid ratio of validation datasets
   @param split_test ratio of test datasets
   """

    if split_train + split_valid + split_test != 1:
        ValueError("summed Split ratios not equal 1! Please adjust the values")

    n = damping.shape[0]
    n_train = int(n * split_train)
    n_valid = int(n * split_valid)
    n_test = int(n * split_test)

    if n_train + n_valid + n_test != n:
        n_test = n - n_train - n_valid

    damping_train, damping_valid, damping_test = tf.split(
        damping, [n_train, n_valid, n_test], 0)
    data = {
        "train_damping": damping_train,
        "valid_damping": damping_valid,
        "test_damping": damping_test
    }

    return data


def normalize(data):
    """! # normalize input data with:
                               (value     − min_value)
                    data =  _________________________
                               (max_value − min_value)

    alternative is tf.linalg.normalize(
        data, ord='euclidean', axis=None, name=None)

    @param data  dataset

    """
    # return tf.divide(
    #         tf.subtract(
    #             data,
    #             tf.reduce_min(data)
    #         ),
    #         tf.subtract(
    #             tf.reduce_max(data),
    #             tf.reduce_min(data)
    #         )
    #     )
    return tf.linalg.normalize(data, ord='euclidean', axis=None, name=None)[0]


if __name__ == "__main__":
    # TODO: Save contact matrix depending on the damping.
    # In the actual state it might be enough to save the regular one and the damping

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')

    input_width = 5
    label_width = 30
    num_runs = 1000
    generate_data(num_runs, path_data, input_width,
                  label_width, normalize_labels=True)
