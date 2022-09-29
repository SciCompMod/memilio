from cProfile import label
import json
import pickle
from memilio.simulation import UncertainContactMatrix, ContactMatrix, Damping, ContactMatrixGroup
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
from sklearn.preprocessing import FunctionTransformer
import seaborn as sns  # plot after normalization


def run_secir_groups_simulation(days, t1, t2, t3, populations):
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

    baseline = getBaselineMatrix()
    minimum = getMinimumMatrix()

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = minimum

    damping1 = np.ones((num_groups, num_groups)
                       ) * np.float16(random.uniform(0, 0.6))

    damping2 = np.ones((num_groups, num_groups)
                       ) * np.float16(random.uniform(0, 0.6))

    damping3 = np.ones((num_groups, num_groups)
                       ) * np.float16(random.uniform(0, 0.6))

    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=(damping1), t=t1, level=0, type=0))

    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=(damping2), t=t2, level=0, type=0))

    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=(damping3), t=t3, level=0, type=0))

    correct_damping_matrix1 = model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
        t1+1)
    correct_damping_matrix2 = model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
        t2+1)
    correct_damping_matrix3 = model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
        t3+1)

    # Apply mathematical constraints to parameters
    model.apply_constraints()

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

    return dataset, correct_damping_matrix1, correct_damping_matrix2, correct_damping_matrix3


def getBaselineMatrix():

    baseline_contact_matrix0 = os.path.join(
        "./data/contacts/baseline_home.txt")
    baseline_contact_matrix1 = os.path.join(
        "./data/contacts/baseline_school_pf_eig.txt")
    baseline_contact_matrix2 = os.path.join(
        "./data/contacts/baseline_work.txt")
    baseline_contact_matrix3 = os.path.join(
        "./data/contacts/baseline_other.txt")

    baseline = np.loadtxt(baseline_contact_matrix0) \
        + np.loadtxt(baseline_contact_matrix1) + \
        np.loadtxt(baseline_contact_matrix2) + \
        np.loadtxt(baseline_contact_matrix3)

    return baseline


def getMinimumMatrix():

    minimum_contact_matrix0 = os.path.join(
        "./data/contacts/minimum_home.txt")
    minimum_contact_matrix1 = os.path.join(
        "./data/contacts/minimum_school_pf_eig.txt")
    minimum_contact_matrix2 = os.path.join(
        "./data/contacts/minimum_work.txt")
    minimum_contact_matrix3 = os.path.join(
        "./data/contacts/minimum_other.txt")

    minimum = np.loadtxt(minimum_contact_matrix0) \
        + np.loadtxt(minimum_contact_matrix1) + \
        np.loadtxt(minimum_contact_matrix2) + \
        np.loadtxt(minimum_contact_matrix3)

    return minimum


def generate_data(num_runs, path,
                  normalize_labels=False, save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often

   @param num_runs Number of times, the function run_secir_simulation is called.
   @param path Path, where the datasets are stored.
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """
    data = {
        "inputs": [],
        "labels": [],
        "contact_matrix": []
    }

    t1 = 5    # maybe it would be better to put them as a paramenter of the generate data function
    t2 = 20
    t3 = 40
    days = 60
    #days = t3+30
    input_width = 5
    label_width = days-input_width

    # get population sizes
    population = get_population()

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)

    for _ in range(num_runs):

        data_run, damping_matrix1, damping_matrix2, damping_matrix3 = run_secir_groups_simulation(
            days, t1, t2, t3,  population[random.randint(0, len(population) - 1)])

        data["contact_matrix"].append(
            [damping_matrix1, damping_matrix2, damping_matrix3])

        inputs = data_run[:input_width]
        data["inputs"].append(inputs)
        # data["labels"].append(data_run[input_width:])
        data['labels'].append(data_run[input_width:days])
        bar.next()

    bar.finish()

    if save_data:

        transformer = FunctionTransformer(np.log1p, validate=True)
        inputs = np.asarray(data['inputs']).transpose(2, 0, 1).reshape(48, -1)
        scaled_inputs = transformer.transform(inputs)
        scaled_inputs = scaled_inputs.transpose().reshape(num_runs, input_width, 48)
        scaled_inputs_list = scaled_inputs.tolist()

        labels = np.asarray(data['labels']).transpose(2, 0, 1).reshape(48, -1)
        scaled_labels = transformer.transform(labels)
        scaled_labels = scaled_labels.transpose().reshape(num_runs, label_width, 48)
        scaled_labels_list = scaled_labels.tolist()

        # cast dfs to tensors
        data['inputs'] = tf.stack(scaled_inputs_list)
        data['labels'] = tf.stack(scaled_labels_list)

        # check if data directory exists. If necessary create it.
        if not os.path.isdir(path):
            os.mkdir(path)

        # save dict to json file
        with open(os.path.join(path, 'data_secir_age_groups.pickle'), 'wb') as f:
            pickle.dump(data, f)


def get_population(path="/home/schm_a45/Documents/Code/memilio/memilio/data/pydata/Germany/county_population_dim401.json"):

    with open(path) as f:
        data = json.load(f)
    population = []
    for data_entry in data:
        population_county = []
        population_county.append(
            data_entry['<3 years'] + data_entry['3-5 years'] / 2)
        population_county.append(data_entry['6-14 years'])
        population_county.append(
            data_entry['15-17 years'] + data_entry['18-24 years'] +
            data_entry['25-29 years'] + data_entry['30-39 years'] / 2)
        population_county.append(
            data_entry['30-39 years'] / 2 + data_entry['40-49 years'] +
            data_entry['50-64 years'] * 2 / 3)
        population_county.append(
            data_entry['65-74 years'] + data_entry['>74 years'] * 0.2 +
            data_entry['50-64 years'] * 1 / 3)
        population_county.append(
            data_entry['>74 years'] * 0.8)

        population.append(population_county)
    return population


if __name__ == "__main__":

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data_simulation')

    num_runs = 200
    generate_data(num_runs, path_data, normalize_labels=True)
