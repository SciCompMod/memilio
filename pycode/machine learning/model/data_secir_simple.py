from cProfile import label
from memilio.simulation import UncertainContactMatrix, ContactMatrix, Damping
from memilio.simulation.secir import SecirModel, simulate, AgeGroup, Index_InfectionState, SecirSimulation
from memilio.simulation.secir import InfectionState as State
import numpy as np
import pandas as pd
import pickle
from datetime import date
from math import ceil
import random
import os
from progress.bar import Bar # pip install progess
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf
import matplotlib.pyplot as plt
import memilio.simulation as mio

import seaborn as sns # plot after normalization

def run_secir_simulation(days):
    """! Slightly modified example of Secir Simple.

   Here, the initial values are choosen randomly, so the model is no longer deterministic.
   Necessary to create the dataset for later training process.
   This method is called within a loop in the function 'generate_data', which also sets the number of runs.

   """
    mio.set_log_level(mio.LogLevel.Off)

    # Define Comartment names
    compartments = ['Susceptible', 'Exposed', 'Carrier',
                    'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead']

    # Define random population number in intervall (0,500000)
    populations = [500000 * random.random()]

    # since we defined the number of days, we already know the dimension of our dataset.
    # We aim to save the compartment population for each day.
    # i.e. for days = 35 and len(compartments = 8) -> 35 x  8 array
    data = np.zeros((days,len(compartments)))
    dataset = []

    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1
    num_groups = 1

    # Initialize Parameters
    model = SecirModel(1)

    A0 = AgeGroup(0)


    # Define random values for damping coef and start date
    # damping_factor = 1
    # damping_date = ceil(random.random() * days) # scale with days since random is in Interval (0,1)

    # save damping factor in vector. -1 for index correction
    # data[damping_date-1:,-1] = damping_factor

    # Set parameters

    # Compartment transition duration
    model.parameters.IncubationTime[A0] = 5.2  # R_2^(-1)+R_3^(-1)
    model.parameters.InfectiousTimeMild[A0] = 6.  # 4-14  (=R4^(-1))
    # 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
    model.parameters.SerialInterval[A0] = 4.2
    model.parameters.HospitalizedToHomeTime[A0] = 12.  # 7-16 (=R5^(-1))
    model.parameters.HomeToHospitalizedTime[A0] = 5.  # 2.5-7 (=R6^(-1))
    model.parameters.HospitalizedToICUTime[A0] = 2.  # 1-3.5 (=R7^(-1))
    model.parameters.ICUToHomeTime[A0] = 8.  # 5-16 (=R8^(-1))
    model.parameters.ICUToDeathTime[A0] = 5.  # 3.5-7 (=R5^(-1))

    # Initial number of people in each compartment with random numbers 
    model.populations[A0, Index_InfectionState(State.Exposed)] = 250 * random.random()
    model.populations[A0, Index_InfectionState(State.Carrier)] = 120 * random.random()
    model.populations[A0, Index_InfectionState(State.Infected)] = 120 * random.random()
    model.populations[A0, Index_InfectionState(State.Hospitalized)] = 50 * random.random()
    model.populations[A0, Index_InfectionState(State.ICU)] = 30 * random.random()
    model.populations[A0, Index_InfectionState(State.Recovered)] = 30 * random.random()
    model.populations[A0, Index_InfectionState(State.Dead)] = 0
    model.populations.set_difference_from_total(
        (A0, Index_InfectionState(State.Susceptible)), populations[0])

    # Compartment transition propabilities
    model.parameters.RelativeCarrierInfectability[A0] = 0.67
    model.parameters.InfectionProbabilityFromContact[A0] = 1.0
    model.parameters.AsymptomaticCasesPerInfectious[A0] = 0.09  # 0.01-0.16
    model.parameters.RiskOfInfectionFromSymptomatic[A0] = 0.25  # 0.05-0.5
    model.parameters.HospitalizedCasesPerInfectious[A0] = 0.2  # 0.1-0.35
    model.parameters.ICUCasesPerHospitalized[A0] = 0.25  # 0.15-0.4
    model.parameters.DeathsPerICU[A0] = 0.3  # 0.15-0.77
    # twice the value of RiskOfInfectionFromSymptomatic
    model.parameters.MaxRiskOfInfectionFromSymptomatic[A0] = 0.5

    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

    # model.parameters.ContactPatterns.cont_freq_mat[0] = ContactMatrix(np.r_[0.5])
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (num_groups, num_groups)) * 1
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
        (num_groups, num_groups)) * 0
    # model.parameters.ContactPatterns.cont_freq_mat.add_damping(
    #     Damping(coeffs=np.r_[damping_factor], t=damping_date, level=0, type=0))

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    data[0] = model.populations.get_compartments()
    # dataset.append(model.populations.get_compartments())
    for day in range(1, days):
        result = simulate(0, day, dt, model)
        data[day] = result.get_last_value()[:]
        # dataset.append(val)

    # when directly using a list instead of an array, we get problems wutg references
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
        "inputs" : [],
        "labels" : []
    }

    days = input_width + label_width

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)
    for _ in range(0, num_runs):
        data_run = run_secir_simulation(days)
        data["inputs"].append(data_run[:input_width])
        data["labels"].append(data_run[input_width:])
        bar.next()
    
    bar.finish()

    if save_data:

        # cast dfs to tensors
        data["inputs"] = tf.stack(data["inputs"])
        data["labels"] = tf.stack(data["labels"])

        # normalize data 
        data["inputs"] = normalize(data["inputs"])
        if normalize_labels:
            data["labels"] = normalize(data["labels"])

        # data = splitdata(data["inputs"], data["labels"])

        # check if data directory exists. If necessary create it.
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

    inputs_train, inputs_valid, inputs_test = tf.split(inputs, [n_train, n_valid, n_test], 0)
    labels_train, labels_valid, labels_test = tf.split(labels, [n_train, n_valid, n_test], 0)

    data = {
        "train_inputs" : inputs_train,
        "train_labels" :labels_train,
        "valid_inputs" : inputs_valid,
        "valid_labels" :labels_valid,
        "test_inputs" : inputs_test,
        "test_labels" :labels_test
    }

    return data

def normalize(data):
    """! # normalize input data with:
                               (value     − min_value)
                    data =  _________________________
                               (max_value − min_value)
                    
    alternative is tf.linalg.normalize(data, ord='euclidean', axis=None, name=None)

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
    path_data = os.path.join(os.path.dirname(os.path.realpath(path)), 'data')

    input_width = 5
    label_width = 30
    num_runs = 500
    generate_data(num_runs, path_data, input_width, label_width, normalize_labels=True)
    
    
