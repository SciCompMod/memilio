import time 
import os
import numpy as np
from datetime import date
import copy
import pickle
import random 
import pandas as pd 
from progress.bar import Bar
import tensorflow as tf
import matplotlib.pyplot as plt
import pandas as pd
import json 

import spektral
import scipy.sparse as sp
import tensorflow as tf

from keras.layers import Dense
from keras.metrics import mean_absolute_percentage_error
from keras.models import Model

from spektral.data import  MixedLoader
from spektral.layers import ARMAConv
from spektral.transforms.normalize_adj import NormalizeAdj
from spektral.utils.convolution import gcn_filter, normalized_laplacian, rescale_laplacian, normalized_adjacency

from memilio.simulation.secir import InfectionState
from memilio.simulation import (ContactMatrix, Damping, LogLevel,
                                UncertainContactMatrix, set_log_level, AgeGroup)
from memilio.simulation.secir import ( Index_InfectionState,
                                      InfectionState, Model, Simulation,
                                      interpolate_simulation_result, simulate)


iterations = 1
runs = 1 

########## GNN no damp days ######################
# set parameters for best model
layer_name='ARMAConv'
number_of_layers = 3
number_of_channels = 128
layer = ARMAConv

# open commuter data
numer_of_nodes = 400 
# path = os.path.dirname(os.path.realpath(__file__))
# path_data = os.path.join(path,
#                          'data')


# commuter_file = open(os.path.join(
#     path_data, 'commuter_migration_scaled.txt'), 'rb')
# commuter_data = pd.read_csv(commuter_file, sep=" ", header=None)
commuter_data = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/data/commuter_migration_scaled.txt', sep=" ", header=None)
sub_matrix = commuter_data.iloc[:numer_of_nodes, 0:numer_of_nodes]

adjacency_matrix = np.asarray(sub_matrix) # non binary adjacency matrix 



weight_names = ['best_weights_ARMAConv_30days_nodamp_w1.pickle', 'best_weights_ARMAConv_60days_nodamp_w.pickle', 
               'best_weights_ARMAConv_120days_nodamp_w.pickle']
filenames_input = ['inputs_30days_w.pickle', 'inputs_60days_w.pickle', 
             'inputs_120days_w.pickle']

filenames_labels = ['labels_30days_w.pickle', 'labels_60days_w.pickle',
             'labels_120days_w.pickle']
#days = [30,60,90,120]
days = [30,60,120]

df_GNN_days = pd.DataFrame(columns = ['days', 'time_per_run'])
for weight_name, filename_i, filename_l in zip(weight_names, filenames_input, filenames_labels):

    # Define the model architecture
    class Net(tf.keras.Model):
        def __init__(self):
            super().__init__()

            self.conv1 = layer(number_of_channels, activation='elu')
            self.conv2 = layer(number_of_channels, activation='elu')
            self.conv3 = layer(number_of_channels, activation='elu')
            self.dense = tf.keras.layers.Dense(data.n_labels, activation="linear")

        def call(self, inputs):
            x, a = inputs
            a = np.asarray(a)

            x = self.conv1([x, a])
            x = self.conv2([x, a])
            x = self.conv2([x, a])

            output = self.dense(x)

            return output

    # Load weights
    weight_path = os.path.join('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_model_weights', weight_name)
    with open(weight_path, "rb") as fp:
        b = pickle.load(fp)
        best_weights = b    

    # Load test data
    file_i = open(os.path.join('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_inputs', filename_i), 'rb')
    test_inputs = pickle.load(file_i)
    file_i.close()

    file_l = open(os.path.join('/home/schm_a45/Documents/Code/memilio/memilio/pycode/machine-learning/model/GNN/saved_labels', filename_l), 'rb')
    test_labels = pickle.load(file_l)
    file_l.close()

    node_features = test_inputs[0]

    class MyDataset(spektral.data.dataset.Dataset):
        def read(self):
            self.a = spektral.utils.normalized_laplacian(spektral.utils.rescale_laplacian(spektral.utils.normalized_laplacian(adjacency_matrix)))
            return [spektral.data.Graph(x=x, y=y) for x, y in zip(node_features, test_labels)]

    data = MyDataset(transforms=spektral.transforms.NormalizeAdj())
    batch_size = node_features.shape[0]

    model = Net()

    # Build the model before setting weights
    #model.build((batch_size, test_inputs))  # Replace ... with the input shape
    model(test_inputs, training=False)

    model.set_weights(best_weights)

    times = []
    for i in range(iterations):
        for node_features, a in zip(test_inputs[0][:runs], test_inputs[1]):
            #reshaped_input = input_data[None, :]
            start = time.time()
            pred = model(tuple([test_inputs[0], test_inputs[1]]), training=False)
            end = time.time()
            times.append(end - start)


    print('GNN : Mean time per rum for ' + str(days) + 'days: ' + str(np.asarray(times).mean()))
    df_GNN_days.loc[len(df_GNN_days)]=[days, np.asarray(times).mean()]




####### GNN with ODE ##############
from memilio.simulation import (ContactMatrix, Damping, LogLevel,
                                UncertainContactMatrix, set_log_level, AgeGroup)
from memilio.simulation.secir import ( Index_InfectionState,
                                      InfectionState, Model, Simulation,
                                      interpolate_simulation_result, simulate)

def get_population(path="/home/schm_a45/Documents/Code/memilio/memilio/data/pydata/Germany/county_population.json"):
    
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

def getBaselineMatrix():
    """! loads the baselinematrix
    """

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

def remove_confirmed_compartments(dataset_entries, num_groups):
    new_dataset_entries = []
    for i in dataset_entries : 
      dataset_entries_reshaped  = i.reshape([num_groups, int(np.asarray(dataset_entries).shape[1]/num_groups) ])
      sum_inf_no_symp = np.sum(dataset_entries_reshaped [:, [2, 3]], axis=1)
      sum_inf_symp = np.sum(dataset_entries_reshaped [:, [4, 5]], axis=1)
      dataset_entries_reshaped[:, 2] = sum_inf_no_symp
      dataset_entries_reshaped[:, 4] = sum_inf_symp
      new_dataset_entries.append(np.delete(dataset_entries_reshaped , [3, 5], axis=1).flatten())
    return new_dataset_entries


def run_secir_groups_simulation(days, populations):

    """
    Runs the c++ secir model using mulitple age groups
    and plots the results
    @param days number of days to be simulated. 
    @param t1 damping day. 
    @param populations array of population data. 
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

    model = Model(num_groups)
    # Initialize Parameters

    # set parameters
    for i in range(num_groups):

        # Compartment transition duration
        model.parameters.IncubationTime[AgeGroup(i)] = 5.2
        model.parameters.TimeInfectedSymptoms[AgeGroup(i)] = 6.
        model.parameters.SerialInterval[AgeGroup(i)] = 4.2
        model.parameters.TimeInfectedSevere[AgeGroup(i)] = 12.
        model.parameters.TimeInfectedCritical[AgeGroup(i)] = 8.

        # Initial number of people in each compartment with random numbers
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Exposed)] = random.uniform(
            0.00025, 0.005) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedNoSymptoms)] = random.uniform(
            0.0001, 0.0035) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSymptoms)] = random.uniform(
            0.00007, 0.001) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSevere)] = random.uniform(
            0.00003, 0.0006) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedCritical)] = random.uniform(
            0.00001, 0.00002) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Recovered)] = random.uniform(
            0.002, 0.08) * populations[i]
        model.populations[AgeGroup(i),
                          Index_InfectionState(InfectionState.Dead)] = random.uniform(
            0, 0.0003) * populations[i]
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

    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

    # set contact rates and emulate some mitigations
    # set contact frequency matrix
    # Load baseline and minimum contact matrix and assign them to the model
    path = '/unsecured_KP/schm_a45/memilio/memilio/data/contacts/baseline_matrix'
    with open(path, "rb") as fp:
       baseline = pickle.load(fp)

    #baseline = getBaselineMatrix()
    #minimum = getMinimumMatrix()

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
        (num_groups, num_groups)) * 0
    #model.parameters.ContactPatterns.cont_freq_mat[0].minimum = minimum

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days-1, dt, model)
    # Interpolate simulation result on days time scale
    result = interpolate_simulation_result(result)

    # Using an array instead of a list to avoid problems with references
    result_array = result.as_ndarray()
    dataset = []
    # Omit first column, as the time points are not of interest here.
    dataset_entries = copy.deepcopy(result_array[1:, :].transpose())

    dataset_entires_without_confirmed = remove_confirmed_compartments(dataset_entries, num_groups)
    return dataset_entires_without_confirmed

label_widths = [25,55,85,115]
df_ODE_days = pd.DataFrame(columns = ['days', 'time_per_run'])
for lw in label_widths:
        input_width = 5 
        label_width = lw
        num_runs = runs 
        data = {
                "inputs": [],
                "labels": []
        }
        days = input_width + label_width - 1

        path_population = os.path.abspath(
                r"data//pydata//Germany//county_population.json")

        # Load population data
        population = get_population(path_population)

        times = []
        for i in range(iterations):
                

                # show progess in terminal for longer runs
                # Due to the random structure, theres currently no need to shuffle the data

                for _ in range(0, num_runs):
                        pop = population[random.randint(0, len(population) - 1)]
                        start = time.time()
                        data_run = run_secir_groups_simulation(
                                        days, pop)

                        end = time.time()
                        times.append(end - start)
        #print(np.asarray(times).mean())
        print('GNN Days : Mean time per rum for ' + str(input_width + label_width) + 'days: ' + str(np.asarray(times).mean()))
        df_ODE_days.loc[len(df_ODE_days)]=[days+1, np.asarray(times).mean()]
print(df_ODE_days)