
import json
import pickle
from memilio.simulation import Damping
#from memilio.simulation.secir import Model, simulate, AgeGroup
from memilio.simulation.secir import InfectionState as State
import numpy as np
from datetime import date
import random
import copy
import os
from progress.bar import Bar  # pip install progess
import tensorflow as tf
from datetime import date
from sklearn.preprocessing import FunctionTransformer
from random import randrange

from memilio.simulation import (ContactMatrix, Damping, LogLevel,
                                UncertainContactMatrix, set_log_level, AgeGroup)
from memilio.simulation.secir import ( Index_InfectionState,
                                      InfectionState, Model, Simulation,
                                      interpolate_simulation_result, simulate)




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

def run_secir_groups_simulation(days, populations, damping_days):
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

    model = Model(len(populations))
    # Initialize Parameters

    # set parameters
    for i in range(num_groups):
        # Compartment transition duration

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
            0.00001, 0.0002) * populations[i]
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

    baseline = getBaselineMatrix()
    #minimum = getMinimumMatrix()

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline
    #model.parameters.ContactPatterns.cont_freq_mat[0].minimum = minimum
    #set minum matrix to 0 
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
        (num_groups, num_groups)) * 0

    # Generate a damping matrix and assign it to the model
    damped_matrices = []
    damping_coeff = []
    for day in damping_days:
        damping = np.ones((num_groups, num_groups)
                      ) * np.float16(random.uniform(0, 0.5))

        model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
            coeffs=(damping), t=day, level=0, type=0))

        damped_matrices.append(model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
            day+1))
        damping_coeff.append(damping[0][0])

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)

    # Interpolate simulation result on days time scale
    result = interpolate_simulation_result(result)

    # Using an array instead of a list to avoid problems with references
    result_array = result.as_ndarray()

    # Omit first column, as the time points are not of interest here.
    dataset_entries = copy.deepcopy(result_array[1:, :].transpose())
    dataset_entires_without_confirmed = remove_confirmed_compartments(dataset_entries, num_groups)
    return dataset_entires_without_confirmed, damped_matrices, damping_days, damping_coeff





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


def getMinimumMatrix():
    """! loads the minimum matrix
    """

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



def generate_data(
        num_runs, path, path_population,  input_width, label_width,number_of_populations, number_of_dampings,  file_name, normalize = True, save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often
   @param num_runs Number of times, the function run_secir_simulation is called.
   @param path Path, where the datasets are stored.
   @param input_width number of time steps used for model input.
   @param label_width number of time steps (days) used as model output/label.  
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """
    all_data = {"inputs": [],
                "labels": [],
                "damping_coeff": [],
                "damping_day": [], 
                "damped_matrix": []}

    num_groups = 6
        # The number of days is the same as the sum of input and label width.
    # Since the first day of the input is day 0, we still need to subtract 1.
    days = input_width + label_width - 1



    damping_days = generate_dampings_withshadowdamp(number_of_dampings = number_of_dampings, days= days, min_distance=7, min_damping_day=input_width, n_runs = num_runs)


    for i in range(number_of_populations):
        # days = damping_days[-1]+20
        data = {"inputs": [],
                "labels": [],
                "damping_coeff": [],
                "damping_day": [],
                "damped_matrix": []}

                # Load population data
        population = get_population(path_population)

        # show progess in terminal for longer runs
        # Due to the random structure, theres currently no need to shuffle the data
        bar = Bar('Number of Runs done', max=num_runs)

        for j in range(num_runs):
            

            data_run, damped_contact_matrix, damping_days_s, damping_factor = run_secir_groups_simulation(
                 days,  population[i],  damping_days[j])

    

            # for i in damping_matrix:
            #    data["contact_matrix"].append(i)

            inputs = data_run[:input_width]
            data["inputs"].append(inputs)

            data["labels"].append(data_run[input_width:])
            data['damping_coeff'].append(damping_factor)
            data['damping_day'].append(damping_days_s)
            data['damped_matrix'].append(damped_contact_matrix)
            bar.next()

        bar.finish()
        print('Number of populations done:', i)

        if save_data:

            transformer = FunctionTransformer(np.log1p, validate=True)
            inputs = np.asarray(
                data['inputs']).transpose(
                2, 0, 1).reshape(
                48, -1)
            scaled_inputs = transformer.transform(inputs)
            scaled_inputs = scaled_inputs.transpose().reshape(num_runs, input_width, 48)
            scaled_inputs_list = scaled_inputs.tolist()

            labels = np.asarray(
                data['labels']).transpose(
                2, 0, 1).reshape(
                48, -1)
            scaled_labels = transformer.transform(labels)
            scaled_labels = scaled_labels.transpose().reshape(
                num_runs, np.asarray(data['labels']).shape[1], 48)
            scaled_labels_list = scaled_labels.tolist()

            data['inputs'] = tf.stack(scaled_inputs_list)
            data['labels'] = tf.stack(scaled_labels_list)

            all_data['inputs'].append(data['inputs'])
            all_data['labels'].append(data['labels'])
            all_data['damping_coeff'].append(data['damping_coeff'])
            all_data['damping_day'].append(data['damping_day'])
            all_data['damped_matrix'].append(data['damped_matrix'])
            # check if data directory exists. If necessary create it.
            if not os.path.isdir(path):
                os.mkdir(path)

            # save dict to json file
            with open(os.path.join(path, file_name), 'wb') as f:
                pickle.dump(all_data, f)
    return all_data



def generate_dampings_withshadowdamp(number_of_dampings, days, min_distance, min_damping_day, n_runs):
    
    number_of_dampings = number_of_dampings
    days = days
    min_distance = min_distance 
    min_damping_day = min_damping_day
    number_of_runs= n_runs

    all_dampings = []
    count_runs = 0 
    count_shadow = 0
    while len(all_dampings)<number_of_runs:

        days_list = list(range((min_damping_day),days))
        dampings = []
        if count_shadow <2:
            for i in range(number_of_dampings):
                damp = random.choice(days_list)
                days_before = list(range(damp-(min_distance), damp))
                days_after = list(range(damp, damp+(min_distance+1)))
                dampings.append(damp)
                days_list = [ele for ele in days_list if ele not in (days_before+days_after)]
        else: 
            # chose a forbidden damping 
            damp = random.choice(list(range((0-min_distance),0))+ list(range(days+1, (days+min_distance+1))))
                
            days_before = list(range(damp-(min_distance), damp))
            days_after = list(range(damp, damp+(min_distance+1)))
            days_list = [ele for ele in days_list if ele not in (days_before+days_after)]
            dampings.append(damp)
            for i in range(number_of_dampings):
                
                
                damp = random.choice(days_list)
                days_before = list(range(damp-(min_distance), damp))
                days_after = list(range(damp, damp+(min_distance+1)))
                dampings.append(damp)
                days_list = [ele for ele in days_list if ele not in (days_before+days_after)]
                count_shadow = 0
        
            
        forbidden_damping_values = list(range((0-min_distance),0))+ list(range(days+1, (days+min_distance+1)))
        dampings = [ele for ele in dampings if ele not in forbidden_damping_values]
        count_runs+=1
        count_shadow +=1
        # select first or last five dampings
        if len(dampings) >= number_of_dampings:
            #dampings = random.sample(dampings, 5)
            all_dampings.append(sorted(dampings))
        #     if count_runs % 2 == 0:
        
    return np.asarray(all_dampings)



def get_population(path):
    """! read population data in list from dataset
    @param path Path to the dataset containing the population data
    """

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


def flat_input(input):
    dim = tf.reduce_prod(tf.shape(input)[1:])
    return tf.reshape(input, [-1, dim])


def split_contact_matrices(contact_matrices, split_train=0.7,
                           split_valid=0.2, split_test=0.1):
    """! Split dampings in train, valid and test
   @param contact_matrices contact matrices
   @param labels label dataset
   @param split_train ratio of train datasets
   @param split_valid ratio of validation datasets
   @param split_test ratio of test datasets
   """

    if split_train + split_valid + split_test != 1:
        ValueError("summed Split ratios not equal 1! Please adjust the values")

    n = contact_matrices.shape[0]
    n_train = int(n * split_train)
    n_valid = int(n * split_valid)
    n_test = int(n * split_test)

    if n_train + n_valid + n_test != n:
        n_test = n - n_train - n_valid

    contact_matrices_train, contact_matrices_valid, contact_matrices_test = tf.split(
        contact_matrices, [n_train, n_valid, n_test], 0)
    data = {
        "train": contact_matrices_train,
        "valid": contact_matrices_valid,
        "test": contact_matrices_test
    }

    return data


print('x')
if __name__ == "__main__":

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'data_GNN_400pop_4var_damp_100days_1k_w')
    
    file_name = 'GNN_400pop_damp_w.pickle'
    
    #path_population = os.path.abspath(
    #    r"data//pydata//Germany//county_population.json")
    path_population = '/home/schm_a45/Documents/Code/memilio/memilio/data/pydata/Germany/county_population.json'

    input_width = 5
    label_width = 100
    num_runs = 1000
    number_of_populations = 400
    number_of_dampings = 4
    generate_data(num_runs, path_data, path_population, input_width,
                         label_width, number_of_populations, number_of_dampings, file_name, normalize = True, save_data=True)
    
     
    