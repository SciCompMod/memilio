import json
import pickle
from memilio.simulation import Damping
from memilio.simulation.secir import SecirModel, simulate, AgeGroup
from memilio.simulation.secir import InfectionState as State
import numpy as np
from datetime import date
import random
import os
from progress.bar import Bar  # pip install progess
import tensorflow as tf
from datetime import date
from sklearn.preprocessing import FunctionTransformer


def run_secir_groups_simulation(days, t1, populations):
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

    model = SecirModel(len(populations))
    # Initialize Parameters

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
                       ) * np.float16(random.uniform(0, 0.5))

    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=(damping1), t=t1, level=0, type=0))

    correct_damping_matrix1 = model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
        t1+1)

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

    return dataset, correct_damping_matrix1


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


def generate_data(num_runs, path, input_width, label_width,
                  save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often

   @param num_runs Number of times, the function run_secir_simulation is called.
   @param path Path, where the datasets are stored.
   @param input_width number of time steps used for model input.
   @param label_wodth number of time steps (days) used as model output/label.  
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """
    data = {
        "inputs": [],
        "labels": [],
        "contact_matrix": []
    }

    t1 = 5
    days = t1+label_width

    population = get_population()

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)

    for _ in range(num_runs):

        data_run, damping_matrix1 = run_secir_groups_simulation(
            days, t1, population[random.randint(0, len(population) - 1)])

        data["contact_matrix"].append(
            [damping_matrix1])

        inputs = data_run[:input_width]
        data["inputs"].append(inputs)

        data["labels"].append(data_run[input_width:])

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


# def get_population(path="data\pydata\Germany\county_population_dim401.json"):
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


if __name__ == "__main__":

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data_groups')

    input_width = 5
    label_width = 30
    num_runs = 50000
    generate_data(num_runs, path_data, input_width,
                  label_width)
