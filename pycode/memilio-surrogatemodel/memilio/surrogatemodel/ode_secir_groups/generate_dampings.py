import os
import pickle

import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import tensorflow as tf

import sklearn as sk
from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import train_test_split

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

from memilio.surrogatemodel.ode_secir_groups.data_generation import (
    remove_confirmed_compartments, get_population, getBaselineMatrix)
from memilio.surrogatemodel.ode_secir_groups.data_generation_test import transform_data


def calc_damps_var1(n_days, gamma_pos=0, alpha=-1, p0=0.5, t1_max=-0.3, t1_min=-1.2, t2_max=2, t2_min=-0.5):
    '''
    :param n_days: Number of days per time series
    :param gamma_pos: upper bound for h value
    :param alpha: upper bound for h value, where active damping should be stopped
    :param p0: probability for a change of the damping factor
    :param t1_max: upper end point for size of active damping changes
    :param t1_min: lower end point for size of active damping changes
    :param t2_max: upper end point for size of active damping changes
    :param t2_min: lower end point for size of base damping changes
    :param t3_max: upper end
    '''
    h = 0
    k = 0
    dampings = []
    while k < n_days:
        if h > gamma_pos:
            while h > alpha and k < n_days:
                delta_h = np.random.uniform(t1_min, t1_max)
                h = h + delta_h
                dampings.append(1-np.exp(h))
                k = k+1

        else:
            if np.random.binomial(1, p0):
                delta_h = np.random.uniform(t2_min, t2_max)

            else:
                delta_h = 0
            h = h+delta_h
            dampings.append(1 - np.exp(h))
            k = k+1

    return dampings


def generate_dampings(N, ndays, ndampings, parameters):
    """
    Generating N samples of dampings for a timeseries of length ndays.

    :param N: Number of samples to produce
    :param ndays: Number of days per sample
    :param ndampings: maximal number of dampings 
    :param parameters: Tuple of parameters used to generate the damping, this should be of the form
                        (gamma_pos, alpha, p0, t1_max, t1_min, t2_max, t2_min)
    :returns: dictionary of lists, "damping_days" is a list of length N containing list of days of length ndays,
    "damping_coeff" is a list of length N, where entry is a list of length ndays containing the damping factors.

    """
    # Unpacking the parameters
    gamma_pos, alpha, p0, t1_max, t1_min, t2_max, t2_min = parameters

    # Initializing the returns
    damping_days = []
    damping_coeff = []
    k = ndays//ndampings

    days = [j*k for j in np.arange(ndampings)]

    # Generating the damping coefficients for each sample
    for i in np.arange(N):
        damping_days.append(days)
        damping_coeff.append(calc_damps_var1(
            n_days=ndampings, gamma_pos=gamma_pos, alpha=alpha, p0=p0, t1_max=t1_max, t1_min=t1_min,
            t2_max=t2_max, t2_min=t2_min
        ))

    return {
        "damping_days": damping_days,
        "damping_coeff": damping_coeff
    }


def run_secir_simulation_v1_1(days, damping_days, damping_factors, populations, divi_start):
    """! Uses an ODE SECIR model allowing for asymptomatic infection with 6 different age groups.
    The model is not stratified by region.
    Virus-specific parameters are fixed and initial number of persons in the particular
    infection states are chosen randomly from defined ranges.
    :param days: Describes how many days we simulate within a single run.
    :param damping_days: Days, where a damping is possible
    :param damping_factors: Damping factors associated to the damping days.
    :param populations: List containing the population in each age group.
    :param divi_start: Start value for critical infected.
    :return: List containing the populations in each compartment used to initialize the run.
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

        # Initial number of people in each compartment with random numbers
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Exposed)] = random.uniform(
            0.00025, 0.0005) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedNoSymptoms)] = random.uniform(
            0.0001, 0.00035) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedNoSymptomsConfirmed)] = random.uniform(
            0.0001, 0.00035) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSymptoms)] = random.uniform(
            0.00007, 0.0001) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSymptomsConfirmed)] = 0
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSevere)] = random.uniform(
            0.00003, 0.00006) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedCritical)] = 1/6 * divi_start
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Recovered)] = random.uniform(
            0.002, 0.008) * populations[i]
        model.populations[AgeGroup(i),
                          Index_InfectionState(InfectionState.Dead)] = random.uniform(
            0, 0.0003) * populations[i]
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
    # set baseline to baseline matrix
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
        (num_groups, num_groups)) * 0

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline

    # Generate a damping matrix and assign it to the model
    damped_matrices = []
    damping_coefficients = []

    for i in np.arange(len(damping_days)):
        factor = damping_factors[i]
        # Catching to large values
        if factor > 20000:
            factor = 20000
        if factor < -20000:
            factor = -20000

        day = damping_days[i]
        f = np.float16(factor)
        damping = np.ones((num_groups, num_groups)
                          ) * f
        model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
            coeffs=(damping), t=day, level=0, type=0))

        damped_matrices.append(model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
            day+1))
        damping_coefficients.append(damping)

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

    return dataset_entries.tolist(), damped_matrices, damping_coefficients


def generate_data_v1_1(
        num_runs, parameters, path_out, population, divi_start,  num_days, filename,
        normalize=False, save_data=True):
    """
    Generating data for different population and parameter values. Using ODE-SECIR model without spatial resolution.

    :param num_runs: Number of generated samples
    :param parameters: Tuple containing the parameters used in the generation of the damping factors
    :param path_out: If the generated data should be stored, this is the associated path to the output file
    :param population: List of the age group resolved data
    :param divi_start: Start value for total critical infected
    :param num_days: Length of the simulated time series
    :param filename: Name of the file, which should contain the output
    :param normalize: Boolean, whether a log-normalization should be applied
    :param save_data: Boolean, whether the generated data should be stored or not.
    :returns" A dictionary of "inputs" and "labels"

    """
    num_dampings = 9

    data = {
        "inputs": [],
        "labels": [],
        "contact_matrix": [],
        "damping_day": [],
        "damping_coeff": []
    }

    # The number of days is the same as the sum of input and label width.
    # Since the first day of the input is day 0, we still need to subtract 1.
    days = num_days - 1

    # Load population data
    damping_data = generate_dampings(
        num_runs, num_days, num_dampings, parameters)
    damping_days = damping_data["damping_days"]
    damping_factors = damping_data["damping_coeff"]
    # show progess in terminal for longer runs
    # Due to the random structure, there's currently no need to shuffle the data

    # Setting seed to get the same initializations
    np.random.seed(42)
    for i in range(0, num_runs):

        data_run, _, _ = run_secir_simulation_v1_1(
            days, damping_days=damping_days[i], damping_factors=damping_factors[i], populations=population, divi_start=divi_start)
        data_run = np.nan_to_num(data_run)
        if len(data_run) != 90:
            print("Hääääääääääää??????")
        data['inputs'].append(data_run)
        data['labels'].append(data_run)

    if normalize:
        # logarithmic normalization
        transformer = FunctionTransformer(np.log1p, validate=True)

        # transform inputs and labels
        data['inputs'] = transform_data(data['inputs'], transformer, num_runs)
        data['labels'] = transform_data(data['labels'], transformer, num_runs)
    else:

        data_inputs = tf.convert_to_tensor(data['inputs'])
        data_labels = tf.convert_to_tensor(data['labels'])

    if save_data:
        # check if data directory exists. If necessary, create it.
        if not os.path.isdir(path_out):
            os.mkdir(path_out)

        # save dict to json file
        with open(os.path.join(path_out, filename), 'wb') as f:
            pickle.dump(data, f)
    return {"inputs": data_inputs,
            "labels": data_labels}


def reduce_to_infected(data, num_runs):
    """
    Reduce the age-resolved results to current number of infected.

    :param data: Data, that should be reduced
    :num_runs: Number of samples represented by data
    :returns: A numpy array of shape num_runs*num_days
    """
    infected = []
    infected_critical = []
    for i in np.arange(num_runs):
        d = data[i].numpy()
        d = np.reshape(d, (90, 6, 8))
        d = np.sum(d, axis=1)
        infected_critical.append(d[:, 5])
        infected.append(d[:, 2]+d[:, 3]+d[:, 4] + d[:, 5])
    return (np.asarray(infected), np.asarray(infected_critical))


def calc_min_error(data_prod, real_data):
    """
    Calculating the minimal error of an sample to the given real world data.

    :param data_prod: artificial data samples
    :param real_data: real data sample
    """
    neigh = KNeighborsRegressor(n_neighbors=1)
    neigh.fit(data_prod, data_prod)
    real_data = real_data.reshape(1, -1)
    nearest = neigh.predict(real_data)
    return sk.metrics.mean_absolute_percentage_error(real_data, nearest)


def evaluate_parameters(path_test_data, path_population, path_divi,  path_output):
    """
    Evaluating a set of possible parameter values in a grid search manner to find the optimal one.

    :param path_test_data: Path, containing the case data.
    :param path_population: Path, containing the population data.
    :param path_divi: Path, containing the divi case data
    :param path_output: Path to store the results
    """

    # Setting
    num_runs = 100
    num_days = 90
    n_populations = 250

    # Initializing result array
    results = []

    # Loading case and population data
    with open(os.path.join(path_population, "data_county_population.pickle"), 'rb') as f:
        population_data = pickle.load(f)

    with open(os.path.join(path_test_data, "data_county_cases_90_days.pickle"), "rb") as file:
        case_data = pickle.load(file)

    with open(os.path.join(path_divi, "data_county_divi_cases_90_days.pickle"), "rb") as file:
        divi_data = pickle.load(file)

    print("data loaded")
    # defining possible parameter values
    gamma_pos = [-2, -1, 0, 1, 2]
    delta_alpha = [1, 2, 3]
    p0 = [0.1, 0.2, 0.5, 0.7, 0.9]
    t1 = [-1, -0.5, 0, 0.5]  # , 0.5, 0, -0.5]
    delta_t1 = [0.5, 1, 1.5]  # , 1.5]
    t2 = [0, 0.5, 1, 1.5]
    delta_t2 = [3, 2, 1, 0.5]

    parameters = [
        (gamma_p, gamma_p - delta_a, p, t1_m, t1_m-d_t1, t2_m, t2_m-d_t2) for
        gamma_p in gamma_pos for delta_a in delta_alpha for p in p0 for t1_m in t1 for d_t1 in delta_t1
        for t2_m in t2 for d_t2 in delta_t2
    ]

    bar = Bar("Parameter checked", max=len(parameters))
    for parameter in parameters:
        error_inf = []
        error_crit_inf = []
        for i in np.arange(n_populations):
            data = generate_data_v1_1(
                num_runs=num_runs,
                parameters=parameter,
                population=population_data[i],
                divi_start=divi_data[i][0],
                num_days=num_days,
                path_out="",
                filename="",
                normalize=False,
                save_data=False)
            data_infected, data_critical = reduce_to_infected(
                data["inputs"], num_runs=num_runs)
            case_data_critical = divi_data[i]
            error_inf.append(calc_min_error(data_infected, case_data[i]))
            error_crit_inf.append(calc_min_error(
                data_critical, case_data_critical))
        mean_error_inf = np.mean(error_inf)
        var_error_inf = np.var(error_inf)
        mean_error_crit_inf = np.mean(error_crit_inf)
        var_error_crit_inf = np.var(error_crit_inf)
        res = {
            "gamma": parameter[0],
            "alpha": parameter[1],
            "p0": parameter[2],
            "t1_max": parameter[3],
            "t1_min": parameter[4],
            "t2_max": parameter[5],
            "t2_min": parameter[6],
            "mean_error_inf": mean_error_inf,
            "variance_inf": var_error_inf,
            "mean_error_crit_inf": mean_error_crit_inf,
            "variance_error_crit_inf": var_error_crit_inf
        }
        print(res)
        bar.next()
        results.append(res)
    bar.finish()
    # Saving results
    filename = "data_parameter_optimization.pickle"
    with open(os.path.join(path_output, filename), "wb") as f:
        pickle.dump(results, f)


if __name__ == "__main__":
    path = os.path.dirname(os.path.realpath(__file__))
    path = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')
    evaluate_parameters(path, path, path, path)
