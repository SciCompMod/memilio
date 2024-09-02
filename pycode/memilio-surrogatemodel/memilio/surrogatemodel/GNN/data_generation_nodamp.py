import copy
import os
import pickle
import random
import numpy as np
import pandas as pd

from progress.bar import Bar
from sklearn.preprocessing import FunctionTransformer
from datetime import date

from memilio.simulation import (AgeGroup, LogLevel, set_log_level)
from memilio.simulation.osecir import (Index_InfectionState, interpolate_simulation_result, ParameterStudy,
                                       InfectionState, Model, interpolate_simulation_result)
from memilio.epidata import geoModificationGermany as geoger
import memilio.epidata.getPopulationData as gpd
from .GNN_utils import (getBaselineMatrix, transform_mobility_directory,
                        make_graph, remove_confirmed_compartments, get_population)


def run_secir_groups_simulation(days, populations):
    """! Uses an ODE SECIR model allowing for asymptomatic infection with 6 
        different age groups. The model is not stratified by region. 
        Virus-specific parameters are fixed and initial number of persons 
        in the particular infection states are chosen randomly from defined ranges.
    @param Days Describes how many days we simulate within a single run.
    @param populations List containing the population in each age group.
    @return List containing the populations in each compartment used to initialize 
            the run.
   """

    set_log_level(LogLevel.Off)

    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1

    # get county ids
    countykey_list = geoger.get_county_ids(merge_eisenach=True, zfill=True)

    # Define age Groups
    groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']
    num_groups = len(groups)
    num_regions = len(populations)
    models = []

    # Initialize Parameters
    for region in range(num_regions):
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
            # Numbers are chosen heuristically based on experience
            model.populations[AgeGroup(i), Index_InfectionState(
                InfectionState.Exposed)] = random.uniform(
                0.00025, 0.005) * populations[region][i]
            model.populations[AgeGroup(i), Index_InfectionState(
                InfectionState.InfectedNoSymptoms)] = random.uniform(
                0.0001, 0.0035) * populations[region][i]
            model.populations[AgeGroup(i), Index_InfectionState(
                InfectionState.InfectedNoSymptomsConfirmed)] = 0
            model.populations[AgeGroup(i), Index_InfectionState(
                InfectionState.InfectedSymptoms)] = random.uniform(
                0.00007, 0.001) * populations[region][i]
            model.populations[AgeGroup(i), Index_InfectionState(
                InfectionState.InfectedSymptomsConfirmed)] = 0
            model.populations[AgeGroup(i), Index_InfectionState(
                InfectionState.InfectedSevere)] = random.uniform(
                0.00003, 0.0006) * populations[region][i]
            model.populations[AgeGroup(i), Index_InfectionState(
                InfectionState.InfectedCritical)] = random.uniform(
                0.00001, 0.0002) * populations[region][i]
            model.populations[AgeGroup(i), Index_InfectionState(
                InfectionState.Recovered)] = random.uniform(
                0.002, 0.08) * populations[region][i]
            model.populations[AgeGroup(i),
                              Index_InfectionState(InfectionState.Dead)] = random.uniform(
                0, 0.0003) * populations[region][i]
            model.populations.set_difference_from_group_total_AgeGroup(
                (AgeGroup(i), Index_InfectionState(InfectionState.Susceptible)),
                populations[region][i])

            # Compartment transition propabilities
            model.parameters.RelativeTransmissionNoSymptoms[AgeGroup(i)] = 0.5
            model.parameters.TransmissionProbabilityOnContact[AgeGroup(
                i)] = 0.1
            model.parameters.RecoveredPerInfectedNoSymptoms[AgeGroup(i)] = 0.09
            model.parameters.RiskOfInfectionFromSymptomatic[AgeGroup(i)] = 0.25
            model.parameters.SeverePerInfectedSymptoms[AgeGroup(i)] = 0.2
            model.parameters.CriticalPerSevere[AgeGroup(i)] = 0.25
            model.parameters.DeathsPerCritical[AgeGroup(i)] = 0.3
            model.parameters.MaxRiskOfInfectionFromSymptomatic[AgeGroup(
                i)] = 0.5

        # StartDay is the n-th day of the year
        model.parameters.StartDay = (
            date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

        # Load baseline and minimum contact matrix and assign them to the model
        baseline = getBaselineMatrix()

        model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline
        model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
            (num_groups, num_groups)) * 0

        # Apply mathematical constraints to parameters
        model.apply_constraints()
        models.append(model)

    directory = transform_mobility_directory()
    graph = make_graph(directory, num_regions, countykey_list, models)

    study = ParameterStudy(graph, 0, days, dt=dt, num_runs=1)
    study.run()

    graph_run = study.run()[0]
    results = interpolate_simulation_result(graph_run)

    for result_indx in range(len(results)):
        results[result_indx] = remove_confirmed_compartments(
            np.transpose(results[result_indx].as_ndarray()[1:, :]), num_groups)

    # Omit first column, as the time points are not of interest here.
    dataset_entry = copy.deepcopy(results)

    return dataset_entry


def generate_data(
        num_runs, path, input_width, days, save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often
   @param num_runs Number of times, the function run_secir_simulation is called.
   @param path Path, where the datasets are stored.
   @param input_width number of time steps used for model input.
   @param label_width number of time steps (days) used as model output/label.  
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """

    population = get_population()
    days_sum = days + input_width - 1

    data = {"inputs": [],
            "labels": [],
            }

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)

    for _ in range(num_runs):

        data_run = run_secir_groups_simulation(
            days_sum, population)

        inputs = np.asarray(data_run).transpose(1, 2, 0)[: input_width]
        data["inputs"].append(inputs)

        data["labels"].append(np.asarray(
            data_run).transpose(1, 2, 0)[input_width:])

        bar.next()

    bar.finish()

    if save_data:

        # we use a logistic transformer to reduce the
        # influence of ouliers and to correct the skewness of out dataset
        # as a result we obtaina dataset which is handled better by the NN

        num_groups = int(np.asarray(data['inputs']).shape[2] / 8)
        transformer = FunctionTransformer(np.log1p, validate=True)

        # Scale inputs
        inputs = np.asarray(
            data['inputs']).transpose(2, 0, 1, 3).reshape(num_groups * 8, -1)
        scaled_inputs = transformer.transform(inputs)
        original_shape_input = np.asarray(data['inputs']).shape

        # Reverse the reshape
        reshaped_back = scaled_inputs.reshape(original_shape_input[2],
                                              original_shape_input[0],
                                              original_shape_input[1],
                                              original_shape_input[3])

        # Reverse the transpose
        original_inputs = reshaped_back.transpose(1, 2, 0, 3)
        scaled_inputs = original_inputs.transpose(0, 3, 1, 2)

        # Scale labels
        labels = np.asarray(
            data['labels']).transpose(2, 0, 1, 3).reshape(num_groups * 8, -1)
        scaled_labels = transformer.transform(labels)
        original_shape_labels = np.asarray(data['labels']).shape

        # Reverse the reshape
        reshaped_back = scaled_labels.reshape(original_shape_labels[2],
                                              original_shape_labels[0],
                                              original_shape_labels[1],
                                              original_shape_labels[3])

        # Reverse the transpose
        original_labels = reshaped_back.transpose(1, 2, 0, 3)
        scaled_labels = original_labels.transpose(0, 3, 1, 2)

        all_data = {"inputs": scaled_inputs,
                    "labels": scaled_labels,
                    }

        # check if data directory exists. If necessary create it.
        if not os.path.isdir(path):
            os.mkdir(path)

        # save dict to json file
        with open(os.path.join(path, 'data_secir_age_groups.pickle'), 'wb') as f:
            pickle.dump(all_data, f)

    return data


if __name__ == "__main__":

    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'data_GNN_nodamp')

    input_width = 5
    days = 30
    num_runs = 1000

    generate_data(num_runs, path_data, input_width,
                  days)
