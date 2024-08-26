import copy
import os
import pickle
import random
import json
from datetime import date
import numpy as np
 
from progress.bar import Bar
from sklearn.preprocessing import FunctionTransformer

from memilio.simulation import (AgeGroup, LogLevel, set_log_level)
from memilio.simulation.osecir import (Index_InfectionState, interpolate_simulation_result, ParameterStudy,
                                      InfectionState, Model, ModelGraph, 
                                      interpolate_simulation_result, set_edges)
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import transformMobilityData as tmd
from memilio.epidata import getDataIntoPandasDataFrame as gd


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
                              Index_InfectionState(InfectionState.Dead)] =  random.uniform(
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

    graph = make_graph(num_regions, countykey_list, models)

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


def remove_confirmed_compartments(dataset_entries, num_groups):
    """! The compartments which contain confirmed cases are not needed and are 
        therefore omitted by summarizing the confirmed compartment with the 
        original compartment. 
    @param dataset_entries Array that contains the compartmental data with 
            confirmed compartments. 
    @param num_groups Number of age groups.
    @return Array that contains the compartmental data without confirmed compartments. 
   """
    
    new_dataset_entries = []
    for i in dataset_entries : 
      dataset_entries_reshaped  = i.reshape(
          [num_groups, int(np.asarray(dataset_entries).shape[1]/num_groups)]
          )
      sum_inf_no_symp = np.sum(dataset_entries_reshaped [:, [2, 3]], axis=1)
      sum_inf_symp = np.sum(dataset_entries_reshaped [:, [4, 5]], axis=1)
      dataset_entries_reshaped[:, 2] = sum_inf_no_symp
      dataset_entries_reshaped[:, 4] = sum_inf_symp
      new_dataset_entries.append(
          np.delete(dataset_entries_reshaped , [3, 5], axis=1).flatten()
          )
    return new_dataset_entries


def get_population(path="data/pydata/Germany/county_population.json"):
    """! loads population data 
    @param path Path to population file. 
    @return List with all 400 populations and 6 age groups. 
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

def make_graph(num_regions, countykey_list, models):
    """! 
    @param num_regions Number (int) of counties that should be added to the 
            grap-ODE model. Equals 400 for whole Germany. 
    @param countykey_list List of keys/IDs for each county. 
    @models models List of osecir Model with one model per population. 
    @return graph Graph-ODE model. 
   """
    graph = ModelGraph()
    for i in range(num_regions):
        graph.add_node(int(countykey_list[i]), models[i])

    # get mobility data directory
    arg_dict = gd.cli("commuter_official")

    directory = arg_dict['out_folder'].split('/pydata')[0]
    directory = os.path.join(directory, 'mobility/')  

    # Merge Eisenach and Wartbugkreis in Input Data 
    tmd.updateMobility2022(directory, mobility_file='twitter_scaled_1252')
    tmd.updateMobility2022(directory, mobility_file='commuter_migration_scaled')

    num_locations = 4

    set_edges(os.path.abspath(os.path.join(directory, os.pardir)), 
                            graph, num_locations)
    return graph
    

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

            data["labels"].append(np.asarray(data_run).transpose(1, 2, 0)[input_width:])

            bar.next()

    bar.finish()

    if save_data:
            num_groups = int(np.asarray(data['inputs']).shape[2] / 8)
            transformer = FunctionTransformer(np.log1p, validate=True)

            # Scale inputs
            inputs = np.asarray(
                data['inputs']).transpose(2, 0, 1, 3).reshape(num_groups * 8, -1)
            scaled_inputs = transformer.transform(inputs)
            original_shape_input = np.asarray(data['inputs']).shape
            
            # Step 1: Reverse the reshape
            reshaped_back = scaled_inputs.reshape(original_shape_input[2], 
                                                  original_shape_input[0], 
                                                  original_shape_input[1], 
                                                  original_shape_input[3])

            # Step 2: Reverse the transpose
            original_inputs = reshaped_back.transpose(1, 2, 0, 3)
            scaled_inputs = original_inputs.transpose(0, 3, 1, 2)

            
            # Scale labels
            labels = np.asarray(
                data['labels']).transpose(2, 0, 1, 3).reshape(num_groups * 8, -1)
            scaled_labels = transformer.transform(labels)
            original_shape_labels = np.asarray(data['labels']).shape
            
            # Step 1: Reverse the reshape
            reshaped_back = scaled_labels.reshape(original_shape_labels[2], 
                                                  original_shape_labels[0], 
                                                  original_shape_labels[1], 
                                                  original_shape_labels[3])

            # Step 2: Reverse the transpose
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
    
