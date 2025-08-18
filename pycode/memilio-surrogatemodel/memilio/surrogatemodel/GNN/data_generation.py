#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Agatha Schmidt, Henrik Zunker
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################

import copy
import os
import pickle
import random
import time
import memilio.simulation as mio
import memilio.simulation.osecir as osecir
import numpy as np

from progress.bar import Bar

from memilio.simulation import (AgeGroup, LogLevel, set_log_level, Damping)
from memilio.simulation.osecir import (Index_InfectionState, interpolate_simulation_result, ParameterStudy,
                                       InfectionState, Model, interpolate_simulation_result)

from memilio.surrogatemodel.GNN.GNN_utils import (transform_mobility_directory,
                                                  make_graph, scale_data, getBaselineMatrix, remove_confirmed_compartments)
import memilio.surrogatemodel.utils.dampings as dampings

from enum import Enum


# Enumerate the different locations
class Location(Enum):
    Home = 0
    School = 1
    Work = 2
    Other = 3


# Define the start and the latest end date for the simulation
start_date = mio.Date(2019, 1, 1)
end_date = mio.Date(2021, 12, 31)


def set_covid_parameters(model, num_groups=6):
    """Setting COVID-parameters for the different age groups. 

    :param model: memilio model, whose parameters should be clarified 
    :param num_groups: Number of age groups 
    """

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

    for i, rho, muCR, muHI, muUH, muDU, tc, ti, th, tu in zip(range(num_groups),
                                                              TransmissionProbabilityOnContact, RecoveredPerInfectedNoSymptoms,
                                                              SeverePerInfectedSymptoms, CriticalPerSevere, DeathsPerCritical,
                                                              TimeInfectedNoSymptoms, TimeInfectedSymptoms,
                                                              TimeInfectedSevere, TimeInfectedCritical):
        # Compartment transition duration
        model.parameters.TimeExposed[AgeGroup(i)] = 3.335
        model.parameters.TimeInfectedNoSymptoms[AgeGroup(i)] = tc
        model.parameters.TimeInfectedSymptoms[AgeGroup(i)] = ti
        model.parameters.TimeInfectedSevere[AgeGroup(i)] = th
        model.parameters.TimeInfectedCritical[AgeGroup(i)] = tu

        # Compartment transition probabilities
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
    model.parameters.StartDay = start_date.day_in_year


def set_contact_matrices(model, data_dir, num_groups=6):
    """Setting the contact matrices for a model 

    :param model: memilio ODE-model, whose contact matrices should be modified. 
    :param data_dir: directory, where the contact data is stored (should contain folder "contacts")
    :param num_groups: Number of age groups considered 

    """

    contact_matrices = mio.ContactMatrixGroup(
        len(list(Location)), num_groups)
    locations = ["home", "school_pf_eig", "work", "other"]

    # Loading contact matrices for each location from .txt file
    for i, location in enumerate(locations):
        baseline_file = os.path.join(
            data_dir, "Germany", "contacts", "baseline_" + location + ".txt")

        contact_matrices[i] = mio.ContactMatrix(
            mio.read_mobility_plain(baseline_file),
        )
    model.parameters.ContactPatterns.cont_freq_mat = contact_matrices


def get_graph(num_groups, data_dir, mobility_directory):
    """ Generate the associated graph given the mobility data. 

    :param num_groups: Number of age groups 
    :param data_dir: Directory, where the contact data is stored (should contain folder "contacts")
    :param mobility_directory: Directory containing the mobility data
    """
    # Generating and Initializing the model
    model = Model(num_groups)
    set_covid_parameters(model)
    set_contact_matrices(model, data_dir)

    # Generating the graph
    graph = osecir.ModelGraph()

    # Setting the parameters
    scaling_factor_infected = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
    scaling_factor_icu = 1.0
    tnt_capacity_factor = 7.5 / 100000.

    # Path containing the population data
    data_dir_Germany = os.path.join(data_dir, "Germany")
    pydata_dir = os.path.join(data_dir_Germany, "pydata")

    path_population_data = os.path.join(pydata_dir,
                                        "county_current_population.json")

    # Setting node information based on model parameters
    mio.osecir.set_nodes(
        model.parameters,
        mio.Date(start_date.year,
                 start_date.month, start_date.day),
        mio.Date(end_date.year,
                 end_date.month, end_date.day), pydata_dir,
        path_population_data, True, graph, scaling_factor_infected,
        scaling_factor_icu, tnt_capacity_factor, 0, False)

    # Setting edge information based on the mobility data
    mio.osecir.set_edges(
        mobility_directory, graph, len(Location))

    return graph


def run_secir_groups_simulation1(days, damping_days, damping_factors, graph, num_groups=6, start_date=mio.Date(2020, 6, 1)):
    """ Uses an ODE SECIR model allowing for asymptomatic infection with 6
        different age groups. The model is not stratified by region.
        Virus-specific parameters are fixed and initial number of person
        in the particular infection states are chosen randomly from defined ranges.

    :param days: Number of days simulated within a single run.
    :param damping_days: Days, where a damping is applied
    :param damping_factors: damping factors associated to the damping days 
    :param graph: Graph initialized for the start_date with the population data which
            is sampled during the run.
    :param num_groups: Number of age groups considered in the simulation
    :param start_date: Date, when the simulation starts
    :returns: List containing the populations in each compartment used to initialize
            the run.
   """
    if len(damping_days) != len(damping_factors):
        raise ValueError("Length of damping_days and damping_factors differ!")

    min_date = mio.Date(2020, 6, 1)
    min_date_num = min_date.day_in_year
    start_date_num = start_date.day_in_year - min_date_num

    # Load the ground truth data
    pydata_dir = os.path.join(data_dir, "Germany", "pydata")
    ground_truth_dir = os.path.join(
        pydata_dir, "ground_truth_all_nodes.pickle")
    with open(ground_truth_dir, 'rb') as f:
        ground_truth_all_nodes = pickle.load(f)

    # Initialize model for each node, using the population data and sampling the number of
    # individuals in the different compartments
    for node_indx in range(graph.num_nodes):
        model = graph.get_node(node_indx).property
        data = ground_truth_all_nodes[node_indx][start_date_num]

        # Iterate over the different age groups
        for i in range(num_groups):
            age_group = AgeGroup(i)
            pop_age_group = model.populations.get_group_total_AgeGroup(
                age_group)

            # Generating valid, noisy configuration of the compartments
            valid_configuration = False
            while valid_configuration is False:
                comp_data = data[:, i]
                ratios = np.random.uniform(0.8, 1.2, size=8)
                init_data = comp_data * ratios
                if np.sum(init_data[1:]) < pop_age_group:
                    valid_configuration = True

            # Set the populations for the different compartments
            model.populations[age_group, Index_InfectionState(
                InfectionState.Exposed)] = init_data[1]
            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedNoSymptoms)] = init_data[2]
            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedSymptoms)] = init_data[3]
            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedSevere)] = init_data[4]
            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedCritical)] = init_data[5]
            model.populations[age_group, Index_InfectionState(
                InfectionState.Recovered)] = init_data[6]
            model.populations[age_group, Index_InfectionState(
                InfectionState.Dead)] = init_data[7]
            model.populations.set_difference_from_group_total_AgeGroup((
                age_group, InfectionState.Susceptible), pop_age_group)

        # Introduce the damping information, in general dampings can be local, but till now the code just allows global dampings
        damped_matrices = []
        damping_coefficients = []

        for i in np.arange(len(damping_days)):
            day = damping_days[i]
            factor = damping_factors[i]

            damping = np.ones((num_groups, num_groups)
                              ) * np.float16(factor)
            model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
                coeffs=(damping), t=day, level=0, type=0))
            damped_matrices.append(model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
                day+1))
            damping_coefficients.append(damping)

        # Apply mathematical constraints to parameters
        model.apply_constraints()

        # set model to graph
        graph.get_node(node_indx).property.populations = model.populations

    # Start simulation
    study = ParameterStudy(graph, 0, days, dt=0.5, num_runs=1)
    start_time = time.perf_counter()
    study.run()
    runtime = time.perf_counter() - start_time

    graph_run = study.run()[0]
    results = interpolate_simulation_result(graph_run)

    for result_indx in range(len(results)):
        results[result_indx] = remove_confirmed_compartments(
            np.asarray(results[result_indx]), num_groups)

    dataset_entry = copy.deepcopy(results)

    return dataset_entry, damped_matrices, damping_coefficients, runtime


def run_secir_groups_simulation(days, damping_days, damping_factors, graph, num_groups=6, start_date=mio.Date(2020, 6, 1)):
    """ Uses an ODE SECIR model allowing for asymptomatic infection with 6
        different age groups. The model is not stratified by region.
        Virus-specific parameters are fixed and initial number of person
        in the particular infection states are chosen randomly from defined ranges.

    :param days: Number of days simulated within a single run.
    :param damping_days: Days, where a damping is applied
    :param damping_factors: damping factors associated to the damping days 
    :param graph: Graph initialized for the start_date with the population data which
            is sampled during the run.
    :param num_groups: Number of age groups considered in the simulation
    :param start_date: Date, when the simulation starts
    :returns: List containing the populations in each compartment used to initialize
            the run.
   """
    if len(damping_days) != len(damping_factors):
        raise ValueError("Length of damping_days and damping_factors differ!")

    min_date = mio.Date(2020, 6, 1)
    min_date_num = min_date.day_in_year
    start_date_num = start_date.day_in_year - min_date_num

    # Load the ground truth data
    pydata_dir = os.path.join(data_dir, "Germany", "pydata")
    upper_bound_dir = os.path.join(
        pydata_dir, "ground_truth_upper_bound.pickle")
    lower_bound_dir = os.path.join(
        pydata_dir, "ground_truth_lower_bound.pickle")
    with open(upper_bound_dir, 'rb') as f:
        ground_truth_upper_bound = pickle.load(f)

    with open(lower_bound_dir, 'rb') as f:
        ground_truth_lower_bound = pickle.load(f)

    # Initialize model for each node, using the population data and sampling the number of
    # individuals in the different compartments
    for node_indx in range(graph.num_nodes):
        model = graph.get_node(node_indx).property
        max_data = ground_truth_upper_bound[node_indx]
        min_data = ground_truth_lower_bound[node_indx]

        # Iterate over the different age groups
        for i in range(num_groups):
            age_group = AgeGroup(i)
            pop_age_group = model.populations.get_group_total_AgeGroup(
                age_group)

            # Generating valid, noisy configuration of the compartments
            valid_configuration = False
            while valid_configuration is False:
                init_data = np.asarray([0 for _ in range(8)])
                max_val = max_data[:, i]
                min_val = min_data[:, i]
                for j in range(1, 8):
                    init_data[j] = random.uniform(min_val[j], max_val[j])
                if np.sum(init_data[1:]) < pop_age_group:
                    valid_configuration = True

            # Set the populations for the different compartments
            model.populations[age_group, Index_InfectionState(
                InfectionState.Exposed)] = init_data[1]
            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedNoSymptoms)] = init_data[2]
            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedSymptoms)] = init_data[3]
            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedSevere)] = init_data[4]
            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedCritical)] = init_data[5]
            model.populations[age_group, Index_InfectionState(
                InfectionState.Recovered)] = init_data[6]
            model.populations[age_group, Index_InfectionState(
                InfectionState.Dead)] = init_data[7]
            model.populations.set_difference_from_group_total_AgeGroup((
                age_group, InfectionState.Susceptible), pop_age_group)

        # Introduce the damping information, in general dampings can be local, but till now the code just allows global dampings
        damped_matrices = []
        damping_coefficients = []

        for i in np.arange(len(damping_days)):
            day = damping_days[i]
            factor = damping_factors[i]

            damping = np.ones((num_groups, num_groups)
                              ) * np.float16(factor)
            model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
                coeffs=(damping), t=day, level=0, type=0))
            damped_matrices.append(model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
                day+1))
            damping_coefficients.append(damping)

        # Apply mathematical constraints to parameters
        model.apply_constraints()

        # set model to graph
        graph.get_node(node_indx).property.populations = model.populations

    # Start simulation
    study = ParameterStudy(graph, 0, days, dt=0.5, num_runs=1)
    start_time = time.perf_counter()
    study.run()
    runtime = time.perf_counter() - start_time

    graph_run = study.run()[0]
    results = interpolate_simulation_result(graph_run)

    for result_indx in range(len(results)):
        results[result_indx] = remove_confirmed_compartments(
            np.asarray(results[result_indx]), num_groups)

    dataset_entry = copy.deepcopy(results)

    return dataset_entry, damped_matrices, damping_coefficients, runtime


def generate_data(
        num_runs, data_dir, path, input_width, label_width, save_data=True,
        transform=True, damping_method="classic", max_number_damping=3):
    """ Generate dataset by calling run_secir_simulation (num_runs)-often

    :param num_runs: Number of times, the function run_secir_simulation is called.
    :param data_dir: Directory with all data needed to initialize the models.
    :param path: Path, where the datasets are stored.
    :param input_width: number of time steps used for model input.
    :param label_width: number of time steps (days) used as model output/label.
    :param save_data: Option to deactivate the save of the dataset. Per default True.
    :param damping_method: String specifying the damping method, that should be used. Possible values "classic", "active", "random".
    :param max_number_damping: Maximal number of possible dampings. 
    :returns: Dictionary containing inputs, labels, contact_matrix, damping_days and damping_factors
   """
    set_log_level(mio.LogLevel.Error)
    days = label_width + input_width - 1

    # Preparing output dictionary
    data = {
        "inputs": [],
        "labels": [],
        "contact_matrix": [],
        "damping_days": [],
        "damping_factors": []
    }

    # Setting basic parameter
    num_groups = 6
    mobility_dir = data_dir + "/Germany/mobility/commuter_mobility_2022.txt"
    graph = get_graph(num_groups, data_dir, mobility_dir)
    # Define possible start dates for the simulation
    # start_dates = [
    #    mio.Date(2020, 6, 1),
    #    mio.Date(2020, 7, 1),
    #    mio.Date(2020, 8, 1),
    #    mio.Date(2020, 9, 1),
    #    mio.Date(2020, 10, 1),
    #    mio.Date(2020, 11, 1),
    #    mio.Date(2020, 12, 1)
    # ]

    # show progess in terminal for longer runs
    # Due to the random structure, there is currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)

    times = []
    for i in range(0, num_runs):
        # start_date = start_dates[i % len(start_dates)]
        # Generate random damping days and damping factors
        if max_number_damping > 0:
            damping_days, damping_factors = dampings.generate_dampings(
                days, max_number_damping, method=damping_method, min_distance=2,
                min_damping_day=2)
        else:
            damping_days = []
            damping_factors = []
        # Run simulation
        data_run, damped_matrices, damping_coefficients, t_run = run_secir_groups_simulation(
            days, damping_days, damping_factors, graph, num_groups, start_date)

        times.append(t_run)

        inputs = np.asarray(data_run).transpose(1, 0, 2)[: input_width]
        labels = np.asarray(data_run).transpose(1, 0, 2)[input_width:]

        data["inputs"].append(inputs)
        data["labels"].append(labels)
        data["contact_matrix"].append(np.array(damped_matrices))
        data["damping_factors"].append(damping_coefficients)
        data["damping_days"].append(damping_days)

        bar.next()

    bar.finish()

    print(
        f"For Days = {days}, AVG runtime: {np.mean(times)}s, Median runtime: {np.median(times)}s")

    if save_data:
        if transform:
            inputs, labels = scale_data(data)

        all_data = {"inputs": inputs,
                    "labels": labels,
                    "damping_day": data["damping_days"],
                    "contact_matrix": data["contact_matrix"],
                    "damping_coeff": data["damping_factors"]
                    }

        # check if data directory exists. If necessary create it.
        if not os.path.isdir(path):
            os.mkdir(path)

        # generate the filename
        if num_runs < 1000:
            filename = 'GNN_data_%ddays_%ddampings_' % (
                label_width, max_number_damping) + damping_method+'%d.pickle' % (num_runs)
        else:
            filename = 'GNN_data_%ddays_%ddampings_' % (
                label_width, max_number_damping) + damping_method+'%dk.pickle' % (num_runs//1000)

        # save dict to pickle file
        with open(os.path.join(path, filename), 'wb') as f:
            pickle.dump(all_data, f)

    return data


if __name__ == "__main__":

    path = os.getcwd()
    path_output = os.path.join(os.getcwd(), 'saves')
    # data_dir = os.path.join(os.getcwd(), 'data')
    data_dir = "/localdata1/hege_mn/memilio/data"
    input_width = 5
    number_of_dampings = 0
    num_runs = 100
    label_width_list = [30]

    random.seed(10)
    for label_width in label_width_list:
        generate_data(num_runs=num_runs,
                      data_dir=data_dir,
                      path=path_output,
                      input_width=input_width,
                      label_width=label_width,
                      save_data=True,
                      damping_method="active",
                      max_number_damping=number_of_dampings)
