import copy
import os
import pickle
import random
import numpy as np
import time
import memilio.simulation as mio
import memilio.simulation.osecir as osecir
from datetime import date
from enum import Enum

from progress.bar import Bar
from sklearn.preprocessing import FunctionTransformer

from memilio.simulation import (AgeGroup, LogLevel, set_log_level, Damping)
from memilio.simulation.osecir import (Index_InfectionState, interpolate_simulation_result, ParameterStudy,
                                       InfectionState, Model,
                                       interpolate_simulation_result, ModelGraph, set_edges)
from memilio.epidata import geoModificationGermany as geoger

from memilio.surrogatemodel.GNN.GNN_utils import (transform_mobility_directory,
                                                  make_graph, scale_data)
from memilio.surrogatemodel.utils_surrogatemodel import (
    getBaselineMatrix, remove_confirmed_compartments, get_population, getMinimumMatrix)


class Location(Enum):
    Home = 0
    School = 1
    Work = 2
    Other = 3


start_date = mio.Date(2019, 1, 1)
end_date = mio.Date(2020, 12, 31)


def set_covid_parameters(model, num_groups=6):
    for i in range(num_groups):
        # Compartment transition duration
        model.parameters.TimeExposed[AgeGroup(i)] = 3.2
        model.parameters.TimeInfectedNoSymptoms[AgeGroup(i)] = 2.
        model.parameters.TimeInfectedSymptoms[AgeGroup(i)] = 6.
        model.parameters.TimeInfectedSevere[AgeGroup(i)] = 12.
        model.parameters.TimeInfectedCritical[AgeGroup(i)] = 8.

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
    model.parameters.StartDay = start_date.day_in_year


def set_contact_matrices(model, data_dir, num_groups=6):
    contact_matrices = mio.ContactMatrixGroup(
        len(list(Location)), num_groups)
    locations = ["home", "school_pf_eig", "work", "other"]

    for i, location in enumerate(locations):
        baseline_file = os.path.join(
            data_dir, "contacts", "baseline_" + location + ".txt")
        minimum_file = os.path.join(
            data_dir, "contacts", "minimum_" + location + ".txt")
        contact_matrices[i] = mio.ContactMatrix(
            mio.read_mobility_plain(baseline_file),
            mio.read_mobility_plain(minimum_file)
        )
    model.parameters.ContactPatterns.cont_freq_mat = contact_matrices


def get_graph(num_groups, data_dir):
    model = Model(num_groups)
    set_covid_parameters(model)
    set_contact_matrices(model, data_dir)

    graph = osecir.ModelGraph()

    scaling_factor_infected = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
    scaling_factor_icu = 1.0
    tnt_capacity_factor = 7.5 / 100000.

    path_population_data = os.path.join(
        data_dir, "pydata", "Germany",
        "county_current_population.json")

    mio.osecir.set_nodes(
        model.parameters,
        mio.Date(start_date.year,
                 start_date.month, start_date.day),
        mio.Date(end_date.year,
                 end_date.month, end_date.day), data_dir,
        path_population_data, True, graph, scaling_factor_infected,
        scaling_factor_icu, tnt_capacity_factor, 0, False)

    mio.osecir.set_edges(
        data_dir, graph, len(Location))

    return graph


def run_secir_groups_simulation(days, graph, dampings, num_groups=6):
    """! Uses an ODE SECIR model allowing for asymptomatic infection 
        with 6 different age groups. The model is not stratified by region. 
        Virus-specific parameters are fixed and initial number of persons 
        in the particular infection states are chosen randomly from defined ranges.
    @param days (int) Describes how many days we simulate within a single run.
    @param Graph Graph initilized for the start_date with the population data which
            is sampled during the run.
    @param damping_day (int) The day when damping is applied.
    @return List containing the populations in each compartment 
            used to initialize the run.
   """
    for node_indx in range(graph.num_nodes):
        model = graph.get_node(node_indx).property

        # Set parameters
        # TODO: Put This in the draw_sample function in the ParameterStudy
        for i in range(num_groups):
            age_group = AgeGroup(i)
            pop_age_group = model.populations.get_group_total_AgeGroup(
                age_group)

            # Initial number of people in each compartment with random numbers
            # Numbers are chosen heuristically based on experience
            model.populations[age_group, Index_InfectionState(InfectionState.Exposed)] = random.uniform(
                0.00025, 0.005) * pop_age_group
            model.populations[age_group, Index_InfectionState(InfectionState.InfectedNoSymptoms)] = random.uniform(
                0.0001, 0.0035) * pop_age_group
            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedNoSymptomsConfirmed)] = 0
            model.populations[age_group, Index_InfectionState(InfectionState.InfectedSymptoms)] = random.uniform(
                0.00007, 0.001) * pop_age_group
            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedSymptomsConfirmed)] = 0
            model.populations[age_group, Index_InfectionState(InfectionState.InfectedSevere)] = random.uniform(
                0.00003, 0.0006) * pop_age_group
            model.populations[age_group, Index_InfectionState(InfectionState.InfectedCritical)] = random.uniform(
                0.00001, 0.0002) * pop_age_group
            model.populations[age_group, Index_InfectionState(InfectionState.Recovered)] = random.uniform(
                0.002, 0.08) * pop_age_group
            model.populations[age_group, Index_InfectionState(InfectionState.Dead)] = random.uniform(
                0, 0.0003) * pop_age_group
            model.populations.set_difference_from_group_total_AgeGroup(
                (age_group, Index_InfectionState(InfectionState.Susceptible)),
                pop_age_group)

        # Apply mathematical constraints to parameters
        model.apply_constraints()

        # Generate a damping matrix and assign it to the model
        # TODO: This can be done outside and is (currently) static for all models
        damped_matrices = []
        damping_coeff = []
        for day in dampings:

            # generat a random damping factor
            damping = np.ones((num_groups, num_groups)
                              ) * np.float16(random.uniform(0, 0.5))

            # add damping to model
            model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
                coeffs=(damping), t=day, level=0, type=0))

            damped_matrices.append(model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
                day+1))
            damping_coeff.append(damping[0][0])

        # Apply mathematical constraints to parameters
        model.apply_constraints()
        # set model to graph
        graph.get_node(node_indx).property.populations = model.populations
        graph.get_node(node_indx).property.parameters = model.parameters

    study = ParameterStudy(graph, 0, days, dt=0.5, num_runs=1)
    start_time = time.time()
    study.run()
    print("Simulation took: ", time.time() - start_time)

    graph_run = study.run()[0]
    results = interpolate_simulation_result(graph_run)

    for result_indx in range(len(results)):
        results[result_indx] = remove_confirmed_compartments(
            np.transpose(results[result_indx].as_ndarray()[1:, :]), num_groups)

    # Omit first column, as the time points are not of interest here.
    dataset_entry = copy.deepcopy(results)

    return dataset_entry, damped_matrices, dampings, damping_coeff


def generate_dampings_withshadowdamp(number_of_dampings, days, min_distance, min_damping_day, n_runs):
    """! Draw damping days while keeping a minimum distance between the 
        damping days. This method aims to create a uniform ditribution of 
        drawn damping days. 
   @param num_of_dampings (int) Number of dampings that have to be drawn.
   @param days (int) Number of days which are simulated (label_width).
   @param min_distance (int) The minimum number of days between two dampings. 
   @param min_damping_day (int) The earliest day of the simualtion where a damping 
        can take place.   
   @param n_runs 8int) Number of simulation runs. 
   """

    all_dampings = []
    count_runs = 0
    count_shadow = 0
    while len(all_dampings) < n_runs:

        days_list = list(range((min_damping_day), days))
        dampings = []
        if count_shadow < 2:
            for i in range(number_of_dampings):

                damp = random.choice(days_list)
                days_before = list(range(damp-(min_distance), damp))
                days_after = list(range(damp, damp+(min_distance+1)))
                dampings.append(damp)
                days_list = [ele for ele in days_list if ele not in (
                    days_before+days_after)]
        else:
            # chose a forbidden damping
            damp = random.choice(
                list(range((0-min_distance), 0)) + list(range(days+1, (days+min_distance+1))))

            days_before = list(range(damp-(min_distance), damp))
            days_after = list(range(damp, damp+(min_distance+1)))
            days_list = [ele for ele in days_list if ele not in (
                days_before+days_after)]
            dampings.append(damp)
            for i in range(number_of_dampings):

                damp = random.choice(days_list)
                days_before = list(range(damp-(min_distance), damp))
                days_after = list(range(damp, damp+(min_distance+1)))
                dampings.append(damp)
                days_list = [ele for ele in days_list if ele not in (
                    days_before+days_after)]
                count_shadow = 0

        forbidden_damping_values = list(
            range((0-min_distance), 0)) + list(range(days+1, (days+min_distance+1)))
        dampings = [
            ele for ele in dampings if ele not in forbidden_damping_values]
        count_runs += 1
        count_shadow += 1
        # select first or last five dampings
        if len(dampings) >= number_of_dampings:
            all_dampings.append(sorted(dampings))

    return np.asarray(all_dampings)


def generate_data(
        num_runs, data_dir, path, input_width, label_width, number_of_dampings, save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often
   @param num_runs Number of times, the function run_secir_simulation is called.
   @param data_dir Directory with all data needed to initialize the models.
   @param path Path, where the datasets are stored.
   @param input_width Number of time steps used for model input.
   @param label_width Number of time steps (days) used as model output/label.  
   @param number_of_dampings (int) The number of contact change points applied to the simulation. 
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """
    set_log_level(mio.LogLevel.Error)
    days_sum = label_width+input_width-1

    num_groups = 6
    graph = get_graph(num_groups, data_dir)

    # generate dampings
    damping_days = generate_dampings_withshadowdamp(
        number_of_dampings=number_of_dampings, days=label_width,
        min_distance=7, min_damping_day=input_width, n_runs=num_runs
    )

    # all data including damping information
    all_data = {"inputs": [],
                "labels": [],
                "damping_coeff": [],
                "damping_day": [],
                "damped_matrix": []}

    # data that needs to be scaled
    data = {"inputs": [],
            "labels": [],
            "damping_coeff": [],
            "damping_day": [],
            "damped_matrix": []}

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)

    model_params = copy.deepcopy(graph.get_node(0).property.parameters)

    for i in range(num_runs):
        params_run = copy.deepcopy(model_params)
        # reset contact matrix in each node
        for node_indx in range(graph.num_nodes):
            graph.get_node(node_indx).property.parameters = params_run
        data_run, damped_contact_matrix, damping_days_s, damping_factor = run_secir_groups_simulation(
            days_sum,  graph,  damping_days[i])

        inputs = np.asarray(data_run).transpose(1, 2, 0)[:input_width]
        data["inputs"].append(inputs)
        data["labels"].append(np.asarray(
            data_run).transpose(1, 2, 0)[input_width:])
        data["damping_coeff"].append(damping_factor)
        data["damping_day"].append(damping_days_s)
        data["damped_matrix"].append(damped_contact_matrix)

        bar.next()

    bar.finish()

    if save_data:

        scaled_inputs, scaled_labels = scale_data(data)

        all_data = {"inputs": scaled_inputs,
                    "labels": scaled_labels,
                    "damping_coeff": data['damping_coeff'],
                    "damping_day": data['damping_day'],
                    "damped_matrix": data['damped_matrix']}

        # check if data directory exists. If necessary create it.
        if not os.path.isdir(path):
            os.mkdir(path)

        # save dict to json file
        with open(os.path.join(path, 'data_secir_age_groups.pickle'), 'wb') as f:
            pickle.dump(all_data, f)
    return data


if __name__ == "__main__":

    input_width = 5
    label_width = 95
    number_of_dampings = 3
    num_runs = 5
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'data_GNN_with_'+str(number_of_dampings)+'_dampings_test')

    data_dir = os.path.join(os.getcwd(), 'data')

    generate_data(num_runs, data_dir, path_data, input_width,
                  label_width, number_of_dampings, save_data=False)
