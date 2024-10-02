import copy
import os
import pickle
import random
import time
import memilio.simulation as mio
import memilio.simulation.osecir as osecir
import numpy as np

from progress.bar import Bar

from datetime import date

from memilio.simulation import (AgeGroup, LogLevel, set_log_level)
from memilio.simulation.osecir import (Index_InfectionState, interpolate_simulation_result, ParameterStudy,
                                       InfectionState, Model, interpolate_simulation_result)
from memilio.epidata import geoModificationGermany as geoger
from memilio.surrogatemodel.GNN.GNN_utils import (transform_mobility_directory,
                                                  make_graph, scale_data)
from memilio.surrogatemodel.utils_surrogatemodel import (
    getBaselineMatrix, remove_confirmed_compartments, get_population)
from enum import Enum


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


def run_secir_groups_simulation(days, graph, num_groups=6):
    """! Uses an ODE SECIR model allowing for asymptomatic infection with 6
        different age groups. The model is not stratified by region.
        Virus-specific parameters are fixed and initial number of persons
        in the particular infection states are chosen randomly from defined ranges.
    @param Days Describes how many days we simulate within a single run.
    @param Graph Graph initilized for the start_date with the population data which
            is sampled during the run.
    @return List containing the populations in each compartment used to initialize
            the run.
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

        # set model to graph
        graph.get_node(node_indx).property.populations = model.populations

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

    return dataset_entry


def generate_data(
        num_runs, data_dir, path, input_width, days, save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often
   @param num_runs Number of times, the function run_secir_simulation is called.
   @param data_dir Directory with all data needed to initialize the models.
   @param path Path, where the datasets are stored.
   @param input_width number of time steps used for model input.
   @param label_width number of time steps (days) used as model output/label.  
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """
    set_log_level(mio.LogLevel.Error)
    days_sum = days + input_width - 1

    data = {"inputs": [],
            "labels": [],
            }

    num_groups = 6
    graph = get_graph(num_groups, data_dir)

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)

    for _ in range(num_runs):

        data_run = run_secir_groups_simulation(
            days_sum, graph)

        inputs = np.asarray(data_run).transpose(1, 2, 0)[: input_width]
        data["inputs"].append(inputs)

        data["labels"].append(np.asarray(
            data_run).transpose(1, 2, 0)[input_width:])

        bar.next()

    bar.finish()

    if save_data:

        scaled_inputs, scaled_labels = scale_data(data)

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
        'data_GNN_nodamp_test')

    data_dir = os.path.join(os.getcwd(), 'data')

    input_width = 5
    days = 30
    num_runs = 1

    generate_data(num_runs, data_dir,  path_data, input_width,
                  days, save_data=True)
