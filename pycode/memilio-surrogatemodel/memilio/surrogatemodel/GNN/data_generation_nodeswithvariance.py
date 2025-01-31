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
                                                  make_graph, scale_data, getBaselineMatrix, remove_confirmed_compartments)
from memilio.surrogatemodel.utils_surrogatemodel import (
    get_population)
from enum import Enum


class Location(Enum):
    Home = 0
    School = 1
    Work = 2
    Other = 3


start_date = mio.Date(2019, 1, 1)
end_date = mio.Date(2020, 12, 31)


def set_covid_parameters(model, num_groups=6):

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
    model.parameters.StartDay = start_date.day_in_year


def set_contact_matrices(model, data_dir, num_groups=6):
    contact_matrices = mio.ContactMatrixGroup(
        len(list(Location)), num_groups)
    locations = ["home", "school_pf_eig", "work", "other"]

    for i, location in enumerate(locations):
        baseline_file = os.path.join(
            data_dir, "contacts", "baseline_" + location + ".txt")

        contact_matrices[i] = mio.ContactMatrix(
            mio.read_mobility_plain(baseline_file),
        )
    model.parameters.ContactPatterns.cont_freq_mat = contact_matrices


def get_graph(num_groups, data_dir, mobility_directory):
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
        os.path.dirname(os.path.realpath(
            mobility_directory)), graph, len(Location))

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

    p_infected = random.uniform(0.0001, 0.05)
    p_exposed = p_infected * random.uniform(0.1, 5)
    p_ins = p_infected*random.uniform(0.1, 5)
    p_hosp = p_infected * random.uniform(0.001, 1)
    p_icu = p_hosp*random.uniform(0.001, 1)
    p_dead = p_icu*random.uniform(0.001, 1)

    sum_randoms = p_infected + p_exposed+p_ins + p_hosp + p_icu + p_dead
    p_recovered = random.uniform(0, (1-sum_randoms))

    for node_indx in range(graph.num_nodes):
        model = graph.get_node(node_indx).property

        # Set parameters
        # TODO: Put This in the draw_sample function in the ParameterStudy
        for i in range(num_groups):
            age_group = AgeGroup(i)
            pop_age_group = model.populations.get_group_total_AgeGroup(
                age_group)

            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedSymptoms)] = pop_age_group * p_infected * random.uniform(0.1, 1)

            model.populations[age_group, Index_InfectionState(
                InfectionState.Exposed)] = model.populations[age_group,
                                                             InfectionState.InfectedSymptoms].value * p_exposed * random.uniform(0.1, 1)

            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedNoSymptoms)] = model.populations[age_group,
                                                                        InfectionState.InfectedSymptoms].value * p_ins * random.uniform(0.1, 1)

            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedSevere)] = model.populations[age_group,
                                                                    InfectionState.InfectedSymptoms].value * p_hosp * random.uniform(0.1, 1)
            model.populations[age_group, Index_InfectionState(
                InfectionState.InfectedCritical)] = model.populations[age_group,
                                                                      InfectionState.InfectedSevere].value * p_icu * random.uniform(0.1, 1)

            model.populations[age_group, Index_InfectionState(
                InfectionState.Dead)] = model.populations[age_group,
                                                          InfectionState.InfectedCritical].value * p_dead * random.uniform(0.1, 1)

            subtotal = (model.populations[age_group, InfectionState.InfectedSymptoms].value
                        + model.populations[age_group,
                                            InfectionState.Exposed].value
                        + model.populations[age_group,
                                            InfectionState.InfectedNoSymptoms].value
                        + model.populations[age_group,
                                            InfectionState.InfectedSevere].value
                        + model.populations[age_group,
                                            InfectionState.InfectedCritical].value
                        + model.populations[age_group,
                                            InfectionState.Dead].value
                        )

            # if subtotal > pop_age_group:
            #    print('Subtotal is larger than population!')
            # print('The remaining population before Recovered is: '+ str(pop_age_group-subtotal))
            # model.populations[age_group, Index_InfectionState(
            #    InfectionState.Recovered)] = random.uniform(0, ((1-np.asarray(randoms).sum())*pop_age_group))
            model.populations[age_group, Index_InfectionState(
                InfectionState.Recovered)] = pop_age_group*p_recovered * random.uniform(0.1, 1)

            # if (subtotal+model.populations[age_group, InfectionState.Recovered].value) >= pop_age_group:
            #    print('Subtotal is larger or equal than population!')
            # print(': ' +
            #      str(model.populations.get_group_total_AgeGroup(age_group)))
            # model.populations[age_group, InfectionState.Susceptible] = pop_age_group - (
            #    subtotal + model.populations[age_group, InfectionState.Recovered].value)
            model.populations.set_difference_from_group_total_AgeGroup((
                age_group, InfectionState.Susceptible), pop_age_group)

            # model.populations.get_group_total_AgeGroup(age_group))

            # print('Susceptible is set to: ' + str(model.populations[age_group,
            #                                                        InfectionState.Susceptible].value))

        # Apply mathematical constraints to parameters
        model.apply_constraints()

        # set model to graph
        graph.get_node(node_indx).property.populations = model.populations

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

    return dataset_entry, runtime


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
    mobility_dir = transform_mobility_directory()
    graph = get_graph(num_groups, data_dir, mobility_dir)

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)

    times = []

    for _ in range(num_runs):

        data_run, t_run = run_secir_groups_simulation(
            days_sum, graph)

        times.append(t_run)

        inputs = np.asarray(data_run).transpose(1, 0, 2)[: input_width]
        data["inputs"].append(inputs)

        data["labels"].append(np.asarray(
            data_run).transpose(1, 0, 2)[input_width:])

        bar.next()

    bar.finish()

    # print median and average runtime
    print(
        f"For Days = {days}, AVG runtime: {np.mean(times)}s, Median runtime: {np.median(times)}s")

    if save_data:

        scaled_inputs, scaled_labels = scale_data(data)

        all_data = {"inputs": scaled_inputs,
                    "labels": scaled_labels,
                    }

        # check if data directory exists. If necessary create it.
        if not os.path.isdir(path):
            os.mkdir(path)

        # save dict to json file
        with open(os.path.join(path, 'GNN_data_30days_nodeswithvariance_10k.pickle'), 'wb') as f:
            pickle.dump(all_data, f)

    return data


if __name__ == "__main__":

    path = os.getcwd()
    # path_output = '/localdata1/gnn_paper_2024/data_Iteration2/GNNs'
    path_output = os.path.join(os.getcwd(), 'saves')
    data_dir = os.path.join(os.getcwd(), 'data')

    input_width = 5
    days = 30
    num_runs = 100
    random.seed(10)
    generate_data(num_runs, data_dir,  path_output, input_width,
                  days, save_data=False)
