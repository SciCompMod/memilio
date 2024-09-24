import copy
import os
import pickle
import random
import numpy as np
from datetime import date

from progress.bar import Bar
from sklearn.preprocessing import FunctionTransformer

from memilio.simulation import (AgeGroup, LogLevel, set_log_level, Damping)
from memilio.simulation.osecir import (Index_InfectionState, interpolate_simulation_result, ParameterStudy,
                                       InfectionState, Model,
                                       interpolate_simulation_result)
from memilio.epidata import geoModificationGermany as geoger

from memilio.surrogatemodel.GNN.GNN_utils import (transform_mobility_directory,
                                                  make_graph, scale_data)
from memilio.surrogatemodel.utils_surrogatemodel import (
    getBaselineMatrix, remove_confirmed_compartments, get_population, getMinimumMatrix)


def run_secir_groups_simulation(days, populations, dampings, countykey_list, directory, baseline):
    """! Uses an ODE SECIR model allowing for asymptomatic infection 
        with 6 different age groups. The model is not stratified by region. 
        Virus-specific parameters are fixed and initial number of persons 
        in the particular infection states are chosen randomly from defined ranges.
    @param days (int) Describes how many days we simulate within a single run.
    @param damping_day (int) The day when damping is applied.
    @param populations List containing the population in each age group.
    @return List containing the populations in each compartment 
            used to initialize the run.
   """
    set_log_level(LogLevel.Off)

    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1

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
                (AgeGroup(i), Index_InfectionState(InfectionState.Susceptible)),  populations[region][i])

            # Compartment transition propabilities
            model.parameters.RelativeTransmissionNoSymptoms[AgeGroup(i)] = 0.5
            model.parameters.TransmissionProbabilityOnContact[AgeGroup(
                i)] = 0.1
            model.parameters.RecoveredPerInfectedNoSymptoms[AgeGroup(i)] = 0.09
            model.parameters.RiskOfInfectionFromSymptomatic[AgeGroup(i)] = 0.25
            model.parameters.SeverePerInfectedSymptoms[AgeGroup(i)] = 0.2
            model.parameters.CriticalPerSevere[AgeGroup(i)] = 0.25
            model.parameters.DeathsPerCritical[AgeGroup(i)] = 0.3
            # twice the value of RiskOfInfectionFromSymptomatic
            model.parameters.MaxRiskOfInfectionFromSymptomatic[AgeGroup(
                i)] = 0.5

        # StartDay is the n-th day of the year
        model.parameters.StartDay = (
            date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

        model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline
        model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.ones(
            (num_groups, num_groups)) * 0

        # Generate a damping matrix and assign it to the model
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
        models.append(model)

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
        num_runs, path, input_width, label_width, number_of_dampings, save_data=True):
    """! Generate dataset by calling run_secir_simulation (num_runs)-often
   @param num_runs Number of times, the function run_secir_simulation is called.
   @param path Path, where the datasets are stored.
   @param input_width Number of time steps used for model input.
   @param label_width Number of time steps (days) used as model output/label.  
   @param number_of_dampings (int) The number of contact change points applied to the simulation. 
   @param save_data Option to deactivate the save of the dataset. Per default true.
   """
    population = get_population()
    days_sum = label_width+input_width-1

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

    # get county ids
    countykey_list = geoger.get_county_ids(merge_eisenach=True, zfill=True)
    directory = transform_mobility_directory()

    # Load baseline matrix
    baseline = getBaselineMatrix()

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)

    for i in range(num_runs):

        data_run, damped_contact_matrix, damping_days_s, damping_factor = run_secir_groups_simulation(
            days_sum,  population,  damping_days[i], countykey_list, directory, baseline)

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
    label_width = 30
    number_of_dampings = 1
    num_runs = 2
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(
        os.path.dirname(
            os.path.realpath(os.path.dirname(os.path.realpath(path)))),
        'data_GNN_with_'+str(number_of_dampings)+'_dampings_test')

    generate_data(num_runs, path_data, input_width,
                  label_width, number_of_dampings, save_data=False)
