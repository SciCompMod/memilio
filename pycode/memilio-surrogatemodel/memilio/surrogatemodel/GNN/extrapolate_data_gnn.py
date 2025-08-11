import os
import numpy as np
import pickle
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import memilio.simulation as mio
import tensorflow as tf

from enum import Enum
from memilio.simulation import (AgeGroup, LogLevel, set_log_level, Damping)
from memilio.simulation.osecir import (
    InfectionState, Model, ModelGraph, set_nodes)
from memilio.surrogatemodel.GNN.GNN_utils import remove_confirmed_compartments


class Location(Enum):
    Home = 0
    School = 1
    Work = 2
    Other = 3


def set_covid_parameters(model, start_date, num_groups=6):
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


def extrapolate_data(start_date, num_days, data_dir):
    model = Model(6)
    set_covid_parameters(model, start_date)
    set_contact_matrices(model, data_dir)

    graph = ModelGraph()

    scaling_factor_infected = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    scaling_factor_icu = 1.0
    tnt_capacity_factor = 7.5 / 100000.

    # Path containing the population data
    data_dir_Germany = os.path.join(data_dir, "Germany")
    pydata_dir = os.path.join(data_dir_Germany, "pydata")

    path_population_data = os.path.join(pydata_dir,
                                        "county_current_population.json")

    set_nodes(
        model.parameters,
        start_date, start_date + num_days,
        pydata_dir,
        path_population_data, True, graph, scaling_factor_infected,
        scaling_factor_icu, tnt_capacity_factor, num_days, True)


def read_results_h5(path, group_key='Total'):
    with h5py.File(path, 'r') as f:
        keys = list(f.keys())
        res = {}
        for i, key in enumerate(keys):
            group = f[key]
            total = group[group_key][()]
            # remove confirmed compartments
            sum_inf_no_symp = np.sum(total[:, [2, 3]], axis=1)
            sum_inf_symp = np.sum(total[:, [4, 5]], axis=1)
            total[:, 2] = sum_inf_no_symp
            total[:, 4] = sum_inf_symp
            total = np.delete(total, [3, 5], axis=1)
            res[key] = total
    return res


def get_ground_truth_data(start_date, num_days, data_dir, create_new=False):
    path_rki_h5 = os.path.join(data_dir, "Germany", "pydata", "Results_rki.h5")
    if not os.path.isfile(path_rki_h5) or create_new:
        print("Generating real data from C++ backend...")
        extrapolate_data(start_date, num_days, data_dir)
        print("Data generation complete.")

    # The GNN expects age-stratified input.
    num_age_groups = 6
    all_age_data_list = []
    group_keys = ['Group1', 'Group2', 'Group3', 'Group4', 'Group5', 'Group6']
    for age in range(num_age_groups):
        age_data_dict = read_results_h5(path_rki_h5, group_keys[age])
        age_data_np = np.array(list(age_data_dict.values()))
        all_age_data_list.append(age_data_np)

    # Combine age groups to get shape (num_nodes, timesteps, num_features=48)
    all_age_data = np.stack(all_age_data_list, axis=-1)
    num_nodes, timesteps, _, _ = all_age_data.shape
    ground_truth_all_nodes_np = np.reshape(
        all_age_data.transpose(0, 1, 3, 2), (num_nodes, timesteps, -1))

    return ground_truth_all_nodes_np


def extrapolate_ground_truth_data(data_dir, num_days=360, num_input_days=0):
    start_dates = [
        mio.Date(2020, 6, 1)
        # mio.Date(2020, 7, 1),
        # mio.Date(2020, 8, 1),
        # mio.Date(2020, 9, 1),
        # mio.Date(2020, 10, 1),
        # mio.Date(2020, 11, 1)
        # mio.Date(2020, 12, 1),
    ]

    for start_date in start_dates:
        ground_truth_np = get_ground_truth_data(
            start_date, num_days, data_dir, create_new=True)
    path_rki_h5 = "/localdata1/hege_mn/memilio/data/Germany/pydata/Results_rki.h5"
    num_age_groups = 6
    all_age_data_list = []
    group_keys = ['Group1', 'Group2', 'Group3', 'Group4', 'Group5', 'Group6']
    for age in range(num_age_groups):
        age_data_dict = read_results_h5(path_rki_h5, group_keys[age])
        age_data_np = np.array(list(age_data_dict.values()))
        all_age_data_list.append(age_data_np)
    # Combine age groups to get shape (num_nodes, timesteps, num_features=48)
    all_age_data = np.stack(all_age_data_list, axis=-1)

    # Save the ground truth data
    output_path = "/localdata1/hege_mn/memilio/data/Germany/pydata/ground_truth_all_nodes.pickle"
    with open(output_path, 'wb') as f:
        pd.to_pickle(all_age_data, f)
    print(f"Ground truth data saved to {output_path}")


def main():
    cwd = os.getcwd()
    start_date = mio.Date(2020, 8, 1)
    num_days = 360
    num_input_days = 5
    number_of_channels = 512
    data_dir = os.path.join(cwd, "data")

    # path_weights = '/localdata1/hege_mn/memilio/pycode/memilio-surrogatemodel/memilio/surrogatemodel/GNN/saved_weights_paper/'
    # path_model_weights = os.path.join(
    #    path_weights, "GNN_30days_nodeswithvariance_1k_3Damp_test2.pickle")

    # compare_and_predict(
    #     start_date,
    #     num_days,
    #     data_dir,
    #     path_model_weights,
    #     number_of_channels,
    #     num_input_days
    # )

    extrapolate_ground_truth_data(
        data_dir,
        num_days,
        num_input_days=0
    )


if __name__ == "__main__":
    main()
