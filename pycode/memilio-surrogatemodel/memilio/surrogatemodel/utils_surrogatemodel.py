import numpy as np
import pandas as pd
import os
import json
from memilio.epidata import modifyDataframeSeries as mdfs


def remove_confirmed_compartments(result_array):
    """! Removes the confirmed compartments which are not used in the data generation.
    @param result_array Array containing the simulation results.
    @return Array containing the simulation results without the confirmed compartments.
    """
    num_groups = int(result_array.shape[1] / 10)
    delete_indices = [index for i in range(
        num_groups) for index in (3+10*i, 5+10*i)]
    return np.delete(result_array, delete_indices, axis=1)

def remove_confirmed_compartments_groups(dataset_entries, num_groups):
    """! The compartments which contain confirmed cases are not needed and are 
        therefore omitted by summarizing the confirmed compartment with the 
        original compartment. 
    @param dataset_entries Array that contains the compartmental data with 
            confirmed compartments. 
    @param num_groups Number of age groups.
    @return Array that contains the compartmental data without confirmed compartments. 
   """

    new_dataset_entries = []
    for i in dataset_entries:
        dataset_entries_reshaped = i.reshape(
            [num_groups, int(np.asarray(dataset_entries).shape[1]/num_groups)]
        )
        sum_inf_no_symp = np.sum(dataset_entries_reshaped[:, [2, 3]], axis=1)
        sum_inf_symp = np.sum(dataset_entries_reshaped[:, [4, 5]], axis=1)
        dataset_entries_reshaped[:, 2] = sum_inf_no_symp
        dataset_entries_reshaped[:, 4] = sum_inf_symp
        new_dataset_entries.append(
            np.delete(dataset_entries_reshaped, [3, 5], axis=1).flatten()
        )
    return new_dataset_entries


def getBaselineMatrix():
    """! loads the baselinematrix
    """

    baseline_contact_matrix0 = os.path.join(
        "./memilio/data/contacts/baseline_home.txt")
    baseline_contact_matrix1 = os.path.join(
        "./memilio/data/contacts/baseline_school_pf_eig.txt")
    baseline_contact_matrix2 = os.path.join(
        "./memilio/data/contacts/baseline_work.txt")
    baseline_contact_matrix3 = os.path.join(
        "./memilio/data/contacts/baseline_other.txt")

    baseline = np.loadtxt(baseline_contact_matrix0) \
        + np.loadtxt(baseline_contact_matrix1) + \
        np.loadtxt(baseline_contact_matrix2) + \
        np.loadtxt(baseline_contact_matrix3)

    return baseline

# def get_population():
#     df_population = pd.read_json(
#         'data/pydata/Germany/county_current_population.json')
#     age_groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80-130']

#     df_population_agegroups = pd.DataFrame(
#         columns=[df_population.columns[0]] + age_groups)
#     for region_id in df_population.iloc[:, 0]:
#         df_population_agegroups.loc[len(df_population_agegroups.index), :] = [int(region_id)] + list(
#             mdfs.fit_age_group_intervals(df_population[df_population.iloc[:, 0] == int(region_id)].iloc[:, 2:], age_groups))

#     population = df_population_agegroups.values.tolist()

#     return population

def interpolate_age_groups(data_entry):
    """! Interpolates the age groups from the population data into the age groups used in the simulation. 
    We assume that the people in the age groups are uniformly distributed.
    @param data_entry Data entry containing the population data.
    @return List containing the population in each age group used in the simulation.
    """
    age_groups = {
        "A00-A04": data_entry['<3 years'] + data_entry['3-5 years'] * 2 / 3,
        "A05-A14": data_entry['3-5 years'] * 1 / 3 + data_entry['6-14 years'],
        "A15-A34": data_entry['15-17 years'] + data_entry['18-24 years'] + data_entry['25-29 years'] + data_entry['30-39 years'] * 1 / 2,
        "A35-A59": data_entry['30-39 years'] * 1 / 2 + data_entry['40-49 years'] + data_entry['50-64 years'] * 2 / 3,
        "A60-A79": data_entry['50-64 years'] * 1 / 3 + data_entry['65-74 years'] + data_entry['>74 years'] * 1 / 5,
        "A80+": data_entry['>74 years'] * 4 / 5
    }
    return [age_groups[key] for key in age_groups]


def get_population(path):
    """! read population data in list from dataset
    @param path Path to the dataset containing the population data
    """

    with open(path) as f:
        data = json.load(f)
    population = []
    for data_entry in data:
        population.append(interpolate_age_groups(data_entry))
    return population

