import numpy as np
import pandas as pd
import os
from memilio.epidata import modifyDataframeSeries as mdfs


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


def getMinimumMatrix():
    """! loads the minimum matrix
    """

    minimum_contact_matrix0 = os.path.join(
        "./data/contacts/minimum_home.txt")
    minimum_contact_matrix1 = os.path.join(
        "./data/contacts/minimum_school_pf_eig.txt")
    minimum_contact_matrix2 = os.path.join(
        "./data/contacts/minimum_work.txt")
    minimum_contact_matrix3 = os.path.join(
        "./data/contacts/minimum_other.txt")

    minimum = np.loadtxt(minimum_contact_matrix0) \
        + np.loadtxt(minimum_contact_matrix1) + \
        np.loadtxt(minimum_contact_matrix2) + \
        np.loadtxt(minimum_contact_matrix3)

    return minimum


def get_population():
    df_population = pd.read_json(
        'data/pydata/Germany/county_population.json')
    age_groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80-130']

    df_population_agegroups = pd.DataFrame(
        columns=[df_population.columns[0]] + age_groups)
    for region_id in df_population.iloc[:, 0]:
        df_population_agegroups.loc[len(df_population_agegroups.index), :] = [int(region_id)] + list(
            mdfs.fit_age_group_intervals(df_population[df_population.iloc[:, 0] == int(region_id)].iloc[:, 2:], age_groups))

    population = df_population_agegroups.values.tolist()

    return population
