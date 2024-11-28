import numpy as np
import pandas as pd
import os
import json
from memilio.epidata import modifyDataframeSeries as mdfs
import random


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


def generate_dampings_withshadowdamp(number_of_dampings, days, min_distance, min_damping_day, n_runs):

    # the idea is to draw dampings with a minimum distance while traying to keep
    # the distribution of damping days uniformly. We create a list of all possible days,
    # draw one damping day and delete all days before and after the damping that
    # are within the range of the min_distance. To ensure that the the data is not biased,
    # we include days outside the usual range. A day x in the middle of the list can
    # be removed from the list by a drawn day before and after x. A day in the beggining
    # of the list can be removed only by drawn days y , y>x. This leads to the effect that
    # the first and last days are chosen more often. By drawing days ouside of the allowed range
    # (forbidden dampings) which are removed after, we ensure that also the days atthe beginning and
    # end of the list can be removed from the list because of the minimum distance.
    number_of_dampings = number_of_dampings
    days = days
    min_distance = min_distance
    min_damping_day = min_damping_day
    number_of_runs = n_runs

    all_dampings = []
    count_runs = 0
    count_shadow = 0
    while len(all_dampings) < number_of_runs:
        # Reset the days list and dampings for each run
        days_list = list(range(min_damping_day, days))
        dampings = []

        if count_shadow < 2:
            for _ in range(number_of_dampings):
                if len(days_list) > 0:
                    damp = random.choice(days_list)
                    days_before = list(range(damp - min_distance, damp))
                    days_after = list(range(damp, damp + min_distance + 1))
                    dampings.append(damp)
                    days_list = [ele for ele in days_list if ele not in (
                        days_before + days_after)]
                else:
                    # Restart the process when days_list is empty
                    break
            else:
                # Exit loop only if dampings were successfully drawn
                forbidden_damping_values = list(
                    range(0 - min_distance, 0)) + list(range(days + 1, days + min_distance + 1))
                dampings = [
                    ele for ele in dampings if ele not in forbidden_damping_values]
                if len(dampings) >= number_of_dampings:
                    all_dampings.append(sorted(dampings))
                continue
        else:
            # Generate forbidden damping
            damp = random.choice(
                list(range(0 - min_distance, 0)) +
                list(range(days + 1, days + min_distance + 1))
            )
            dampings.append(damp)
            for _ in range(number_of_dampings):
                if len(days_list) > 0:
                    damp = random.choice(days_list)
                    days_before = list(range(damp - min_distance, damp))
                    days_after = list(range(damp, damp + min_distance + 1))
                    dampings.append(damp)
                    days_list = [ele for ele in days_list if ele not in (
                        days_before + days_after)]
                else:
                    # Restart the process when days_list is empty
                    break
            else:
                # Reset shadow count only if dampings were successfully drawn
                count_shadow = 0
                forbidden_damping_values = list(
                    range(0 - min_distance, 0)) + list(range(days + 1, days + min_distance + 1))
                dampings = [
                    ele for ele in dampings if ele not in forbidden_damping_values]
                if len(dampings) >= number_of_dampings:
                    all_dampings.append(sorted(dampings))
                continue

        # Restart process if any issue occurred
        count_runs += 1
        count_shadow += 1

    return all_dampings
