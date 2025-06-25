#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Manuel Heger, Henrik Zunker
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

import random
import numpy as np


def generate_dampings(days, max_number_damping, method, min_distance=2,
                      min_damping_day=2):
    """
    Producing dampings for timeseries of length days according to the used method. 

    :param days: Number of days per time series 
    :param max_number_damping: maximal possible number of active dampings, actual number of active dampings can be smaller. 
    :param method: Method used to generate the dampings, possible values "classic", "active", "random"
    :param min_distance: Minimal distance between two dampings 
    :min_damping_day: First day, where a damping can be applied 
    :returns: two lists containing the damping days and the associated damping factors
    """
    if method == "classic":
        damp_days, damp_factors = dampings_classic(
            days, max_number_damping, min_distance, min_damping_day)
    elif method == "active":
        damp_days, damp_factors = dampings_active(days, max_number_damping)
    elif method == "random":
        damp_days, damp_factors = dampings_random(
            days, max_number_damping, min_distance, min_damping_day)
    else:
        raise ValueError(
            "The method argument has to be one of the following: 'classic', 'active' or 'random'.")

    return damp_days, damp_factors


# Active Damping
def dampings_active(days, max_number_damping):
    """"
    Generating list of damping days and corresponding damping factors using the active method. 

    The damping days are created with equal distance on the interval [1, days-3].

    :param days: Number of simulated days 
    :param max_number_damping: Maximal number of dampings per simulation
    :returns: list of days and damping factors
    """
    if max_number_damping > days-3:
        raise ValueError(
            "Number of dampings must be smaller than total number of days!"
        )

    # Setting parameters
    gamma_pos = 0
    alpha = -1
    p0 = 0.5
    t1_max = -0.3
    t1_min = -1.2
    t2_max = 2
    t2_min = -0.5

    # Defining possible damping days
    distance_between_days = np.floor((days-3)/max_number_damping)
    damp_days = [distance_between_days*(i+1)
                 for i in np.arange(max_number_damping)]

    # Generate damping factors
    dampings = calc_factors_active(
        max_number_damping, gamma_pos, alpha, p0, t1_max, t1_min, t2_max, t2_min)

    return damp_days, dampings


def calc_factors_active(n_ddays, gamma_pos=0, alpha=-1, p0=0.5, t1_max=-0.3, t1_min=-1.2, t2_max=2, t2_min=-0.5):
    '''
    Producing damping factors using active damping method. 

    The idea is the following: Damping factors are produced randomly until a threshold value is achieved. 
    In this case the factors are reduced stepwise until a moderate level is reached. 

    :param n_ddays: Number of damping days in time series
    :param gamma_pos: upper bound for h value
    :param alpha: upper bound for h value, where active damping should be stopped
    :param p0: probability for a change of the damping factor
    :param t1_max: upper end point for size of active damping changes
    :param t1_min: lower end point for size of active damping changes
    :param t2_max: upper end point for size of active damping changes
    :param t2_min: lower end point for size of base damping changes
    :param t3_max: upper end
    '''
    h = 0
    k = 0
    dampings = []

    while k < n_ddays:
        # If the threshold value is reached, active damping is started
        if h > gamma_pos:
            # active reducing the damping factor
            while h > alpha and k < n_ddays:
                delta_h = np.random.uniform(t1_min, t1_max)
                h = h + delta_h
                dampings.append(1-np.exp(h))
                k = k+1
        # otherwise changes of the damping factor are generated randomly
        else:
            # Whether or not a non-trivial change of the damping factor is applied
            if np.random.binomial(1, p0):
                delta_h = np.random.uniform(t2_min, t2_max)

            else:
                delta_h = 0
            h = h+delta_h
            dampings.append(1 - np.exp(h))
            k = k+1

    return dampings


# Classic Damping
def dampings_classic(days, max_number_damping,  min_distance=2,
                     min_damping_day=2):
    """
    Generate the damping days using shadow damping and picking days uniformly with a given minimal distance. 

    The corresponding factors are drawn uniformly from the interval (0,0.5)

    :param days: Number of days simulated per run. 
    :param max_number_damping: Number of damping days generated. 
    :param min_distance: Minimal distance between two dampings 
    :param min_damping_day: First day, where a damping can be applied 
    :returns: Two lists of length max_number_dampings containing the days and the factors.
    """
    # Generating damping days
    if min_distance*max_number_damping+min_damping_day > days:
        raise ValueError("Invalid input: It's not possible to generate this number of damping"
                         "in the desired time interval.")
    damp_days = generate_dampings_withshadowdamp(
        max_number_damping, days, min_distance, min_damping_day, 1)[0]

    # Generating damping factors
    damp_factors = np.random.uniform(0, 0.5, max_number_damping).tolist()

    return damp_days, damp_factors


def generate_dampings_withshadowdamp(number_of_dampings, days, min_distance, min_damping_day, n_runs):
    """
    Sampling the damping days according to the established method. 

    The idea is to draw dampings with a minimum distance while traying to keep
    the distribution of damping days uniformly. We create a list of all possible days,
    draw one damping day and delete all days before and after the damping that
    are within the range of the min_distance. To ensure that the the data is not biased,
    we include days outside the usual range. A day x in the middle of the list can
    be removed from the list by a drawn day before and after x. A day in the beggining
    of the list can be removed only by drawn days y , y>x. This leads to the effect that
    the first and last days are chosen more often. By drawing days ouside of the allowed range
    (forbidden dampings) which are removed after, we ensure that also the days at the beginning and
    end of the list can be removed from the list because of the minimum distance.

    :param number_of_dampings: Number of damping days per run 
    :param days: Total number of days per run 
    :param min_distance: Minimal distance between two damping days
    :param min_damping_day: First day when a damping can be applied 
    :n_runs: Number of runs for which damping days should be generated 
    :returns: list of list of damping days. 
    """

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


# Random Damping

def dampings_random(days, max_number_damping, min_damping_day=2,
                    min_distance_damping_day=2):
    """
    Generate dampings according to an easy random rule. 

    The days are drawn using geometrical distributed waiting times and fixed minimal distance between two damping days. 
    The factors are drwan uniformly on the interval (0,0.5)


    :param days: Number of days simulated per run. 
    :param max_number_damping: Number of damping days generated. 
    :param min_distance: Minimal distance between two dampings 
    :param min_damping_day: First day, where a damping can be applied 
    :returns: Two lists of length max_number_dampings containing the days and the factors.
    """
    # Setting parameters

    # Calculating the expected distance between two dampings
    distance_between_days = np.floor((days-min_damping_day)/max_number_damping)
    # Reducing due to minimal distance restriction
    reduced_distance = distance_between_days - min_distance_damping_day

    if reduced_distance <= 0:
        raise ValueError("Invalid input: It's not possible to generate this number of damping"
                         "in the desired time interval.")

    # Try till one admissible configuration of waiting times is produced
    running = True
    while running:
        dist = np.random.geometric(1/reduced_distance, max_number_damping)
        if np.sum(dist) + min_damping_day + max_number_damping*min_distance_damping_day < days:
            running = False

    # Reconstructing the days using the waiting times
    ddays = []
    day = min_damping_day
    for k in np.arange(len(dist)):
        day = day + dist[k]
        ddays.append(day)
        day = day + min_distance_damping_day

    # Generating the associated damping factors
    damping_factors = np.random.uniform(0, 0.5, max_number_damping)

    return ddays, damping_factors
