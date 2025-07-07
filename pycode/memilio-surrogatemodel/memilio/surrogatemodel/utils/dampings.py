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


def calc_dist_days(days, min_day, n_dampings, min_distance=1):
    """
    Calculating distance between two dampings if there are n_dampings on the interval
    (min_day, days)

    :param days: Total number of days 
    :min_day: First day on which a damping can be applied 
    :n_dampings: Number of dampings
    :min_distance: Lower bound for the distance between two dampings
    """
    res = np.floor((days-min_day)/n_dampings)

    if res < min_distance:
        raise ValueError(
            "It's not possible to arrange this number of dampings in the desired interval with the given minimal distance.")

    return res


def generate_dampings(days, number_dampings, method, min_distance=2,
                      min_damping_day=2):
    """
    Producing dampings for timeseries of length days according to the used method. 

    :param days: Number of days per time series 
    :param number_dampings: Number of days on which damping can occur.
    :param method: Method used to generate the dampings, possible values "classic", "active", "random"
    :param min_distance: Minimal distance between two dampings 
    :min_damping_day: First day, where a damping can be applied 
    :returns: two lists containing the damping days and the associated damping factors
    """
    if method == "classic":
        damp_days, damp_factors = dampings_classic(
            days, number_dampings, min_distance, min_damping_day)
    elif method == "active":
        damp_days, damp_factors = dampings_active(
            days, number_dampings, min_damping_day)
    elif method == "random":
        damp_days, damp_factors = dampings_random(
            days, number_dampings, min_distance, min_damping_day)
    else:
        raise ValueError(
            "The method argument has to be one of the following: 'classic', 'active' or 'random'.")

    return damp_days, damp_factors


# Active Damping
def dampings_active(days, number_dampings, min_damping_day):
    """"
    Generating list of damping days and corresponding damping factors using the active method. 

    The damping days are created with equal distance on the interval [1, days-3].

    :param days: Number of simulated days 
    :param number_dampings: Maximal number of dampings per simulation
    :param min_damping_day: First day, where a damping can be applied 
    :returns: list of days and damping factors
    """

    # Setting parameters
    gamma_pos = -2
    alpha = -4
    p0 = 0.4
    t1_max = -1
    t1_min = -2.5
    t2_max = 0.95
    t2_min = -0.25

    # Defining possible damping days
    distance_between_days = calc_dist_days(
        days, min_damping_day, number_dampings)
    damp_days = [min_damping_day + distance_between_days*(i+1)
                 for i in np.arange(number_dampings)]

    # Generate damping factors
    dampings = calc_factors_active(
        number_dampings, gamma_pos, alpha, p0, t1_max, t1_min, t2_max, t2_min)

    return damp_days, dampings


def calc_factors_active(n_ddays, gamma_pos=0, alpha=-1, p0=0.5, t1_max=-0.3, t1_min=-1.2, t2_max=2, t2_min=-0.5):
    '''
    Producing damping factors using active damping method. 

    The idea is the following: Damping factors are produced randomly until a threshold value is achieved. 
    In this case the factors are reduced stepwise until a moderate level is reached. 
    The new contact matrix is always calculated with resepect to the initial matrix according 
    to the following rule: 

             M = exp(h)*M_0

    :param n_ddays: Number of damping days in time series
    :param gamma_pos: upper bound for h value
    :param alpha: upper bound for h value, where active damping should be stopped
    :param p0: probability for a change of the damping factor
    :param t1_max: upper end point for size of active damping changes
    :param t1_min: lower end point for size of active damping changes
    :param t2_max: upper end point for size of active damping changes
    :param t2_min: lower end point for size of base damping changes
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
def dampings_classic(days, number_dampings,  min_distance=2,
                     min_damping_day=2):
    """
    Generate the damping days using shadow damping and picking days uniformly with a given minimal distance. 

    The idea behind shadow damping is the following: Days are picked randomly in the interval between min_damping_day and days. 
    After picking one day a fixed number of days before and after the chosen day is blocked. 
    The procedure is repeated till number_damping many days are chosen. To overcome the problem of higher probability at the boundary, 
    the interval is artificially increased, artificial days are not counted. 
    The corresponding factors are drawn uniformly from the interval (0,0.5)

    :param days: Number of days simulated per run. 
    :param number_dampings: Number of damping days generated. 
    :param min_distance: Minimal distance between two dampings 
    :param min_damping_day: First day, where a damping can be applied 
    :returns: Two lists of length number_dampingss containing the days and the factors.
    """
    # Checking, if the given parameters are compatible
    calc_dist_days(days, min_damping_day, number_dampings, min_distance)

    # Generating damping days
    damp_days = generate_dampings_withshadowdamp(
        number_dampings, days, min_distance, min_damping_day, 1)[0]

    # Generating damping factors
    damp_factors = np.random.uniform(0, 0.5, number_dampings).tolist()

    return damp_days, damp_factors


def generate_dampings_withshadowdamp(number_of_dampings, days, min_distance, min_damping_day, number_of_runs):
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
    :number_of_runs: Number of runs for which damping days should be generated 
    :returns: list of list of damping days. 
    """

    all_dampings = []
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
                if len(dampings) == number_of_dampings:
                    all_dampings.append(sorted(dampings))
                continue

        # Restart process if any issue occurred
        count_shadow += 1

    return all_dampings


# Random Damping

def dampings_random(days, number_dampings, min_damping_day=2,
                    min_distance_damping_day=2):
    """
    Generate random damping days according to the following rule. 

    The days are drawn using geometrical distributed waiting times and a fixed minimal distance betweem two 
    damping days. The first damping can occure at min_damping_day. The associated damping factors are drawn uniformly 
    between 0 and 0.5.

    :param days: Number of days simulated per run. 
    :param number_dampings: Number of damping days generated. 
    :param min_distance: Minimal distance between two dampings 
    :param min_damping_day: First day, where a damping can be applied 
    :returns: Two lists of length number_dampingss containing the days and the factors.
    """
    # Setting parameters

    # Calculating the expected distance between two dampings
    distance_between_days = calc_dist_days(
        days, min_damping_day, number_dampings, min_distance_damping_day)
    # Reducing due to minimal distance restriction
    reduced_distance = distance_between_days - min_distance_damping_day

    # Try till one admissible configuration of waiting times is produced
    running = True
    while running:
        dist = np.random.geometric(1/reduced_distance, number_dampings)
        if np.sum(dist) + min_damping_day + number_dampings*min_distance_damping_day < days:
            running = False

    # Reconstructing the days using the waiting times
    ddays = []
    day = min_damping_day
    for k in np.arange(len(dist)):
        day = day + dist[k]
        ddays.append(day)
        day = day + min_distance_damping_day

    # Generating the associated damping factors
    damping_factors = np.random.uniform(0, 0.5, number_dampings)

    return ddays, damping_factors
