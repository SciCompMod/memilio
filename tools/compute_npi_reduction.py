#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn, David Kerkmann
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
"""
@file cumpute_npi_reduction.py

@brief Computes the overall reduction in contacts by enforcing a NPI.
"""

import sys
import numpy as np
import os
import pandas as pd
from memilio.epidata import defaultDict as dd


def find_county(county):
    """"! Attempts to find a given county by name or ID and returns both if successful.

    @param county Name or ID of county.
    @return Name and ID of county.
    """

    if county == 'Germany':
        return county, 0
    else:
        # try to find county ID
        county_id = dd.invert_dict(dd.County).get(county, False)
        # if county was not found, try to check if ID was passed instead
        if not county_id:
            county_id = county
            county = dd.County.get(county_id)
            if county:
                print('Found county ' + county + ' with input ID ' + str(county_id) + '.')
            else:
                raise ValueError('County ' + str(county_id) + ' was not found in County list. Aborting.')
    return county, county_id


def get_data_path():
    """"! Returns data path of memilio repository. Currently only works if repository is called 'memilio'.

    @return Data path of memilio repository.
    """

    data_path = os.getcwd().split('memilio')[0] + 'memilio/data/'
    #print("My current directory is : " + data_path)
    return data_path


def compute_population(county_id, filename):
    """"! Computes current population across age groups from given county.

    @param county_id ID of county.
    @param filename Name of file which has the age-resolved data.
    @return Array containing the population for each age group in the given county.
    """

    # TODO: We currently don't have age-resolved data for each county.
    #       Once we have this data, find correct file and adapt extraction to that file format.

    data_path = get_data_path()

    df = pd.read_json(data_path + '/Germany/' + filename)
    
    if county_id == 0:
        population = np.array([
        len(df[df['Age_RKI'] == 'A00-A04']),
        len(df[df['Age_RKI'] == 'A05-A14']),
        len(df[df['Age_RKI'] == 'A15-A34']),
        len(df[df['Age_RKI'] == 'A35-A59']),
        len(df[df['Age_RKI'] == 'A60-A79']),
        len(df[df['Age_RKI'] == 'A80+'])  ])

        # provide population of age groups as used in contact matrices
        # population = np.array([3961376,7429883,19117865,28919134,18057318,5681135])
        
    else:
        df_county = df[df['ID_County'] == county_id]
        population = np.array([
        len(df_county[df_county['Age_RKI'] == 'A00-A04']),
        len(df_county[df_county['Age_RKI'] == 'A05-A14']),
        len(df_county[df_county['Age_RKI'] == 'A15-A34']),
        len(df_county[df_county['Age_RKI'] == 'A35-A59']),
        len(df_county[df_county['Age_RKI'] == 'A60-A79']),
        len(df_county[df_county['Age_RKI'] == 'A80+'])  ])
    
    #print(df_county[(df_county['Age_RKI'] != 'A05-A14') & (df_county['Age_RKI'] != 'A15-A34') & (df_county['Age_RKI'] != 'A35-A59') & (df_county['Age_RKI'] != 'A60-A79') & (df_county['Age_RKI'] != 'A80+') & (df_county['Age_RKI'] != 'A00-A04')].head())
    #print(sum(population))
    # Some ages are unknown in this specific list. This will hopefully be fixed when the correct list is used.
    
    return population


def compute_base_contacts(population, use_minimum):
    """"! Computes base contacts for a given population.

    @param population Array containing the population for each age group.
    @param use_minimum Boolean if the minimum contact matrix should be used. If false, minimum contacts are set to 0.
    @return Arrays containing the total amount of contacts for each age group at baseline and at minimum contact state.
    """

    data_path = get_data_path()

    # names for contact location data files
    home = 'home'
    school = 'school_pf_eig'
    work = 'work'
    other = 'other'

    # list of names
    contact_locs = [home, school, work, other]

    # read contact matrices
    # assuming baseline contact matrices are constant across the whole population
    cm_bl = np.zeros((4,6,6))
    cm_bl_min = np.zeros((4,6,6))

    for k, loc in enumerate(contact_locs):
        cm_bl[k] = np.loadtxt(data_path + '/contacts/baseline_' + loc + '.txt')
        if use_minimum:
            cm_bl_min[k] = np.loadtxt(data_path + '/contacts/minimum_' + loc + '.txt')
        # else: already initialized to zeros

    # multiply matrices from left and sum up in each category
    contacts_total_base = np.sum(population@cm_bl,1)
    contacts_total_min = np.sum(population@cm_bl_min,1)

    return contacts_total_base, contacts_total_min


def apply_NPI(NPI, contacts_total_base, contacts_total_min):
    """"! Applies given NPI and computes mean reduction of contacts for a given number of contacts.

    @param NPI NPI that represents minimum and maximum contact reductions. For now, this is not used but hardcoded.
    @param contacts_total_base Array containing the total amount of contacts for each age group at baseline.
    @param contacts_total_min Array containing the total amount of contacts for each age group at minimum contact state. 
    @return Factor of mean reduction in contacts between 0 and 1. 1: 100% reduction, 0: 0% reduction.
    """

    # TODO: think about how to read in NPI. For now, it is hardcoded.

    # define intervals of contact reduction for different levels
    reduc_factors_min = np.array([[0.5, 0.3, 0.6, 0.6], [0.0, 0.25, 0.25, 0.25]])
    reduc_factors_max = np.array([[0.7, 0.5, 0.8, 0.8], [0.0, 0.35, 0.35, 0.35]])
    reduc_factors_mean = (reduc_factors_min + reduc_factors_max)/2

    # compute total reduction factor for each level
    # these are the factors for baseline; the factors for minimum contacts are [1 - this value]
    reduc_factors_mean_total = np.prod(1 - reduc_factors_mean, 0)

    # compute total mean reduction in contacts (that is 1/2*min + 1/2*max)
    contacts_total_reduced_mean = reduc_factors_mean_total*contacts_total_base + (1-reduc_factors_mean_total)*contacts_total_min

    reduc_total = 1 - np.sum(contacts_total_reduced_mean) / np.sum(contacts_total_base)
    return reduc_total


def compute_npi_reduction(NPI = 0, county = 'Germany', filename = 'cases_all_county_age.json', use_minimum = False):
    """"! Main function that computes the reduction from a NPI in a given county and prints the result.

    @param NPI NPI that represents minimum and maximum contact reductions. For now, this is not used but hardcoded.
    @param county [Default: Germany] County in which the NPI is applied. 
    @param filename [Default: 'cases_all_county_age.json'] Name of file which has the age-resolved data.
    @use_minimum [Default: False] Boolean if the minimum contact matrix should be used. If false, minimum contacts are set to 0.
    @return Factor of mean reduction in contacts between 0 and 1. 1: 100% reduction, 0: 0% reduction.
    """

    # TODO: Find correct file, adapt 'compute_population' to it.
    county, county_id = find_county(county)
    
    print('NPI reduction for citizens in ' + county + ':')
    population = compute_population(county_id, filename)

    contacts_total_base, contacts_total_min = compute_base_contacts(population, use_minimum)

    # print out basic info
    print('Number of baseline in home: '
          + str(np.round(100*contacts_total_base[0] / sum(contacts_total_base), 2)) + '%, school: '
          + str(np.round(100*contacts_total_base[1] / sum(contacts_total_base), 2)) + '%, work: '
          + str(np.round(100*contacts_total_base[2] / sum(contacts_total_base), 2)) + '%, other: '
          + str(np.round(100*contacts_total_base[3] / sum(contacts_total_base), 2)) + '%')

    reduc_total = apply_NPI(NPI, contacts_total_base, contacts_total_min)

    print('Total contact reduction (including protection effects): '
        + str(np.round(100*reduc_total,2)) + '%')

    return reduc_total

    # OLD CODE (kept to have fast access to other computations (min/max) and outputs):
    # make three copies and compute for min, mean and max
    #contacts_total_reduced = [contacts_total_base.copy(
    #), contacts_total_base.copy(), contacts_total_base.copy()]
    #if reduc_factors_min.shape == reduc_factors_max.shape:
    #    for i in range(len(contacts_total_base)):
    #        reduc = [1, 1, 1]
    #        for j in range(reduc_factors_min.shape[0]):
    #            # compute combined (all levels) reduction factors
    #            reduc[0] *= (1 - reduc_factors_min[j][i]) # min
    #            reduc[1] *= (1 - 0.5*(reduc_factors_min[j][i] + reduc_factors_max[j][i])) # mean
    #            reduc[2] *= (1 - reduc_factors_max[j][i]) # max
    #        for k in range(len(reduc)):
    #            contacts_total_reduced[k][i] = contacts_total_reduced[k][i] - (1 - reduc[k]) * (contacts_total_reduced[k][i] - contacts_total_min[i])
    #
    #    print(contacts_total_reduced)
    #
    #    print('Summed contact reduction (including protection effects): ' + \
    #          str(np.round(100*(1 - sum(contacts_total_reduced[1]) / sum(contacts_total_base)), 2))
    #          + '% (' + str(np.round(100*(1 - sum(contacts_total_reduced[0]) / sum(contacts_total_base)), 2)) 
    #          + ' - ' + str(np.round(100*(1 - sum(contacts_total_reduced[2]) / sum(contacts_total_base)), 2)) + '%)')
    #
    #    print('Average number of contact reduction in home: '
    #        + str(np.round(100*(1 - contacts_total_reduced[1][0] / contacts_total_base[0]), 2)) + '%, school: '
    #        + str(np.round(100*(1 - contacts_total_reduced[1][1] / contacts_total_base[1]), 2)) + '%, work: '
    #        + str(np.round(100*(1 - contacts_total_reduced[1][2] / contacts_total_base[2]), 2)) + '%, other: '
    #        + str(np.round(100*(1 - contacts_total_reduced[1][3] / contacts_total_base[3]), 2)) + '%')
    #
    #    print('Number of unprotected contacts after reduction in home: '
    #        + str(np.round(100*contacts_total_reduced[1][0] / sum(contacts_total_reduced), 2)) + '%, school: '
    #        + str(np.round(100*contacts_total_reduced[1][1] / sum(contacts_total_reduced), 2)) + '%, work: '
    #        + str(np.round(100*contacts_total_reduced[1][2] / sum(contacts_total_reduced), 2)) + '%, other: '
    #        + str(np.round(100*contacts_total_reduced[1][3] / sum(contacts_total_reduced), 2)) + '%')
    #else:
    #    print('Error in reduction factor dimension. Please provide equally '
    #        'sized arrays for minimum and maximum.')

if __name__ == "__main__":
    compute_npi_reduction(*sys.argv[1:])
