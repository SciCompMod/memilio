#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn, Wadim Koslow, Daniel Abele
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
import argparse
import pandas as pd
import numpy as np

from memilio.simulation import abm
from memilio.simulation.abm import AgeGroup
from memilio.simulation.abm import VaccinationState
from memilio.simulation.abm import VirusVariant
from memilio.simulation.abm import History
from memilio.simulation.abm import Infection

# class used to map input areas to abm locations


class LocationMapping:

    def __init__(self):
        self.inputId = None
        self.modelId = []

# infection parameters are dependent on age group and vaccination state


def set_infection_parameters():
    infection_params = abm.GlobalInfectionParameters()

    # AgeGroup 0-4
    infection_params.IncubationPeriod[VirusVariant.Wildtype, AgeGroup.Age0to4,
                                      VaccinationState.Unvaccinated] = 3
    infection_params.InfectedNoSymptomsToSymptoms[VirusVariant.Wildtype, AgeGroup.Age0to4,
                                                  VaccinationState.Unvaccinated] = 2.2
    infection_params.InfectedNoSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age0to4,
                                                   VaccinationState.Unvaccinated] = 9.2
    infection_params.InfectedSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age0to4,
                                                 VaccinationState.Unvaccinated] = 7
    infection_params.InfectedSymptomsToSevere[VirusVariant.Wildtype, AgeGroup.Age0to4,
                                              VaccinationState.Unvaccinated] = 10.5
    infection_params.SevereToRecovered[VirusVariant.Wildtype, AgeGroup.Age0to4,
                                       VaccinationState.Unvaccinated] = 5
    infection_params.SevereToCritical[VirusVariant.Wildtype, AgeGroup.Age0to4,
                                      VaccinationState.Unvaccinated] = 5
    infection_params.CriticalToRecovered[VirusVariant.Wildtype, AgeGroup.Age0to4,
                                         VaccinationState.Unvaccinated] = 7
    infection_params.CriticalToDead[VirusVariant.Wildtype, AgeGroup.Age0to4,
                                    VaccinationState.Unvaccinated] = 6
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age0to4,
                                      VaccinationState.Unvaccinated, 8.1, 8.1, 2.0, 2.0, -0.17, -0.17)   
    abm.set_infectivity_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age0to4, -7.0, -7.0, 1.0, 1.0)

    # AgeGroup 5-14
    infection_params.IncubationPeriod[VirusVariant.Wildtype, AgeGroup.Age5to14,
                                      VaccinationState.Unvaccinated] = 3
    infection_params.InfectedNoSymptomsToSymptoms[VirusVariant.Wildtype, AgeGroup.Age5to14,
                                                  VaccinationState.Unvaccinated] = 2.2
    infection_params.InfectedNoSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age5to14,
                                                   VaccinationState.Unvaccinated] = 9.2
    infection_params.InfectedSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age5to14,
                                                 VaccinationState.Unvaccinated] = 7
    infection_params.InfectedSymptomsToSevere[VirusVariant.Wildtype, AgeGroup.Age5to14,
                                              VaccinationState.Unvaccinated] = 10.5
    infection_params.SevereToRecovered[VirusVariant.Wildtype, AgeGroup.Age5to14,
                                       VaccinationState.Unvaccinated] = 5
    infection_params.SevereToCritical[VirusVariant.Wildtype, AgeGroup.Age5to14,
                                      VaccinationState.Unvaccinated] = 5
    infection_params.CriticalToRecovered[VirusVariant.Wildtype, AgeGroup.Age5to14,
                                         VaccinationState.Unvaccinated] = 7
    infection_params.CriticalToDead[VirusVariant.Wildtype, AgeGroup.Age5to14,
                                    VaccinationState.Unvaccinated] = 6
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age5to14,
                                      VaccinationState.Unvaccinated, 8.1, 8.1, 2.0, 2.0, -0.17, -0.17)   
    abm.set_infectivity_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age5to14, -7.0, -7.0, 1.0, 1.0)

    # AgeGroup 15-34
    infection_params.IncubationPeriod[VirusVariant.Wildtype, AgeGroup.Age15to34,
                                      VaccinationState.Unvaccinated] = 3
    infection_params.InfectedNoSymptomsToSymptoms[VirusVariant.Wildtype, AgeGroup.Age15to34,
                                                  VaccinationState.Unvaccinated] = 2.2
    infection_params.InfectedNoSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age15to34,
                                                   VaccinationState.Unvaccinated] = 9.2
    infection_params.InfectedSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age15to34,
                                                 VaccinationState.Unvaccinated] = 7
    infection_params.InfectedSymptomsToSevere[VirusVariant.Wildtype, AgeGroup.Age15to34,
                                              VaccinationState.Unvaccinated] = 10.5
    infection_params.SevereToRecovered[VirusVariant.Wildtype, AgeGroup.Age15to34,
                                       VaccinationState.Unvaccinated] = 6
    infection_params.SevereToCritical[VirusVariant.Wildtype, AgeGroup.Age15to34,
                                      VaccinationState.Unvaccinated] = 5
    infection_params.CriticalToRecovered[VirusVariant.Wildtype, AgeGroup.Age15to34,
                                         VaccinationState.Unvaccinated] = 7
    infection_params.CriticalToDead[VirusVariant.Wildtype, AgeGroup.Age15to34,
                                    VaccinationState.Unvaccinated] = 6
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age15to34,
                                      VaccinationState.Unvaccinated, 8.1, 8.1, 2.0, 2.0, -0.17, -0.17)   
    abm.set_infectivity_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age15to34, -7.0, -7.0, 1.0, 1.0)

    # AgeGroup 35-59
    infection_params.IncubationPeriod[VirusVariant.Wildtype, AgeGroup.Age35to59,
                                      VaccinationState.Unvaccinated] = 3
    infection_params.InfectedNoSymptomsToSymptoms[VirusVariant.Wildtype, AgeGroup.Age35to59,
                                                  VaccinationState.Unvaccinated] = 2.2
    infection_params.InfectedNoSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age35to59,
                                                   VaccinationState.Unvaccinated] = 9.2
    infection_params.InfectedSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age35to59,
                                                 VaccinationState.Unvaccinated] = 7
    infection_params.InfectedSymptomsToSevere[VirusVariant.Wildtype, AgeGroup.Age35to59,
                                              VaccinationState.Unvaccinated] = 6
    infection_params.SevereToRecovered[VirusVariant.Wildtype, AgeGroup.Age35to59,
                                       VaccinationState.Unvaccinated] = 8
    infection_params.SevereToCritical[VirusVariant.Wildtype, AgeGroup.Age35to59,
                                      VaccinationState.Unvaccinated] = 5
    infection_params.CriticalToRecovered[VirusVariant.Wildtype, AgeGroup.Age35to59,
                                         VaccinationState.Unvaccinated] = 17.5
    infection_params.CriticalToDead[VirusVariant.Wildtype, AgeGroup.Age35to59,
                                    VaccinationState.Unvaccinated] = 16.5
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age35to59,
                                      VaccinationState.Unvaccinated, 8.1, 8.1, 2.0, 2.0, -0.17, -0.17)   
    abm.set_infectivity_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age35to59, -7.0, -7.0, 1.0, 1.0)

    # AgeGroup 60-79
    infection_params.IncubationPeriod[VirusVariant.Wildtype, AgeGroup.Age60to79,
                                      VaccinationState.Unvaccinated] = 3
    infection_params.InfectedNoSymptomsToSymptoms[VirusVariant.Wildtype, AgeGroup.Age60to79,
                                                  VaccinationState.Unvaccinated] = 2.2
    infection_params.InfectedNoSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age60to79,
                                                   VaccinationState.Unvaccinated] = 9.2
    infection_params.InfectedSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age60to79,
                                                 VaccinationState.Unvaccinated] = 7.0
    infection_params.InfectedSymptomsToSevere[VirusVariant.Wildtype, AgeGroup.Age60to79,
                                              VaccinationState.Unvaccinated] = 6.0
    infection_params.SevereToRecovered[VirusVariant.Wildtype, AgeGroup.Age60to79,
                                       VaccinationState.Unvaccinated] = 10.0
    infection_params.SevereToCritical[VirusVariant.Wildtype, AgeGroup.Age60to79,
                                      VaccinationState.Unvaccinated] = 5
    infection_params.CriticalToRecovered[VirusVariant.Wildtype, AgeGroup.Age60to79,
                                         VaccinationState.Unvaccinated] = 17.5
    infection_params.CriticalToDead[VirusVariant.Wildtype, AgeGroup.Age60to79,
                                    VaccinationState.Unvaccinated] = 16.5
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age60to79,
                                      VaccinationState.Unvaccinated, 8.1, 8.1, 2.0, 2.0, -0.17, -0.17)    
    abm.set_infectivity_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age60to70, -7.0, -7.0, 1.0, 1.0)
    

    # AgeGroup 80+
    infection_params.IncubationPeriod[VirusVariant.Wildtype, AgeGroup.Age80plus,
                                      VaccinationState.Unvaccinated] = 3
    infection_params.InfectedNoSymptomsToSymptoms[VirusVariant.Wildtype, AgeGroup.Age80plus,
                                                  VaccinationState.Unvaccinated] = 2.2
    infection_params.InfectedNoSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age80plus,
                                                   VaccinationState.Unvaccinated] = 9.2
    infection_params.InfectedSymptomsToRecovered[VirusVariant.Wildtype, AgeGroup.Age80plus,
                                                 VaccinationState.Unvaccinated] = 7.0
    infection_params.InfectedSymptomsToSevere[VirusVariant.Wildtype, AgeGroup.Age80plus,
                                              VaccinationState.Unvaccinated] = 6.0
    infection_params.SevereToRecovered[VirusVariant.Wildtype, AgeGroup.Age80plus,
                                       VaccinationState.Unvaccinated] = 15.0
    infection_params.SevereToCritical[VirusVariant.Wildtype, AgeGroup.Age80plus,
                                      VaccinationState.Unvaccinated] = 5
    infection_params.CriticalToRecovered[VirusVariant.Wildtype, AgeGroup.Age80plus,
                                         VaccinationState.Unvaccinated] = 12.5
    infection_params.CriticalToDead[VirusVariant.Wildtype, AgeGroup.Age80plus,
                                    VaccinationState.Unvaccinated] = 11    
    abm.set_viral_load_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age80plus,
                                      VaccinationState.Unvaccinated, 8.1, 8.1, 2.0, 2.0, -0.17, -0.17)    
    abm.set_infectivity_parameters(infection_params, VirusVariant.Wildtype, AgeGroup.Age80plus, -7.0, -7.0, 1.0, 1.0)

    return infection_params


def read_txt(path):
    # read input file and save it in a pd.Dataframe
    return pd.read_csv(path, sep='\t')


def make_one_person_households(number_of_households):
    # one-person household member
    one_person_household_member = abm.HouseholdMember()
    # set weights for household member
    one_person_household_member.set_age_weight(AgeGroup.Age15to34, 4364)
    one_person_household_member.set_age_weight(AgeGroup.Age35to59, 7283)
    one_person_household_member.set_age_weight(AgeGroup.Age60to79, 4100)
    one_person_household_member.set_age_weight(AgeGroup.Age80plus, 1800)

    # create one-person household group
    household_group = abm.HouseholdGroup()
    for hh in range(int(number_of_households)):
        household = abm.Household()
        household.add_members(one_person_household_member, 1)
        household_group.add_households(household, 1)

    return household_group


def make_multiple_person_households(household_size,
                                    number_of_two_parent_households, number_of_one_parent_households,
                                    number_of_other_households):

    # members for multiple person households
    child = abm.HouseholdMember()
    child.set_age_weight(AgeGroup.Age0to4, 1)
    child.set_age_weight(AgeGroup.Age5to14, 1)

    parent = abm.HouseholdMember()
    parent.set_age_weight(AgeGroup.Age15to34, 2)
    parent.set_age_weight(AgeGroup.Age35to59, 2)
    parent.set_age_weight(AgeGroup.Age60to79, 1)

    other = abm.HouseholdMember()
    other.set_age_weight(AgeGroup.Age0to4, 5000)
    other.set_age_weight(AgeGroup.Age5to14, 6000)
    other.set_age_weight(AgeGroup.Age15to34, 14943)
    other.set_age_weight(AgeGroup.Age35to59, 22259)
    other.set_age_weight(AgeGroup.Age60to79, 11998)
    other.set_age_weight(AgeGroup.Age80plus, 5038)

    household_group = abm.HouseholdGroup()
    # add two parent households
    hh_two_parents = abm.Household()
    hh_two_parents.add_members(child, household_size - 2)
    hh_two_parents.add_members(parent, 2)
    household_group.add_households(
        hh_two_parents, number_of_two_parent_households)
    # add one parent households
    hh_one_parent = abm.Household()
    hh_one_parent.add_members(child, household_size - 1)
    hh_one_parent.add_members(parent, 1)
    household_group.add_households(
        hh_one_parent, number_of_one_parent_households)
    # add other households
    hh_other = abm.Household()
    hh_other.add_members(other, household_size)
    household_group.add_households(hh_other, number_of_other_households)

    return household_group


def add_households(world, distribution, num_inhabitants):
    household_sizes = np.zeros(len(distribution))
    locationIds = []
    new_index = len(world.locations)
    # distribute inhabintants to households
    while num_inhabitants > 0:
        size = np.random.choice(
            np.arange(0, len(distribution)), p=distribution)
        household_sizes[size] += 1
        num_inhabitants -= (size + 1)

    # one-person household
    one_person_household_group = make_one_person_households(household_sizes[0])
    abm.add_household_group_to_world(world, one_person_household_group)

    # two-person households
    two_person_two_parents = int(0.85 * household_sizes[1])
    two_person_one_parent = int(0.13 * household_sizes[1])
    two_person_other = int(
        household_sizes[1] - two_person_two_parents - two_person_one_parent)
    two_person_household_group = make_multiple_person_households(2, two_person_two_parents,
                                                                 two_person_one_parent, two_person_other)
    abm.add_household_group_to_world(world, two_person_household_group)

    # three-person households
    three_person_two_parents = int(0.83 * household_sizes[2])
    three_person_one_parent = int(0.13 * household_sizes[2])
    three_person_other = int(
        household_sizes[2] - three_person_two_parents - three_person_one_parent)
    three_person_household_group = make_multiple_person_households(3, three_person_two_parents,
                                                                   three_person_one_parent, three_person_other)
    abm.add_household_group_to_world(world, three_person_household_group)

    # four-person households
    four_person_two_parents = int(0.93 * household_sizes[3])
    four_person_one_parent = int(0.03 * household_sizes[3])
    four_person_other = int(
        household_sizes[3] - four_person_two_parents - four_person_one_parent)
    four_person_household_group = make_multiple_person_households(4, four_person_two_parents,
                                                                  four_person_one_parent, four_person_other)
    abm.add_household_group_to_world(world, four_person_household_group)

    # five-person households
    five_person_two_parents = int(0.88 * household_sizes[4])
    five_person_one_parent = int(0.05 * household_sizes[4])
    five_person_other = int(
        household_sizes[4] - five_person_two_parents - five_person_one_parent)
    five_person_household_group = make_multiple_person_households(5, five_person_two_parents,
                                                                  five_person_one_parent, five_person_other)
    abm.add_household_group_to_world(world, five_person_household_group)

    for hh in range(int(sum(household_sizes))):
        id = abm.LocationId(new_index, abm.LocationType.Home)
        locationIds.append(id)
        new_index += 1
    return locationIds


def insert_locations_to_map(mapping, inputId, locationIds):
    map = LocationMapping()
    map.inputId = inputId
    for id in locationIds:
        type = str(int(id.type))
        index = str(id.index)
        if int(id.type) < 10:
            type = "0" + type
        if int(id.index) < 10:
            index = "0" + index
        map.modelId.append(type+index)
    mapping.append(map)
    return mapping


def create_locations_from_input(world, input_areas, household_distribution):
    # set seeds to have fixed locations for given input
    np.random.seed(0)
    # map input area ids to corresponding abm location ids
    mapping = []
    # bools to make sure the world has a school and a hospital
    has_school = False
    has_hospital = False
    for index, area in input_areas.iterrows():
        locationIds = []
        if ('residential' in area.type):
            # area 'residential' corresponds to location type 'Home'
            locationIds = add_households(
                world, household_distribution, area.inhabitants)
        elif (area.type == 'recreational'):
            # area 'recreational' corresponds to location type 'SocialEvent'
            location = world.add_location(abm.LocationType.SocialEvent)
            # set maximum contacts and capacity for social events
            world.get_individualized_location(
                location).infection_parameters.MaximumContacts = 30.
            world.get_individualized_location(location).set_capacity(30, 40, 0)
            locationIds.append(location)
        elif (area.type == 'shopping_business'):
            # area 'shopping_business' corresponds to location type School, Hospital, BasicsShops, Work
            if (not has_school):
                # if world does not have a school yet, a school is added
                location = world.add_location(abm.LocationType.School)
                # set maximum contacts and capacity for school
                world.get_individualized_location(
                    location).infection_parameters.MaximumContacts = 40.
                world.get_individualized_location(
                    location).set_capacity(500, 2000, 0)
                locationIds.append(location)
                has_school = True
            elif (not has_hospital):
                # if world does not have a hospital yet, a hospital and a icu is added
                locHosp = world.add_location(abm.LocationType.Hospital)
                # set maximum contacts and capacity for hospital
                world.get_individualized_location(
                    locHosp).infection_parameters.MaximumContacts = 5.
                world.get_individualized_location(
                    locHosp).set_capacity(300, 10000, 0)
                locICU = world.add_location(abm.LocationType.ICU)
                # set maximum contacts and capacity for icu
                world.get_individualized_location(
                    locICU).infection_parameters.MaximumContacts = 5.
                world.get_individualized_location(
                    locICU).set_capacity(30, 1000, 0)
                locationIds.append(locHosp)
                locationIds.append(locICU)
                has_hospital = True
            else:
                # when hospital and school has been added, the area 'shopping_business' is either
                # transformed to location type BasicsShop or location type Work with same probability
                type = np.random.choice(np.arange(0, 2), p=[0.5, 0.5])
                if (type):
                    location = world.add_location(abm.LocationType.BasicsShop)
                    # set maximum contacts and capacity for basics shops
                    world.get_individualized_location(
                        location).infection_parameters.MaximumContacts = 20.
                    world.get_individualized_location(
                        location).set_capacity(100, 1000, 0)
                    locationIds.append(location)
                else:
                    location = world.add_location(abm.LocationType.Work)
                    # set maximum contacts and capacity for work
                    world.get_individualized_location(
                        location).infection_parameters.MaximumContacts = 40.
                    world.get_individualized_location(
                        location).set_capacity(300, 2000, 0)
                    locationIds.append(location)
        elif (area.type == 'university'):
            # area 'university' corresponds to location type 'Work'
            location = world.add_location(abm.LocationType.Work)
            # set maximum contacts and capacity for work
            world.get_individualized_location(
                location).infection_parameters.MaximumContacts = 50.
            world.get_individualized_location(
                location).set_capacity(200, 4000, 0)
            locationIds.append(location)
        elif (area.type == 'mixed'):
            # area 'mixed' corresponds either to location type 'Work' or location type 'Home' with same probability
            type = np.random.choice(np.arange(0, 2), p=[0.5, 0.5])
            if (type):
                location = world.add_location(abm.LocationType.Work)
                # set maximum contacts and capacity for work
                world.get_individualized_location(
                    location).infection_parameters.MaximumContacts = 40.
                world.get_individualized_location(
                    location).set_capacity(100, 2000, 0)
                locationIds.append(location)
            else:
                locationIds = add_households(
                    world, household_distribution, area.inhabitants)
        insert_locations_to_map(mapping, str(area.id), locationIds)
    return mapping


def assign_infection_states(world, t0, exposed_pct, infected_no_symptoms_pct, infected_symptoms_pct,
                            infected_severe_pct, infected_critical_pct, recovered_pct):
    susceptible_pct = 1 - exposed_pct - infected_no_symptoms_pct - \
        infected_symptoms_pct - infected_severe_pct - \
        infected_critical_pct - recovered_pct
    for person in world.persons:
        # draw infection state from distribution for every agent
        infection_state = np.random.choice(np.arange(0, int(abm.InfectionState.Count)),
                                           p=[susceptible_pct, exposed_pct, infected_no_symptoms_pct,
                                               infected_symptoms_pct, infected_severe_pct, infected_critical_pct, recovered_pct, 0.0])
        if (abm.InfectionState(infection_state) != abm.InfectionState.Susceptible):
            person.add_new_infection(Infection(VirusVariant.Wildtype, person.age,
                                     world.infection_parameters, t0, abm.InfectionState(infection_state), False))


def find_all_locations_of_type(world, type):
    locations = []
    for loc in world.locations:
        if (loc.type == type):
            locations.append(abm.LocationId(loc.index, type))
    return locations


def assign_locations(world):
    # get locations from world
    schools = find_all_locations_of_type(world, abm.LocationType.School)
    school_weights = [(1/len(schools)) for i in range(len(schools))]
    hospitals = find_all_locations_of_type(world, abm.LocationType.Hospital)
    hospital_weights = [(1/len(hospitals)) for i in range(len(hospitals))]
    icus = find_all_locations_of_type(world, abm.LocationType.ICU)
    icu_weights = [(1/len(icus)) for i in range(len(icus))]
    workplaces = find_all_locations_of_type(world, abm.LocationType.Work)
    workplace_weights = [(1/len(workplaces)) for i in range(len(workplaces))]
    basic_shops = find_all_locations_of_type(
        world, abm.LocationType.BasicsShop)
    shop_weights = [(1/len(basic_shops)) for i in range(len(basic_shops))]
    social_events = find_all_locations_of_type(
        world, abm.LocationType.SocialEvent)
    event_weights = [(1/len(social_events)) for i in range(len(social_events))]

    # assign locations to agents
    for person in world.persons:
        shop = np.random.choice(np.arange(0, len(basic_shops)), p=shop_weights)
        person.set_assigned_location(abm.LocationId(
            basic_shops[int(shop)].index, basic_shops[int(shop)].type))

        hospital = np.random.choice(
            np.arange(0, len(hospitals)), p=hospital_weights)
        person.set_assigned_location(abm.LocationId(
            hospitals[int(hospital)].index, hospitals[int(hospital)].type))

        icu = np.random.choice(np.arange(0, len(icus)), p=icu_weights)
        person.set_assigned_location(abm.LocationId(
            icus[int(icu)].index, icus[int(icu)].type))

        event = np.random.choice(
            np.arange(0, len(social_events)), p=event_weights)
        person.set_assigned_location(abm.LocationId(
            social_events[int(event)].index, social_events[int(event)].type))

        # assign school to agents between 5 and 14 years
        if (person.age == AgeGroup.Age5to14):
            school = np.random.choice(
                np.arange(0, len(schools)), p=school_weights)
            person.set_assigned_location(abm.LocationId(
                schools[int(school)].index, schools[int(school)].type))
        # assign work to agents between 15 and 59
        if (person.age == AgeGroup.Age15to34 or person.age == AgeGroup.Age35to59):
            work = np.random.choice(
                np.arange(0, len(workplaces)), p=workplace_weights)
            person.set_assigned_location(abm.LocationId(
                workplaces[int(work)].index, workplaces[int(work)].type))


def convert_loc_id_to_string(loc_id):
    type = str(int(loc_id[0]))
    index = str(loc_id[1])
    if int(int(loc_id[0])) < 10:
        type = "0" + type
    if int(loc_id[1]) < 10:
        index = "0" + index

    return type + index


def get_agents_per_location(loc_id, agents):
    agents_per_loc = []
    for a in agents:
        if (int(a[0].type) == int(loc_id[0]) and a[0].index == loc_id[1]):
            agents_per_loc.append(a)
    return agents_per_loc


def write_results_to_file(path, log):
    location_ids = log[1][0]
    time_points = log[0]
    agents = log[2]
    with open(path, 'w') as f:
        for location_id_index in range(len(location_ids)):
            line = convert_loc_id_to_string(
                location_ids[location_id_index]) + " " + str(len(time_points))
            for t in range(len(time_points)):
                agents_per_loc = get_agents_per_location(
                    location_ids[location_id_index], agents[t])
                line += " " + str(time_points[t]) + \
                    " " + str(len(agents_per_loc))
                for a in agents_per_loc:
                    if (a[2].days > 1000):
                        time_since_transmission = -1.0
                    else:
                        time_since_transmission = a[2].hours
                    line += " " + str(a[1]) + " " + \
                        str(time_since_transmission)
            f.write(line)
            f.write('\n')
    f.close()


def write_location_mapping_to_file(path, mapping):
    with open(path, 'w') as f:
        for id in mapping:
            line = id.inputId + " "
            for modelId in id.modelId:
                line += modelId + " "
            f.write(line)
            f.write('\n')
        f.close()


def set_sim_result_at_start(sim):
    for location in sim.world.locations:
        result = sim.result.get_last_value()
        result += location.population.get_last_value()


def run_abm_simulation():

    input_path = 'C:/Users/bick_ju/Documents/INSIDe/Demonstrator/INSIDeDemonstrator/INSIDe_Demonstrator_AreaList_modified.txt'
    output_path = 'C:/Users/bick_ju/Documents/INSIDe/Demonstrator/INSIDeDemonstrator/output/'
    # set seed for fixed model initialization (locations and initial infection states)
    np.random.seed(0)
    # starting time point
    t0 = abm.TimePoint(0)
    # end time point: simulation will run 14 days
    tmax = t0 + abm.days(14)
    # distribution used to distribute the residential areas to one-, two-, three-, four- and five-person households
    household_distribution = [0.4084, 0.3375, 0.1199, 0.0965, 0.0377]
    # read txt file with inputs
    areas = read_txt(input_path)
    # create simulation
    sim = abm.Simulation(t0)
    # set infection parameters
    sim.world.infection_parameters = set_infection_parameters()
    # as input areas do fit one-to-one to abm location types they have there has to be a mapping
    mapping = create_locations_from_input(
        sim.world, areas, household_distribution)
    # assign initial infection states according to distribution
    assign_infection_states(sim.world, t0, 0.005, 0.001,
                            0.001, 0.0001, 0.0, 0.0)
    # assign locations to agents
    assign_locations(sim.world)
    # output object
    history = History()
    # advance simulation until tmax
    sim.advance(tmax, history)
    # results collected during the simulation
    log = history.log
    # write simulation results to txt file
    write_results_to_file(output_path + 'output.txt', log)
    # write location mapping to txt file
    write_location_mapping_to_file(
        output_path + 'location_mapping.txt', mapping)

    # print compartment values to csv
    # only used for validation purposes
    with open(output_path + 'console_output.csv', 'w') as f:
        f.write("t S E C I I_s I_c R D\n")
        for t in range(sim.result.get_num_time_points()):
            line = str(sim.result.get_time(t)) + " "
            comps = sim.result.get_value(t)
            for c in range(len(comps)):
                line += str(comps[c]) + " "
            f.write(line + "\n")
    print('done')


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'abm demonstrator',
        description='Example demonstrating the agent-based model for a synthetic population.')
    args = arg_parser.parse_args()
    run_abm_simulation(**args.__dict__)
