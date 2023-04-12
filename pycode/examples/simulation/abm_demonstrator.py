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


class LocationMapping:

    def __init__(self):
        self.inputId = None
        self.modelId = []


def set_infection_parameters():
    infection_params = abm.GlobalInfectionParameters()

    # 0-4
    infection_params.IncubationPeriod[AgeGroup.Age0to4,
                                      VaccinationState.Unvaccinated] = 4
    infection_params.SusceptibleToExposedByCarrier[AgeGroup.Age0to4,
                                                   VaccinationState.Unvaccinated] = 0.05
    infection_params.SusceptibleToExposedByInfected[AgeGroup.Age0to4,
                                                    VaccinationState.Unvaccinated] = 0.05
    infection_params.CarrierToInfected[AgeGroup.Age0to4,
                                       VaccinationState.Unvaccinated] = 0.276
    infection_params.CarrierToRecovered[AgeGroup.Age0to4,
                                        VaccinationState.Unvaccinated] = 0.092
    infection_params.InfectedToRecovered[AgeGroup.Age0to4,
                                         VaccinationState.Unvaccinated] = 0.142
    infection_params.InfectedToSevere[AgeGroup.Age0to4,
                                      VaccinationState.Unvaccinated] = 0.001
    infection_params.SevereToRecovered[AgeGroup.Age0to4,
                                       VaccinationState.Unvaccinated] = 0.186
    infection_params.SevereToCritical[AgeGroup.Age0to4,
                                      VaccinationState.Unvaccinated] = 0.015
    infection_params.CriticalToRecovered[AgeGroup.Age0to4,
                                         VaccinationState.Unvaccinated] = 0.143
    infection_params.CriticalToDead[AgeGroup.Age0to4,
                                    VaccinationState.Unvaccinated] = 0.001
    infection_params.RecoveredToSusceptible[AgeGroup.Age0to4,
                                            VaccinationState.Unvaccinated] = 0

    # 5-14
    infection_params.IncubationPeriod[AgeGroup.Age5to14,
                                      VaccinationState.Unvaccinated] = 4
    infection_params.SusceptibleToExposedByCarrier[AgeGroup.Age5to14,
                                                   VaccinationState.Unvaccinated] = 0.1
    infection_params.SusceptibleToExposedByInfected[AgeGroup.Age5to14,
                                                    VaccinationState.Unvaccinated] = 0.1
    infection_params.CarrierToInfected[AgeGroup.Age5to14,
                                       VaccinationState.Unvaccinated] = 0.276
    infection_params.CarrierToRecovered[AgeGroup.Age5to14,
                                        VaccinationState.Unvaccinated] = 0.092
    infection_params.InfectedToRecovered[AgeGroup.Age5to14,
                                         VaccinationState.Unvaccinated] = 0.142
    infection_params.InfectedToSevere[AgeGroup.Age5to14,
                                      VaccinationState.Unvaccinated] = 0.001
    infection_params.SevereToRecovered[AgeGroup.Age5to14,
                                       VaccinationState.Unvaccinated] = 0.186
    infection_params.SevereToCritical[AgeGroup.Age5to14,
                                      VaccinationState.Unvaccinated] = 0.015
    infection_params.CriticalToRecovered[AgeGroup.Age5to14,
                                         VaccinationState.Unvaccinated] = 0.143
    infection_params.CriticalToDead[AgeGroup.Age5to14,
                                    VaccinationState.Unvaccinated] = 0.001
    infection_params.RecoveredToSusceptible[AgeGroup.Age5to14,
                                            VaccinationState.Unvaccinated] = 0.

    # 15-34
    infection_params.IncubationPeriod[AgeGroup.Age15to34,
                                      VaccinationState.Unvaccinated] = 4
    infection_params.SusceptibleToExposedByCarrier[AgeGroup.Age15to34,
                                                   VaccinationState.Unvaccinated] = 0.13
    infection_params.SusceptibleToExposedByInfected[AgeGroup.Age15to34,
                                                    VaccinationState.Unvaccinated] = 0.13
    infection_params.CarrierToInfected[AgeGroup.Age15to34,
                                       VaccinationState.Unvaccinated] = 0.315
    infection_params.CarrierToRecovered[AgeGroup.Age15to34,
                                        VaccinationState.Unvaccinated] = 0.079
    infection_params.InfectedToRecovered[AgeGroup.Age15to34,
                                         VaccinationState.Unvaccinated] = 0.139
    infection_params.InfectedToSevere[AgeGroup.Age15to34,
                                      VaccinationState.Unvaccinated] = 0.003
    infection_params.SevereToRecovered[AgeGroup.Age15to34,
                                       VaccinationState.Unvaccinated] = 0.157
    infection_params.SevereToCritical[AgeGroup.Age15to34,
                                      VaccinationState.Unvaccinated] = 0.013
    infection_params.CriticalToRecovered[AgeGroup.Age15to34,
                                         VaccinationState.Unvaccinated] = 0.126
    infection_params.CriticalToDead[AgeGroup.Age15to34,
                                    VaccinationState.Unvaccinated] = 0.021
    infection_params.RecoveredToSusceptible[AgeGroup.Age15to34,
                                            VaccinationState.Unvaccinated] = 0.

    # 35-59
    infection_params.IncubationPeriod[AgeGroup.Age35to59,
                                      VaccinationState.Unvaccinated] = 4
    infection_params.SusceptibleToExposedByCarrier[AgeGroup.Age35to59,
                                                   VaccinationState.Unvaccinated] = 0.11
    infection_params.SusceptibleToExposedByInfected[AgeGroup.Age35to59,
                                                    VaccinationState.Unvaccinated] = 0.11
    infection_params.CarrierToInfected[AgeGroup.Age35to59,
                                       VaccinationState.Unvaccinated] = 0.315
    infection_params.CarrierToRecovered[AgeGroup.Age35to59,
                                        VaccinationState.Unvaccinated] = 0.079
    infection_params.InfectedToRecovered[AgeGroup.Age35to59,
                                         VaccinationState.Unvaccinated] = 0.136
    infection_params.InfectedToSevere[AgeGroup.Age35to59,
                                      VaccinationState.Unvaccinated] = 0.009
    infection_params.SevereToRecovered[AgeGroup.Age35to59,
                                       VaccinationState.Unvaccinated] = 0.113
    infection_params.SevereToCritical[AgeGroup.Age35to59,
                                      VaccinationState.Unvaccinated] = 0.02
    infection_params.CriticalToRecovered[AgeGroup.Age35to59,
                                         VaccinationState.Unvaccinated] = 0.05
    infection_params.CriticalToDead[AgeGroup.Age35to59,
                                    VaccinationState.Unvaccinated] = 0.008
    infection_params.RecoveredToSusceptible[AgeGroup.Age35to59,
                                            VaccinationState.Unvaccinated] = 0.

    # 60-79
    infection_params.IncubationPeriod[AgeGroup.Age60to79,
                                      VaccinationState.Unvaccinated] = 4
    infection_params.SusceptibleToExposedByCarrier[AgeGroup.Age60to79,
                                                   VaccinationState.Unvaccinated] = 0.04
    infection_params.SusceptibleToExposedByInfected[AgeGroup.Age60to79,
                                                    VaccinationState.Unvaccinated] = 0.04
    infection_params.CarrierToInfected[AgeGroup.Age60to79,
                                       VaccinationState.Unvaccinated] = 0.315
    infection_params.CarrierToRecovered[AgeGroup.Age60to79,
                                        VaccinationState.Unvaccinated] = 0.079
    infection_params.InfectedToRecovered[AgeGroup.Age60to79,
                                         VaccinationState.Unvaccinated] = 0.123
    infection_params.InfectedToSevere[AgeGroup.Age60to79,
                                      VaccinationState.Unvaccinated] = 0.024
    infection_params.SevereToRecovered[AgeGroup.Age60to79,
                                       VaccinationState.Unvaccinated] = 0.083
    infection_params.SevereToCritical[AgeGroup.Age60to79,
                                      VaccinationState.Unvaccinated] = 0.035
    infection_params.CriticalToRecovered[AgeGroup.Age60to79,
                                         VaccinationState.Unvaccinated] = 0.035
    infection_params.CriticalToDead[AgeGroup.Age60to79,
                                    VaccinationState.Unvaccinated] = 0.023
    infection_params.RecoveredToSusceptible[AgeGroup.Age60to79,
                                            VaccinationState.Unvaccinated] = 0.

    # 80+
    infection_params.IncubationPeriod[AgeGroup.Age80plus,
                                      VaccinationState.Unvaccinated] = 4
    infection_params.SusceptibleToExposedByCarrier[AgeGroup.Age80plus,
                                                   VaccinationState.Unvaccinated] = 0.07
    infection_params.SusceptibleToExposedByInfected[AgeGroup.Age80plus,
                                                    VaccinationState.Unvaccinated] = 0.07
    infection_params.CarrierToInfected[AgeGroup.Age80plus,
                                       VaccinationState.Unvaccinated] = 0.315
    infection_params.CarrierToRecovered[AgeGroup.Age80plus,
                                        VaccinationState.Unvaccinated] = 0.079
    infection_params.InfectedToRecovered[AgeGroup.Age80plus,
                                         VaccinationState.Unvaccinated] = 0.115
    infection_params.InfectedToSevere[AgeGroup.Age80plus,
                                      VaccinationState.Unvaccinated] = 0.033
    infection_params.SevereToRecovered[AgeGroup.Age80plus,
                                       VaccinationState.Unvaccinated] = 0.055
    infection_params.SevereToCritical[AgeGroup.Age80plus,
                                      VaccinationState.Unvaccinated] = 0.036
    infection_params.CriticalToRecovered[AgeGroup.Age80plus,
                                         VaccinationState.Unvaccinated] = 0.035
    infection_params.CriticalToDead[AgeGroup.Age80plus,
                                    VaccinationState.Unvaccinated] = 0.052
    infection_params.RecoveredToSusceptible[AgeGroup.Age80plus,
                                            VaccinationState.Unvaccinated] = 0.

    return infection_params


def read_txt(path):
    input = pd.read_csv(path, sep='\t')
    return input


def make_one_person_households(number_of_households):
    # one-person household member
    one_person_household_member = abm.HouseholdMember()
    one_person_household_member.set_age_weight(AgeGroup.Age15to34, 5)
    one_person_household_member.set_age_weight(AgeGroup.Age35to59, 6)
    one_person_household_member.set_age_weight(AgeGroup.Age60to79, 4)
    one_person_household_member.set_age_weight(AgeGroup.Age80plus, 2)

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
    other.set_age_weight(AgeGroup.Age0to4, 1)
    other.set_age_weight(AgeGroup.Age5to14, 2)
    other.set_age_weight(AgeGroup.Age15to34, 3)
    other.set_age_weight(AgeGroup.Age35to59, 3)
    other.set_age_weight(AgeGroup.Age60to79, 2)
    other.set_age_weight(AgeGroup.Age80plus, 2)

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
    new_index = len(world.locations[int(abm.LocationType.Home)])
    while num_inhabitants > 0:
        size = np.random.choice(
            np.arange(0, len(distribution)), p=distribution)
        household_sizes[size] += 1
        num_inhabitants -= (size + 1)

    # one-person household
    one_person_household_group = make_one_person_households(household_sizes[0])
    abm.add_household_group_to_world(world, one_person_household_group)

    # two-person households
    two_person_two_parents = int(0.4 * household_sizes[1])
    two_person_one_parent = int(0.4 * household_sizes[1])
    two_person_other = int(
        household_sizes[1] - two_person_two_parents - two_person_one_parent)
    two_person_household_group = make_multiple_person_households(2, two_person_two_parents,
                                                                 two_person_one_parent, two_person_other)
    abm.add_household_group_to_world(world, two_person_household_group)

    # three-person households
    three_person_two_parents = int(0.4 * household_sizes[2])
    three_person_one_parent = int(0.4 * household_sizes[2])
    three_person_other = int(
        household_sizes[2] - three_person_two_parents - three_person_one_parent)
    three_person_household_group = make_multiple_person_households(3, three_person_two_parents,
                                                                   three_person_one_parent, three_person_other)
    abm.add_household_group_to_world(world, three_person_household_group)

    # four-person households
    four_person_two_parents = int(0.5 * household_sizes[3])
    four_person_one_parent = int(0.2 * household_sizes[3])
    four_person_other = int(
        household_sizes[3] - four_person_two_parents - four_person_one_parent)
    four_person_household_group = make_multiple_person_households(4, four_person_two_parents,
                                                                  four_person_one_parent, four_person_other)
    abm.add_household_group_to_world(world, four_person_household_group)

    # five-person households
    five_person_two_parents = int(0.6 * household_sizes[4])
    five_person_one_parent = int(0.1 * household_sizes[4])
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


def create_locations_from_input(world, input_areas):
    mapping = []
    has_school = False
    has_hospital = False
    print(input_areas)
    for index, area in input_areas.iterrows():
        locationIds = []
        if ('residential' in area.type):
            locationIds = add_households(
                world, [0.2, 0.2, 0.2, 0.2, 0.2], area.inhabitants)
        elif (area.type == 'recreational'):
            location = world.add_location(abm.LocationType.SocialEvent)
            locationIds.append(location)
            #world.get_individualized_location(location).infection_parameters.MaximumContacts = 30.
            #world.get_individualized_location(location).set_capacity(30, 40)
        elif (area.type == 'shopping_business'):
            if (not has_school):
                location = world.add_location(abm.LocationType.School)
                locationIds.append(location)
                has_school = True
            elif (not has_hospital):
                locHosp = world.add_location(abm.LocationType.Hospital)
                locICU = world.add_location(abm.LocationType.ICU)
                locationIds.append(locHosp)
                locationIds.append(locICU)
                has_hospital = True
            else:
                type = np.random.choice(np.arange(0, 2), p=[0.5, 0.5])
                if (type):
                    location = world.add_location(abm.LocationType.BasicsShop)
                    locationIds.append(location)
                else:
                    location = world.add_location(abm.LocationType.Work)
                    locationIds.append(location)
        elif (area.type == 'university'):
            location = world.add_location(abm.LocationType.Work)
            locationIds.append(location)
        elif (area.type == 'mixed'):
            type = np.random.choice(np.arange(0, 2), p=[0.5, 0.5])
            if (type):
                location = world.add_location(abm.LocationType.Work)
                locationIds.append(location)
            else:
                locationIds = add_households(
                    world, [0.2, 0.2, 0.2, 0.2, 0.2], area.inhabitants)
        insert_locations_to_map(mapping, area.id, locationIds)


def assign_infection_states(world, exposed_pct, infected_no_symptoms_pct, infected_symptoms_pct,
                            infected_severe_pct, infected_critical_pct, recovered_infected_no_symptoms,
                            recovered_infected):
    susceptible_pct = 1 - exposed_pct - infected_no_symptoms_pct - infected_symptoms_pct - \
        infected_severe_pct - infected_critical_pct - \
        recovered_infected_no_symptoms - recovered_infected
    for person in world.persons:
        infection_state = np.random.choice(np.arange(0, int(abm.InfectionState.Count)),
                                           p=[susceptible_pct, exposed_pct, infected_no_symptoms_pct,
                                               infected_symptoms_pct, infected_severe_pct, infected_critical_pct, recovered_infected_no_symptoms,
                                               recovered_infected, 0.0])
        person.infection_state = abm.InfectionState(infection_state)


def assign_locations(world):
    schools = world.locations[abm.LocationType.School]
    school_weights = [(1/len(schools)) for i in range(len(schools))]

    hospitals = world.locations[abm.LocationType.Hospital]
    hospital_weights = [(1/len(hospitals)) for i in range(len(hospitals))]

    icus = world.locations[abm.LocationType.ICU]
    icu_weights = [(1/len(icus)) for i in range(len(icus))]

    workplaces = world.locations[abm.LocationType.Work]
    workplace_weights = [(1/len(workplaces)) for i in range(len(workplaces))]

    basic_shops = world.locations[abm.LocationType.BasicsShop]
    shop_weights = [(1/len(basic_shops)) for i in range(len(basic_shops))]

    social_events = world.locations[abm.LocationType.SocialEvent]
    event_weights = [(1/len(social_events)) for i in range(len(social_events))]

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

        if (person.age == AgeGroup.Age5to14):
            school = np.random.choice(
                np.arange(0, len(schools)), p=school_weights)
            person.set_assigned_location(abm.LocationId(
                schools[int(school)].index, schools[int(school)].type))

        if (person.age == AgeGroup.Age15to34 or person.age == AgeGroup.Age35to59):
            work = np.random.choice(
                np.arange(0, len(workplaces)), p=workplace_weights)
            person.set_assigned_location(abm.LocationId(
                workplaces[int(work)].index, workplaces[int(work)].type))


def run_abm_simulation():

    LocationIds = []

    infection_params = set_infection_parameters()
    world = abm.World(infection_params)

    t0 = abm.TimePoint(0)
    tmax = t0 + abm.days(14)

    areas = read_txt(
        '/home/bick_ju/Documents/INSIDeDemonstrator/INSIDe_Demonstrator_AreaList.txt')
    create_locations_from_input(world, areas)
    assign_infection_states(world, 0.005, 0.001, 0.001, 0.0001, 0.0, 0.0, 0.0)
    assign_locations(world)
    sim = abm.Simulation(t0)
    sim.advance(tmax)
    print('done')


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'abm demonstrator',
        description='Example demonstrating the agent-based model for a synthetic population.')
    args = arg_parser.parse_args()
    run_abm_simulation(**args.__dict__)
