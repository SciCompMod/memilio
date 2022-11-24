#############################################################################
# Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
#
# Authors: Daniel Abele
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

import unittest

import numpy as np

import memilio.simulation as mio
import memilio.simulation.abm as abm


class TestAbm(unittest.TestCase):
    def test_world(self):
        t0 = abm.TimePoint(0)
        sim = abm.Simulation(t0)
        world = sim.world
        self.assertEqual(len(world.persons), 0)
        self.assertEqual(sum(map(len, world.locations)), 0)
        self.assertEqual(len(sim.result), 1)

    def test_locations(self):
        t0 = abm.TimePoint(0)
        sim = abm.Simulation(t0)
        world = sim.world

        home_id = world.add_location(abm.LocationType.Home)
        social_event_id = world.add_location(abm.LocationType.SocialEvent)
        self.assertEqual(sum(map(len, world.locations)), 2)

        home = world.locations[home_id.type][home_id.index]
        self.assertEqual(home.type, abm.LocationType.Home)

        home.infection_parameters.MaximumContacts = 10
        self.assertEqual(home.infection_parameters.MaximumContacts, 10)

        testing_ages = [abm.AgeGroup.Age0to4]
        testing_locations = [abm.LocationType.Home]
        testing_inf_states = []
        testing_crit = [abm.TestingCriteria(
            testing_ages, testing_locations, testing_inf_states)]
        testing_scheme = abm.TestingScheme(testing_crit, abm.days(
            1), t0, t0 + abm.days(1), abm.AntigenTest(), 1.0)
        # initially false, will only active once simulation starts
        self.assertEqual(testing_scheme.active, False)

    def test_persons(self):
        t0 = abm.TimePoint(0)
        sim = abm.Simulation(t0)
        world = sim.world

        home_id = world.add_location(abm.LocationType.Home)
        social_event_id = world.add_location(abm.LocationType.SocialEvent)

        p1 = world.add_person(
            home_id, abm.InfectionState.Carrier, abm.AgeGroup.Age15to34)
        p2 = world.add_person(
            social_event_id, abm.InfectionState.Recovered_Infected, abm.AgeGroup.Age80plus)

        # check persons
        self.assertEqual(len(world.persons), 2)
        self.assertEqual(p1.age, abm.AgeGroup.Age15to34)
        self.assertEqual(p1.location_id, home_id)
        self.assertEqual(p2.infection_state,
                         abm.InfectionState.Recovered_Infected)
        self.assertEqual(world.persons[0], p1)
        self.assertEqual(world.persons[1], p2)

    def test_simulation(self):
        t0 = abm.TimePoint(0)
        sim = abm.Simulation(t0)
        world = sim.world

        # add some locations and persons
        for type in abm.LocationType.values():
            world.add_location(type)
        home_id = abm.LocationId(0, abm.LocationType.Home)
        social_event_id = abm.LocationId(0, abm.LocationType.SocialEvent)
        work_id = abm.LocationId(0, abm.LocationType.Work)
        p1 = world.add_person(
            home_id, abm.InfectionState.Infected, abm.AgeGroup.Age0to4)
        p2 = world.add_person(
            home_id, abm.InfectionState.Recovered_Carrier, abm.AgeGroup.Age15to34)
        for type in abm.LocationType.values():
            p1.set_assigned_location(abm.LocationId(0, type))
            p2.set_assigned_location(abm.LocationId(0, type))

        # parameters so that the infected person doesn't randomly change state and gets tested reliably
        # DUE TO THE CURRENT IMPLEMENTATION OF DIFFERENT TEST TYPES, THIS IS NOT POSSIBLE, NEEDS TO BE CHANGED IN THE FUTURE
        social_event = world.locations[social_event_id.type][social_event_id.index]
        #social_event.testing_scheme = abm.TestingScheme(abm.days(1), 1.0)
        #world.testing_parameters.AntigenTest = abm.TestParameters(1, 1)
        world.infection_parameters.InfectedToSevere[abm.AgeGroup.Age0to4,
                                                    abm.VaccinationState.Unvaccinated] = 0.0
        world.infection_parameters.InfectedToRecovered[abm.AgeGroup.Age0to4,
                                                       abm.VaccinationState.Unvaccinated] = 0.0

        # trips
        trip_list = abm.TripList()
        trip_list.add_trip(abm.Trip(0, abm.TimePoint(
            0) + abm.hours(8), social_event_id, home_id))
        trip_list.add_trip(abm.Trip(1, abm.TimePoint(0) +
                           abm.hours(8), work_id, home_id))
        world.trip_list = trip_list
        world.use_migration_rules = False
        self.assertEqual(world.trip_list.num_trips, 2)

        # run
        t1 = t0 + abm.days(1)
        sim.advance(t1)
        self.assertEqual(sim.result.get_num_time_points(), 25)

        # check effect of trips
        # self.assertEqual(p1.location_id, home_id) #person 1 is tested when goging to social event
        #self.assertEqual(p1.is_in_quarantine, True)
        # self.assertEqual(p2.location_id, work_id) #person 2 goes to work


if __name__ == '__main__':
    unittest.main()
