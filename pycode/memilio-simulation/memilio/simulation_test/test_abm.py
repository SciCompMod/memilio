#############################################################################
# Copyright (C) 2020-2024 MEmilio
#
# Authors: Daniel Abele, Khoa Nguyen
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

global num_age_groups
num_age_groups = 6


class TestAbm(unittest.TestCase):
    def test_model(self):
        t0 = abm.TimePoint(0)
        sim = abm.Simulation(t0, num_age_groups)
        model = sim.model
        self.assertEqual(len(model.persons), 0)
        self.assertEqual(len(model.locations), 1)

    def test_locations(self):
        t0 = abm.TimePoint(0)
        sim = abm.Simulation(t0, num_age_groups)
        model = sim.model

        home_id = model.add_location(abm.LocationType.Home)
        social_event_id = model.add_location(abm.LocationType.SocialEvent)
        self.assertEqual(len(model.locations), 3)

        home = model.locations[home_id.index]
        self.assertEqual(home.type, abm.LocationType.Home)

        testing_ages = [mio.AgeGroup(0)]

        home.infection_parameters.MaximumContacts = 10
        self.assertEqual(home.infection_parameters.MaximumContacts, 10)

        testing_inf_states = []
        testing_crit = abm.TestingCriteria(
            testing_ages, testing_inf_states)
        testing_scheme = abm.TestingScheme(testing_crit, abm.days(
            1), t0, t0 + abm.days(1), abm.AntigenTest(), 1.0)
        # initially false, will only active once simulation starts
        self.assertEqual(testing_scheme.active, False)

    def test_persons(self):
        t0 = abm.TimePoint(0)
        sim = abm.Simulation(t0, 6)
        model = sim.model

        home_id = model.add_location(abm.LocationType.Home)
        social_event_id = model.add_location(abm.LocationType.SocialEvent)

        p1 = model.add_person(
            home_id, mio.AgeGroup(2))
        p2 = model.add_person(
            social_event_id, mio.AgeGroup(5))

        # check persons
        self.assertEqual(len(model.persons), 2)
        self.assertEqual(p1.age, mio.AgeGroup(2))
        self.assertEqual(p1.location.index, 1)
        self.assertEqual(model.persons[0], p1)
        self.assertEqual(model.persons[1], p2)

    def test_simulation(self):
        t0 = abm.TimePoint(0)
        sim = abm.Simulation(t0, num_age_groups)
        model = sim.model

        # add some locations and persons
        for type in abm.LocationType.values():
            model.add_location(type)
        home_id = abm.LocationId(0, abm.LocationType.Home)
        social_event_id = abm.LocationId(0, abm.LocationType.SocialEvent)
        work_id = abm.LocationId(0, abm.LocationType.Work)
        p1 = model.add_person(
            home_id, mio.AgeGroup(0))
        p2 = model.add_person(
            home_id, mio.AgeGroup(2))
        for type in abm.LocationType.values():
            p1.set_assigned_location(abm.LocationId(0, type))
            p2.set_assigned_location(abm.LocationId(0, type))

        social_event = model.locations[social_event_id.index]

        model.parameters.InfectedSymptomsToSevere[abm.VirusVariant.Wildtype, mio.AgeGroup(
            0)] = 0.0
        model.parameters.InfectedSymptomsToRecovered[abm.VirusVariant.Wildtype, mio.AgeGroup(
            0)] = 0.0

        # trips
        trip_list = abm.TripList()
        trip_list.add_trip(abm.Trip(0, abm.TimePoint(
            0) + abm.hours(8), social_event_id, home_id))
        trip_list.add_trip(abm.Trip(1, abm.TimePoint(0) +
                           abm.hours(8), work_id, home_id))
        model.trip_list = trip_list
        model.use_movement_rules = False
        self.assertEqual(model.trip_list.num_trips(), 2)

        # vaccination
        vaccine = abm.Vaccination(
            abm.ExposureType.GenericVaccine, abm.TimePoint(0))
        self.assertEqual(vaccine.exposure_type,
                         abm.ExposureType.GenericVaccine)

        # run
        t1 = t0 + abm.days(1)
        sim.advance(t1)


if __name__ == '__main__':
    unittest.main()
