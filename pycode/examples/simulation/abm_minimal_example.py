#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Carlotta Gerstein
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

from memilio.simulation import AgeGroup
from memilio.simulation.abm import Model, VirusVariant
import memilio.simulation as mio

import numpy as np
from typing import Tuple

num_age_groups = 4

model = Model(num_age_groups)

# Set parameters

for age_group in range(num_age_groups):
    model.parameters.TimeExposedToNoSymptoms[VirusVariant.Wildtype, mio.AgeGroup(age_group)] = mio.AbstractParameterDistribution(mio.ParameterDistributionLogNormal(
        4., 1.))

    model.parameters.AgeGroupGotoSchool[AgeGroup(age_group)] = False
    model.parameters.AgeGroupGotoWork[AgeGroup(age_group)] = False

model.parameters.AgeGroupGotoSchool[AgeGroup(1)] = True
model.parameters.AgeGroupGotoWork[AgeGroup(2)] = True
model.parameters.AgeGroupGotoWork[AgeGroup(3)] = True

model.parameters.check_constraints()

# Set populations

n_households = 10

child = mio.abm.HouseholdMember(num_age_groups)
child.age_weights[AgeGroup(0)] = 1.
child.age_weights[AgeGroup(1)] = 1.

parent = mio.abm.HouseholdMember(num_age_groups)
parent.age_weights[AgeGroup(2)] = 1.
parent.age_weights[AgeGroup(3)] = 1.

twoPersonHousehold_group = mio.abm.HouseholdGroup()
twoPersonHousehold_full = mio.abm.Household()
twoPersonHousehold_full.add_members(child, 1)
twoPersonHousehold_full.add_members(parent, 1)
twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households)
mio.abm.add_household_group_to_model(model, twoPersonHousehold_group)

threePersonHousehold_group = mio.abm.HouseholdGroup()
threePersonHousehold_full = mio.abm.Household()
threePersonHousehold_full.add_members(child, 1)
threePersonHousehold_full.add_members(parent, 2)
threePersonHousehold_group.add_households(
    threePersonHousehold_full, n_households)
mio.abm.add_household_group_to_model(model, threePersonHousehold_group)

# Set locations

event = model.add_location(mio.abm.LocationType.SocialEvent)
model.get_location(event).infection_parameters.MaximumContacts = 5

hospital = model.add_location(mio.abm.LocationType.Hospital)
model.get_location(hospital).infection_parameters.MaximumContacts = 5
icu = model.add_location(mio.abm.LocationType.ICU)
model.get_location(icu).infection_parameters.MaximumContacts = 5

shop = model.add_location(mio.abm.LocationType.BasicsShop)
model.get_location(shop).infection_parameters.MaximumContacts = 20

school = model.add_location(mio.abm.LocationType.School)
model.get_location(school).infection_parameters.MaximumContacts = 20

work = model.add_location(mio.abm.LocationType.Work)
model.get_location(work).infection_parameters.MaximumContacts = 20

model.parameters.AerosolTransmissionRates[VirusVariant.Wildtype] = 10

contacts = np.zeros((num_age_groups, num_age_groups))
contacts[2, 3] = 10

model.get_location(
    work).infection_parameters.ContactRates.baseline = contacts

# Testing Schemes

validity_period = mio.abm.days(1)
probability = 0.5
start_date = mio.abm.TimePoint(0)
end_date = mio.abm.TimePoint(0) + mio.abm.days(10)
test_type = mio.abm.TestType.Antigen
test_parameters = model.parameters.TestData[test_type]

testing_criteria_work = mio.abm.TestingCriteria()
testing_scheme_work = mio.abm.TestingScheme(
    testing_criteria_work, validity_period, start_date, end_date, test_parameters, probability)

model.testing_strategy.add_scheme(
    mio.abm.LocationType.Work, testing_scheme_work)

# Seed infections

infection_distribution = [0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.0, 0.0]
rng = np.random.default_rng()
for person in model.persons:
    infection_state = mio.abm.InfectionState(rng.choice(
        len(infection_distribution), p=infection_distribution))
    prng = mio.abm.PersonalRandomNumberGenerator(person)

    if infection_state != mio.abm.InfectionState.Susceptible:
        person.add_new_infection(prng, mio.abm.VirusVariant.Wildtype,
                                 person.age, model.parameters, start_date, infection_state)

for person in model.persons:
    id = person.id

    model.assign_location(id, event)
    model.assign_location(id, shop)

    model.assign_location(id, hospital)
    model.assign_location(id, icu)

    if person.age == AgeGroup(1):
        model.assign_location(id, school)

    if person.age == AgeGroup(2) or person.age == AgeGroup(3):
        model.assign_location(id, work)

t_lockdown = mio.abm.TimePoint(0) + mio.abm.days(10)
mio.abm.close_social_events(t_lockdown, 0.9, model.parameters)

t0 = mio.abm.TimePoint(0)
tmax = t0 + mio.abm.days(10)
sim = mio.abm.Simulation(t0, model)


history = mio.abm.History(mio.TimeSeries(num_age_groups))

sim.advance(tmax, history)

for person in sim.model.persons:
    # print("start_date: ", person.get_infection_state(start_date))
    if (person.get_infection_state(start_date) == mio.abm.InfectionState.Susceptible):

        print("end_date: ", person.get_infection_state(end_date))

history.get_log().print_table()
