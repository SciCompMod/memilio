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
import memilio.simulation.abm as abm
import memilio.simulation as mio

import numpy as np
import random

num_age_groups = 4

model = abm.Model(num_age_groups)

# Set parameters

for age_group in range(num_age_groups):
    model.parameters.TimeExposedToNoSymptoms[abm.VirusVariant.Wildtype, AgeGroup(age_group)] = mio.AbstractParameterDistribution(mio.ParameterDistributionLogNormal(
        4., 1.))

model.parameters.AgeGroupGotoSchool[AgeGroup(1)] = True
model.parameters.AgeGroupGotoWork[AgeGroup(2)] = True
model.parameters.AgeGroupGotoWork[AgeGroup(3)] = True

for age in range(num_age_groups):
    model.parameters.InfectionProtectionFactor[abm.ProtectionType.GenericVaccine, AgeGroup(
        age), abm.VirusVariant.Wildtype] = mio.TimeSeriesFunctor(
        [[0, 0.0], [14, 0.67], [180, 0.4]])

    model.parameters.SeverityProtectionFactor[abm.ProtectionType.GenericVaccine, AgeGroup(
        age), abm.VirusVariant.Wildtype] = mio.TimeSeriesFunctor(
        [[0, 0.0], [14, 0.85], [180, 0.7]])

model.parameters.check_constraints()

# Set populations

n_households = 10

child = abm.HouseholdMember(num_age_groups)
child.age_weights[AgeGroup(0)] = 1.
child.age_weights[AgeGroup(1)] = 1.

parent = abm.HouseholdMember(num_age_groups)
parent.age_weights[AgeGroup(2)] = 1.
parent.age_weights[AgeGroup(3)] = 1.

twoPersonHousehold_group = abm.HouseholdGroup()
twoPersonHousehold_full = abm.Household()
twoPersonHousehold_full.add_members(child, 1)
twoPersonHousehold_full.add_members(parent, 1)
twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households)
abm.add_household_group_to_model(model, twoPersonHousehold_group)

threePersonHousehold_group = abm.HouseholdGroup()
threePersonHousehold_full = abm.Household()
threePersonHousehold_full.add_members(child, 1)
threePersonHousehold_full.add_members(parent, 2)
threePersonHousehold_group.add_households(
    threePersonHousehold_full, n_households)
abm.add_household_group_to_model(model, threePersonHousehold_group)

# Set locations

event = model.add_location(abm.LocationType.SocialEvent)
model.get_location(event).infection_parameters.MaximumContacts = 5

hospital = model.add_location(abm.LocationType.Hospital)
model.get_location(hospital).infection_parameters.MaximumContacts = 5
icu = model.add_location(abm.LocationType.ICU)
model.get_location(icu).infection_parameters.MaximumContacts = 5

shop = model.add_location(abm.LocationType.BasicsShop)
model.get_location(shop).infection_parameters.MaximumContacts = 20

school = model.add_location(abm.LocationType.School)
model.get_location(school).infection_parameters.MaximumContacts = 20

work = model.add_location(abm.LocationType.Work)
model.get_location(work).infection_parameters.MaximumContacts = 20

model.parameters.AerosolTransmissionRates[abm.VirusVariant.Wildtype] = 10

contacts = np.zeros((num_age_groups, num_age_groups))
contacts[2, 3] = 10

model.get_location(
    work).infection_parameters.ContactRates.baseline = contacts

# Testing Schemes

validity_period = abm.days(1)
probability = 0.5
start_date = abm.TimePoint(0)
end_date = abm.TimePoint(0) + abm.days(10)
test_type = abm.TestType.Antigen
test_parameters = model.parameters.TestData[test_type]

testing_criteria_work = abm.TestingCriteria()
testing_scheme_work = abm.TestingScheme(
    testing_criteria_work, validity_period, start_date, end_date, test_parameters, probability)

model.testing_strategy.add_scheme(
    abm.LocationType.Work, testing_scheme_work)

# Seed infections

infection_distribution = [0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.0, 0.0]
for person in model.persons:
    infection_state = abm.InfectionState(
        mio.DiscreteDistribution.get_instance()(model.rng, infection_distribution))
    prng = abm.PersonalRandomNumberGenerator(model.rng, person)
    if infection_state != abm.InfectionState.Susceptible:
        person.add_new_infection(mio.abm.Infection(
            prng, abm.VirusVariant.Wildtype, person.age, model.parameters, start_date, infection_state))

# Assign locations

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

# Vaccinations

vacc_rate = 0.7
vaccination_priority = [AgeGroup(3), AgeGroup(2), AgeGroup(1)]
vaccination_time = start_date - abm.days(20)

persons_by_age = [[] for _ in range(num_age_groups)]
for idx, person in enumerate(model.persons):
    persons_by_age[person.age.get()].append(idx)

for age in vaccination_priority:
    indices = persons_by_age[age.get()]

    random.shuffle(indices)

    temp = vacc_rate * len(indices)
    n_to_vaccinate = int(np.round(vacc_rate * len(indices)))

    count = 0
    for i in range(n_to_vaccinate):
        person = model.persons[indices[i]]
        if person.get_infection_state(vaccination_time) == abm.InfectionState.Susceptible:
            person.add_new_vaccination(
                abm.ProtectionType.GenericVaccine, vaccination_time)

# Simulate

t_lockdown = start_date + abm.days(10)
abm.close_social_events(t_lockdown, 0.9, model.parameters)

t0 = start_date
tmax = t0 + abm.days(10)
sim = abm.Simulation(t0, model)


history = abm.TimeSeriesWriterLogInfectionStateHistory(
    mio.TimeSeries(len(abm.InfectionState.values())))

sim.advance(tmax, history)

history.get_log().print_table()
