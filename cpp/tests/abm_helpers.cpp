/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: David Kerkmann, Khoa Nguyen
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "abm_helpers.h"
#include "abm/person.h"
#include "memilio/utils/random_number_generator.h"

mio::abm::Person<double> make_test_person(mio::abm::Location<double>& location, mio::AgeGroup age,
                                          mio::abm::InfectionState infection_state, mio::abm::TimePoint t,
                                          mio::abm::Parameters<double> params)
{
    assert(age.get() < params.get_num_groups());
    auto rng           = mio::RandomNumberGenerator();
    mio::abm::Person p = mio::abm::Person(rng, location, age);
    if (infection_state != mio::abm::InfectionState::Susceptible) {
        auto rng_p = mio::abm::Person<double>::RandomNumberGenerator(rng, p);
        p.add_new_infection(
            mio::abm::Infection(rng_p, static_cast<mio::abm::VirusVariant>(0), age, params, t, infection_state));
    }
    return p;
}

mio::abm::Person<double>& add_test_person(mio::abm::World<double>& world, mio::abm::LocationId loc_id,
                                          mio::AgeGroup age, mio::abm::InfectionState infection_state,
                                          mio::abm::TimePoint t)
{
    assert(age.get() < world.parameters.get_num_groups());
    mio::abm::Person<double>& p = world.add_person(loc_id, age);
    if (infection_state != mio::abm::InfectionState::Susceptible) {
        auto rng_p = mio::abm::Person<double>::RandomNumberGenerator(world.get_rng(), p);
        p.add_new_infection(mio::abm::Infection(rng_p, static_cast<mio::abm::VirusVariant>(0), age, world.parameters, t,
                                                infection_state));
    }
    return p;
}
