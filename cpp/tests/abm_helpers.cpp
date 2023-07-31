/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: David Kerkmann
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

mio::abm::Person make_test_person(mio::abm::Location& location, mio::abm::AgeGroup age,
                                  mio::abm::InfectionState infection_state, mio::abm::TimePoint t,
                                  mio::abm::GlobalInfectionParameters params)
{
    mio::abm::Person p = mio::abm::Person(location, age);
    if (infection_state != mio::abm::InfectionState::Susceptible) {
        p.add_new_infection(
            mio::abm::Infection(static_cast<mio::abm::VirusVariant>(0), age, params, t, infection_state));
    }
    return p;
}

mio::abm::Person& add_test_person(mio::abm::World& world, mio::abm::LocationId loc_id, mio::abm::AgeGroup age,
                                  mio::abm::InfectionState infection_state, mio::abm::TimePoint t)
{
    mio::abm::Person& p = world.add_person(loc_id, age);
    if (infection_state != mio::abm::InfectionState::Susceptible) {
        p.add_new_infection(mio::abm::Infection(static_cast<mio::abm::VirusVariant>(0), age,
                                                world.get_global_infection_parameters(), t, infection_state));
    }
    return p;
}
