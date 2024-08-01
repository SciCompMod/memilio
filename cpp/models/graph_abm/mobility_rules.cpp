/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Julia Bicker
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

#include "graph_abm/mobility_rules.h"
#include "abm/location.h"
#include "abm/infection_state.h"

namespace mio
{

abm::LocationType apply_commuting(const abm::Person& person, abm::TimePoint t, const abm::Parameters& params)
{
    abm::LocationType current_loc = person.get_location_type();

    if (current_loc == abm::LocationType::Home && params.get<abm::AgeGroupGotoWork>()[person.get_age()] &&
        t < params.get<abm::LockdownDate>() && t.day_of_week() < 5 && person.goes_to_work(t, params) &&
        !person.is_in_quarantine(t, params)) {
        return abm::LocationType::Work;
    }

    //person is at hospital in non-home model
    if (current_loc == abm::LocationType::Hospital &&
        person.get_location_model_id() != person.get_assigned_location_model_id(abm::LocationType::Home)) {
        //if person is still infected it stays at hospital
        if (person.is_infected(t)) {
            return current_loc;
        }
        //if person has recovered it goes home
        else if (person.get_infection_state(t) == abm::InfectionState::Recovered) {
            return abm::LocationType::Home;
        }
        //if person is dead it is sent to cementary
        else {
            return abm::LocationType::Cemetery;
        }
    }
    //person is at location in Home model (and should not go to work) it stays at that location
    if (person.get_location_model_id() == person.get_assigned_location_model_id(abm::LocationType::Home)) {
        return current_loc;
    }

    // agents are sent home or to work (if none of the above cases occurs) every time this function is called
    // i.e. if it is called too often they will be sent to work multiple times
    return abm::LocationType::Home;
}
} // namespace mio