/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Sascha Korf
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

#ifndef ABM_COMMON_LOGGERS_H
#define ABM_COMMON_LOGGERS_H

#include "memilio/io/history.h"
#include "models/abm/location_type.h"
#include "abm/movement_data.h"
#include "abm/abm.h"

struct LogLocationInformation : mio::LogOnce {
    using Type = std::vector<std::tuple<uint32_t, mio::abm::GeographicalLocation, size_t, int>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        Type location_information{};
        for (auto&& location : sim.get_world().get_locations()) {
            auto n_cells     = location.get_cells().size();
            int loc_capacity = 0;
            for (int i = 0; i < n_cells; i++) {
                loc_capacity += location.get_capacity(i).persons;
            }
            location_information.push_back(
                std::make_tuple(location.get_index(), location.get_geographical_location(), n_cells, loc_capacity));
        }
        return location_information;
    }
};

struct LogPersonInformation : mio::LogOnce {
    using Type = std::vector<std::tuple<uint32_t, uint32_t, mio::abm::AgeGroup>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        Type person_information{};
        for (auto&& person : sim.get_world().get_persons()) {
            person_information.push_back(std::make_tuple(
                person.get_person_id(), sim.get_world().find_location(mio::abm::LocationType::Home, person).get_index(),
                person.get_age()));
        }
        return person_information;
    }
};

mio::abm::movement_data calculate_movement_data(const mio::abm::Simulation& sim)
{
    mio::abm::movement_data movement_data{};
    mio::unused(sim);
    // movement_data.agent_id        = sim.get_world().get_persons().get_person_id();
    // movement_data.from_id         = sim.get_world().find_location(mio::abm::LocationType::Home, person).get_index();
    // movement_data.to_id           = sim.get_world().find_location(mio::abm::LocationType::Home, person).get_index();
    // movement_data.start_time      = sim.get_world().get_time();
    // movement_data.end_time        = sim.get_world().get_time();
    // movement_data.transport_mode  = sim.get_world().get_persons().get_transport_mode();
    // movement_data.activity_type   = sim.get_world().get_persons().get_activity_type();
    // movement_data.infection_state = sim.get_world().get_persons().get_infection_state();
    return movement_data;
}
struct LogMovementData : mio::LogAlways {
    using Type = int;
    static Type log(const mio::abm::Simulation& sim)
    {
        return calculate_movement_data(sim).agent_id;
    }
};

#endif //ABM_COMMON_LOGGERS_H
