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
namespace mio
{
namespace abm
{
struct LogLocationInformation : mio::LogOnceStart {
    using Type = std::vector<std::tuple<uint32_t, mio::abm::GeographicalLocation, size_t, int>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        Type location_information{};
        for (auto&& location : sim.get_world().get_locations()) {
            auto n_cells     = location.get_cells().size();
            int loc_capacity = 0;
            for (int i = 0; i < (int)n_cells; i++) {
                loc_capacity += location.get_capacity(i).persons;
            }
            location_information.push_back(
                std::make_tuple(location.get_index(), location.get_geographical_location(), n_cells, loc_capacity));
        }
        return location_information;
    }
};

struct LogPersonInformation : mio::LogOnceStart {
    using Type = std::vector<std::tuple<uint32_t, uint32_t, mio::AgeGroup>>;
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

struct LogDataForMovement : mio::LogAlways {

    using Type = std::vector<std::tuple<uint32_t, uint32_t, mio::abm::TimePoint, mio::abm::TransportMode,
                                        mio::abm::ActivityType, mio::abm::InfectionState>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        Type movement_data{};
        for (Person p : sim.get_world().get_persons()) {
            movement_data.push_back(std::make_tuple(
                p.get_person_id(), p.get_location().get_index(), sim.get_time(), mio::abm::TransportMode::Unknown,
                mio::abm::ActivityType::UnknownActivity, p.get_infection_state(sim.get_time())));
        }
        return movement_data;
    }
};
} // namespace abm
} // namespace mio
#endif //ABM_COMMON_LOGGERS_H
