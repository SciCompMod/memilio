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

#include "abm/movement_data.h"
#include "abm/abm.h"

struct LogLocationInformation : mio::LogOnce {
    using Type = std::vector<std::tuple<uint32_t, mio::abm::GeographicalLocation>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        Type location_information{};
        for (auto&& location : sim.get_world().get_locations()) {
            location_information.push_back(std::make_tuple(location.get_index(), location.get_geographical_location()));
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

struct LogMovementData : mio::LogAlways {
    using Type = std::vector<mio::abm::movement_data>;
    static Type log(const mio::abm::Simulation& sim)
    {
        return sim.get_world().get_movement_data();
    }
};

#endif //ABM_COMMON_LOGGERS_H