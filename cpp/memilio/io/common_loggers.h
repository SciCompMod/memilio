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

#ifndef COMMON_LOGGERS_H
#define COMMON_LOGGERS_H

#include "memilio/io/history.h"
#include "models/abm/location_type.h"

struct LogLocationInformation : mio::LogOnce {
    using Type = std::vector<std::tuple<uint32_t, mio::abm::GeographicalLocation>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<std::tuple<mio::abm::LocationType, uint32_t>> location_information{};
        for (auto&& location : sim.get_world().get_locations()) {
            location_information.push_back(std::make_tuple(location.get_index(), location.get_geographical_location()));
        }
        return location_information;
    }
};

struct LogPersonInformation : mio::LogOnce {
    using Type = std::vector<std::tuple<uint32_t, uint32_t, AgeGroup>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<std::tuple<uint32_t, unsigned, AgeGroup>> person_information{};
        for (auto&& person : sim.get_world().get_persons()) {
            person_information.push_back(
                std::make_tuple(person.get_person_id(),
                                sim.get_world().find_location(mio::abm::home, person).get_index(), person.get_age()));
        }
        return person_information;
    }
};

struct LogPersonsPerLocationAndInfectionTime : mio::LogAlways {
    using Type = std::vector<std::tuple<mio::abm::LocationId, uint32_t, mio::abm::TimeSpan, mio::abm::InfectionState>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<std::tuple<mio::abm::LocationId, uint32_t, mio::abm::TimeSpan, mio::abm::InfectionState>>
            location_ids_person{};
        for (auto&& person : sim.get_world().get_persons()) {
            location_ids_person.push_back(std::make_tuple(person.get_location().get_id(), person.get_person_id(),
                                                          person.get_time_since_transmission(),
                                                          person.get_infection_state(sim.get_time())));
        }
        return location_ids_person;
    }
};

#endif //COMMON_LOGGERS_H