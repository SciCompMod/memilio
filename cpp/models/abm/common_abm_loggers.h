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
#include "memilio/utils/mioomp.h"
namespace mio
{
namespace abm
{

/**
 * @brief Struct to save spcific movement data of an agent.
 * The Data consists of :
 * 
 */
struct movement_data {
    uint32_t agent_id;
    uint32_t from_id;
    uint32_t to_id;
    mio::abm::TimePoint start_time;
    mio::abm::TimePoint end_time;
    mio::abm::TransportMode transport_mode;
    mio::abm::ActivityType activity_type;
    mio::abm::InfectionState infection_state;
};

mio::abm::ActivityType guess_activity_type(mio::abm::LocationType current_location)
{
    switch (current_location) {
    case mio::abm::LocationType::Home:
        return mio::abm::ActivityType::Home;
    case mio::abm::LocationType::Work:
        return mio::abm::ActivityType::Workplace;
    case mio::abm::LocationType::School:
        return mio::abm::ActivityType::Education;
    case mio::abm::LocationType::SocialEvent:
        return mio::abm::ActivityType::Leisure;
    case mio::abm::LocationType::BasicsShop:
        return mio::abm::ActivityType::Shopping;
    case mio::abm::LocationType::ICU:
        return mio::abm::ActivityType::OtherActivity;
    case mio::abm::LocationType::Hospital:
        return mio::abm::ActivityType::OtherActivity;
    case mio::abm::LocationType::Cemetery:
        return mio::abm::ActivityType::OtherActivity;
    default:
        return mio::abm::ActivityType::UnknownActivity;
    }
}

/**
 * @brief Logger to log the LocationInformation of the simulation.
 */
struct LogLocationInformation : mio::LogOnce {
    using Type = std::vector<std::tuple<uint32_t, mio::abm::LocationType, mio::abm::GeographicalLocation, size_t, int>>;
    /**
     * @brief Log the LocationInformation of the simulation. 
     * @param[in] sim The simulation of the abm.
     * @return A vector of tuples with the LocationInformation, where each tuple contains the following information:
     * -# The index of the location.
     * -# The type of the location.
     * -# The geographical location of the location.
     * -# The number of cells in the location.
     * -# The capacity of the location.
    */
    static Type log(const mio::abm::Simulation& sim)
    {
        Type location_information{};
        for (auto&& location : sim.get_world().get_locations()) {
            auto n_cells     = location.get_cells().size();
            int loc_capacity = 0;
            for (int i = 0; i < (int)n_cells; i++) {
                loc_capacity += location.get_capacity(i).persons;
            }
            location_information.push_back(std::make_tuple(location.get_index(), location.get_type(),
                                                           location.get_geographical_location(), n_cells,
                                                           loc_capacity));
        }
        return location_information;
    }
};

/**
 * @brief Logger to log the Person%s Information in the simulation.
 */
struct LogPersonInformation : mio::LogOnce {
    using Type = std::vector<std::tuple<uint32_t, uint32_t, mio::AgeGroup>>;
    /** 
     * @brief Log the LocationInformation of the simulation. 
     * @param[in] sim The simulation of the abm.
     * @return A vector of tuples with the LocationInformation, where each tuple contains the following information:
     * -# The person id.
     * -# The index of the home location.
     * -# The age group of the person.
    */
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

/**
 * @brief Logger to log Movement Data of the agents in the simulation.
 */
struct LogDataForMovement : mio::LogAlways {
    using Type = std::vector<std::tuple<uint32_t, uint32_t, mio::abm::TimePoint, mio::abm::TransportMode,
                                        mio::abm::ActivityType, mio::abm::InfectionState>>;
    /** 
     * @brief Log the Movement Data of the agents in the simulation.
     * @param[in] sim The simulation of the abm.
     * @return A vector of tuples with the Movement Data, where each tuple contains the following information:
     * -# The person id.
     * -# The index of the location.
     * -# The time point.
     * -# The transport mode.
     * -# The activity type.
     * -# The infection state.
     */
    static Type log(const mio::abm::Simulation& sim)
    {
        Type movement_data{};
        for (Person p : sim.get_world().get_persons()) {
            movement_data.push_back(std::make_tuple(
                p.get_person_id(), p.get_location().get_index(), sim.get_time(), p.get_last_transport_mode(),
                guess_activity_type(p.get_location().get_type()), p.get_infection_state(sim.get_time())));
        }
        return movement_data;
    }
};

/**
* @brief Logger to log the TimeSeries of the number of Person%s in an #InfectionState.
*/
struct LogInfectionState : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    /** 
     * @brief Log the TimeSeries of the number of Person%s in an #InfectionState.
     * @param[in] sim The simulation of the abm.
     * @return A pair of the TimePoint and the TimeSeries of the number of Person%s in an #InfectionState.
     */
    static Type log(const mio::abm::Simulation& sim)
    {

        Eigen::VectorXd sum = Eigen::VectorXd::Zero(Eigen::Index(mio::abm::InfectionState::Count));
        auto curr_time      = sim.get_time();
        PRAGMA_OMP(for)
        for (auto&& location : sim.get_world().get_locations()) {
            for (uint32_t inf_state = 0; inf_state < (int)mio::abm::InfectionState::Count; inf_state++) {
                sum[inf_state] += location.get_subpopulation(curr_time, mio::abm::InfectionState(inf_state));
            }
        }
        return std::make_pair(curr_time, sum);
    }
};

/**
* @brief This is like the DataWriterToMemory, but it only logs time series data.
* @tparam Loggers The loggers that are used to log data. The loggers must return a touple with a TimePoint and a value.
*/
template <class... Loggers>
struct TimeSeriesWriter {
    using Data = std::tuple<mio::TimeSeries<ScalarType>>;
    template <class Logger>
    /**
     * @brief This function adds an entry to the TimeSeries consisting of the TimePoint and the value. The Loggers must return a touple with a TimePoint and a value of return type Eigen::VectorXd.
     * @param[in] t The data from the logger.
     * @param[in,out] data The data tuple.
    */
    static void add_record(const typename Logger::Type& t, Data& data)
    {
        std::get<index_of_type_v<Logger, Loggers...>>(data).add_time_point(t.first.days(), t.second);
    }
};

/**
* @brief This class writes data retrieved from loggers to memory. It can be used as the Writer template parameter for the History class.
* This specialization just saves the difference to the last saved data. Suitable when one wants to save huge data with a few changes.
* @tparam Loggers The loggers that are used to log data.
*/
template <class... Loggers>
struct DataWriterToMemoryDelta {
    using Data = std::tuple<std::vector<typename Loggers::Type>...>;
    template <class Logger>
    /**
     * @brief This function adds an entry to the data vector. 
     * @param[in] t The data from the logger.
     * @param[in,out] data The data tuple.
     * @details The data is only added if it differs from the last entry. For this we use the first entry as a reference for the current position.
    */
    static void add_record(const typename Logger::Type& t, Data& data)
    {

        if (std::get<index_of_type_v<Logger, Loggers...>>(data).size() > 0) {
            typename Logger::Type diff_vector{};
            auto& current_state_vec = std::get<index_of_type_v<Logger, Loggers...>>(data).front();
            for (auto i = 0; i < (int)current_state_vec.size(); i++) {
                if (std::get<1>(t[i]) != std::get<1>(current_state_vec[i])) {
                    std::get<1>(current_state_vec[i]) = std::get<1>(t[i]);
                    diff_vector.push_back(t[i]);
                }
            }
            std::get<index_of_type_v<Logger, Loggers...>>(data).push_back(diff_vector);
        }
        else {
            std::get<index_of_type_v<Logger, Loggers...>>(data).push_back(
                t); // We use the first entry as a reference for the current position.
            std::get<index_of_type_v<Logger, Loggers...>>(data).push_back(t);
        }
    }
};

} // namespace abm
} // namespace mio
#endif //ABM_COMMON_LOGGERS_H
