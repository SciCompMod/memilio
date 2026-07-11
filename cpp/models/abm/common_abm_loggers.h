/*
* Copyright (C) 2020-2026 MEmilio
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

#include "abm/infection_state.h"
#include "abm/person_id.h"
#include "abm/simulation.h"
#include "memilio/io/history.h"
#include "memilio/utils/time_series.h"
#include "models/abm/location_type.h"
#include "abm/mobility_data.h"
#include "memilio/utils/mioomp.h"

namespace mio
{
namespace abm
{

/**
 * @brief Struct to save specific mobility data of an agent.
 * The data consists of:
 *
 */
struct mobility_data {
    uint32_t agent_id;
    uint32_t from_id;
    uint32_t to_id;
    mio::abm::TimePoint start_time;
    mio::abm::TimePoint end_time;
    mio::abm::TransportMode transport_mode;
    mio::abm::ActivityType activity_type;
    mio::abm::InfectionState infection_state;
};

constexpr mio::abm::ActivityType guess_activity_type(mio::abm::LocationType current_location)
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
    using Type = std::vector<
        std::tuple<mio::abm::LocationId, mio::abm::LocationType, mio::geo::GeographicalLocation, size_t, int>>;
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
    static Type log(const mio::abm::Simulation<>& sim)
    {
        Type location_information{};
        for (auto& location : sim.get_model().get_locations()) {
            auto n_cells     = location.get_cells().size();
            int loc_capacity = 0;
            for (int i = 0; i < (int)n_cells; i++) {
                loc_capacity += location.get_capacity(i).persons;
            }
            location_information.push_back(std::make_tuple(
                location.get_id(), location.get_type(), location.get_geographical_location(), n_cells, loc_capacity));
        }
        return location_information;
    }
};

/**
 * @brief Logger to log the Person%s Information in the simulation.
 */
struct LogPersonInformation : mio::LogOnce {
    using Type = std::vector<std::tuple<mio::abm::PersonId, mio::abm::LocationId, mio::AgeGroup>>;
    /**
     * @brief Log the LocationInformation of the simulation.
     * @param[in] sim The simulation of the abm.
     * @return A vector of tuples with the LocationInformation, where each tuple contains the following information:
     * -# The person id.
     * -# The index of the home location.
     * -# The age group of the person.
    */
    static Type log(const mio::abm::Simulation<>& sim)
    {
        Type person_information{};
        person_information.reserve(sim.get_model().get_persons().size());
        for (auto& person : sim.get_model().get_persons()) {
            person_information.push_back(std::make_tuple(
                person.get_id(), sim.get_model().find_location(mio::abm::LocationType::Home, person.get_id()),
                person.get_age()));
        }
        return person_information;
    }
};

/**
 * @brief Logger to log mobility data of the agents in the simulation.
 */
struct LogDataForMobility : mio::LogAlways {
    using Type = std::vector<std::tuple<mio::abm::PersonId, mio::abm::LocationId, mio::abm::TimePoint,
                                        mio::abm::TransportMode, mio::abm::ActivityType, mio::abm::InfectionState>>;
    /**
     * @brief Log the mobility data of the agents in the simulation.
     * @param[in] sim The simulation of the ABM.
     * @return A vector of tuples with the mobility Data, where each tuple contains the following information:
     * -# The person id.
     * -# The index of the location.
     * -# The time point.
     * -# The transport mode.
     * -# The activity type.
     * -# The infection state.
     */
    static Type log(const mio::abm::Simulation<>& sim)
    {
        Type mobility_data{};
        for (Person p : sim.get_model().get_persons()) {
            mobility_data.push_back(
                std::make_tuple(p.get_id(), p.get_location(), sim.get_time(), p.get_last_transport_mode(),
                                guess_activity_type(p.get_location_type()), p.get_infection_state(sim.get_time())));
        }
        return mobility_data;
    }
};

/**
* @brief Logger to log the TimeSeries of the number of Person%s in an #InfectionState.
*/
struct LogInfectionState : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorX<ScalarType>>;
    /**
     * @brief Log the TimeSeries of the number of Person%s in an #InfectionState.
     * @param[in] sim The simulation of the abm.
     * @return A pair of the TimePoint and the TimeSeries of the number of Person%s in an #InfectionState.
     */
    static Type log(const mio::abm::Simulation<>& sim)
    {

        Eigen::VectorX<ScalarType> sum =
            Eigen::VectorX<ScalarType>::Zero(Eigen::Index(mio::abm::InfectionState::Count));
        auto curr_time = sim.get_time();
        PRAGMA_OMP(for)
        for (auto& location : sim.get_model().get_locations()) {
            for (uint32_t inf_state = 0; inf_state < (int)mio::abm::InfectionState::Count; inf_state++) {
                sum[inf_state] += sim.get_model().get_subpopulation(location.get_id(), curr_time,
                                                                    mio::abm::InfectionState(inf_state));
            }
        }
        return std::make_pair(curr_time, sum);
    }
};

/**
 * @brief Iterate over all infected persons, bin them, and accumulate into a (num_bins x num_clusters) histogram.
 *
 * The cluster column is age_index * num_loc_types + location_index, using the person's
 * current location (where a test would occur at this time point). Persons in cemetery
 * are skipped for now. The result is flattened column-major to a single vector for
 * the TimeSeriesWriter: entry (bin, cluster) lives at index cluster * num_bins + bin.
 *
 * @tparam BinFn  Callable (ScalarType viral_load) -> int returning the (unclamped) bin index.
 * @param[in] sim       The ABM simulation.
 * @param[in] num_bins  Number of histogram bins along the value axis.
 * @param[in] bin_fn    Binning function mapping a viral load to a bin index.
 * @return Flattened (num_bins * num_clusters) histogram vector.
 */
template <class BinFn>
Eigen::VectorX<uint32_t> accumulate_infected_histogram(const mio::abm::Simulation<>& sim, int num_bins,
                                                         BinFn bin_fn)
{
    const int num_age_groups = static_cast<int>(sim.get_model().parameters.get_num_groups());
    const int num_loc_types = static_cast<int>(mio::abm::LocationType::Cemetery); // exclude Cemetery
    const int num_clusters  = num_age_groups * num_loc_types;
    const auto curr_time     = sim.get_time();

    Eigen::VectorX<uint32_t> hist = Eigen::VectorX<uint32_t>::Zero(num_bins * num_clusters);

    for (auto& person : sim.get_model().get_persons()) {
        const auto& infection_vector = person.get_infection_vector();
        if (infection_vector.empty()) {
            continue;
        }
        const int loc_idx = static_cast<int>(sim.get_model().get_location(person.get_location()).get_type());
        if (loc_idx >= num_loc_types) {
            continue; // skip Cemetery
        }
        const int age_idx = static_cast<int>(person.get_age().get());
        const int cluster = age_idx * num_loc_types + loc_idx;

        const ScalarType vl = infection_vector[0].get_viral_load(curr_time);
        const int bin       = std::clamp(bin_fn(vl), 0, num_bins - 1);

        hist[cluster * num_bins + bin] += uint32_t(1);
    }
    return hist;
}

/**
 * @brief Logger for the cluster-stratified distribution of viral load among infected persons.
 *
 * Output is a (num_bins x num_clusters) histogram flattened column-major; see
 * accumulate_infected_histogram for the cluster indexing convention.
 */
struct LogViralLoad : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorX<uint32_t>>;
    static Type log(const mio::abm::Simulation<>& sim)
    {
        constexpr int num_bins = 10;
        auto hist = accumulate_infected_histogram(sim, num_bins, [](ScalarType vl) {
            return static_cast<int>(vl);
        });
        return std::make_pair(sim.get_time(), hist);
    }
};

/**
 * @brief Logger for the cluster-stratified distribution of cycle thresholds among infected persons.
 *
 * Output is a (num_bins x num_clusters) histogram flattened column-major; see
 * accumulate_infected_histogram for the cluster indexing convention.
 */
struct LogCycleThreshhold : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorX<uint32_t>>;
    static Type log(const mio::abm::Simulation<>& sim)
    {
        constexpr int max_cycle_threshold = 40;
        constexpr int num_bins            = max_cycle_threshold + 1; // bins 0..40 inclusive
        constexpr double log10_cycle_slope = 3.32192809489; // 1/log10(2)
        auto hist = accumulate_infected_histogram(sim, num_bins, [](ScalarType vl) {
            return static_cast<int>(max_cycle_threshold - vl * log10_cycle_slope);
        });
        return std::make_pair(sim.get_time(), hist);
    }
};





/**
 * @brief Logs CT histograms for one or more (age-group, location-type) clusters per day.
 *
 * Each cluster produces an independent 41-bin histogram (CT 0..40).
 * Uninfected persons (and recovered with zero viral load) are binned at CT=40:
 * 40 = "not detected", 0 = maximum viral load.
 * The sum of each histogram equals the number of persons in that cluster that morning.
 * Fires at a fixed hour each day (default 8), matching LogSchoolCohort's schedule.
 *
 * Output Type: pair<TimePoint, vector<VectorX<uint32_t>>> — one histogram per cluster,
 * in the same order as the clusters passed to the constructor.
 */
struct LogCTCluster {
    using ClusterSpec = std::pair<mio::AgeGroup, mio::abm::LocationType>;
    using Type        = std::pair<mio::abm::TimePoint, std::vector<Eigen::VectorX<uint32_t>>>;

    explicit LogCTCluster(std::vector<ClusterSpec> clusters, int hour = 10)
        : m_clusters(std::move(clusters)), m_hour(hour)
    {
    }

    bool should_log(const mio::abm::Simulation<>& sim)
    {
        return sim.get_time().hour_of_day() == m_hour;
    }

    Type log(const mio::abm::Simulation<>& sim)
    {
        constexpr int num_bins       = 41;
        constexpr int max_ct         = 40;
        constexpr double log10_slope = 3.32192809489; // 1/log10(2)

        const auto t = sim.get_time();
        std::vector<Eigen::VectorX<uint32_t>> hists(m_clusters.size(),
                                                     Eigen::VectorX<uint32_t>::Zero(num_bins));

        for (const auto& person : sim.get_model().get_persons()) {
            const auto age      = person.get_age();
            const auto loc_type = sim.get_model().get_location(person.get_location()).get_type();
            for (size_t c = 0; c < m_clusters.size(); ++c) {
                if (m_clusters[c].first != age || m_clusters[c].second != loc_type)
                    continue;
                const auto& infections = person.get_infection_vector();
                if (infections.empty()) {
                    hists[c][max_ct] += 1;
                } else {
                    const ScalarType vl = infections[0].get_viral_load(t);
                    const int bin = std::clamp(static_cast<int>(max_ct - vl * log10_slope), 0, num_bins - 1);
                    hists[c][bin] += 1;
                }
            }
        }
        return {t, std::move(hists)};
    }

private:
    std::vector<ClusterSpec> m_clusters;
    int                      m_hour;
};

/**
 * @brief Logs CT values for a fixed cohort of school students once per day at 8AM.
 *
 * The cohort is sampled on the first 8AM when students are present at a School location.
 * Each subsequent 8AM the same individuals are re-observed.
 *
 * Encoding (matches LogCTCluster): 0 = max viral load, 40 = not detected (uninfected
 * or below threshold), 255 = unused budget slot (cohort smaller than max_budget).
 */
struct LogSchoolCohort {
    using Type = std::pair<mio::abm::TimePoint, std::vector<uint8_t>>;

    explicit LogSchoolCohort(int max_budget = 100) : m_max_budget(max_budget) {}

    bool should_log(const mio::abm::Simulation<>& sim)
    {
        return sim.get_time().hour_of_day() == 8;
    }

    Type log(const mio::abm::Simulation<>& sim)
    {
        const auto& persons = sim.get_model().get_persons();

        if (m_cohort.empty()) {
            std::vector<size_t> candidates;
            for (size_t i = 0; i < persons.size(); ++i) {
                if (sim.get_model().get_location(persons[i].get_location()).get_type() ==
                    mio::abm::LocationType::School) {
                    candidates.push_back(i);
                }
            }
            std::shuffle(candidates.begin(), candidates.end(), mio::thread_local_rng());
            size_t n = std::min(candidates.size(), size_t(m_max_budget));
            m_cohort.assign(candidates.begin(), candidates.begin() + n);
        }

        auto t = sim.get_time();
        std::vector<uint8_t> cts(m_max_budget, 255); // 255 = unused slot

        for (size_t i = 0; i < m_cohort.size(); ++i) {
            const auto& infections = persons[m_cohort[i]].get_infection_vector();
            if (infections.empty()) {
                cts[i] = static_cast<uint8_t>(s_max_ct); // 40 = not detected
            }
            else {
                ScalarType vl = infections[0].get_viral_load(t);
                int ct        = static_cast<int>(s_max_ct - vl * s_log10_slope);
                cts[i]        = static_cast<uint8_t>(std::clamp(ct, 0, s_max_ct));
            }
        }
        return {t, cts};
    }

private:
    static constexpr int    s_max_ct      = 40;
    static constexpr double s_log10_slope = 3.32192809489; // 1/log10(2)
    int                     m_max_budget;
    std::vector<size_t>     m_cohort;
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
     * @brief This function adds an entry to the TimeSeries consisting of the TimePoint and the value. The Loggers must return a touple with a TimePoint and a value of return type Eigen::VectorX<ScalarType>.
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
