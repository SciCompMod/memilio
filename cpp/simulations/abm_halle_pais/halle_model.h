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
#ifndef MIO_SIMULATIONS_ABM_HALLE_PAIS_HALLE_MODEL_H
#define MIO_SIMULATIONS_ABM_HALLE_PAIS_HALLE_MODEL_H

#include "abm/model.h"
#include "abm/time.h"
#include "memilio/io/io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/random_number_generator.h"

#include <array>
#include <string>
#include <vector>

namespace mio
{
namespace halle
{

/// @brief Number of age groups used throughout this simulation.
constexpr size_t num_age_groups = 6;

/// @brief County id of Halle (Saale), Stadt.
constexpr int halle_county_id = 15002;

/// @brief Index of each output channel in Observable::values.
enum Channel : size_t
{
    channel_deaths     = 0, ///< Cumulative deaths since t0, at the end of the day.
    channel_admissions = 1, ///< New InfectedSevere (hospital admissions) during the day.
    channel_icu        = 2, ///< Persons InfectedCritical (ICU occupancy) at the end of the day.
    channel_reported   = 3, ///< Positive tests (reported infections) during the day.
};

/// @brief Number of output channels, see Channel.
constexpr size_t num_channels = 4;

/// @brief One simulation result: the daily value of every output channel, for the model population.
struct Observable {
    std::vector<double> days{}; ///< Day offset from the start of the simulation, 0 to num_days - 1.
    std::vector<std::array<double, num_channels>> values{}; ///< One entry per day, indexed by Channel.
};

/**
 * @brief Names of the fitted parameters, in the order expected in the parameter vector theta.
 *
 * {viral_shedding_rate, dark_figure, testing_probability_symptomatic, ratio_symptomatic_to_asymptomatic}:
 * - viral_shedding_rate: InfectionRateFromViralShed.
 * - dark_figure: Ratio of all to known infections in the 30 days before t0, >= 1. Scales the recent cohort
 *   infections, which sets the number of infections active at t0.
 * - testing_probability_symptomatic: Probability that a symptomatic Person is tested when entering a school,
 *   workplace, shop or social event.
 * - ratio_symptomatic_to_asymptomatic: The testing probability of all other Persons is the symptomatic one
 *   divided by this ratio, >= 1.
 */
const std::vector<std::string>& parameter_names();

/**
 * @brief Build the Halle ABM.
 *
 * Builds the population and its locations from the Halle population file, sets the contact rates, the
 * fixed infection and PAIS parameters and the testing strategy, and replays the individual histories of the
 * cohort file onto the agents, so that immunity and pre-existing PAIS are in their correct state at @p t0.
 *
 * @param[in] person_file CSV with columns age,home_id,school_id,event_id,shopping_id,work_id.
 * @param[in] contact_dir Directory with the German baseline contact matrices.
 * @param[in] history_file CSV with columns age_group, ih_infection_1_date to ih_infection_4_date and
 *                         ih_vaccine_1_date to ih_vaccine_3_date (YYYY-MM-DD or NA). A sample of the
 *                         population; every agent is given the history of a member drawn from it.
 * @param[in] theta Fitted parameters, see parameter_names().
 * @param[in] start_date Calendar date that @p t0 corresponds to.
 * @param[in] t0 Start of the simulation.
 * @param[in] tmax End of the simulation, until which the testing strategy is active.
 * @param[in] rng Random number generator used for the whole model.
 */
IOResult<abm::Model> make_model(const std::string& person_file, const std::string& contact_dir,
                                const std::string& history_file, const std::vector<double>& theta, Date start_date,
                                abm::TimePoint t0, abm::TimePoint tmax, const RandomNumberGenerator& rng);

/**
 * @brief Build the model and run it once for @p num_days days.
 * @return The daily output channels for the model population, unscaled. Deaths are counted from t0.
 */
IOResult<Observable> run_once(const std::string& person_file, const std::string& contact_dir,
                              const std::string& history_file, const std::vector<double>& theta, Date start_date,
                              int num_days);

} // namespace halle
} // namespace mio

#endif
