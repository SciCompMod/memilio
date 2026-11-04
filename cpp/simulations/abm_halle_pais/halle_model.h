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

/// @brief Number of output channels: cumulative deaths, InfectedSevere, InfectedCritical.
constexpr size_t num_channels = 3;

/// @brief One simulation result: the daily value of every output channel.
struct Observable {
    std::vector<double> days{}; ///< Day offset from the start of the simulation.
    std::vector<std::array<double, num_channels>> values{}; ///< {deaths since t0, severe, critical} per day.
};

/**
 * @brief Build the Halle ABM.
 *
 * Builds the population and its locations from the Halle population file, sets the contact rates and the
 * fixed infection and PAIS parameters, and replays the individual histories of the cohort file onto the
 * agents, so that immunity and pre-existing PAIS are in their correct state at @p t0.
 *
 * @param[in] person_file CSV with columns age,home_id,school_id,event_id,shopping_id,work_id.
 * @param[in] contact_dir Directory with the German baseline contact matrices.
 * @param[in] history_file CSV with columns age_group, ih_infection_1_date to ih_infection_4_date and
 *                         ih_vaccine_1_date to ih_vaccine_3_date (YYYY-MM-DD or NA).
 * @param[in] theta Fitted parameters. theta[0] is the viral shedding rate; further entries are ignored.
 * @param[in] start_date Calendar date that @p t0 corresponds to.
 * @param[in] t0 Start of the simulation.
 * @param[in] rng Random number generator used for the whole model.
 */
IOResult<abm::Model> make_model(const std::string& person_file, const std::string& contact_dir,
                                const std::string& history_file, const std::vector<double>& theta, Date start_date,
                                abm::TimePoint t0, const RandomNumberGenerator& rng);

/**
 * @brief Build the model and run it once for @p num_days days.
 * @return The daily output channels, with deaths rebased onto the start of the simulation.
 */
IOResult<Observable> run_once(const std::string& person_file, const std::string& contact_dir,
                              const std::string& history_file, const std::vector<double>& theta, Date start_date,
                              int num_days);

} // namespace halle
} // namespace mio

#endif
