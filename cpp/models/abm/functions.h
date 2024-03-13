/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Rene Schmieding
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

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

// TODO: find a meaningfull header name

#include "abm/location.h"
#include "abm/person.h"

#include <vector>

namespace mio
{
namespace abm
{

/** 
 * @brief Add the contribution of a person to the local exposure rates.
 * @param[in, out] local_air_exposure Exposure rates by aerosols for the local population.
 * @param[in, out] local_contact_exposure Exposure by rates contacts for the local population.
 * @param[in] person A person from the local population.
 * @param[in] location The person's current location.
 * @param[in] t Current Simulation time.
 * @param[in] dt Length of the current Simulation time step.
 */
void add_exposure_contribution(AirExposureRates& local_air_exposure, ContactExposureRates& local_contact_exposure,
                               const Person& person, const Location& location, const TimePoint t, const TimeSpan dt);

/** 
 * @brief Let a Person interact with the population at its current Location, possibly getting infected.
 * @param[in, out] rng PersonalRandomNumberGenerator for this Person.
 * @param[in, out] person The person to interact with the local population.
 * @param[in] location The person's current location.
 * @param[in] local_air_exposure Precomputed exposure rates by aerosols for the local population.
 * @param[in] local_contact_exposure Precomputed exposure by rates contacts for the local population.
 * @param[in] t Current Simulation time.
 * @param[in] dt Length of the current Simulation time step.
 * @param[in] global_parameters Parameters of the World.
 */
void interact(PersonalRandomNumberGenerator& personal_rng, Person& person, const Location& location,
              const AirExposureRates& local_air_exposure, const ContactExposureRates& local_contact_exposure,
              const TimePoint t, const TimeSpan dt, const Parameters& global_parameters);

// interact, but it computes the correct exposures for you
void interact_testing(PersonalRandomNumberGenerator& personal_rng, Person& person, const Location& location,
                      const std::vector<Person>& local_population, const TimePoint t, const TimeSpan dt,
                      const Parameters& global_parameters);

/**
 * @brief Move a person to another location.
 * If the person already is at the destination, neither mode nor cells are set.
 * @param[in, out] person The person to be moved.
 * @param[in] destination The destination to move to.
 * @param[in] cells The cells within the destination the person should be in.
 * @param[in] mode The transport mode the person uses to move.
 * @return Returns false if the person already is at the given destination, true otherwise.
 */
bool migrate(Person& person, const Location& destination, const std::vector<uint32_t>& cells = {0},
             const TransportMode mode = TransportMode::Unknown);

} // namespace abm
} // namespace mio

#endif
