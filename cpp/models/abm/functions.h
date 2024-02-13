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

// add the contribution of person to the exposure rates at location
void add_exposure_contribution(AirExposureRates& local_air_exposure, ContactExposureRates& local_contact_exposure,
                               const Person& person, const Location& location, const TimePoint t, const TimeSpan dt);

// let a person interact with a location for and at some time
void interact(PersonalRandomNumberGenerator& personal_rng, Person& person, const Location& location,
              const AirExposureRates& local_air_exposure, const ContactExposureRates& local_contact_exposure,
              const TimePoint t, const TimeSpan dt, const Parameters& global_parameters);

// interact, but it computes the correct exposures for you
void interact(PersonalRandomNumberGenerator& personal_rng, Person& person, const Location& location,
              const std::vector<Person>& local_population, const TimePoint t, const TimeSpan dt,
              const Parameters& global_parameters);

// move a person to another location. returns false if the person was at the target location already, true otherwise
bool migrate(Person& person, const Location& destination, const std::vector<uint32_t>& cells = {0},
             const TransportMode mode = TransportMode::Unknown);
} // namespace abm

} // namespace mio

#endif