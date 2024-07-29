/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Julia Bicker
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

#ifndef MIO_ABM_GRAPH_MOBILITY_RULES_H
#define MIO_ABM_GRAPH_MOBILITY_RULES_H

#include "abm/location_type.h"
#include "abm/time.h"
#include "abm/parameters.h"
#include "abm/person.h"

namespace mio
{

/**
 * @brief Once a day commuters go to work in another node.
 * @param[in] person Person the rule is applies to
 * @param[in] t Current time point
 * @param[in] params Parameters of person's Home world
 * @return LocationType the person is going to
 */
abm::LocationType apply_commuting(const abm::Person& person, abm::TimePoint t, const abm::Parameters& params);

} // namespace mio

#endif //MIO_ABM_GRAPH_MOBILITY_RULES_H