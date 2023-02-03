/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth
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
#ifndef EPI_ABM_LOCKDOWN_RULES_H
#define EPI_ABM_LOCKDOWN_RULES_H

#include "abm/time.h"
#include "abm/location_type.h"
#include "abm/person.h"
#include "abm/parameters.h"

#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/contact_matrix.h"

namespace mio
{
namespace abm
{

/**
 * @file LockdownRules implements non phamarceutical interventions via Dampings.
 * For interventions, Person%s are randomly divided into groups, e.g. one group works at home and the other group still
 * goes to work.
 * The probability with which a Person belongs to a certain group is time-dependent. This change
 * in probabilty is implemented by using Dampings.
 */

/**
 * @brief Person%s who are in home office are staying at home instead of going to work.
 * @param[in] t_begin Begin of the intervention.
 * @param[in] p Percentage of Person%s that work in home office.
 * @param[in, out] params Migration parameters that include Damping.
 */
void set_home_office(TimePoint t_begin, double p, MigrationParameters& params);

/**
 * @brief If schools are closed, students stay at home instead of going to school.
 * @param[in] t_begin Begin of the intervention.
 * @param[in] p Percentage of Person%s that are homeschooled.
 * @param[in,out] params Migration parameters that include Damping.
 */
void set_school_closure(TimePoint t_begin, double p, MigrationParameters& params);

/** 
 * @brief During lockdown Person%s join social events less often.
 * Whether a Person joins a social event is a random event (exponentially distributed).
 * The Damping changes the parameter of the exponential distribution, where a Damping of 0 corresponds to no Damping
 * and a Damping of 1 means that no social events are happening.
 * @param[in] t_begin Begin of the intervention.
 * @param[in] p Damping between 0 and 1 that changes the parameter of the exponential distribution.
 * @param[in,out] params Migration parameters that include Damping.
 */
void close_social_events(TimePoint t_begin, double p, MigrationParameters& params);

} // namespace abm
} //namespace mio

#endif // LOCKDOWN_RULES_H
