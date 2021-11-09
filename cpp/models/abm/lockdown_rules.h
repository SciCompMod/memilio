/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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

/**
 * LockdownRules implements non phamarceutical interventions via dampings.
 * For interventions, people are randomly divided into groups, e.g. one group works at home and the other group still goes to work.
 * The probability with which a person belongs to a certain group is time dependet. This change
 * in probabilty is implemented by using dampings.
 */

    
/**
 * Persons who are in home office are staying at home instead of going to work.
 * @param t_begin begin of the intervention
 * @param p percentage of people that work in home office
 * @param params migration parameters that include damping
 */
void set_home_office(TimePoint t_begin, double p, AbmMigrationParameters& params);
 
/**
 * If schools are closed, students stay at home instead of going to school.
 * @param t_begin begin of the intervention
 * @param p percentage of people that are homeschooled
 * @param params migration parameters that include damping
 */
void set_school_closure(TimePoint t_begin, double p, AbmMigrationParameters& params);



/** 
 * During lockdown people join social events less often.
 * If a person joins a social event is a random event (exponentially distributed).
 * The damping changes the parameter of the exponential distribution, where a damping of 0 corresponds to no damping
 * and a damping of 1 means that no social events are happening.
 * @param t_begin begin of the intervention
 * @param p damping between 0 and 1 that changes the parameter of the exponential distribution
 * @param params migration parameters that include damping
 */
void close_social_events(TimePoint t_begin, double p, AbmMigrationParameters& params);



} //namespace mio

#endif // LOCKDOWN_RULES_H