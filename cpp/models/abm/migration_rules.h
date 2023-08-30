/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Daniel Abele, Majid Abedi, Elisabeth Kluth, Khoa Nguyen
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
#ifndef EPI_ABM_MIGRATION_RULES_H
#define EPI_ABM_MIGRATION_RULES_H

#include "abm/location_type.h"
#include "abm/parameters.h"
#include "abm/time.h"

namespace mio
{
namespace abm
{

class Person;

/**
 * @name Rules for migration between Location%s.
 * @param[in] p Person the rule is applied to.
 * @param[in] t Current time.
 * @param[in] dt Length of the time step.
 * @param[in] params Migration parameters.
 * @return Location that the Person migrates to if the rule is applied, the current Location of the person 
 * if the rule is not applied because of age, time, etc.
 * 
 * @{
 */
/**
 * @brief Completely random migration to any other Location.
 */
LocationType random_migration(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * @brief School age children go to school in the morning and return later in the day.
 */
LocationType go_to_school(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/** 
 * @brief Adults may go shopping in their free time.
 */
LocationType go_to_shop(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * @brief Person%s might go to social events.
 */
LocationType go_to_event(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * @brief Adults go to work in the morning and return later in the day.
 */
LocationType go_to_work(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * @brief Person%s who are in quarantine should go home.
 */
LocationType go_to_quarantine(const Person& person, TimePoint /*t*/, TimeSpan /*dt*/,
                              const MigrationParameters& /*params*/);

/**
 * @brief Infected Person%s may be hospitalized.
 */
LocationType go_to_hospital(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * @brief Person%s in the hospital may be put in intensive care.
 */
LocationType go_to_icu(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * @brief Person%s in the hospital/icu return home when they recover.
 */
LocationType return_home_when_recovered(const Person& person, TimePoint t, TimeSpan dt,
                                        const MigrationParameters& params);

/**
 * @brief Person%s in the icu go to cemetery when they are dead.
 */
LocationType get_buried(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params);
/**@}*/

} // namespace abm
} // namespace mio

#endif //EPI_ABM_MIGRATION_RULES_H
