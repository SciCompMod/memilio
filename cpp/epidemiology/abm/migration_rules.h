#ifndef EPI_ABM_MIGRATION_RULES_H
#define EPI_ABM_MIGRATION_RULES_H

#include "epidemiology/abm/location_type.h"
#include "epidemiology/abm/parameters.h"
#include "epidemiology/abm/time.h"

namespace epi
{

class Person;

/**
 * @name rules for migration between locations.
 * @param p person the rule is applied to.
 * @param t current time.
 * @param dt length of the time step.
 * @param paras migration parameters.
 * @return location that the person migrates to if the rule is applied, the current location of the person 
 * if the rule is not applied because of age, time, etc.
 * 
 * @{
 */
/**
 * completely random migration to any other location.
 */
LocationType random_migration(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * school age children go to school in the morning and return later in the day.
 */
LocationType go_to_school(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * working age adults go to work in the morning and return later in the day.
 */
LocationType go_to_work(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * people go to the shop outside work/school except on sunday.
 */
LocationType go_to_shop(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * people go to social events outside work/school.
 */
LocationType go_to_event(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * infected people may be hospitalized.
 */
LocationType go_to_hospital(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * people in the hospital may be put in intensive care.
 */
LocationType go_to_icu(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

/**
 * people in the hospital/icu return home when they recover.
 */
LocationType return_home_when_recovered(const Person& person, TimePoint t, TimeSpan dt,
                                        const MigrationParameters& params);
/**@}*/

} // namespace epi

#endif //EPI_ABM_MIGRATION_RULES_H
