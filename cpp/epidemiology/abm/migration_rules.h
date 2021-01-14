#ifndef EPI_ABM_MIGRATION_RULES_H
#define EPI_ABM_MIGRATION_RULES_H

#include "epidemiology/abm/location_type.h"
#include "epidemiology/abm/parameters.h"
#include "epidemiology/abm/time.h"

namespace epi
{

class Person;

LocationType random_migration(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);
LocationType go_to_school(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);
LocationType go_to_work(const Person& p, TimePoint t, TimeSpan dt, const MigrationParameters& params);

} // namespace epi

#endif //EPI_ABM_MIGRATION_RULES_H
