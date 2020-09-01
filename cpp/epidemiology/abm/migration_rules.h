#ifndef EPI_ABM_MIGRATION_RULES_H
#define EPI_ABM_MIGRATION_RULES_H

#include "epidemiology/abm/location_type.h"
#include "epidemiology/abm/parameters.h"

namespace epi
{

class Person;

LocationType random_migration(const Person& p, double t, double dt, const MigrationParameters& params);

} // namespace epi

#endif //EPI_ABM_MIGRATION_RULES_H