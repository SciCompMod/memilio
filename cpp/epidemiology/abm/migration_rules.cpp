#include "epidemiology/abm/migration_rules.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/abm/random_number_generator.h"
#include "epidemiology/abm/location.h"

#include <random>

namespace epi
{

LocationType random_migration(const Person& person, double t, double dt, const MigrationParameters& params)
{
    if (t < params.lockdown_date) {
        auto u = -std::log(std::uniform_real_distribution<double>()(thread_local_rng()));
        if (u < dt) {
            std::vector<LocationType> allowed_locations = {LocationType::Home, LocationType::School,
                                                           LocationType::Work};
            allowed_locations.erase(
                std::remove(allowed_locations.begin(), allowed_locations.end(), person.get_location().get_type()),
                allowed_locations.end());
            auto random_idx =
                std::uniform_int_distribution<size_t>(0, allowed_locations.size() - 1)(thread_local_rng());
            return allowed_locations[random_idx];
        }
    }
    return person.get_location().get_type();
}

int day_of_week(double t)
{
    //is assumed that the simulation starts from monday
    // Monday: 0
    // Sunday: 6
    int days = int(t) % 24;
    return days % 7;
}

int hour_of_day(double t)
{
    //is assumed that the simulation starts from 00:00
    return (int(t) % 24);
}

LocationType go_to_school(const Person& person, double t, double dt, const MigrationParameters& params)
{
    if (t < params.lockdown_date) {
        int dow = day_of_week(t);
        int hod = hour_of_day(t);
        //TODO: Improve it by changing the specific time to a range of times
        if (dow < 5 and hod == 8)
            return epi::LocationType::School;
    }
    return person.get_location().get_type();
}

} // namespace epi
