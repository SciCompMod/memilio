#include "epidemiology/abm/migration_rules.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/abm/random_number_generator.h"
#include "epidemiology/abm/location.h"

#include <random>

namespace epi
{

LocationType random_migration(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params)
{
    if (t < params.lockdown_date) {
        auto u = -std::log(std::uniform_real_distribution<double>()(thread_local_rng()));
        if (u < dt.days()) {
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

LocationType go_to_school(const Person& person, TimePoint t, TimeSpan /*dt*/, const MigrationParameters& params)
{
    if (t < params.lockdown_date && t.day_of_week() < 5 && t.hour_of_day() >= 8 &&
        person.get_age() == AbmAgeGroup::Age5to14 && person.get_location().get_type() != epi::LocationType::School) {
        return epi::LocationType::School;
    }
    return person.get_location().get_type();
}

LocationType go_to_work(const Person& person, TimePoint t, TimeSpan /*dt*/, const MigrationParameters& params)
{
    if (t < params.lockdown_date &&
        (person.get_age() == AbmAgeGroup::Age15to34 || person.get_age() == AbmAgeGroup::Age35to59) &&
        t.day_of_week() < 5 && t.hour_of_day() >= 8 && person.get_location().get_type() != epi::LocationType::Work) {
        return epi::LocationType::Work;
    }
    return person.get_location().get_type();
}

} // namespace epi
