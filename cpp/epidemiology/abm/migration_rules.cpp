#include "epidemiology/abm/migration_rules.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/abm/random_events.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/utils/random_number_generator.h"
#include "epidemiology/abm/location_type.h"

#include <random>

namespace epi
{

LocationType random_migration(const Person& person, TimePoint t, TimeSpan dt, const AbmMigrationParameters& params)
{
    auto current_loc = person.get_location_id().type;
    auto make_transition = [current_loc](auto l) {
        return std::make_pair(l, l == current_loc ? 0. : 1.);
    };
    if (t < params.get<LockdownDate>()) {
        return random_transition(current_loc, dt,
                                 {make_transition(LocationType::Work), make_transition(LocationType::Home),
                                  make_transition(LocationType::School), make_transition(LocationType::SocialEvent),
                                  make_transition(LocationType::BasicsShop)});
    }
    return current_loc;
}

LocationType go_to_school(const Person& person, TimePoint t, TimeSpan /*dt*/, const AbmMigrationParameters& params)
{
    auto current_loc = person.get_location_id().type;
    if (current_loc == LocationType::Home && t < params.get<LockdownDate>() && t.day_of_week() < 5 &&
        t.hour_of_day() >= 8 && person.get_age() == AbmAgeGroup::Age5to14 && person.goes_to_school(t, params)
        && !person.is_in_quarantine()) {
        return epi::LocationType::School;
    }
    //return home
    if (current_loc == epi::LocationType::School && t.hour_of_day() >= 15) {
        return epi::LocationType::Home;
    }
    return current_loc;
}

LocationType go_to_work(const Person& person, TimePoint t, TimeSpan /*dt*/, const AbmMigrationParameters& params)
{
    auto current_loc = person.get_location_id().type;

    if (current_loc == LocationType::Home && t < params.get<LockdownDate>() &&
        (person.get_age() == AbmAgeGroup::Age15to34 || person.get_age() == AbmAgeGroup::Age35to59) &&
        t.day_of_week() < 5 && t.hour_of_day() >= 8 && person.goes_to_work(t, params) && !person.is_in_quarantine()) {
        return epi::LocationType::Work;
    }
    //return home
    if (current_loc == epi::LocationType::Work && t.hour_of_day() >= 17) {
        return epi::LocationType::Home;
    }
    return current_loc;
}

LocationType go_to_shop(const Person& person, TimePoint t, TimeSpan dt, const AbmMigrationParameters& params)
{
    auto current_loc = person.get_location_id().type;
    //leave
    if (t.day_of_week() < 6 && t.hour_of_day() > 7 && t.hour_of_day() < 22 && current_loc == LocationType::Home
        && !person.is_in_quarantine()) {
        return random_transition(current_loc, dt,
                                 {{LocationType::BasicsShop, params.get<BasicShoppingRate>()[person.get_age()]}});
    }

    //return home
    if (current_loc == LocationType::BasicsShop && person.get_time_at_location() >= hours(1)) {
        return LocationType::Home;
    }

    return current_loc;
}

LocationType go_to_event(const Person& person, TimePoint t, TimeSpan dt, const AbmMigrationParameters& params)
{
    auto current_loc = person.get_location_id().type;
    //leave
    if (current_loc == LocationType::Home && t < params.get<LockdownDate>() &&
        ((t.day_of_week() <= 4 && t.hour_of_day() >= 19) || (t.day_of_week() >= 5 && t.hour_of_day() >= 10))
         && !person.is_in_quarantine()) {
        return random_transition(current_loc, dt,
                                 {{LocationType::SocialEvent, params.get<SocialEventRate>().get_matrix_at(t.days())[(size_t)person.get_age()]}});
    }

    //return home
    if (current_loc == LocationType::SocialEvent && t.hour_of_day() >= 20 &&
        person.get_time_at_location() >= hours(2)) {
        return LocationType::Home;
    }

    return current_loc;
}

LocationType go_to_hospital(const Person& person, TimePoint /*t*/, TimeSpan /*dt*/, const AbmMigrationParameters& /*params*/)
{
    auto current_loc = person.get_location_id().type;
    if (person.get_infection_state() == InfectionState::Infected_Severe) {
        return LocationType::Hospital;
    }
    return current_loc;
}

LocationType go_to_icu(const Person& person, TimePoint /*t*/, TimeSpan /*dt*/, const AbmMigrationParameters& /*params*/)
{
    auto current_loc = person.get_location_id().type;
    if (person.get_infection_state() == InfectionState::Infected_Critical) {
        return LocationType::ICU;
    }
    return current_loc;
}

LocationType return_home_when_recovered(const Person& person, TimePoint /*t*/, TimeSpan /*dt*/,
                                        const AbmMigrationParameters& /*params*/)
{
    auto current_loc = person.get_location_id().type;
    if ((current_loc == LocationType::Hospital || current_loc == LocationType::ICU) &&
        person.get_infection_state() == InfectionState::Recovered_Infected) {
        return LocationType::Home;
    }
    return current_loc;
}

} // namespace epi
