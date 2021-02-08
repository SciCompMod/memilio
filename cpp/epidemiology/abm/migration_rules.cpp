#include "epidemiology/abm/migration_rules.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/abm/random_number_generator.h"
#include "epidemiology/abm/random_events.h"
#include "epidemiology/abm/location.h"

#include <random>

namespace epi
{

LocationType random_migration(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params)
{
    auto current = person.get_location().get_type();
    auto make_transition = [current](auto l) {
        return std::make_pair(l, l == current ? 0. : 1.);
    };
    if (t < params.lockdown_date) {
        return random_transition(current, dt,
                                 {make_transition(LocationType::Work), make_transition(LocationType::Home),
                                  make_transition(LocationType::School), make_transition(LocationType::SocialEvent),
                                  make_transition(LocationType::BasicsShop)});
    }
    return current;
}

LocationType go_to_school(const Person& person, TimePoint t, TimeSpan /*dt*/, const MigrationParameters& params)
{
    if (person.get_location().get_type() == LocationType::Home && t < params.lockdown_date && t.day_of_week() < 5 &&
        t.hour_of_day() >= 8 && person.get_age() == AbmAgeGroup::Age5to14 &&
        person.get_location().get_type() != epi::LocationType::School) {
        return epi::LocationType::School;
    }
    //return home
    if (person.get_location().get_type() == epi::LocationType::School && t.hour_of_day() >= 15) {
        return epi::LocationType::Home;
    }
    return person.get_location().get_type();
}

LocationType go_to_work(const Person& person, TimePoint t, TimeSpan /*dt*/, const MigrationParameters& params)
{
    if (person.get_location().get_type() == LocationType::Home && t < params.lockdown_date &&
        (person.get_age() == AbmAgeGroup::Age15to34 || person.get_age() == AbmAgeGroup::Age35to59) &&
        t.day_of_week() < 5 && t.hour_of_day() >= 8 && person.get_location().get_type() != epi::LocationType::Work) {
        return epi::LocationType::Work;
    }
    //return home
    if (person.get_location().get_type() == epi::LocationType::Work && t.hour_of_day() >= 17) {
        return epi::LocationType::Home;
    }
    return person.get_location().get_type();
}

LocationType go_to_shop(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params)
{
    auto current_loc = person.get_location().get_type();
    //leave
    if (t.day_of_week() < 6 && t.hour_of_day() > 7 && t.hour_of_day() < 22 && current_loc == LocationType::Home) {
        return random_transition(current_loc, dt,
                                 {{LocationType::BasicsShop, params.basic_shopping_rate[person.get_age()]}});
    }

    //return home
    if (current_loc == LocationType::BasicsShop && person.get_time_at_location() >= hours(1)) {
        return LocationType::Home;
    }

    return current_loc;
}

LocationType go_to_event(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params)
{
    auto current_loc = person.get_location().get_type();
    //leave
    if (current_loc == LocationType::Home && t < params.lockdown_date &&
        ((t.day_of_week() <= 4 && t.hour_of_day() >= 19) || (t.day_of_week() >= 5 && t.hour_of_day() >= 10))) {
        return random_transition(current_loc, dt,
                                 {{LocationType::SocialEvent, params.social_event_rate[person.get_age()]}});
    }

    //return home
    if (current_loc == LocationType::SocialEvent && t.hour_of_day() > 20 && person.get_time_at_location() >= hours(2)) {
        return LocationType::Home;
    }

    return current_loc;
}

LocationType go_to_hospital(const Person& person, TimePoint /*t*/, TimeSpan dt, const MigrationParameters& params)
{
    auto current_loc = person.get_location().get_type();
    if (person.get_infection_state() == InfectionState::Infected_Detected) {
        return random_transition(current_loc, dt,
                                 {{LocationType::Hospital, params.hospitalization_rate[person.get_age()]}});
    }
    return current_loc;
}

LocationType go_to_icu(const Person& person, TimePoint /*t*/, TimeSpan dt, const MigrationParameters& params)
{
    auto current_loc = person.get_location().get_type();
    if (current_loc == LocationType::Hospital) {
        return random_transition(current_loc, dt,
                                 {{LocationType::IntensiveCareUnit, params.icu_rate[person.get_age()]}});
    }
    return current_loc;
}

LocationType return_home_when_recovered(const Person& person, TimePoint /*t*/, TimeSpan /*dt*/,
                                        const MigrationParameters& /*params*/)
{
    auto current_loc = person.get_location().get_type();
    if ((current_loc == LocationType::Hospital || current_loc == LocationType::IntensiveCareUnit) &&
        person.get_infection_state() == InfectionState::Recovered_Infected) {
        return LocationType::Home;
    }
    return current_loc;
}

} // namespace epi
