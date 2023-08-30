/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)   
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
#include "abm/migration_rules.h"
#include "abm/person.h"
#include "abm/location.h"
#include "abm/random_events.h"
#include "abm/location.h"
#include "memilio/utils/random_number_generator.h"
#include "abm/location_type.h"

#include <random>

namespace mio
{
namespace abm
{

LocationType random_migration(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params)
{
    auto current_loc     = person.get_location().get_type();
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

LocationType go_to_school(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params)
{
    auto current_loc = person.get_location().get_type();

    if (current_loc == LocationType::Home && t < params.get<LockdownDate>() && t.day_of_week() < 5 &&
        person.get_go_to_school_time(params) >= t.time_since_midnight() &&
        person.get_go_to_school_time(params) < t.time_since_midnight() + dt && person.get_age() == AgeGroup::Age5to14 &&
        person.goes_to_school(t, params) && !person.is_in_quarantine()) {
        return LocationType::School;
    }
    //return home
    if (current_loc == LocationType::School && t.hour_of_day() >= 15) {
        return LocationType::Home;
    }
    return current_loc;
}

LocationType go_to_work(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params)
{
    auto current_loc = person.get_location().get_type();

    if (current_loc == LocationType::Home && t < params.get<LockdownDate>() &&
        (person.get_age() == AgeGroup::Age15to34 || person.get_age() == AgeGroup::Age35to59) && t.day_of_week() < 5 &&
        t.time_since_midnight() + dt > person.get_go_to_work_time(params) &&
        t.time_since_midnight() <= person.get_go_to_work_time(params) && person.goes_to_work(t, params) &&
        !person.is_in_quarantine()) {
        return LocationType::Work;
    }
    //return home
    if (current_loc == LocationType::Work && t.hour_of_day() >= 17) {
        return LocationType::Home;
    }
    return current_loc;
}

LocationType go_to_shop(const Person& person, TimePoint t, TimeSpan dt, const MigrationParameters& params)
{
    auto current_loc = person.get_location().get_type();
    //leave
    if (t.day_of_week() < 6 && t.hour_of_day() > 7 && t.hour_of_day() < 22 && current_loc == LocationType::Home &&
        !person.is_in_quarantine()) {
        return random_transition(current_loc, dt,
                                 {{LocationType::BasicsShop, params.get<BasicShoppingRate>()[person.get_age()]}});
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
    if (current_loc == LocationType::Home && t < params.get<LockdownDate>() &&
        ((t.day_of_week() <= 4 && t.hour_of_day() >= 19) || (t.day_of_week() >= 5 && t.hour_of_day() >= 10)) &&
        !person.is_in_quarantine()) {
        return random_transition(current_loc, dt,
                                 {{LocationType::SocialEvent,
                                   params.get<SocialEventRate>().get_matrix_at(t.days())[(size_t)person.get_age()]}});
    }

    //return home
    if (current_loc == LocationType::SocialEvent && t.hour_of_day() >= 20 &&
        person.get_time_at_location() >= hours(2)) {
        return LocationType::Home;
    }

    return current_loc;
}

LocationType go_to_quarantine(const Person& person, TimePoint /*t*/, TimeSpan /*dt*/,
                              const MigrationParameters& /*params*/)
{
    auto current_loc = person.get_location().get_type();
    if (person.is_in_quarantine() && current_loc != LocationType::Hospital && current_loc != LocationType::ICU) {
        return LocationType::Home;
    }
    return current_loc;
}

LocationType go_to_hospital(const Person& person, const TimePoint t, TimeSpan /*dt*/,
                            const MigrationParameters& /*params*/)
{
    auto current_loc = person.get_location().get_type();
    if (person.get_infection_state(t) == InfectionState::InfectedSevere) {
        return LocationType::Hospital;
    }
    return current_loc;
}

LocationType go_to_icu(const Person& person, const TimePoint t, TimeSpan /*dt*/, const MigrationParameters& /*params*/)
{
    auto current_loc = person.get_location().get_type();
    if (person.get_infection_state(t) == InfectionState::InfectedCritical) {
        return LocationType::ICU;
    }
    return current_loc;
}

LocationType return_home_when_recovered(const Person& person, const TimePoint t, TimeSpan /*dt*/,
                                        const MigrationParameters& /*params*/)
{
    auto current_loc = person.get_location().get_type();
    if ((current_loc == LocationType::Hospital || current_loc == LocationType::ICU) &&
        person.get_infection_state(t) == InfectionState::Recovered) {
        return LocationType::Home;
    }
    return current_loc;
}

LocationType get_buried(const Person& person, const TimePoint t, TimeSpan /*dt*/,
                            const MigrationParameters& /*params*/)
{
    auto current_loc = person.get_location().get_type();
    if (person.get_infection_state(t) == InfectionState::Dead) {
        return LocationType::Cemetery;
    }
    return current_loc;
}

} // namespace abm
} // namespace mio
