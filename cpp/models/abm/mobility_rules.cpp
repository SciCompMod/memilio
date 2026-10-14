/* 
* Copyright (C) 2020-2026 MEmilio   
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
#include "abm/mobility_rules.h"
#include "abm/activity_type.h"
#include "abm/person.h"
#include "abm/random_events.h"
#include "abm/location_type.h"
#include "abm/parameters.h"

namespace mio
{
namespace abm
{

ActivityType random_mobility(PersonalRandomNumberGenerator& rng, const Person& person, TimePoint t, TimeSpan dt,
                             const Parameters& params)
{
    auto current_loc     = person.get_activity_type();
    auto make_transition = [current_loc](auto l) {
        return std::make_pair(l, l == current_loc ? 0. : 1.);
    };
    if (t < params.get<LockdownDate>()) {
        return random_transition(rng, current_loc, dt,
                                 {make_transition(ActivityType::Work), make_transition(ActivityType::Home),
                                  make_transition(ActivityType::School), make_transition(ActivityType::Recreation),
                                  make_transition(ActivityType::BasicsShop)});
    }
    return current_loc;
}

ActivityType go_to_school(PersonalRandomNumberGenerator& /*rng*/, const Person& person, TimePoint t, TimeSpan dt,
                          const Parameters& params)
{
    auto current_loc = person.get_activity_type();

    if (current_loc == ActivityType::Home && t < params.get<LockdownDate>() && t.day_of_week() < 5 &&
        person.get_go_to_school_time(params) >= t.time_since_midnight() &&
        person.get_go_to_school_time(params) < t.time_since_midnight() + dt &&
        params.get<AgeGroupGotoSchool>()[person.get_age()] && person.goes_to_school(t, params) &&
        !person.is_in_quarantine(t, params)) {
        return ActivityType::School;
    }
    //return home
    if (current_loc == ActivityType::School && t.hour_of_day() >= 15) {
        return ActivityType::Home;
    }
    return current_loc;
}

ActivityType go_to_work(PersonalRandomNumberGenerator& /*rng*/, const Person& person, TimePoint t, TimeSpan dt,
                        const Parameters& params)
{
    auto current_loc = person.get_activity_type();

    if (current_loc == ActivityType::Home && t < params.get<LockdownDate>() &&
        params.get<AgeGroupGotoWork>()[person.get_age()] && t.day_of_week() < 5 &&
        t.time_since_midnight() + dt > person.get_go_to_work_time(params) &&
        t.time_since_midnight() <= person.get_go_to_work_time(params) && person.goes_to_work(t, params) &&
        !person.is_in_quarantine(t, params)) {
        return ActivityType::Work;
    }
    //return home
    if (current_loc == ActivityType::Work && t.hour_of_day() >= 17) {
        return ActivityType::Home;
    }
    return current_loc;
}

ActivityType go_to_shop(PersonalRandomNumberGenerator& rng, const Person& person, TimePoint t, TimeSpan dt,
                        const Parameters& params)
{
    auto current_loc = person.get_activity_type();
    //leave
    if (t.day_of_week() < 6 && t.hour_of_day() > 7 && t.hour_of_day() < 22 && current_loc == ActivityType::Home &&
        !person.is_in_quarantine(t, params)) {
        return random_transition(rng, current_loc, dt,
                                 {{ActivityType::BasicsShop, params.get<BasicShoppingRate>()[person.get_age()]}});
    }

    //return home
    if (current_loc == ActivityType::BasicsShop && person.get_time_at_location() >= hours(1)) {
        return ActivityType::Home;
    }

    return current_loc;
}

ActivityType go_to_recreation(PersonalRandomNumberGenerator& rng, const Person& person, TimePoint t, TimeSpan dt,
                              const Parameters& params)
{
    auto current_loc = person.get_activity_type();
    //leave
    if (current_loc == ActivityType::Home && t < params.get<LockdownDate>() &&
        ((t.day_of_week() <= 4 && t.hour_of_day() >= 19 && t.hour_of_day() < 22) ||
         (t.day_of_week() >= 5 && t.hour_of_day() >= 10 && t.hour_of_day() < 22)) &&
        !person.is_in_quarantine(t, params)) {
        return random_transition(
            rng, current_loc, dt,
            {{ActivityType::Recreation, params.get<SocialEventRate>().get_matrix_at(
                                            SimulationTime<ScalarType>(t.days()))[(size_t)person.get_age()]}});
    }

    //return home
    if (current_loc == ActivityType::Recreation && t.hour_of_day() >= 20 && person.get_time_at_location() >= hours(2)) {
        return ActivityType::Home;
    }

    return current_loc;
}

ActivityType go_to_quarantine(PersonalRandomNumberGenerator& /*rng*/, const Person& person, TimePoint t,
                              TimeSpan /*dt*/, const Parameters& params)
{
    auto current_loc = person.get_activity_type();
    if (person.is_in_quarantine(t, params) && current_loc != ActivityType::Hospital &&
        current_loc != ActivityType::ICU) {
        return ActivityType::Home;
    }
    return current_loc;
}

ActivityType go_to_hospital(PersonalRandomNumberGenerator& /*rng*/, const Person& person, const TimePoint t,
                            TimeSpan /*dt*/, const Parameters& /*params*/)
{
    auto current_loc = person.get_activity_type();
    if (person.get_infection_state(t) == InfectionState::InfectedSevere) {
        return ActivityType::Hospital;
    }
    return current_loc;
}

ActivityType go_to_icu(PersonalRandomNumberGenerator& /*rng*/, const Person& person, const TimePoint t, TimeSpan /*dt*/,
                       const Parameters& /*params*/)
{
    auto current_loc = person.get_activity_type();
    if (person.get_infection_state(t) == InfectionState::InfectedCritical) {
        return ActivityType::ICU;
    }
    return current_loc;
}

ActivityType return_home_when_recovered(PersonalRandomNumberGenerator& /*rng*/, const Person& person, const TimePoint t,
                                        TimeSpan /*dt*/, const Parameters& /*params*/)
{
    auto current_loc = person.get_activity_type();
    if ((current_loc == ActivityType::Hospital || current_loc == ActivityType::ICU) &&
        person.get_infection_state(t) == InfectionState::Recovered) {
        return ActivityType::Home;
    }
    return current_loc;
}

ActivityType get_buried(PersonalRandomNumberGenerator& /*rng*/, const Person& person, const TimePoint t,
                        TimeSpan /*dt*/, const Parameters& /*params*/)
{
    auto current_loc = person.get_activity_type();
    if (person.get_infection_state(t) == InfectionState::Dead) {
        return ActivityType::Cemetery;
    }
    return current_loc;
}

} // namespace abm
} // namespace mio
