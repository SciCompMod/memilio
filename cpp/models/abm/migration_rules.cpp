/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Daniel Abele, Majid Abedi, Elisabeth Kluth
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

LocationType random_migration(const Person& person, TimePoint t, TimeSpan dt, const AbmMigrationParameters& params)
{
    auto current_loc     = person.get_location_id().type;
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

LocationType go_to_quarantine(const Person& person, TimePoint /*t*/, TimeSpan /*dt*/,
                              const AbmMigrationParameters& /*params*/)
{
    auto current_loc = person.get_location_id().type;
    if (person.is_in_quarantine()) {
        return LocationType::Home;
    }
    return current_loc;
}

LocationType go_to_hospital(const Person& person, TimePoint /*t*/, TimeSpan /*dt*/,
                            const AbmMigrationParameters& /*params*/)
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

} // namespace mio
