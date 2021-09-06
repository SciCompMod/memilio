/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/world.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/utils/random_number_generator.h"

namespace epi
{

Person::Person(LocationId id, InfectionState infection_state, VaccinationState vaccination_state, AbmAgeGroup age, const GlobalInfectionParameters& global_params)
    : m_location_id(id)
    , m_assigned_locations((uint32_t)LocationType::Count, INVALID_LOCATION_INDEX)
    , m_infection_state(infection_state)
    , m_vaccination_state(vaccination_state)
    , m_quarantine(false)
    , m_age(age)
    , m_time_at_location(std::numeric_limits<int>::max())
{
    m_random_workgroup = UniformDistribution<double>::get_instance()();
    m_random_schoolgroup = UniformDistribution<double>::get_instance()();
    if (m_infection_state == InfectionState::Infected_Detected) {
        m_quarantine = true;
    }
    if (m_infection_state == InfectionState::Exposed){
        m_time_until_carrier = hours(UniformIntDistribution<int>::get_instance()(0, int(global_params.get<IncubationPeriod>()[{m_age,m_vaccination_state}] * 24)));
    }
}


Person::Person(LocationId id, InfectionState infection_state, AbmAgeGroup age, const GlobalInfectionParameters& global_params)
    :  Person(id, infection_state, VaccinationState::Unvaccinated, age, global_params)
{
}

Person::Person(Location& location, InfectionState infection_state, AbmAgeGroup age, const GlobalInfectionParameters& global_params)
    : Person({location.get_index(), location.get_type()}, infection_state, VaccinationState::Unvaccinated, age, global_params)
{
}

void Person::interact(TimeSpan dt, const GlobalInfectionParameters& global_infection_params, Location& loc)
{
    auto infection_state     = m_infection_state;
    auto new_infection_state = infection_state;

    if (infection_state == InfectionState::Exposed) {
        if (m_time_until_carrier <= TimeSpan(0)) {
            new_infection_state = InfectionState::Carrier;
        }
        m_time_until_carrier -= dt;
    }
    else {
        new_infection_state = loc.interact(*this, dt, global_infection_params);
        if (new_infection_state == InfectionState::Exposed) {
            m_time_until_carrier = hours(int(global_infection_params.get<IncubationPeriod>()[{this->m_age, this->m_vaccination_state}] * 24));
        }
    }

    if (new_infection_state == InfectionState::Infected_Detected || new_infection_state == InfectionState::Infected_Severe || new_infection_state == InfectionState::Infected_Critical) {
        m_quarantine = true;
    }
    else {
        m_quarantine = false;
    }

    m_infection_state = new_infection_state;
    if (infection_state != new_infection_state) {
        loc.changed_state(*this, infection_state);
    }

    m_time_at_location += dt;
}

void Person::migrate_to(Location& loc_old, Location& loc_new)
{
    if (&loc_old!= &loc_new) {
        loc_old.remove_person(*this);
        m_location_id = {loc_new.get_index(), loc_new.get_type()};
        loc_new.add_person(*this);
        m_time_at_location = epi::TimeSpan(0);
    }
}

void Person::set_assigned_location(Location& location)
{
    m_assigned_locations[(uint32_t)location.get_type()] = location.get_index();
}

void Person::set_assigned_location (LocationId id)
{
    m_assigned_locations[(uint32_t)id.type] = id.index;
}

uint32_t Person::get_assigned_location_index (LocationType type) const
{
    return m_assigned_locations[(uint32_t)type];
}

bool Person::goes_to_work(TimePoint t, const AbmMigrationParameters& params) const
{
    return m_random_workgroup < params.get<WorkRatio>().get_matrix_at(t.days())[0];
}

bool Person::goes_to_school(TimePoint t, const AbmMigrationParameters& params) const
{
    return m_random_schoolgroup < params.get<SchoolRatio>().get_matrix_at(t.days())[0];
}
} // namespace epi
