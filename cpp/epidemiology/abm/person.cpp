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

Person::Person(LocationId id, InfectionProperties infection_properties, AbmAgeGroup age,
               const GlobalInfectionParameters& global_params)
    : m_location_id(id)
    , m_assigned_locations((uint32_t)LocationType::Count, INVALID_LOCATION_INDEX)
    , m_state(infection_properties.state)
    , m_quarantine(false)
    , m_age(age)
    , m_time_at_location(std::numeric_limits<int>::max() / 2) //avoid overflow on next steps
    , m_time_since_negative_test(std::numeric_limits<int>::max() / 2)
{
    m_random_workgroup   = UniformDistribution<double>::get_instance()();
    m_random_schoolgroup = UniformDistribution<double>::get_instance()();
    if (infection_properties.state == InfectionState::Infected && infection_properties.detected) {
        m_quarantine = true;
    }
    if (infection_properties.state == InfectionState::Exposed) {
        m_time_until_carrier = hours(
            UniformIntDistribution<int>::get_instance()(0, int(global_params.get<IncubationPeriod>()[m_age] * 24)));
    }
}

Person::Person(Location& location, InfectionProperties infection_properties, AbmAgeGroup age,
               const GlobalInfectionParameters& global_params)
    : Person({location.get_index(), location.get_type()}, infection_properties, age, global_params)
{
}

void Person::interact(TimeSpan dt, const GlobalInfectionParameters& global_infection_params, Location& loc,
                      const GlobalTestingParameters& global_testing_params)
{
    auto state     = m_state;
    auto new_state = state;

    if (state == InfectionState::Exposed) {
        if (m_time_until_carrier <= TimeSpan(0)) {
            new_state = InfectionState::Carrier;
        }
        m_time_until_carrier -= dt;
    }
    else {
        new_state = loc.interact(*this, dt, global_infection_params);
        if (new_state == InfectionState::Exposed) {
            m_time_until_carrier = hours(int(global_infection_params.get<IncubationPeriod>()[{this->m_age}] * 24));
        }
    }

    if (new_state == InfectionState::Infected_Severe || new_state == InfectionState::Infected_Critical) {
        m_quarantine = true;
    }
    else if (new_state == InfectionState::Infected) {
        double rand = UniformDistribution<double>::get_instance()();
        if (rand < global_infection_params.get<TestWhileInfected>()[this->m_age] * dt.days()) {
            this->get_tested(global_testing_params.get<AntigenTest>());
        }
    }
    else {
        m_quarantine = false;
    }

    m_state = new_state;
    if (state != new_state) {
        loc.changed_state(*this, state);
    }

    m_time_at_location += dt;
}

void Person::migrate_to(Location& loc_old, Location& loc_new)
{
    if (&loc_old != &loc_new) {
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

void Person::set_assigned_location(LocationId id)
{
    m_assigned_locations[(uint32_t)id.type] = id.index;
}

uint32_t Person::get_assigned_location_index(LocationType type) const
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

bool Person::get_tested(const TestParameters& params)
{
    double random = UniformDistribution<double>::get_instance()();
    if (m_state == InfectionState::Carrier || m_state == InfectionState::Infected ||
        m_state == InfectionState::Infected_Severe || m_state == InfectionState::Infected_Critical) {
        if (random < params.sensitivity) {
            m_quarantine = true;
            return true;
        }
        else {
            m_quarantine               = false;
            m_time_since_negative_test = days(0);
            return false;
        }
    }
    else {
        if (random < params.specificity) {
            m_quarantine               = false;
            m_time_since_negative_test = days(0);
            return false;
        }
        else {
            return true;
        }
    }
}
} // namespace epi
