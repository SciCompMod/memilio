/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth
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
#include "abm/person.h"
#include "abm/world.h"
#include "abm/location.h"
#include "memilio/utils/random_number_generator.h"

namespace mio
{

Person::Person(LocationId id, InfectionProperties infection_properties, AbmAgeGroup age,
               const GlobalInfectionParameters& global_params, VaccinationState vaccination_state, uint32_t person_id)
    : m_location_id(id)
    , m_assigned_locations((uint32_t)LocationType::Count, INVALID_LOCATION_INDEX)
    , m_infection_state(infection_properties.state)
    , m_vaccination_state(vaccination_state)
    , m_quarantine(false)
    , m_age(age)
    , m_time_at_location(std::numeric_limits<int>::max() / 2) //avoid overflow on next steps
    , m_time_since_negative_test(std::numeric_limits<int>::max() / 2)
    , m_person_id(person_id)
{
    m_random_workgroup        = UniformDistribution<double>::get_instance()();
    m_random_schoolgroup      = UniformDistribution<double>::get_instance()();
    m_random_goto_work_hour   = UniformDistribution<double>::get_instance()();
    m_random_goto_school_hour = UniformDistribution<double>::get_instance()();

    if (infection_properties.state == InfectionState::Infected && infection_properties.detected) {
        m_quarantine = true;
    }
    if (infection_properties.state == InfectionState::Exposed) {
        m_time_until_carrier = hours(UniformIntDistribution<int>::get_instance()(
            0, int(global_params.get<IncubationPeriod>()[{m_age, m_vaccination_state}] * 24)));
    }
}

Person::Person(Location& location, InfectionProperties infection_properties, AbmAgeGroup age,
               const GlobalInfectionParameters& global_params, VaccinationState vaccination_state, uint32_t person_id)
    : Person({location.get_index(), location.get_type()}, infection_properties, age, global_params, vaccination_state,
             person_id)
{
}

void Person::interact(TimeSpan dt, const GlobalInfectionParameters& global_infection_params, Location& loc,
                      const GlobalTestingParameters& global_testing_params)
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
            m_time_until_carrier = hours(
                int(global_infection_params.get<IncubationPeriod>()[{this->m_age, this->m_vaccination_state}] * 24));
        }
    }

    if (new_infection_state == InfectionState::Infected_Severe ||
        new_infection_state == InfectionState::Infected_Critical) {
        m_quarantine = true;
    }
    else if (new_infection_state == InfectionState::Infected) {
        double rand = UniformDistribution<double>::get_instance()();
        if (rand < global_infection_params.get<TestWhileInfected>()[this->m_age] * dt.days()) {
            this->get_tested(global_testing_params.get<AntigenTest>());
        }
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

void Person::migrate_to(Location& loc_old, Location& loc_new, const std::vector<uint32_t>& cells)
{
    if (&loc_old != &loc_new) {
        loc_old.remove_person(*this);
        m_location_id = {loc_new.get_index(), loc_new.get_type()};
        m_cells       = cells;
        loc_new.add_person(*this);
        m_time_at_location = mio::TimeSpan(0);
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

void Person::set_infection_state(InfectionState inf_state)
{
    m_infection_state = inf_state;
}

uint32_t Person::get_assigned_location_index(LocationType type) const
{
    return m_assigned_locations[(uint32_t)type];
}

bool Person::goes_to_work(TimePoint t, const AbmMigrationParameters& params) const
{
    return m_random_workgroup < params.get<WorkRatio>().get_matrix_at(t.days())[0];
}

TimeSpan Person::get_go_to_work_time(const AbmMigrationParameters& params) const
{
    TimeSpan minimum_goto_work_time = params.get<GotoWorkTimeMinimum>()[m_age];
    TimeSpan maximum_goto_work_time = params.get<GotoWorkTimeMaximum>()[m_age];
    int timeSlots                   = (maximum_goto_work_time.seconds() - minimum_goto_work_time.seconds());
    int seconds_after_minimum       = int(timeSlots * m_random_goto_work_hour);
    return minimum_goto_work_time + mio::seconds(seconds_after_minimum);
}

TimeSpan Person::get_go_to_school_time(const AbmMigrationParameters& params) const
{
    TimeSpan minimum_goto_school_time = params.get<GotoSchoolTimeMinimum>()[m_age];
    TimeSpan maximum_goto_school_time = params.get<GotoSchoolTimeMaximum>()[m_age];
    int timeSlots                     = (maximum_goto_school_time.seconds() - minimum_goto_school_time.seconds());
    int seconds_after_minimum         = int(timeSlots * m_random_goto_school_hour);
    return minimum_goto_school_time + mio::seconds(seconds_after_minimum);
}

bool Person::goes_to_school(TimePoint t, const AbmMigrationParameters& params) const
{
    return m_random_schoolgroup < params.get<SchoolRatio>().get_matrix_at(t.days())[0];
}

bool Person::get_tested(const TestParameters& params)
{
    double random = UniformDistribution<double>::get_instance()();
    if (m_infection_state == InfectionState::Carrier || m_infection_state == InfectionState::Infected ||
        m_infection_state == InfectionState::Infected_Severe ||
        m_infection_state == InfectionState::Infected_Critical) {
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

uint32_t Person::get_person_id()
{
    return m_person_id;
}

void Person::set_cells(const std::vector<uint32_t>& cells)
{
    m_cells = cells;
}

std::vector<uint32_t>& Person::get_cells()
{
    return m_cells;
}

const std::vector<uint32_t>& Person::get_cells() const
{
    return m_cells;
}
} // namespace mio
