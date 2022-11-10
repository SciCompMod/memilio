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
namespace abm
{

Person::Person(const LocationId& id, const AgeGroup& age, const Infection& infection,
               const VaccinationState& vaccination_state, const uint32_t person_id)
    : Person(id, age, vaccination_state, person_id)
{
    m_infections = std::vector<Infection>{infection};
    if (is_infected() && infection.is_detected()) {
        m_quarantine = true;
    }
}

Person::Person(const Location& location, const AgeGroup& age, const Infection& infection,
               const VaccinationState& vaccination_state, const uint32_t person_id)
    : Person({location.get_index(), location.get_type()}, age, infection, vaccination_state, person_id)
{
}

Person::Person(const LocationId& id, const AgeGroup& age, const VaccinationState& vaccination_state,
               const uint32_t person_id)
    : m_location_id(id)
    , m_assigned_locations((uint32_t)LocationType::Count, INVALID_LOCATION_INDEX)
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
}

Person::Person(const Location& location, const AgeGroup& age, const VaccinationState& vaccination_state,
               const uint32_t person_id)
    : Person({location.get_index(), location.get_type()}, age, vaccination_state, person_id)
{
}

void Person::interact(const TimePoint& t, const TimeSpan& dt,
                      const GlobalInfectionParameters& global_infection_parameters,
                      CustomIndexArray<std::shared_ptr<Virus>, VirusVariant> virus_variants, Location& loc)
{
    auto current_infection_state = get_infection_state();
    if (current_infection_state == InfectionState::Susceptible) { // Susceptible
        Infection infection = loc.interact(*this, t, dt, virus_variants);
        if (m_infections.empty()) {
            m_infections.push_back(infection);
        }
    }

    //        m_time_until_carrier = hours(
    //            int(global_infection_params.get<IncubationPeriod>()[{this->m_age, this->m_vaccination_state}] * 24));

    auto new_infection_state = get_infection_state(t);
    if (current_infection_state != new_infection_state) {
        loc.changed_state(*this, current_infection_state, t);
    }

    m_time_at_location += dt;
}

void Person::migrate_to(Location& loc_old, Location& loc_new, const TimePoint& t, const std::vector<uint32_t>& cells)
{
    if (&loc_old != &loc_new) {
        loc_old.remove_person(*this, t);
        m_location_id = {loc_new.get_index(), loc_new.get_type()};
        m_cells       = cells;
        loc_new.add_person(*this, t);
        m_time_at_location = TimeSpan(0);
    }
}

bool Person::is_infected() const
{
    if (m_infections.empty()) {
        return false;
    }
    // subject to change if Recovered is removed
    if (m_infections.back().get_infection_state() == InfectionState::Susceptible ||
        m_infections.back().get_infection_state() == InfectionState::Recovered_Carrier ||
        m_infections.back().get_infection_state() == InfectionState::Recovered_Infected) {
        return false;
    }
    return true;
}

const InfectionState& Person::get_infection_state() const
{
    if (m_infections.empty()) {
        return std::move(InfectionState::Susceptible);
    }
    else {
        return m_infections.back().get_infection_state();
    }
}

const InfectionState& Person::get_infection_state(const TimePoint& t) const
{
    if (m_infections.empty()) {
        return std::move(InfectionState::Susceptible);
    }
    else {
        return m_infections.back().get_infection_state(t);
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

bool Person::goes_to_work(TimePoint t, const MigrationParameters& params) const
{
    return m_random_workgroup < params.get<WorkRatio>().get_matrix_at(t.days())[0];
}

TimeSpan Person::get_go_to_work_time(const MigrationParameters& params) const
{
    TimeSpan minimum_goto_work_time = params.get<GotoWorkTimeMinimum>()[m_age];
    TimeSpan maximum_goto_work_time = params.get<GotoWorkTimeMaximum>()[m_age];
    int timeSlots                   = (maximum_goto_work_time.seconds() - minimum_goto_work_time.seconds());
    int seconds_after_minimum       = int(timeSlots * m_random_goto_work_hour);
    return minimum_goto_work_time + seconds(seconds_after_minimum);
}

TimeSpan Person::get_go_to_school_time(const MigrationParameters& params) const
{
    TimeSpan minimum_goto_school_time = params.get<GotoSchoolTimeMinimum>()[m_age];
    TimeSpan maximum_goto_school_time = params.get<GotoSchoolTimeMaximum>()[m_age];
    int timeSlots                     = (maximum_goto_school_time.seconds() - minimum_goto_school_time.seconds());
    int seconds_after_minimum         = int(timeSlots * m_random_goto_school_hour);
    return minimum_goto_school_time + seconds(seconds_after_minimum);
}

bool Person::goes_to_school(TimePoint t, const MigrationParameters& params) const
{
    return m_random_schoolgroup < params.get<SchoolRatio>().get_matrix_at(t.days())[0];
}

bool Person::get_tested(const TestParameters& params)
{
    double random = UniformDistribution<double>::get_instance()();
    if (is_infected()) {
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

std::vector<uint32_t>& Person::get_cells()
{
    return m_cells;
}

const std::vector<uint32_t>& Person::get_cells() const
{
    return m_cells;
}

} // namespace abm
} // namespace mio
