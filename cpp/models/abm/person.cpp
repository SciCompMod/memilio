/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann
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
#include "abm/location_type.h"
#include "abm/mask_type.h"
#include "abm/parameters.h"
#include "abm/world.h"
#include "abm/location.h"
#include "memilio/utils/random_number_generator.h"
#include <vector>

namespace mio
{
namespace abm
{

Person::Person(Location& location, AgeGroup age, uint32_t person_id)
    : m_location(&location)
    , m_assigned_locations((uint32_t)LocationType::Count, INVALID_LOCATION_INDEX)
    , m_quarantine(false)
    , m_age(age)
    , m_time_at_location(0)
    , m_time_since_negative_test(std::numeric_limits<int>::max() / 2)
    , m_mask(Mask(MaskType::Community))
    , m_wears_mask(false)
    , m_mask_compliance((uint32_t)LocationType::Count, 0.)
    , m_person_id(person_id)
    , m_cells{0}
{
    m_random_workgroup        = UniformDistribution<double>::get_instance()();
    m_random_schoolgroup      = UniformDistribution<double>::get_instance()();
    m_random_goto_work_hour   = UniformDistribution<double>::get_instance()();
    m_random_goto_school_hour = UniformDistribution<double>::get_instance()();
}

void Person::interact(TimePoint t, TimeSpan dt, const GlobalInfectionParameters& params)
{
    if (get_infection_state(t) == InfectionState::Susceptible) { // Susceptible
        m_location->interact(*this, t, dt, params);
    }
    m_time_at_location += dt;
}

void Person::migrate_to(Location& loc_new, const std::vector<uint32_t>& cells)
{
    if (*m_location != loc_new) {
        m_location->remove_person(*this);
        m_location = &loc_new;
        m_cells    = cells;
        loc_new.add_person(*this, cells);
        m_time_at_location = TimeSpan(0);
    }
}

bool Person::is_infected(TimePoint t) const
{
    if (m_infections.empty()) {
        return false;
    }
    // subject to change if Recovered is removed
    if (m_infections.back().get_infection_state(t) == InfectionState::Susceptible ||
        m_infections.back().get_infection_state(t) == InfectionState::Recovered) {
        return false;
    }
    return true;
}

InfectionState Person::get_infection_state(TimePoint t) const
{
    if (m_infections.empty()) {
        return InfectionState::Susceptible;
    }
    else {
        return m_infections.back().get_infection_state(t);
    }
}

void Person::add_new_infection(Infection&& inf)
{
    m_infections.push_back(std::move(inf));
}

Location& Person::get_location()
{
    return *m_location;
}

const Location& Person::get_location() const
{
    return *m_location;
}

void Person::set_assigned_location(Location& location)
{
    /* TODO: This is not safe if the location is not the same as added in the world, e.g. the index is wrong. We need to check this.
    * For now only use it like this:  auto home_id   = world.add_location(mio::abm::LocationType::Home);
    *                                 person.set_assigned_location(home);
    */
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

void Person::detect_infection(TimePoint t)
{
    if (is_infected(t)) {
        m_infections.back().set_detected();
        m_quarantine = true;
    }
}

void Person::remove_quarantine()
{
    m_quarantine = false;
}

bool Person::get_tested(TimePoint t, const TestParameters& params)
{
    ScalarType random = UniformDistribution<double>::get_instance()();
    if (is_infected(t)) {
        // true positive
        if (random < params.sensitivity) {
            m_quarantine = true;
            return true;
        }
        // false negative
        else {
            m_quarantine               = false;
            m_time_since_negative_test = days(0);
            return false;
        }
    }
    else {
        // true negative
        if (random < params.specificity) {
            m_quarantine               = false;
            m_time_since_negative_test = days(0);
            return false;
        }
        // false positive
        else {
            m_quarantine = true;
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

ScalarType Person::get_mask_protective_factor(const GlobalInfectionParameters& params) const
{
    if (m_wears_mask == false) {
        return 0.;
    }
    else {
        return params.get<MaskProtection>()[m_mask.get_type()];
    }
}

bool Person::apply_mask_intervention(const Location& target)
{
    if (target.get_npi_active() == false) {
        m_wears_mask = false;
        if (get_mask_compliance(target.get_type()) > 0.) {
            // draw if the person wears a mask even if not required
            ScalarType wear_mask = UniformDistribution<double>::get_instance()();
            if (wear_mask < get_mask_compliance(target.get_type())) {
                m_wears_mask = true;
            }
        }
    }
    else {
        m_wears_mask = true;
        if (get_mask_compliance(target.get_type()) < 0.) {
            // draw if a person refuses to wear the required mask
            ScalarType wear_mask = UniformDistribution<double>::get_instance()(-1., 0.);
            if (wear_mask > get_mask_compliance(target.get_type())) {
                m_wears_mask = false;
            }
            return false;
        }
        if (m_wears_mask == true) {

            if (static_cast<int>(m_mask.get_type()) < static_cast<int>(target.get_required_mask())) {
                m_mask.change_mask(target.get_required_mask());
            }
        }
    }
    return true;
}

} // namespace abm
} // namespace mio
