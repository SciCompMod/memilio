/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann, Khoa Nguyen
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
#include "abm/time.h"
#include <vector>

namespace mio
{
namespace abm
{

Person::Person(mio::RandomNumberGenerator& rng, Location& location, AgeGroup age, uint32_t person_id)
    : m_location(&location)
    , m_assigned_locations((uint32_t)LocationType::Count, INVALID_LOCATION_INDEX)
    , m_quarantine_start(TimePoint(-(std::numeric_limits<int>::max() / 2)))
    , m_age(age)
    , m_time_at_location(0)
    , m_time_of_last_test(TimePoint(-(std::numeric_limits<int>::max() / 2)))
    , m_mask(Mask(MaskType::None, TimePoint(-(std::numeric_limits<int>::max() / 2))))
    , m_compliance((uint32_t)InterventionType::Count, 1.)
    , m_person_id(person_id)
    , m_cells{0}
    , m_last_transport_mode(TransportMode::Unknown)
{
    m_random_workgroup        = UniformDistribution<double>::get_instance()(rng);
    m_random_schoolgroup      = UniformDistribution<double>::get_instance()(rng);
    m_random_goto_work_hour   = UniformDistribution<double>::get_instance()(rng);
    m_random_goto_school_hour = UniformDistribution<double>::get_instance()(rng);
}

Person Person::copy_person(Location& location)
{
    Person copied_person     = Person(*this);
    copied_person.m_location = &location;
    location.add_person(*this);
    return copied_person;
}

void Person::interact(RandomNumberGenerator& rng, TimePoint t, TimeSpan dt, const Parameters& params)
{
    if (get_infection_state(t) == InfectionState::Susceptible) { // Susceptible
        m_location->interact(rng, *this, t, dt, params);
    }
    m_time_at_location += dt;
}

void Person::migrate_to(Location& loc_new, mio::abm::TransportMode transport_mode, const std::vector<uint32_t>& cells)
{
    if (*m_location != loc_new) {
        m_location->remove_person(*this);
        m_location = &loc_new;
        m_cells    = cells;
        loc_new.add_person(*this, cells);
        m_time_at_location    = TimeSpan(0);
        m_last_transport_mode = transport_mode;
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

const Infection& Person::get_infection() const
{
    return m_infections.back();
}

Infection& Person::get_infection()
{
    return m_infections.back();
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

bool Person::goes_to_work(TimePoint t, const Parameters& params) const
{
    return m_random_workgroup < params.get<WorkRatio>().get_matrix_at(t.days())[0];
}

TimeSpan Person::get_go_to_work_time(const Parameters& params) const
{
    TimeSpan minimum_goto_work_time = params.get<GotoWorkTimeMinimum>()[m_age];
    TimeSpan maximum_goto_work_time = params.get<GotoWorkTimeMaximum>()[m_age];
    int timeSlots                   = (maximum_goto_work_time.seconds() - minimum_goto_work_time.seconds());
    int seconds_after_minimum       = int(timeSlots * m_random_goto_work_hour);
    return minimum_goto_work_time + seconds(seconds_after_minimum);
}

TimeSpan Person::get_go_to_school_time(const Parameters& params) const
{
    TimeSpan minimum_goto_school_time = params.get<GotoSchoolTimeMinimum>()[m_age];
    TimeSpan maximum_goto_school_time = params.get<GotoSchoolTimeMaximum>()[m_age];
    int timeSlots                     = (maximum_goto_school_time.seconds() - minimum_goto_school_time.seconds());
    int seconds_after_minimum         = int(timeSlots * m_random_goto_school_hour);
    return minimum_goto_school_time + seconds(seconds_after_minimum);
}

bool Person::goes_to_school(TimePoint t, const Parameters& params) const
{
    return m_random_schoolgroup < params.get<SchoolRatio>().get_matrix_at(t.days())[0];
}

void Person::remove_quarantine()
{
    m_quarantine_start = TimePoint(-(std::numeric_limits<int>::max() / 2));
}

bool Person::get_tested(RandomNumberGenerator& rng, TimePoint t, const TestParameters& params)
{
    ScalarType random   = UniformDistribution<double>::get_instance()(rng);
    m_time_of_last_test = t;
    if (is_infected(t)) {
        // true positive
        if (random < params.sensitivity) {
            // If the Person complies to isolation, start the quarantine.
            if (is_compliant(rng, InterventionType::Isolation)) {
                m_quarantine_start = t;
            }
            m_infections.back().set_detected();
            return true;
        }
        // false negative
        else {
            return false;
        }
    }
    else {
        // true negative
        if (random < params.specificity) {
            return false;
        }
        // false positive
        else {
            // If the Person complies to isolation, start the quarantine.
            if (is_compliant(rng, InterventionType::Isolation)) {
                m_quarantine_start = t;
            }
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

ScalarType Person::get_mask_protective_factor(const Parameters& params) const
{
    return params.get<MaskProtection>()[m_mask.get_type()];
}

bool Person::is_compliant(RandomNumberGenerator& rng, InterventionType intervention) const
{
    ScalarType compliance_check = UniformDistribution<double>::get_instance()(rng);
    return compliance_check <= get_compliance(intervention);
}

std::pair<ExposureType, TimePoint> Person::get_latest_protection() const
{
    ExposureType latest_exposure_type = ExposureType::NoProtection;
    TimePoint infection_time          = TimePoint(0);
    if (!m_infections.empty()) {
        latest_exposure_type = ExposureType::NaturalInfection;
        infection_time       = m_infections.back().get_start_date();
    }
    if (!m_vaccinations.empty() && infection_time.days() <= m_vaccinations.back().time.days()) {
        latest_exposure_type = m_vaccinations.back().exposure_type;
        infection_time       = m_vaccinations.back().time;
    }
    return std::make_pair(latest_exposure_type, infection_time);
}

ScalarType Person::get_protection_factor(TimePoint t, VirusVariant virus, const Parameters& params) const
{
    auto latest_protection = get_latest_protection();
    // If there is no previous protection or vaccination, return 0.
    if (latest_protection.first == ExposureType::NoProtection) {
        return 0;
    }
    return params.get<InfectionProtectionFactor>()[{latest_protection.first, m_age, virus}](
        t.days() - latest_protection.second.days());
}

void Person::set_mask(MaskType type, TimePoint t)
{
    m_mask.change_mask(type, t);
}

} // namespace abm
} // namespace mio
