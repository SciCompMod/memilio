/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "abm/infection.h"
#include "abm/location.h"
#include "abm/person_id.h"
#include "memilio/utils/random_number_generator.h"
#include <cstdint>
#include <vector>

namespace mio
{
namespace abm
{

Person::Person(mio::RandomNumberGenerator& rng, LocationType location_type, LocationId location_id,
               int location_model_id, AgeGroup age, PersonId person_id)
    : m_location(location_id)
    , m_location_type(location_type)
    , m_location_model_id(location_model_id)
    , m_assigned_locations((uint32_t)LocationType::Count, LocationId::invalid_id())
    , m_home_isolation_start(TimePoint(-(std::numeric_limits<int>::max() / 2)))
    , m_age(age)
    , m_time_at_location(0)
    , m_mask(Mask(MaskType::None, TimePoint(-(std::numeric_limits<int>::max() / 2))))
    , m_compliance((uint32_t)InterventionType::Count, 1.)
    , m_cells{0}
    , m_last_transport_mode(TransportMode::Unknown)
    , m_test_results({TestType::Count}, TestResult())
    , m_assigned_location_model_ids((int)LocationType::Count)
    , m_person_id(person_id)
    , m_rng_index(static_cast<uint32_t>(person_id.get()))
{
    m_random_workgroup        = UniformDistribution<double>::get_instance()(rng);
    m_random_schoolgroup      = UniformDistribution<double>::get_instance()(rng);
    m_random_goto_work_hour   = UniformDistribution<double>::get_instance()(rng);
    m_random_goto_school_hour = UniformDistribution<double>::get_instance()(rng);
}

Person::Person(const Person& other, PersonId person_id)
    : Person(other)
{
    m_person_id = person_id;
    m_rng_index = static_cast<uint32_t>(person_id.get());
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

LocationId Person::get_location() const
{
    return m_location;
}

void Person::set_location(LocationType type, LocationId id, int model_id)
{
    m_location          = id;
    m_location_type     = type;
    m_location_model_id = model_id;
    m_time_at_location  = TimeSpan(0);
}

const Infection& Person::get_infection() const
{
    return m_infections.back();
}

Infection& Person::get_infection()
{
    return m_infections.back();
}

void Person::set_assigned_location(LocationType type, LocationId id, int model_id = 0)
{
    m_assigned_locations[static_cast<uint32_t>(type)]          = id;
    m_assigned_location_model_ids[static_cast<uint32_t>(type)] = model_id;
}

LocationId Person::get_assigned_location(LocationType type) const
{
    return m_assigned_locations[static_cast<uint32_t>(type)];
}

int Person::get_assigned_location_model_id(LocationType type) const
{
    return m_assigned_location_model_ids[static_cast<uint32_t>(type)];
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
    m_home_isolation_start = TimePoint(-(std::numeric_limits<int>::max() / 2));
}

bool Person::get_tested(PersonalRandomNumberGenerator& rng, TimePoint t, const TestParameters& params)
{
    ScalarType random = UniformDistribution<double>::get_instance()(rng);
    if (is_infected(t)) {
        // true positive
        if (random < params.sensitivity) {
            // If the Person complies to isolation, start the quarantine.
            if (is_compliant(rng, InterventionType::Isolation)) {
                m_home_isolation_start = t;
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
                m_home_isolation_start = t;
            }
            return true;
        }
    }
}

PersonId Person::get_id() const
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

bool Person::is_compliant(PersonalRandomNumberGenerator& rng, InterventionType intervention) const
{
    ScalarType compliance_check = UniformDistribution<double>::get_instance()(rng);
    return compliance_check <= get_compliance(intervention);
}

ProtectionEvent Person::get_latest_protection() const
{
    ProtectionType latest_protection_type = ProtectionType::NoProtection;
    TimePoint infection_time              = TimePoint(0);
    if (!m_infections.empty()) {
        latest_protection_type = ProtectionType::NaturalInfection;
        infection_time         = m_infections.back().get_start_date();
    }
    if (!m_vaccinations.empty() && infection_time.days() <= m_vaccinations.back().time.days()) {
        latest_protection_type = m_vaccinations.back().type;
        infection_time         = m_vaccinations.back().time;
    }
    return ProtectionEvent{latest_protection_type, infection_time};
}

ScalarType Person::get_protection_factor(TimePoint t, VirusVariant virus, const Parameters& params) const
{
    auto latest_protection = get_latest_protection();
    // If there is no previous protection or vaccination, return 0.
    if (latest_protection.type == ProtectionType::NoProtection) {
        return 0;
    }
    return params.get<InfectionProtectionFactor>()[{latest_protection.type, m_age, virus}](
        t.days() - latest_protection.time.days());
}

void Person::set_mask(MaskType type, TimePoint t)
{
    m_mask.change_mask(type, t);
}

void Person::add_test_result(TimePoint t, TestType type, bool result)
{
    // Remove outdated test results or replace the old result of the same type
    m_test_results[{type}] = {t, result};
}

TestResult Person::get_test_result(TestType type) const
{
    return m_test_results[{type}];
}

} // namespace abm
} // namespace mio
