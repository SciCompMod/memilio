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
#include "abm/infection.h"
#include "abm/location.h"
#include "memilio/utils/random_number_generator.h"
#include <vector>

namespace mio
{
namespace abm
{

Person::Person(mio::RandomNumberGenerator& rng, LocationType location_type, LocationId location_id, AgeGroup age,
               PersonId person_id)
    : m_location(location_id)
    , m_location_type(location_type)
    , m_assigned_locations((uint32_t)LocationType::Count, LocationId::invalid_id())
    , m_quarantine_start(TimePoint(-(std::numeric_limits<int>::max() / 2)))
    , m_age(age)
    , m_time_at_location(0)
    , m_time_of_last_test(TimePoint(-(std::numeric_limits<int>::max() / 2)))
    , m_mask(Mask(MaskType::Community))
    , m_wears_mask(false)
    , m_mask_compliance((uint32_t)LocationType::Count, 0.)
    , m_person_id(person_id)
    , m_cells{0}
    , m_last_transport_mode(TransportMode::Unknown)
{
    m_random_workgroup        = UniformDistribution<double>::get_instance()(rng);
    m_random_schoolgroup      = UniformDistribution<double>::get_instance()(rng);
    m_random_goto_work_hour   = UniformDistribution<double>::get_instance()(rng);
    m_random_goto_school_hour = UniformDistribution<double>::get_instance()(rng);
}

Person::Person(const Person& other, PersonId id)
    : Person(other)
{
    m_person_id = id;
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

void Person::set_location(LocationType type, LocationId id)
{
    m_location         = id;
    m_location_type    = type;
    m_time_at_location = TimeSpan(0);
}

const Infection& Person::get_infection() const
{
    return m_infections.back();
}

Infection& Person::get_infection()
{
    return m_infections.back();
}

void Person::set_assigned_location(LocationType type, LocationId id)
{
    m_assigned_locations[static_cast<uint32_t>(type)] = id;
}

LocationId Person::get_assigned_location(LocationType type) const
{
    return m_assigned_locations[static_cast<uint32_t>(type)];
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

bool Person::get_tested(PersonalRandomNumberGenerator& rng, TimePoint t, const TestParameters& params)
{
    ScalarType random   = UniformDistribution<double>::get_instance()(rng);
    m_time_of_last_test = t;
    if (is_infected(t)) {
        // true positive
        if (random < params.sensitivity) {
            m_quarantine_start = t;
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
            m_quarantine_start = t;
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
    if (m_wears_mask == false) {
        return 0.;
    }
    else {
        return params.get<MaskProtection>()[m_mask.get_type()];
    }
}

bool Person::apply_mask_intervention(PersonalRandomNumberGenerator& rng, const Location& target)
{
    if (target.get_npi_active() == false) {
        m_wears_mask = false;
        if (get_mask_compliance(target.get_type()) > 0.) {
            // draw if the person wears a mask even if not required
            ScalarType wear_mask = UniformDistribution<double>::get_instance()(rng);
            if (wear_mask < get_mask_compliance(target.get_type())) {
                m_wears_mask = true;
            }
        }
    }
    else {
        m_wears_mask = true;
        if (get_mask_compliance(target.get_type()) < 0.) {
            // draw if a person refuses to wear the required mask
            ScalarType wear_mask = UniformDistribution<double>::get_instance()(rng, -1., 0.);
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

} // namespace abm
} // namespace mio
