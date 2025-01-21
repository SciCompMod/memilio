/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: David Kerkmann, Sascha Korf, Khoa Nguyen
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

#include "abm/infection.h"
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <math.h>
#include <stdexcept>

namespace mio
{
namespace abm
{

Infection::Infection(PersonalRandomNumberGenerator& rng, VirusVariant virus, AgeGroup age, const Parameters& params,
                     TimePoint init_date, InfectionState init_state, ProtectionEvent latest_protection, bool detected)
    : m_virus_variant(virus)
    , m_detected(detected)
{
    assert(age.get() < params.get_num_groups());
    m_viral_load.start_date = draw_infection_course(rng, age, params, init_date, init_state, latest_protection);

    auto vl_params                    = params.get<ViralLoadDistributions>()[{virus, age}];
    ScalarType high_viral_load_factor = 1;
    if (latest_protection.type != ProtectionType::NoProtection) {
        high_viral_load_factor -= params.get<HighViralLoadProtectionFactor>()[{latest_protection.type, age, virus}](
            init_date.days() - latest_protection.time.days());
    }
    m_viral_load.peak = vl_params.viral_load_peak.get_distribution_instance()(rng, vl_params.viral_load_peak.params) *
                        high_viral_load_factor;
    m_time_is_infected = m_infection_course.back().first - m_infection_course[0].first;
    // the third entry in the infection course is the time an agents switches to InfectedSymptoms or recovers after non-symptomatic
    TimePoint t_peak      = m_infection_course[2].first;
    m_viral_load.incline  = std::abs(m_viral_load.peak / (t_peak - m_infection_course[0].first).days());
    m_viral_load.end_date = m_viral_load.start_date + m_time_is_infected;
    if (m_infection_course.size() < 4) {
        t_peak            = m_viral_load.start_date + TimeSpan(int(0.5 * m_time_is_infected.seconds()));
        m_viral_load.peak = m_viral_load.incline * (t_peak - m_viral_load.start_date).days();
    }
    m_viral_load.decline = -m_viral_load.peak / (m_viral_load.end_date - t_peak).days();

    auto inf_params = params.get<InfectivityDistributions>()[{virus, age}];
    m_log_norm_alpha =
        inf_params.infectivity_alpha.get_distribution_instance()(rng, inf_params.infectivity_alpha.params);
    m_log_norm_beta = inf_params.infectivity_beta.get_distribution_instance()(rng, inf_params.infectivity_beta.params);

    // // print infection course
    // for (size_t i = size_t(0); i < m_infection_course.size(); ++i) {
    //     std::cout << static_cast<int>(m_infection_course[i].second) << ": " << m_infection_course[i].first.days();
    // }
    // std::cout << "\n";
    // for (auto t = m_viral_load.start_date; t <= m_viral_load.end_date + days(3); t += days(1)) {
    //     std::cout << "t: " << t.days() << " i(t): " << get_infectivity(t) << std::endl;
    // }
}

ScalarType Infection::get_viral_load(TimePoint t) const
{
    if (t >= m_viral_load.start_date && t <= m_viral_load.end_date) {
        if (t.days() <= m_viral_load.start_date.days() + m_viral_load.peak / m_viral_load.incline) {
            return m_viral_load.incline * (t - m_viral_load.start_date).days();
        }
        else {
            return m_viral_load.peak + m_viral_load.decline * (t.days() - m_viral_load.peak / m_viral_load.incline -
                                                               m_viral_load.start_date.days());
        }
    }
    else {
        return 0.;
    }
}

ScalarType Infection::get_infectivity(TimePoint t) const
{
    auto time_shift = TimeSpan(int(0.6 * (m_infection_course[1].first - m_infection_course[0].first).seconds()));
    if (m_viral_load.start_date + time_shift >= t)
        return 0;
    ScalarType infectivity = 1 / (1 + exp(-(m_log_norm_alpha + m_log_norm_beta * get_viral_load(t - time_shift))));
    if ((infectivity < 0) ||
        (infectivity > (1.0 / (1 + exp(-(m_log_norm_alpha + m_log_norm_beta * m_viral_load.peak)))))) {
        log_error("Error: Infectivity smaller zero or above peak value.");
    }
    return 0.75 * infectivity;
}

VirusVariant Infection::get_virus_variant() const
{
    return m_virus_variant;
}

InfectionState Infection::get_infection_state(TimePoint t) const
{
    if (t < m_infection_course[0].first)
        return InfectionState::Susceptible;

    return (*std::prev(std::upper_bound(m_infection_course.begin(), m_infection_course.end(), t,
                                        [](const TimePoint& s, std::pair<TimePoint, InfectionState> state) {
                                            return state.first > s;
                                        })))
        .second;
}

TimePoint Infection::get_infection_start() const
{
    return (*std::find_if(m_infection_course.begin(), m_infection_course.end(),
                          [](const std::pair<TimePoint, InfectionState>& inf) {
                              return (inf.second == InfectionState::Exposed);
                          }))
        .first;
}

TimePoint Infection::get_infection_end() const
{
    return (*std::find_if(m_infection_course.begin(), m_infection_course.end(),
                          [](const std::pair<TimePoint, InfectionState>& inf) {
                              return (inf.second == InfectionState::Recovered || inf.second == InfectionState::Dead);
                          }))
        .first;
}

void Infection::set_detected()
{
    m_detected = true;
}

bool Infection::is_detected() const
{
    return m_detected;
}

TimePoint Infection::get_start_date() const
{
    return m_viral_load.start_date;
}

TimeSpan Infection::get_time_in_state(InfectionState state)
{
    auto pos = std::find_if(m_infection_course.begin(), m_infection_course.end(),
                            [state](const std::pair<TimePoint, InfectionState>& inf) {
                                return (inf.second == state);
                            });
    // infection state is not part of infection course
    if (pos == m_infection_course.end()) {
        return TimeSpan(0);
    }
    return ((pos + 1)->first - pos->first);
}

TimePoint Infection::draw_infection_course(PersonalRandomNumberGenerator& rng, AgeGroup age, const Parameters& params,
                                           TimePoint init_date, InfectionState init_state,
                                           ProtectionEvent latest_protection)
{
    assert(age.get() < params.get_num_groups());
    TimePoint start_date = draw_infection_course_backward(rng, age, params, init_date, init_state);
    draw_infection_course_forward(rng, age, params, init_date, init_state, latest_protection);
    return start_date;
}

void Infection::draw_infection_course_forward(PersonalRandomNumberGenerator& rng, AgeGroup age,
                                              const Parameters& params, TimePoint init_date, InfectionState start_state,
                                              ProtectionEvent latest_protection)
{
    assert(age.get() < params.get_num_groups());
    auto t = init_date;
    TimeSpan time_period{}; // time period for current infection state
    auto time_in_state = params.get<IncubationPeriod>()[{
        m_virus_variant, age}]; // time distribution parameters for current infection state
    InfectionState next_state{start_state}; // next state to enter
    m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, next_state));
    auto& uniform_dist = UniformDistribution<double>::get_instance();
    ScalarType v; // random draws
    while ((next_state != InfectionState::Recovered && next_state != InfectionState::Dead)) {
        switch (next_state) {
        case InfectionState::Exposed:
            // roll out how long until infected without symptoms
            time_in_state = params.get<IncubationPeriod>()[{m_virus_variant, age}];
            time_period   = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
            next_state    = InfectionState::InfectedNoSymptoms;
            break;
        case InfectionState::InfectedNoSymptoms:
            // roll out next infection step
            v = uniform_dist(rng);
            if (v < params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}]) {
                time_in_state = params.get<TimeInfectedNoSymptomsToSymptoms>()[{m_virus_variant, age}];
                time_period   = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                next_state    = InfectionState::InfectedSymptoms;
            }
            else {
                time_in_state = params.get<TimeInfectedNoSymptomsToRecovered>()[{m_virus_variant, age}];
                time_period   = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                next_state    = InfectionState::Recovered;
            }

            break;
        case InfectionState::InfectedSymptoms:
            // roll out next infection step
            {
                ScalarType severity_protection_factor = 0.5;
                if (latest_protection.type != ProtectionType::NoProtection) {
                    severity_protection_factor =
                        params.get<SeverityProtectionFactor>()[{latest_protection.type, age, m_virus_variant}](
                            t.days() - latest_protection.time.days());
                }
                // roll out next infection step
                v = uniform_dist(rng);
                if (v < (1 - severity_protection_factor) *
                            params.get<SeverePerInfectedSymptoms>()[{m_virus_variant, age}]) {
                    time_in_state = params.get<TimeInfectedSymptomsToSevere>()[{m_virus_variant, age}];
                    time_period   = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                    next_state    = InfectionState::InfectedSevere;
                }
                else {
                    time_in_state = params.get<TimeInfectedSymptomsToRecovered>()[{m_virus_variant, age}];
                    time_period   = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                    next_state    = InfectionState::Recovered;
                }
                break;
            }
        case InfectionState::InfectedSevere:
            // roll out next infection step
            v = uniform_dist(rng);
            if (v < params.get<CriticalPerInfectedSevere>()[{m_virus_variant, age}]) {
                time_in_state = params.get<TimeInfectedSevereToCritical>()[{m_virus_variant, age}];
                time_period   = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                next_state    = InfectionState::InfectedCritical;
            }
            else {
                time_in_state = params.get<TimeInfectedSevereToRecovered>()[{m_virus_variant, age}];
                time_period   = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                next_state    = InfectionState::Recovered;
            }
            break;
        case InfectionState::InfectedCritical:
            // roll out next infection step
            v = uniform_dist(rng);
            if (v < params.get<DeathsPerInfectedCritical>()[{m_virus_variant, age}]) {
                time_in_state = params.get<TimeInfectedCriticalToDead>()[{m_virus_variant, age}];
                time_period   = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                next_state    = InfectionState::Dead;
            }
            else {
                time_in_state = params.get<TimeInfectedCriticalToRecovered>()[{m_virus_variant, age}];
                time_period   = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                next_state    = InfectionState::Recovered;
            }
            break;
        default:
            break;
        }
        t = t + time_period;
        m_infection_course.push_back({t, next_state});
    }
}

TimePoint Infection::draw_infection_course_backward(PersonalRandomNumberGenerator& rng, AgeGroup age,
                                                    const Parameters& params, TimePoint init_date,
                                                    InfectionState init_state)
{
    auto start_date = init_date;
    TimeSpan time_period{}; // time period for current infection state
    auto time_in_state = params.get<IncubationPeriod>()[{
        m_virus_variant, age}]; // time distribution parameters for current infection state
    InfectionState previous_state{init_state}; // next state to enter
    auto& uniform_dist = UniformDistribution<double>::get_instance();
    ScalarType v; // random draws

    while ((previous_state != InfectionState::Exposed)) {
        switch (previous_state) {

        case InfectionState::InfectedNoSymptoms: {
            time_in_state  = params.get<IncubationPeriod>()[{m_virus_variant, age}];
            time_period    = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
            previous_state = InfectionState::Exposed;
        } break;

        case InfectionState::InfectedSymptoms: {
            time_in_state  = params.get<TimeInfectedNoSymptomsToSymptoms>()[{m_virus_variant, age}];
            time_period    = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
            previous_state = InfectionState::InfectedNoSymptoms;
        } break;

        case InfectionState::InfectedSevere: {
            time_in_state  = params.get<TimeInfectedSymptomsToSevere>()[{m_virus_variant, age}];
            time_period    = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
            previous_state = InfectionState::InfectedSymptoms;
        } break;

        case InfectionState::InfectedCritical: {
            time_in_state  = params.get<TimeInfectedSevereToCritical>()[{m_virus_variant, age}];
            time_period    = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
            previous_state = InfectionState::InfectedSevere;
        } break;

        case InfectionState::Recovered: {
            // roll out next infection step
            v = uniform_dist(rng);
            // compute correct probabilities while factoring out the chance to die
            auto p_death = params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}] *
                           params.get<SeverePerInfectedSymptoms>()[{m_virus_variant, age}] *
                           params.get<CriticalPerInfectedSevere>()[{m_virus_variant, age}] *
                           params.get<DeathsPerInfectedCritical>()[{m_virus_variant, age}];
            if (v < (1 - params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}]) / (1 - p_death)) {
                time_in_state  = params.get<TimeInfectedNoSymptomsToRecovered>()[{m_virus_variant, age}];
                time_period    = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                previous_state = InfectionState::InfectedNoSymptoms;
            }
            else if (v < (1 - params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}] *
                                  (1 - params.get<SeverePerInfectedSymptoms>()[{m_virus_variant, age}])) /
                             (1 - p_death)) {
                time_in_state  = params.get<TimeInfectedSymptomsToRecovered>()[{m_virus_variant, age}];
                time_period    = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                previous_state = InfectionState::InfectedSymptoms;
            }
            else if (v < (1 - params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}] *
                                  params.get<SeverePerInfectedSymptoms>()[{m_virus_variant, age}] *
                                  (1 - params.get<CriticalPerInfectedSevere>()[{m_virus_variant, age}])) /
                             (1 - p_death)) {
                time_in_state  = params.get<TimeInfectedSevereToRecovered>()[{m_virus_variant, age}];
                time_period    = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                previous_state = InfectionState::InfectedSevere;
            }
            else {
                time_in_state  = params.get<TimeInfectedCriticalToRecovered>()[{m_virus_variant, age}];
                time_period    = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                previous_state = InfectionState::InfectedCritical;
            }
        } break;

        case InfectionState::Dead: {
            time_in_state  = params.get<TimeInfectedCriticalToDead>()[{m_virus_variant, age}];
            time_period    = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
            previous_state = InfectionState::InfectedCritical;
        } break;

        default:
            break;
        }
        start_date = start_date - time_period;
        m_infection_course.insert(m_infection_course.begin(), {start_date, previous_state});
    }
    return start_date;
}

} // namespace abm
} // namespace mio
