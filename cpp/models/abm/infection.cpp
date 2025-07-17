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
#include "abm/parameters.h"
#include "memilio/utils/compiler_diagnostics.h"
#include <math.h>

namespace mio
{
namespace abm
{

void Infection::initialize_viral_load(PersonalRandomNumberGenerator& rng, VirusVariant virus, AgeGroup age,
                                      const Parameters& params, TimePoint init_date, ProtectionEvent latest_protection)
{
    auto vl_params                    = params.get<ViralLoadDistributions>()[{virus, age}];
    ScalarType high_viral_load_factor = 1;

    if (latest_protection.type != ProtectionType::NoProtection) {
        high_viral_load_factor -= params.get<HighViralLoadProtectionFactor>()[{latest_protection.type, age, virus}](
            init_date.days() - latest_protection.time.days());
    }

    m_viral_load.peak    = vl_params.viral_load_peak.get(rng) * high_viral_load_factor;
    m_viral_load.incline = vl_params.viral_load_incline.get(rng);
    m_viral_load.decline = vl_params.viral_load_decline.get(rng);
    m_viral_load.end_date =
        m_viral_load.start_date +
        days(int(m_viral_load.peak / m_viral_load.incline - m_viral_load.peak / m_viral_load.decline));
}

void Infection::initialize_infectivity(PersonalRandomNumberGenerator& rng, VirusVariant virus, AgeGroup age,
                                       const Parameters& params)
{
    auto inf_params  = params.get<InfectivityDistributions>()[{virus, age}];
    m_log_norm_alpha = inf_params.infectivity_alpha.get(rng);
    m_log_norm_beta  = inf_params.infectivity_beta.get(rng);
}

void Infection::initialize_virus_shed_factor(PersonalRandomNumberGenerator& rng, VirusVariant virus, AgeGroup age,
                                             const Parameters& params)
{
    auto shedfactor_param          = params.get<VirusShedFactor>()[{virus, age}];
    m_individual_virus_shed_factor = shedfactor_param.get(rng);
}

Infection::Infection(PersonalRandomNumberGenerator& rng, VirusVariant virus, AgeGroup age, const Parameters& params,
                     TimePoint init_date, InfectionState init_state, ProtectionEvent latest_protection, bool detected)
    : m_virus_variant(virus)
    , m_detected(detected)
{
    assert(age.get() < params.get_num_groups());
    m_viral_load.start_date = draw_infection_course(rng, age, params, init_date, init_state, latest_protection);

    initialize_viral_load(rng, virus, age, params, init_date, latest_protection);
    initialize_infectivity(rng, virus, age, params);
    initialize_virus_shed_factor(rng, virus, age, params);
}

ScalarType Infection::get_viral_load(TimePoint t) const
{
    if (t < m_viral_load.start_date || t > m_viral_load.end_date) {
        return 0.0;
    }

    ScalarType time_from_start = (t - m_viral_load.start_date).days();
    ScalarType peak_time       = m_viral_load.peak / m_viral_load.incline;

    if (time_from_start <= peak_time) {
        return m_viral_load.incline * time_from_start;
    }
    else {
        return m_viral_load.peak + m_viral_load.decline * (time_from_start - peak_time);
    }
}

ScalarType Infection::get_infectivity(TimePoint t) const
{
    if (m_viral_load.start_date >= t || get_infection_state(t) == InfectionState::Exposed) {
        return 0;
    }
    return m_individual_virus_shed_factor / (1 + exp(-(m_log_norm_alpha + m_log_norm_beta * get_viral_load(t))));
}

VirusVariant Infection::get_virus_variant() const
{
    return m_virus_variant;
}

InfectionState Infection::get_infection_state(TimePoint t) const
{
    if (t < m_infection_course[0].first) {
        return InfectionState::Susceptible;
    }

    auto it = std::upper_bound(m_infection_course.begin(), m_infection_course.end(), t,
                               [](const TimePoint& s, const std::pair<TimePoint, InfectionState>& state) {
                                   return state.first > s;
                               });
    return std::prev(it)->second;
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

TimePoint Infection::draw_infection_course(PersonalRandomNumberGenerator& rng, AgeGroup age, const Parameters& params,
                                           TimePoint init_date, InfectionState init_state,
                                           ProtectionEvent latest_protection)
{
    assert(age.get() < params.get_num_groups());
    TimePoint start_date = draw_infection_course_backward(rng, age, params, init_date, init_state);
    draw_infection_course_forward(rng, age, params, init_date, init_state, latest_protection);
    return start_date;
}

StateTransition Infection::get_forward_transition(PersonalRandomNumberGenerator& rng, AgeGroup age,
                                                  const Parameters& params, InfectionState current_state,
                                                  TimePoint current_time, ProtectionEvent latest_protection) const
{
    auto& uniform_dist = UniformDistribution<double>::get_instance();
    StateTransition transition{current_state, current_state, TimeSpan{}, 0.0};

    switch (current_state) {
    case InfectionState::Exposed:
        transition.to_state = InfectionState::InfectedNoSymptoms;
        transition.duration = days(params.get<TimeExposedToNoSymptoms>()[{m_virus_variant, age}].get(rng));
        break;

    case InfectionState::InfectedNoSymptoms: {
        ScalarType p = uniform_dist(rng);
        if (p < params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}]) {
            transition.to_state = InfectionState::InfectedSymptoms;
            transition.duration = days(params.get<TimeInfectedNoSymptomsToSymptoms>()[{m_virus_variant, age}].get(rng));
        }
        else {
            transition.to_state = InfectionState::Recovered;
            transition.duration =
                days(params.get<TimeInfectedNoSymptomsToRecovered>()[{m_virus_variant, age}].get(rng));
        }
        break;
    }

    case InfectionState::InfectedSymptoms: {
        ScalarType p = uniform_dist(rng);
        ScalarType severity_protection_factor =
            get_severity_protection_factor(params, latest_protection, age, current_time);
        ScalarType severe_probability =
            (1 - severity_protection_factor) * params.get<SeverePerInfectedSymptoms>()[{m_virus_variant, age}];

        if (p < severe_probability) {
            transition.to_state = InfectionState::InfectedSevere;
            transition.duration = days(params.get<TimeInfectedSymptomsToSevere>()[{m_virus_variant, age}].get(rng));
        }
        else {
            transition.to_state = InfectionState::Recovered;
            transition.duration = days(params.get<TimeInfectedSymptomsToRecovered>()[{m_virus_variant, age}].get(rng));
        }
        break;
    }

    case InfectionState::InfectedSevere: {
        ScalarType p             = uniform_dist(rng);
        ScalarType critical_prob = params.get<CriticalPerInfectedSevere>()[{m_virus_variant, age}];
        ScalarType death_prob    = params.get<DeathsPerInfectedSevere>()[{m_virus_variant, age}];

        if (p < critical_prob) {
            transition.to_state = InfectionState::InfectedCritical;
            transition.duration = days(params.get<TimeInfectedSevereToCritical>()[{m_virus_variant, age}].get(rng));
        }
        else if (p < critical_prob + death_prob) {
            transition.to_state = InfectionState::Dead;
            transition.duration = days(params.get<TimeInfectedSevereToDead>()[{m_virus_variant, age}].get(rng));
        }
        else {
            transition.to_state = InfectionState::Recovered;
            transition.duration = days(params.get<TimeInfectedSevereToRecovered>()[{m_virus_variant, age}].get(rng));
        }
        break;
    }

    case InfectionState::InfectedCritical: {
        ScalarType p = uniform_dist(rng);
        if (p < params.get<DeathsPerInfectedCritical>()[{m_virus_variant, age}]) {
            transition.to_state = InfectionState::Dead;
            transition.duration = days(params.get<TimeInfectedCriticalToDead>()[{m_virus_variant, age}].get(rng));
        }
        else {
            transition.to_state = InfectionState::Recovered;
            transition.duration = days(params.get<TimeInfectedCriticalToRecovered>()[{m_virus_variant, age}].get(rng));
        }
        break;
    }

    default:
        break;
    }

    return transition;
}

StateTransition Infection::get_backward_transition(PersonalRandomNumberGenerator& rng, AgeGroup age,
                                                   const Parameters& params, InfectionState current_state) const
{
    StateTransition transition{current_state, current_state, TimeSpan{}, 0.0};

    switch (current_state) {
    case InfectionState::InfectedNoSymptoms:
        transition.to_state = InfectionState::Exposed;
        transition.duration = days(params.get<TimeExposedToNoSymptoms>()[{m_virus_variant, age}].get(rng));
        break;

    case InfectionState::InfectedSymptoms:
        transition.to_state = InfectionState::InfectedNoSymptoms;
        transition.duration = days(params.get<TimeInfectedNoSymptomsToSymptoms>()[{m_virus_variant, age}].get(rng));
        break;

    case InfectionState::InfectedSevere:
        transition.to_state = InfectionState::InfectedSymptoms;
        transition.duration = days(params.get<TimeInfectedSymptomsToSevere>()[{m_virus_variant, age}].get(rng));
        break;

    case InfectionState::InfectedCritical:
        transition.to_state = InfectionState::InfectedSevere;
        transition.duration = days(params.get<TimeInfectedSevereToCritical>()[{m_virus_variant, age}].get(rng));
        break;

    case InfectionState::Recovered:
        transition = get_recovered_backward_transition(rng, age, params);
        break;

    case InfectionState::Dead:
        transition = get_dead_backward_transition(rng, age, params);
        break;

    default:
        break;
    }

    return transition;
}

StateTransition Infection::get_recovered_backward_transition(PersonalRandomNumberGenerator& rng, AgeGroup age,
                                                             const Parameters& params) const
{
    auto& uniform_dist = UniformDistribution<double>::get_instance();
    ScalarType p       = uniform_dist(rng);

    // Compute death probability to factor it out
    ScalarType p_death = calculate_death_probability(age, params);

    ScalarType symptoms_prob = params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}];
    ScalarType severe_prob   = params.get<SeverePerInfectedSymptoms>()[{m_virus_variant, age}];
    ScalarType critical_prob = params.get<CriticalPerInfectedSevere>()[{m_virus_variant, age}];

    StateTransition transition{InfectionState::Recovered, InfectionState::InfectedNoSymptoms, TimeSpan{}, 0.0};

    if (p > symptoms_prob / (1 - p_death)) {
        transition.to_state = InfectionState::InfectedNoSymptoms;
        transition.duration = days(params.get<TimeInfectedNoSymptomsToRecovered>()[{m_virus_variant, age}].get(rng));
    }
    else if (p > symptoms_prob * severe_prob / (1 - p_death)) {
        transition.to_state = InfectionState::InfectedSymptoms;
        transition.duration = days(params.get<TimeInfectedSymptomsToRecovered>()[{m_virus_variant, age}].get(rng));
    }
    else if (p > symptoms_prob * severe_prob * critical_prob / (1 - p_death)) {
        transition.to_state = InfectionState::InfectedSevere;
        transition.duration = days(params.get<TimeInfectedSevereToRecovered>()[{m_virus_variant, age}].get(rng));
    }
    else {
        transition.to_state = InfectionState::InfectedCritical;
        transition.duration = days(params.get<TimeInfectedCriticalToRecovered>()[{m_virus_variant, age}].get(rng));
    }

    return transition;
}

StateTransition Infection::get_dead_backward_transition(PersonalRandomNumberGenerator& rng, AgeGroup age,
                                                        const Parameters& params) const
{
    auto& uniform_dist = UniformDistribution<double>::get_instance();
    ScalarType p       = uniform_dist(rng);

    StateTransition transition{InfectionState::Dead, InfectionState::InfectedSevere, TimeSpan{}, 0.0};

    if (p < params.get<DeathsPerInfectedSevere>()[{m_virus_variant, age}]) {
        transition.to_state = InfectionState::InfectedSevere;
        transition.duration = days(params.get<TimeInfectedSevereToDead>()[{m_virus_variant, age}].get(rng));
    }
    else {
        transition.to_state = InfectionState::InfectedCritical;
        transition.duration = days(params.get<TimeInfectedCriticalToDead>()[{m_virus_variant, age}].get(rng));
    }

    return transition;
}

ScalarType Infection::calculate_death_probability(AgeGroup age, const Parameters& params) const
{
    return params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}] *
           params.get<SeverePerInfectedSymptoms>()[{m_virus_variant, age}] *
           params.get<CriticalPerInfectedSevere>()[{m_virus_variant, age}] *
           params.get<DeathsPerInfectedCritical>()[{m_virus_variant, age}];
}

ScalarType Infection::get_severity_protection_factor(const Parameters& params, ProtectionEvent latest_protection,
                                                     AgeGroup age, TimePoint current_time) const
{
    if (latest_protection.type == ProtectionType::NoProtection) {
        return 0.0;
    }
    return params.get<SeverityProtectionFactor>()[{latest_protection.type, age, m_virus_variant}](
        current_time.days() - latest_protection.time.days());
}

void Infection::draw_infection_course_forward(PersonalRandomNumberGenerator& rng, AgeGroup age,
                                              const Parameters& params, TimePoint init_date, InfectionState start_state,
                                              ProtectionEvent latest_protection)
{
    assert(age.get() < params.get_num_groups());

    TimePoint current_time       = init_date;
    InfectionState current_state = start_state;
    m_infection_course.push_back({current_time, current_state});

    while (current_state != InfectionState::Recovered && current_state != InfectionState::Dead) {
        StateTransition transition =
            get_forward_transition(rng, age, params, current_state, current_time, latest_protection);
        current_time += transition.duration;
        current_state = transition.to_state;
        m_infection_course.push_back({current_time, current_state});
    }
}

TimePoint Infection::draw_infection_course_backward(PersonalRandomNumberGenerator& rng, AgeGroup age,
                                                    const Parameters& params, TimePoint init_date,
                                                    InfectionState init_state)
{
    assert(age.get() < params.get_num_groups());

    TimePoint current_time       = init_date;
    InfectionState current_state = init_state;

    while (current_state != InfectionState::Exposed) {
        StateTransition transition = get_backward_transition(rng, age, params, current_state);
        current_time -= transition.duration;
        current_state = transition.to_state;
        m_infection_course.insert(m_infection_course.begin(), {current_time, current_state});
    }

    return current_time;
}

} // namespace abm
} // namespace mio
