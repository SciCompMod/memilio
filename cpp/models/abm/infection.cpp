/* 
* Copyright (C) 2020-2024 MEmilio
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
#include <math.h>

namespace mio
{
namespace abm
{

Infection::Infection(Person::RandomNumberGenerator& rng, VirusVariant virus, AgeGroup age, const Parameters& params,
                     TimePoint init_date, InfectionState init_state, std::pair<ExposureType, TimePoint> latest_exposure,
                     bool detected)
    : m_virus_variant(virus)
    , m_detected(detected)
{
    assert(age.get() < params.get_num_groups());
    m_viral_load.start_date = draw_infection_course(rng, age, params, init_date, init_state, latest_exposure);

    auto vl_params                    = params.get<ViralLoadDistributions>()[{virus, age}];
    ScalarType high_viral_load_factor = 1;
    if (latest_exposure.first != ExposureType::NoProtection) {
        high_viral_load_factor -=
            params.get<HighViralLoadProtectionFactor>()(init_date.days() - latest_exposure.second.days());
    }
    m_viral_load.peak = vl_params.viral_load_peak.get_distribution_instance()(rng, vl_params.viral_load_peak.params) *
                        high_viral_load_factor;
    m_viral_load.incline =
        vl_params.viral_load_incline.get_distribution_instance()(rng, vl_params.viral_load_incline.params);
    m_viral_load.decline =
        vl_params.viral_load_decline.get_distribution_instance()(rng, vl_params.viral_load_decline.params);
    m_viral_load.end_date =
        m_viral_load.start_date +
        days(int(m_viral_load.peak / m_viral_load.incline - m_viral_load.peak / m_viral_load.decline));

    auto inf_params = params.get<InfectivityDistributions>()[{virus, age}];
    m_log_norm_alpha =
        inf_params.infectivity_alpha.get_distribution_instance()(rng, inf_params.infectivity_alpha.params);
    m_log_norm_beta = inf_params.infectivity_beta.get_distribution_instance()(rng, inf_params.infectivity_beta.params);

    auto shedfactor_param          = params.get<VirusShedFactor>()[{virus, age}];
    m_individual_virus_shed_factor = shedfactor_param.get_distribution_instance()(rng, shedfactor_param.params);
}

ScalarType Infection::get_viral_load(TimePoint t) const
{
    if (t.days() <= m_viral_load.start_date.days() + m_viral_load.peak / m_viral_load.incline) {
        return m_viral_load.incline * (t - m_viral_load.start_date).days();
    }
    else {
        return m_viral_load.peak + m_viral_load.decline * (t.days() - m_viral_load.peak / m_viral_load.incline -
                                                           m_viral_load.start_date.days());
    }
}

ScalarType Infection::get_viral_shed(TimePoint t) const
{
    if (m_viral_load.start_date >= t || get_infection_state(t) == InfectionState::Exposed)
        return 0;
    return m_individual_virus_shed_factor / (1 + exp(-(m_log_norm_alpha + m_log_norm_beta * get_viral_load(t))));
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

TimePoint Infection::draw_infection_course(Person::RandomNumberGenerator& rng, AgeGroup age, const Parameters& params,
                                           TimePoint init_date, InfectionState init_state,
                                           std::pair<ExposureType, TimePoint> latest_protection)
{
    assert(age.get() < params.get_num_groups());
    TimePoint start_of_init_state =
        draw_infection_course_forward(rng, age, params, init_date, init_state, latest_protection);
    TimePoint start_date = draw_infection_course_backward(rng, age, params, start_of_init_state, init_state);

    return start_date;
}

TimePoint Infection::draw_infection_course_forward(Person::RandomNumberGenerator& rng, AgeGroup age,
                                                   const Parameters& params, TimePoint init_date,
                                                   InfectionState start_state,
                                                   std::pair<ExposureType, TimePoint> latest_exposure)
{
    assert(age.get() < params.get_num_groups());
    auto t = init_date;
    TimeSpan time_period{}; // time period for current infection state
    auto time_in_state = params.get<IncubationPeriod>()[{
        m_virus_variant, age}]; // time distribution parameters for current infection state
    InfectionState next_state{start_state}; // next state to enter

    auto& uniform_dist            = UniformDistribution<double>::get_instance();
    ScalarType p                  = 1; // uniform random draws from [0, 1]
    bool init_state               = false;
    TimePoint start_of_init_state = init_date;

    // helper lambda for init_state
    auto determine_time_period = [&]() {
        if (!init_state) {
            // randomly initialize the infetion when it is not initialized at the time of transmission
            p = uniform_dist(rng); // a random amount of time that has passed since begin of the initial InfectionState
            auto time_draw      = time_in_state.get_distribution_instance()(rng, time_in_state.params);
            time_period         = days(p * time_draw);
            start_of_init_state = t - days((1 - p) * time_draw);
            init_state          = true;
        }
        else {
            time_period = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
        }
    };

    // randomly draw the rest of the InfectionState tree
    while ((next_state != InfectionState::Recovered && next_state != InfectionState::Dead)) {
        switch (next_state) {

        case InfectionState::Susceptible: {
            next_state  = InfectionState::Exposed;
            time_period = mio::abm::hours(0);
            init_state  = true;
        } break;
        case InfectionState::Exposed: {
            // roll out how long until infected without symptoms
            time_in_state = params.get<IncubationPeriod>()[{m_virus_variant, age}];
            determine_time_period();
            next_state = InfectionState::InfectedNoSymptoms;
        } break;
        case InfectionState::InfectedNoSymptoms: {
            // roll out next infection step

            p = uniform_dist(rng);
            if (p < params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}]) {
                time_in_state = params.get<TimeInfectedNoSymptomsToSymptoms>()[{m_virus_variant, age}];
                determine_time_period();
                next_state = InfectionState::InfectedSymptoms;
            }
            else {
                time_in_state = params.get<TimeInfectedNoSymptomsToRecovered>()[{m_virus_variant, age}];
                determine_time_period();
                next_state = InfectionState::Recovered;
            }
        } break;
        case InfectionState::InfectedSymptoms: {
            // roll out next infection step

            ScalarType severity_protection_factor = 0.;
            p                                     = uniform_dist(rng);
            if (latest_exposure.first != ExposureType::NoProtection) {
                severity_protection_factor =
                    params.get<SeverityProtectionFactor>()[{latest_exposure.first, age, m_virus_variant}](
                        t.days() - latest_exposure.second.days());
            }
            if (p <
                (1 - severity_protection_factor) * params.get<SeverePerInfectedSymptoms>()[{m_virus_variant, age}]) {
                time_in_state = params.get<TimeInfectedSymptomsToSevere>()[{m_virus_variant, age}];
                determine_time_period();
                next_state = InfectionState::InfectedSevere;
            }
            else {
                time_in_state = params.get<TimeInfectedSymptomsToRecovered>()[{m_virus_variant, age}];
                determine_time_period();
                next_state = InfectionState::Recovered;
            }
        } break;

        case InfectionState::InfectedSevere: {
            // roll out next infection step

            p = uniform_dist(rng);
            if (p < params.get<CriticalPerInfectedSevere>()[{m_virus_variant, age}]) {
                time_in_state = params.get<TimeInfectedSevereToCritical>()[{m_virus_variant, age}];
                determine_time_period();
                next_state = InfectionState::InfectedCritical;
            }
            else {
                time_in_state = params.get<TimeInfectedSevereToRecovered>()[{m_virus_variant, age}];
                determine_time_period();
                next_state = InfectionState::Recovered;
            }
        } break;

        case InfectionState::InfectedCritical: {
            // roll out next infection step

            p = uniform_dist(rng);
            if (p < params.get<DeathsPerInfectedCritical>()[{m_virus_variant, age}]) {
                time_in_state = params.get<TimeInfectedCriticalToDead>()[{m_virus_variant, age}];
                determine_time_period();
                next_state = InfectionState::Dead;
            }
            else {
                time_in_state = params.get<TimeInfectedCriticalToRecovered>()[{m_virus_variant, age}];
                determine_time_period();
                next_state = InfectionState::Recovered;
            }
        } break;

        default:
            break;
        }
        t = t + time_period;
        m_infection_course.push_back({t, next_state});
    }
    return start_of_init_state;
}

TimePoint Infection::draw_infection_course_backward(Person::RandomNumberGenerator& rng, AgeGroup age,
                                                    const Parameters& params, TimePoint init_date,
                                                    InfectionState init_state)
{
    assert(age.get() < params.get_num_groups());
    auto t = init_date;
    TimeSpan time_period{}; // time period for current infection state
    auto time_in_state = params.get<IncubationPeriod>()[{
        m_virus_variant, age}]; // time distribution parameters for current infection state
    InfectionState previous_state{init_state}; // previous state to enter
    auto& uniform_dist = UniformDistribution<double>::get_instance();
    ScalarType p; // uniform random draws from [0, 1]

    if (previous_state != InfectionState::Susceptible) {
        m_infection_course.insert(m_infection_course.begin(), {t, previous_state});
    }

    while ((previous_state != InfectionState::Susceptible && previous_state != InfectionState::Exposed)) {
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
            p = uniform_dist(rng);
            // compute correct probabilities while factoring out the chance to die
            auto p_death = params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}] *
                           params.get<SeverePerInfectedSymptoms>()[{m_virus_variant, age}] *
                           params.get<CriticalPerInfectedSevere>()[{m_virus_variant, age}] *
                           params.get<DeathsPerInfectedCritical>()[{m_virus_variant, age}];
            if (p > (1 - params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}]) / (1 - p_death)) {
                time_in_state  = params.get<TimeInfectedNoSymptomsToRecovered>()[{m_virus_variant, age}];
                time_period    = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                previous_state = InfectionState::InfectedNoSymptoms;
            }
            else if (p > params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}] *
                             (1 - params.get<SeverePerInfectedSymptoms>()[{m_virus_variant, age}]) / (1 - p_death)) {
                time_in_state  = params.get<TimeInfectedSymptomsToRecovered>()[{m_virus_variant, age}];
                time_period    = days(time_in_state.get_distribution_instance()(rng, time_in_state.params));
                previous_state = InfectionState::InfectedSymptoms;
            }
            else if (p > params.get<SymptomsPerInfectedNoSymptoms>()[{m_virus_variant, age}] *
                             params.get<SeverePerInfectedSymptoms>()[{m_virus_variant, age}] *
                             (1 - params.get<CriticalPerInfectedSevere>()[{m_virus_variant, age}]) / (1 - p_death)) {
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
        t = t - time_period;
        m_infection_course.insert(m_infection_course.begin(), {t, previous_state});
    }
    return t;
}

} // namespace abm
} // namespace mio
