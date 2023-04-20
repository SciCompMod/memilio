/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: David Kerkmann, Sascha Korf
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

Infection::Infection(VirusVariant virus, AgeGroup age, const GlobalInfectionParameters& params, TimePoint init_date,
                     InfectionState init_state, bool detected)
    : m_virus_variant(virus)
    , m_detected(detected)
{
    m_viral_load.start_date = draw_infection_course(age, params, init_date, init_state);

    auto vl_params    = params.get<ViralLoadDistributions>()[{virus, age}]; // TODO: change vaccination state
    m_viral_load.peak = vl_params.viral_load_peak.get_distribution_instance()(vl_params.viral_load_peak.params);
    m_viral_load.incline =
        vl_params.viral_load_incline.get_distribution_instance()(vl_params.viral_load_incline.params);
    m_viral_load.decline =
        vl_params.viral_load_decline.get_distribution_instance()(vl_params.viral_load_decline.params);
    m_viral_load.end_date =
        m_viral_load.start_date +
        TimeSpan(int(m_viral_load.peak / m_viral_load.incline - m_viral_load.peak / m_viral_load.decline));

    m_viral_load.end_date =
        m_viral_load.start_date +
        TimeSpan(int(m_viral_load.peak / m_viral_load.incline - m_viral_load.peak / m_viral_load.decline));

    auto inf_params  = params.get<InfectivityDistributions>()[{virus, age}];
    m_log_norm_alpha = inf_params.infectivity_alpha.get_distribution_instance()(inf_params.infectivity_alpha.params);
    m_log_norm_beta  = inf_params.infectivity_beta.get_distribution_instance()(inf_params.infectivity_beta.params);
}

ScalarType Infection::get_viral_load(TimePoint t) const
{
    if (t >= m_viral_load.start_date && t <= m_viral_load.end_date) {
        if (t.seconds() <= m_viral_load.start_date.seconds() + m_viral_load.peak / m_viral_load.incline) {
            return m_viral_load.incline * (t - m_viral_load.start_date).seconds();
        }
        else {
            return m_viral_load.peak + m_viral_load.decline * (t.seconds() - m_viral_load.peak / m_viral_load.incline -
                                                               m_viral_load.start_date.seconds());
        }
    }
    else {
        return 0.;
    }
}

ScalarType Infection::get_infectivity(TimePoint t) const
{
    if (m_viral_load.start_date >= t)
        return 0;
    return 1 / (1 + exp(-(m_log_norm_alpha + m_log_norm_beta * get_viral_load(t))));
}

const VirusVariant& Infection::get_virus_variant() const
{
    return m_virus_variant;
}

const InfectionState& Infection::get_infection_state(TimePoint t) const
{
    return (*std::prev(std::lower_bound(m_infection_course.begin(), m_infection_course.end(), t,
                                        [](std::pair<TimePoint, InfectionState> state, const TimePoint& s) {
                                            return state.first <= s;
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

TimePoint Infection::draw_infection_course(AgeGroup age, const GlobalInfectionParameters& params, TimePoint init_date,
                                           InfectionState init_state)
{
    TimePoint start_date = draw_infection_course_backward(age, params, init_date, init_state);
    draw_infection_course_forward(age, params, init_date, init_state);
    return start_date;
}

void Infection::draw_infection_course_forward(AgeGroup age, const GlobalInfectionParameters& params,
                                              TimePoint init_date, InfectionState start_state)
{
    auto t = init_date;
    TimeSpan time_period{}; // time period for current infection state
    InfectionState next_state{start_state}; // next state to enter
    m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, next_state));
    auto uniform_dist = UniformDistribution<double>::get_instance();
    ScalarType v; // random draws
    while ((next_state != InfectionState::Recovered_Infected && next_state != InfectionState::Recovered_Carrier &&
            next_state != InfectionState::Dead)) {
        switch (next_state) {
        case InfectionState::Exposed:
            // roll out how long until carrier
            time_period = TimeSpan((int)params.get<IncubationPeriod>()[{m_virus_variant, age}]); // subject to change
            next_state  = InfectionState::Carrier;
            break;
        case InfectionState::Carrier:
            // roll out next infection step
            v = uniform_dist();
            if (v < 0.5) { // TODO: subject to change
                time_period =
                    TimeSpan((int)params.get<CarrierToInfected>()[{m_virus_variant, age}]); // TODO: subject to change
                next_state = InfectionState::Infected;
            }
            else {
                time_period =
                    TimeSpan((int)params.get<CarrierToRecovered>()[{m_virus_variant, age}]); // TODO: subject to change
                next_state = InfectionState::Recovered_Carrier;
            }

            break;
        case InfectionState::Infected:
            // roll out next infection step
            v = uniform_dist();
            if (v < 0.5) { // TODO: subject to change
                time_period =
                    TimeSpan((int)params.get<InfectedToSevere>()[{m_virus_variant, age}]); // TODO: subject to change
                next_state = InfectionState::Infected_Severe;
            }
            else {
                time_period =
                    TimeSpan((int)params.get<InfectedToRecovered>()[{m_virus_variant, age}]); // TODO: subject to change
                next_state = InfectionState::Recovered_Infected;
            }
            break;
        case InfectionState::Infected_Severe:
            // roll out next infection step
            v = uniform_dist();
            if (v < 0.5) { // TODO: subject to change
                time_period =
                    TimeSpan((int)params.get<SevereToCritical>()[{m_virus_variant, age}]); // TODO: subject to change
                next_state = InfectionState::Infected_Critical;
            }
            else {
                time_period =
                    TimeSpan((int)params.get<SevereToRecovered>()[{m_virus_variant, age}]); // TODO: subject to change
                next_state = InfectionState::Recovered_Infected;
            }
            break;
        case InfectionState::Infected_Critical:
            // roll out next infection step
            v = uniform_dist();
            if (v < 0.5) { // TODO: subject to change
                time_period =
                    TimeSpan((int)params.get<CriticalToDead>()[{m_virus_variant, age}]); // TODO: subject to change
                next_state = InfectionState::Dead;
            }
            else {
                time_period =
                    TimeSpan((int)params.get<CriticalToRecovered>()[{m_virus_variant, age}]); // TODO: subject to change
                next_state = InfectionState::Recovered_Infected;
            }
            break;
        default:
            break;
        }
        t = t + time_period;
        m_infection_course.push_back({t, next_state});
    }
}

TimePoint Infection::draw_infection_course_backward(AgeGroup age, const GlobalInfectionParameters& params,
                                                    TimePoint init_date, InfectionState init_state)
{
    assert(init_state != InfectionState::Dead && "Cannot initialize dead person.");

    auto start_date = init_date;
    TimeSpan time_period{}; // time period for current infection state
    InfectionState previous_state{init_state}; // next state to enter
    auto uniform_dist = UniformDistribution<double>::get_instance();
    ScalarType v; // random draws
    while ((previous_state != InfectionState::Exposed)) {
        switch (previous_state) {

        case InfectionState::Carrier:
            time_period =
                TimeSpan((int)params.get<IncubationPeriod>()[{m_virus_variant, age}]); // TODO: subject to change
            previous_state = InfectionState::Exposed;
            break;

        case InfectionState::Infected:
            time_period =
                TimeSpan((int)params.get<CarrierToInfected>()[{m_virus_variant, age}]); // TODO: subject to change
            previous_state = InfectionState::Carrier;
            break;

        case InfectionState::Infected_Severe:
            time_period =
                TimeSpan((int)params.get<InfectedToSevere>()[{m_virus_variant, age}]); // TODO: subject to change
            previous_state = InfectionState::Infected;
            break;

        case InfectionState::Infected_Critical:
            time_period =
                TimeSpan((int)params.get<SevereToCritical>()[{m_virus_variant, age}]); // TODO: subject to change
            previous_state = InfectionState::Infected_Severe;
            break;

        case InfectionState::Recovered_Carrier:
            time_period =
                TimeSpan((int)params.get<CarrierToRecovered>()[{m_virus_variant, age}]); // TODO: subject to change
            previous_state = InfectionState::Carrier;
            break;

        case InfectionState::Recovered_Infected:
            // roll out next infection step
            v = uniform_dist();
            if (v < 1 / 3) { // TODO: subject to change
                time_period =
                    TimeSpan((int)params.get<InfectedToRecovered>()[{m_virus_variant, age}]); // TODO: subject to change
                previous_state = InfectionState::Infected;
            }
            else if (v < 2 / 3) {
                time_period =
                    TimeSpan((int)params.get<SevereToRecovered>()[{m_virus_variant, age}]); // TODO: subject to change
                previous_state = InfectionState::Infected_Severe;
            }
            else {
                time_period =
                    TimeSpan((int)params.get<CriticalToRecovered>()[{m_virus_variant, age}]); // TODO: subject to change
                previous_state = InfectionState::Infected_Critical;
            }
            break;

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
