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

ViralLoad::ViralLoad(VirusVariant virus, AgeGroup age, TimePoint start_date, GlobalInfectionParameters& params)
    : m_start_date(start_date)
{
    draw_viral_load(virus, age, VaccinationState::Unvaccinated,
                    params); // to be changed once the immunity level is implemented
}

void ViralLoad::draw_viral_load(VirusVariant virus, AgeGroup age, VaccinationState vaccination_state,
                                GlobalInfectionParameters& params)
{
    auto draws = params.get<ViralLoadParameters>()[{virus, age, vaccination_state}].draw_samples();
    m_peak     = draws[0];
    m_incline  = draws[1];
    m_decline  = draws[2];
    m_end_date = m_start_date + TimeSpan(int(m_peak / m_incline - m_peak / m_decline));
}

ScalarType ViralLoad::get_viral_load(TimePoint t) const
{
    if (t >= m_start_date && t <= m_end_date) {
        if (t.seconds() <= m_start_date.seconds() + m_peak / m_incline) {
            return m_incline * (t - m_start_date).seconds();
        }
        else {
            return m_peak + m_decline * (t.seconds() - m_peak / m_incline - m_start_date.seconds());
        }
    }
    else {
        return 0.;
    }
}

Infection::Infection(VirusVariant virus, AgeGroup age, GlobalInfectionParameters& params, TimePoint start_date,
                     bool detected)
    : m_virus_variant(virus)
    , m_viral_load(virus, age, start_date, params)
    , m_detected(detected)
{
    draw_infection_course(age, params, start_date);

    auto draws       = params.get<InfectivityParameters>()[{virus, age}].draw_samples();
    m_log_norm_alpha = draws[0];
    m_log_norm_beta  = draws[1];
}

ScalarType Infection::get_infectivity(TimePoint t) const
{
    return 1 / (1 + exp(-(m_log_norm_alpha + m_log_norm_beta * m_viral_load.get_viral_load(t))));
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

void Infection::draw_infection_course(AgeGroup age, const GlobalInfectionParameters& params, TimePoint start_date)
{
    auto t = start_date;
    TimeSpan time_period{}; // time period for current infection state
    InfectionState next_state{}; // next state to enter
    m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, InfectionState::Exposed));
    auto uniform_dist = UniformDistribution<double>::get_instance();
    ScalarType v; // random draws
    while ((m_infection_course.back().second != InfectionState::Recovered_Infected &&
            m_infection_course.back().second != InfectionState::Recovered_Carrier &&
            m_infection_course.back().second != InfectionState::Dead)) {
        switch (m_infection_course.back().second) {
        case InfectionState::Exposed:
            // roll out how long until carrier
            time_period = TimeSpan(params.get<IncubationPeriod>()[{
                m_virus_variant, age, VaccinationState::Unvaccinated}]); // subject to change
            next_state  = InfectionState::Carrier;
            break;
        case InfectionState::Carrier:
            // roll out next infection step
            v = uniform_dist();
            if (v < 0.5) { // subject to change
                time_period = TimeSpan(params.get<CarrierToInfected>()[{
                    m_virus_variant, age, VaccinationState::Unvaccinated}]); // subject to change
                next_state  = InfectionState::Infected;
            }
            else {
                time_period = TimeSpan(params.get<CarrierToRecovered>()[{
                    m_virus_variant, age, VaccinationState::Unvaccinated}]); // subject to change
                next_state  = InfectionState::Recovered_Carrier;
            }

            break;
        case InfectionState::Infected:
            // roll out next infection step
            v = uniform_dist();
            if (v < 0.5) { // subject to change
                time_period = TimeSpan(params.get<InfectedToSevere>()[{
                    m_virus_variant, age, VaccinationState::Unvaccinated}]); // subject to change
                next_state  = InfectionState::Infected_Severe;
            }
            else {
                time_period = TimeSpan(params.get<InfectedToRecovered>()[{
                    m_virus_variant, age, VaccinationState::Unvaccinated}]); // subject to change
                next_state  = InfectionState::Recovered_Infected;
            }
            break;
        case InfectionState::Infected_Severe:
            // roll out next infection step
            v = uniform_dist();
            if (v < 0.5) { // subject to change
                time_period = TimeSpan(params.get<SevereToCritical>()[{
                    m_virus_variant, age, VaccinationState::Unvaccinated}]); // subject to change
                next_state  = InfectionState::Infected_Critical;
            }
            else {
                time_period = TimeSpan(params.get<SevereToRecovered>()[{
                    m_virus_variant, age, VaccinationState::Unvaccinated}]); // subject to change
                next_state  = InfectionState::Recovered_Infected;
            }
            break;
        case InfectionState::Infected_Critical:
            // roll out next infection step
            v = uniform_dist();
            if (v < 0.5) { // subject to change
                time_period = TimeSpan(params.get<CriticalToDead>()[{
                    m_virus_variant, age, VaccinationState::Unvaccinated}]); // subject to change
                next_state  = InfectionState::Dead;
            }
            else {
                time_period = TimeSpan(params.get<CriticalToRecovered>()[{
                    m_virus_variant, age, VaccinationState::Unvaccinated}]); // subject to change
                next_state  = InfectionState::Recovered_Infected;
            }
            break;
        default:
            break;
        }
        t = t + time_period;
        m_infection_course.push_back({t, next_state});
    }
}

} // namespace abm
} // namespace mio
