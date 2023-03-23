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
    m_peak     = draws[1];
    m_incline  = draws[2];
    m_decline  = draws[3];
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
                     InfectionState start_state, bool detected)
    : m_virus_variant(virus)
    , m_viral_load(virus, age, start_date, params)
    , m_detected(detected)
{
    auto draws       = params.get<InfectivityParameters>()[{virus, age}].draw_samples();
    m_log_norm_alpha = draws[1];
    m_log_norm_beta  = draws[2];

    draw_infection_course(age, params, start_date, start_state);
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

void Infection::draw_infection_course(AgeGroup age, const GlobalInfectionParameters& params, TimePoint start_date,
                                      InfectionState start_state)
{
    auto t = start_date;
    m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, start_state));
    while ((m_infection_course.back().second != InfectionState::Recovered_Infected &&
            m_infection_course.back().second != InfectionState::Recovered_Carrier &&
            m_infection_course.back().second != InfectionState::Dead)) {
        switch (m_infection_course.back().second) {
        case InfectionState::Exposed:
            // roll out how long until carrier
            t = t + TimeSpan(params.get<IncubationPeriod>()[{m_virus_variant, age,
                                                             VaccinationState::Unvaccinated}]); // subject to change
            m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, InfectionState::Carrier));
            break;
        case InfectionState::Carrier:
            // roll out if and how long until infected
            t = t + TimeSpan(params.get<CarrierToInfected>()[{m_virus_variant, age,
                                                              VaccinationState::Unvaccinated}]); // subject to change
            m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, InfectionState::Infected));
            break;
        case InfectionState::Infected:
            // roll out next infection step
            t = t + TimeSpan(params.get<InfectedToSevere>()[{m_virus_variant, age,
                                                             VaccinationState::Unvaccinated}]); // subject to change
            m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, InfectionState::Infected_Severe));
            break;
        case InfectionState::Infected_Severe:
            // roll out next infection step
            t = t + TimeSpan(params.get<SevereToCritical>()[{m_virus_variant, age,
                                                             VaccinationState::Unvaccinated}]); // subject to change
            m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, InfectionState::Infected_Critical));
            break;
        case InfectionState::Infected_Critical:
            // roll out next infection step
            t = t + TimeSpan(params.get<CriticalToDead>()[{m_virus_variant, age,
                                                           VaccinationState::Unvaccinated}]); // subject to change
            m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, InfectionState::Dead));
            break;
        default:
            break;
        }
    }
}

} // namespace abm
} // namespace mio
