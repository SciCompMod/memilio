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

TimePoint mio::abm::ViralLoad::determine_end_date(TimePoint start_date)
{
    return start_date + TimeSpan(int(m_peak / m_incline - m_peak / m_decline));
}

double mio::abm::Infection::get_infectivity(const TimePoint t) const
{
    return 1 / (1 + exp(-(m_log_norm_alpha + m_log_norm_beta * get_viral_load(t))));
}

double mio::abm::Infection::get_viral_load(const TimePoint t) const
{
    if (t >= m_start_date && t <= m_end_date) {
        return m_viral_load.get_viral_load(t, m_start_date);
    }
    else {
        return 0;
    }
}

// double mio::abm::ViralLoad::get_viral_load(const TimePoint t, const TimePoint start_date) const
// {
//     if (t <= m_peak / m_incline + start_date) {
//         return m_incline * (t - start_date);
//     }
//     else {
//         return m_peak + m_decline * (t - m_peak / m_incline - start_date);
//     }
// }
// void mio::abm::ViralLoad::draw_viral_load()
// {
//     // These numbers are subject to change, They are going to be based on distributions backed from data.
//     m_peak    = 5.0;
//     m_incline = 1.0;
//     m_decline = -1.0;
// }

// InfectionState mio::abm::Infection::get_infection_state(const TimePoint t) const
// {
//     return (*std::prev(std::lower_bound(m_infection_course.begin(), m_infection_course.end(), t))).second;
// }

// void mio::abm::Infection::draw_infection_course(InfectionState start_state)
// {
//     m_infection_course.push_back(std::pair<TimePoint, InfectionState>(m_start_date, start_state));
//     while ((m_infection_course.back().second != InfectionState::Recovered_Infected ||
//             m_infection_course.back().second != InfectionState::Recovered_Carrier ||
//             m_infection_course.back().second != InfectionState::Dead)) {
//         switch (m_infection_course.back().second) {
//         case InfectionState::Susceptible:

//             break;
//         case InfectionState::Infected:
//             break;
//         case InfectionState::Carrier:
//             break;
//         case InfectionState::Infected_Severe:
//             break;
//         case InfectionState::Infected_Critical:
//             break;
//         default:
//             break;
//         }
//     }
//}