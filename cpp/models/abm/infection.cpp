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

ViralLoad::ViralLoad(const TimePoint& start_date)
    : m_start_date(start_date)
{
    draw_viral_load();
}

void ViralLoad::draw_viral_load()
{
    // These numbers are subject to change, They are going to be based on distributions backed from data.
    m_peak     = 5.0;
    m_incline  = 1.0;
    m_decline  = -1.0;
    m_end_date = m_start_date + TimeSpan(int(m_peak / m_incline - m_peak / m_decline));
}

double ViralLoad::get_viral_load(const TimePoint& t) const
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

Infection::Infection(std::shared_ptr<Virus> virus, const TimePoint& start_date, const InfectionState& start_state,
                     const bool detected)
    : m_virus(virus)
    , m_viral_load(start_date)
    , m_detected(detected)
{
    draw_infection_course(start_date, start_state);
};

double Infection::get_viral_load(const TimePoint& t) const
{
    return m_viral_load.get_viral_load(t);
}

double Infection::get_infectivity(const TimePoint& t) const
{
    return 1 / (1 + exp(-(m_log_norm_alpha + m_log_norm_beta * get_viral_load(t))));
}

const Virus Infection::get_virus_type() const
{
    return *m_virus.get();
}

const InfectionState& Infection::get_infection_state(const TimePoint& t) const
{
    return (*std::prev(std::lower_bound(m_infection_course.begin(), m_infection_course.end(), t,
                                        [](std::pair<TimePoint, InfectionState> state, const TimePoint& t) {
                                            return state.first <= t;
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

void Infection::draw_infection_course(const TimePoint& start_date, const InfectionState& start_state)
{
    auto t = start_date;
    m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, start_state));
    while ((m_infection_course.back().second != InfectionState::Recovered_Infected ||
            m_infection_course.back().second != InfectionState::Recovered_Carrier ||
            m_infection_course.back().second != InfectionState::Dead)) {
        switch (m_infection_course.back().second) {
        case InfectionState::Exposed:
            // roll out how long until carrier
            t = t + days(5); // subject to change
            m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, InfectionState::Carrier));
            break;
        case InfectionState::Carrier:
            // roll out if and how long until infected
            t = t + days(4); // subject to change
            m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, InfectionState::Infected));
            break;
        case InfectionState::Infected:
            // roll out next infection step
            t = t + days(3); // subject to change
            m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, InfectionState::Infected_Severe));
            break;
        case InfectionState::Infected_Severe:
            // roll out next infection step
            t = t + days(2); // subject to change
            m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, InfectionState::Infected_Critical));
            break;
        case InfectionState::Infected_Critical:
            // roll out next infection step
            t = t + days(1); // subject to change
            m_infection_course.push_back(std::pair<TimePoint, InfectionState>(t, InfectionState::Dead));
            break;
        default:
            break;
        }
    }
}

} // namespace abm
} // namespace mio
