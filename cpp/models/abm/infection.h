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
#ifndef EPI_ABM_INFECTION_H
#define EPI_ABM_INFECTION_H

#include "abm/time.h"
#include <vector>
#include "abm/state.h"
#include <memory>

namespace mio
{
namespace abm
{

class Virus
{
public:
    Virus(double infectivity, double mortality, double recovery)
        : m_infectivity(infectivity)
        , m_mortality(mortality)
        , m_recovery(recovery)
    {
    }

    double get_infectivity() const
    {
        return m_infectivity;
    }

    double get_mortality() const
    {
        return m_mortality;
    }

    double get_recovery() const
    {
        return m_recovery;
    }

private:
    double m_infectivity;
    double m_mortality;
    double m_recovery;
};

class ViralLoad
{
public:
    ViralLoad()
    {
        draw_viral_load();
    }
    void draw_viral_load();
    TimePoint determine_end_date(TimePoint start_date)
    {
        return start_date + TimeSpan(int(m_peak / m_incline - m_peak / m_decline));
    };
    double get_viral_load(const TimePoint t, const TimePoint start_date) const;

private:
    double m_peak;
    double m_incline;
    double m_decline; // always negative
};

class Infection
{

public:
    /**
     * create an infection for a single person.
     * @param start_date starting date of the infection
     * @param virus_type virus type of the infection
     */
    Infection(TimePoint start_date, std::shared_ptr<Virus> virus)
        : m_start_date(start_date)
        , m_virus(virus)
        , m_viral_load()
    {
        m_end_date = m_viral_load.determine_end_date(m_start_date);
    };

    /**
     * get viral load at a given time
     * @param t time point of the querry
     * @return viral load at given time point
     */
    double get_viral_load(const TimePoint t) const;

    /**
     * get infectivity at a given time
     * @param t time point of the querry
     * @return infectivity at given time point
     */
    double get_infectivity(const TimePoint t) const;

    /**
     * get virus type
     * @return virus type of the infection
     */
    VirusType get_virus_type() const;

    InfectionState get_infection_state(const TimePoint t) const;

    void draw_infection_course(InfectionState start_state = InfectionState::Susceptible);

private:
    /**
     * determine viral load course and infection course
     */

    std::shared_ptr<Virus> m_virus;
    double m_log_norm_alpha, m_log_norm_beta; // have to ask for distribution/parametrization of the infectivity
    std::vector<std::pair<mio::abm::TimePoint, mio::abm::InfectionState>> m_infection_course;
    ViralLoad m_viral_load;
    TimePoint m_start_date;
    TimePoint m_end_date;
    bool m_detected = false;
};

} // namespace abm
} // namespace mio

#endif
