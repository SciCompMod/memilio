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
#include <boost/optional.hpp>

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
    ViralLoad();
    void draw_viral_load();
    TimePoint determine_end_date(const TimePoint& start_date);
    double get_viral_load(const TimePoint& t, const TimePoint& start_date) const;

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
     * @param virus virus type of the infection
     * @param start_date starting date of the infection
     * @param start_state starting infection state of the person, default susceptible. Only for initialization.
     */
    Infection(std::shared_ptr<Virus> virus, const TimePoint& start_date, const InfectionState& start_state = InfectionState::Exposed, const bool detected = false);

    /**
     * get viral load at a given time
     * @param t time point of the querry
     * @return viral load at given time point
     * Computed depending on indiviudal hat function depending on viral load parameters
     * corresponding to https://www.science.org/doi/full/10.1126/science.abi5273
     *
     */
    double get_viral_load(const TimePoint& t) const;

    /**
     * get infectivity at a given time
     * @param t time point of the querry
     * @return infectivity at given time point.
     * Computed depending on current viral load and individual invlogit function of each person
     * corresponding to https://www.science.org/doi/full/10.1126/science.abi5273
     */
    double get_infectivity(const TimePoint& t) const;

    /**
     * get virus type
     * @return virus type of the infection
     */
    const Virus get_virus_type() const;

    const InfectionState& get_infection_state() const;
    
    const InfectionState& get_infection_state(const TimePoint& t) const;
    
    InfectionState& get_infection_state(const TimePoint& t);
    
    void update_infection_state(const TimePoint& t);
    
    void set_detected();
    bool is_detected() const;

private:
    /**
     * determine viral load course and infection course
     */
    void draw_infection_course(const InfectionState& start_state = InfectionState::Exposed);

    std::shared_ptr<Virus> m_virus;
    std::vector<std::pair<mio::abm::TimePoint, mio::abm::InfectionState>> m_infection_course; // start date of each infection state
    InfectionState* current_infection_state;
    ViralLoad m_viral_load;
    TimePoint m_start_date;
    TimePoint m_end_date;
    double m_log_norm_alpha, m_log_norm_beta; // have to ask for distribution/parametrization of the infectivity
    bool m_detected;
};

} // namespace abm
} // namespace mio

#endif
