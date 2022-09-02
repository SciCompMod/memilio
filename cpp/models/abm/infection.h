/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: David Kerkmann
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
    ViralLoad(double peak, double increase, double decrease)
        : m_peak(peak)
        , m_increase_slope(increase)
        , m_decline_slope(decrease)
    {
    }

private:
    double m_peak;
    double m_increase_slope;
    double m_decline_slope;
}

class Infection
{

public:
    /**
     * create an infection for a single person.
     * @param start_date starting date of the infection
     * @param virus_type virus type of the infection
     */
    Infection(TimePoint start_date, Virus virus);

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

private:
    /**
     * determine viral load course and infection course
     */
    void draw_infection_course();
    //peak = gauss::get_instance();
    void determine_end_date();

    Virus m_virus;
    double m_infectivity_parameter; // have to ask for distribution/parametrization of the infectivity
    std::vector<std::pair<mio::abm::TimePoint, mio::abm::InfectionState>> m_infection_course;

    TimePoint m_start_date;
    TimePoint m_end_date;
    bool m_detected = false;
};

} // namespace abm
} // namespace mio

#endif
