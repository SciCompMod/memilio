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
#include "abm/infection_state.h"
#include "abm/virus_variant.h"
#include "abm/parameters.h"
#include "memilio/utils/uncertain_value.h"

#include <vector>
#include <memory>

namespace mio
{
namespace abm
{

class ViralLoad
{
public:
    /**
     * @brief Constructor for ViralLoad.
     * @param[in] start_day TimePoint of construction.
     * @param[in] params Global infection parameters.
     */
    ViralLoad(TimePoint start_day, const GlobalInfectionParameters& params);

    /**
     * @brief Draws the viral load of the infection from a set of distributions.
     * @param[in] params Global infection parameters.
     */
    void draw_viral_load(const GlobalInfectionParameters& params);

    /**
     * @brief Gets the viral load of the infection at a given TimePoint.
     * @param[in] params TimePoint of querry.
     */
    ScalarType get_viral_load(TimePoint t) const;

private:
    TimePoint m_start_date;
    TimePoint m_end_date;
    UncertainValue m_peak;
    UncertainValue m_incline;
    UncertainValue m_decline; // always negative
};

class Infection
{

public:
    /**
     * @brief Create an infection for a single person.
     * @param virus virus type of the infection
     * @param start_date starting date of the infection
     * @param start_state starting infection state of the person, default susceptible. Only for initialization.
     */
    explicit Infection(VirusVariant virus, const GlobalInfectionParameters& params, TimePoint start_date,
                       InfectionState start_state = InfectionState::Exposed, bool detected = false);

    /**
     * @brief Get infectivity at a given time.
     * @param[in] t time point of the querry
     * @return infectivity at given time point.
     * Computed depending on current viral load and individual invlogit function of each person
     * corresponding to https://www.science.org/doi/full/10.1126/science.abi5273
     */
    ScalarType get_infectivity(TimePoint t) const;

    /**
     * @brief: Get virus type.
     * @return Virus type of the infection.
     */
    const VirusVariant& get_virus_variant() const;

    /**
     * @brief Get the infection state of the infection.
     * @param[in] t TimePoint of the querry.
     * @return InfectionState at the given TimePoint.
     */
    const InfectionState& get_infection_state(TimePoint t) const;

    /**
     * @brief Set the infection to detected.
    */
    void set_detected();

    /**
     * @returns Get the detected state.
    */
    bool is_detected() const;

private:
    /**
     * @brief Determine viral load course and infection course.
     * @param[in] start_date Start date of the Infection.
     * @param[in] params Global infection parameters.
     * @param[in] start_state [Default: InfectionState::Exposed] Start state of the Infection.
     */
    void draw_infection_course(TimePoint start_date, const GlobalInfectionParameters& params,
                               InfectionState start_state = InfectionState::Exposed);

    std::vector<std::pair<TimePoint, InfectionState>> m_infection_course; // start date of each infection state
    VirusVariant m_virus_variant;
    ViralLoad m_viral_load;
    ScalarType m_log_norm_alpha, m_log_norm_beta; // have to ask for distribution/parametrization of the infectivity
    bool m_detected;
};

} // namespace abm
} // namespace mio

#endif
