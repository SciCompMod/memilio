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
     * @param[in] virus VirusVariant to determine the ViralLoad course.
     * @param[in] age AgeGroup to determine the ViralLoad course.
     * @param[in] start_day TimePoint of construction.
     * @param[in,out] params Global infection parameters.
     */
    explicit ViralLoad(VirusVariant virus, AgeGroup age, TimePoint start_day, GlobalInfectionParameters& params);

    /**
     * @brief Gets the viral load of the infection at a given TimePoint.
     * @param[in] t TimePoint of querry.
     */
    ScalarType get_viral_load(TimePoint t) const;

private:
    /**
     * @brief Draws the ViralLoad of the Infection from a set of distributions in the global parameters.
     * @param[in] virus VirusVariant to determine the ViralLoad course.
     * @param[in] age AgeGroup to determine the ViralLoad course.
     * @param[in] params Global infection parameters.
     */
    void draw_viral_load(VirusVariant virus, AgeGroup age, GlobalInfectionParameters& params);

    TimePoint m_start_date;
    TimePoint m_end_date;
    ScalarType m_peak;
    ScalarType m_incline;
    ScalarType m_decline; // always negative
};

class Infection
{

public:
    /**
     * @brief Create an Infection for a single Person.
     * @param[in] virus Virus type of the Infection.
     * @param[in] age AgeGroup to determine the ViralLoad course.
     * @param[in] start_date Starting date of the Infection.
     * @param[in] start_state [Default: InfectionState::Susceptible] Starting InfectionState of the Person. Only for initialization.
     * @param[in] detected [Default: false] If the Infection is detected.
     */
    explicit Infection(VirusVariant virus, AgeGroup age,
                       GlobalInfectionParameters& params, TimePoint start_date,
                       InfectionState start_state = InfectionState::Exposed, bool detected = false);

    /**
     * @brief Get infectivity at a given time.
     * @param[in] t TimePoint of the querry.
     * @return Infectivity at given time point.
     * Computed depending on current viral load and individual invlogit function of each person
     * corresponding to https://www.science.org/doi/full/10.1126/science.abi5273
     */
    ScalarType get_infectivity(TimePoint t) const;

    /**
     * @brief: Get VirusVariant.
     * @return VirusVariant of the Infection.
     */
    const VirusVariant& get_virus_variant() const;

    /**
     * @brief Get the InfectionState of the Infection.
     * @param[in] t TimePoint of the querry.
     * @return InfectionState at the given TimePoint.
     */
    const InfectionState& get_infection_state(TimePoint t) const;

    /**
     * @brief Set the Infection to detected.
    */
    void set_detected();

    /**
     * @returns Get the detected state.
    */
    bool is_detected() const;

private:
    /**
     * @brief Determine ViralLoad course and Infection course.
     * @param[in] start_date Start date of the Infection.
     * @param[in] params GlobalInfectionParameters.
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
