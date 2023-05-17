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

#include <vector>

namespace mio
{
namespace abm
{

/**
 * Models the ViralLoad for an Infection, modelled on a log_10 scale.
*/
struct ViralLoad {
    TimePoint start_date;
    TimePoint end_date;
    ScalarType peak;
    ScalarType incline;
    ScalarType decline; // always negative
};

class Infection
{

public:
    /**
     * @brief Create an Infection for a single Person.
     * @param[in] virus Virus type of the Infection.
     * @param[in] age AgeGroup to determine the ViralLoad course.
     * @param[in] params Parameters of the simulation.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state [Default: InfectionState::Exposed] InfectionState at time of initializing the the Infection.
     * @param[in] detected [Default: false] If the Infection is detected.
     */
    Infection(VirusVariant virus, AgeGroup age, const Parameters& params, TimePoint start_date,
              InfectionState start_state = InfectionState::Exposed, bool detected = false);

    /**
     * @brief Gets the viral load of the infection at a given TimePoint.
     * @param[in] t TimePoint of querry.
     */
    ScalarType get_viral_load(TimePoint t) const;

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
    VirusVariant get_virus_variant() const;

    /**
     * @brief Get the InfectionState of the Infection.
     * @param[in] t TimePoint of the querry.
     * @return InfectionState at the given TimePoint.
     */
    InfectionState get_infection_state(TimePoint t) const;

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
     * @brief Determine ViralLoad course and Infection course based on init_state. Calls draw_infection_course_backward for all InfectionState%s prior and draw_infection_course_forward for all subsequent InfectionState%s.
     * @param[in] age AgeGroup of the person.
     * @param[in] params Parameters of the simulation.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state InfectionState at time of initializing the the Infection.
     *
     */
    TimePoint draw_infection_course(AgeGroup age, const Parameters& params, TimePoint init_date,
                                    InfectionState init_state);

    /**
     * @brief Determine ViralLoad course and Infection course prior to the given start_state.
     * @param[in] age AgeGroup of the person.
     * @param[in] params Parameters of the simulation.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state InfectionState at time of initializing the the Infection.
     */
    void draw_infection_course_forward(AgeGroup age, const Parameters& params, TimePoint init_date,
                                       InfectionState init_state);

    /**
     * @brief Determine ViralLoad course and Infection course subsequent to the given start_state.
     * @param[in] age AgeGroup of the person.
     * @param[in] params Parameters of the simulation.
     * @param[in] init_date Date of initializing the Infection.
     * @param[in] init_state InfectionState at time of initializing the the Infection.
     */
    TimePoint draw_infection_course_backward(AgeGroup age, const Parameters& params, TimePoint init_date,
                                             InfectionState init_state);

    std::vector<std::pair<TimePoint, InfectionState>> m_infection_course; // start date of each infection state
    VirusVariant m_virus_variant;
    ViralLoad m_viral_load;
    ScalarType m_log_norm_alpha, m_log_norm_beta; // have to ask for distribution/parametrization of the infectivity
    bool m_detected;
};

} // namespace abm
} // namespace mio

#endif
