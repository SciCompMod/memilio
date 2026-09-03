/*
* Copyright (C) 2020-2026 MEmilio
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
#ifndef MIO_ABM_PAIS_H
#define MIO_ABM_PAIS_H

#include "abm/pais_state.h"
#include "memilio/io/default_serialize.h"
#include "abm/time.h"
#include "abm/personal_rng.h"
#include "abm/infection.h"
#include "abm/parameters.h"

#include <vector>

namespace mio
{
namespace abm
{

class Person; // forward declaration to avoid circular dependency between Person and PAIS

/**
 * @brief Represents a PAIS (Post-Acute Infection Syndrome).
 */
struct PAIS {
    std::vector<std::pair<TimePoint, PAISState>> severity{}; ///< Time series of the severity of the PAIS.

    /**
     * @brief Update the severity of the PAIS.
     * @param[in] params The Parameters of the Simulation.
     * @param[in] rng Personal RandomNumberGenerator.
     * @param[in] t TimePoint of querry. Usually the current time of the Simulation.
     * @param[in] dt The time step size of the Simulation.
     */
    void update_severity(const Parameters& params, PersonalRandomNumberGenerator& rng, TimePoint t, TimeSpan dt);

    /**
     * @brief Initialize or refresh the PAIS status based on a new infection.
     * @param[in] params The Parameters of the Simulation.
     * @param[in] p The Person with the new Infection.
     * @param[in] rng Personal RandomNumberGenerator.
     * @param[in] inf The new Infection.
     * @param[in] t TimePoint of querry. Usually the current time of the Simulation.
     */
    void init_or_refresh(const Parameters& params, Person& p, PersonalRandomNumberGenerator& rng, const Infection& inf,
                         TimePoint t);

    /**
     * @brief Add a new severity state to the PAIS based on the highest InfectionState of the infection.
      * If the highest InfectionState is InfectedSevere or InfectedCritical, the new severity state is Severe, otherwise it is Medium.
      * The new severity state is only added if it is different from the current severity state and if the last update was before t.
      * @param[in] t TimePoint of querry. Usually the current time of the Simulation.
      * @param[in] highest_state The highest InfectionState of the infection.
     */
    void add_new_severity(TimePoint t, std::pair<TimePoint, InfectionState> highest_state);

    /**
     * @brief Get the severity of the PAIS at a given time.
     * @param[in] t TimePoint of querry. Usually the current time of the Simulation.
     * @return The severity of the PAIS at time t.
     */
    PAISState get_severity(TimePoint t) const;

    auto default_serialize()
    {
        return Members("PAIS").add("severity", severity);
    }
};

} // namespace abm
} // namespace mio

#endif
