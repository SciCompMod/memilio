/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: René Schmieding, Julia Bicker
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

#ifndef MIO_D_ABM_MODEL_H
#define MIO_D_ABM_MODEL_H

namespace mio
{

namespace dabm
{

/**
 * @brief Wrap an implementation of a diffusive ABM so it can be run by the d_abm::Simulation.
 * Uses the CRTP. See comments on using statements for expected function signatures.
 * @tparam Implementation A class implementing all functions and types marked with the using keyword in Model.
 */
template <class Implementation>
class Model : public Implementation
{
public:
    /// Use the constructors defined by the Implementation
    using Implementation::Implementation;

    /**
     * @brief Set the status of an agent.
     * Expected signature: `void adopt(Agent&, const Status&)`
     */
    using Implementation::adopt;

    /**
     * @brief Calculate the current adoption rate of an agent from its status to the given one.
     * Expected signature: `ScalarType adoption_rate(const Agent&, const Status&)`
     */
    using Implementation::adoption_rate;

    /**
     * @brief Change the Position of an Agent, depending on its state, the current time and step size.
     * Expected signature: `void move(const ScalarType, const ScalarType, Agent&)`
     * The first argument is time, the second step size.
     */
    using Implementation::move;

    /**
     * @brief Get the Implementations RNG.
     * Expected signature: `mio::RandomNumberGenerator& get_rng()`
     */
    using Implementation::get_rng;

    /**
     * @brief Aggregate the population by their Status for the simulation result.
     * Expected signature: `Eigen::VectorX<ScalarType> time_point()`
     */
    using Implementation::time_point;

    /// @brief The status of an agent.
    using Status = typename Implementation::Status;
    /// @brief An agent is expected to contain at least a status and a position.
    using Agent = typename Implementation::Agent;

    /// @brief Empty function for compatability with MEmilio.
    inline constexpr void check_constraints() const
    {
    }
};

} // namespace dabm
} // namespace mio

#endif
