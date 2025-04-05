/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Khoa Nguyen
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
#ifndef MIO_ABM_SIMULATION_H
#define MIO_ABM_SIMULATION_H

#include "abm/model.h"
#include "abm/time.h"
#include "memilio/io/history.h"

namespace mio
{
namespace abm
{

/**
 * @brief Run the Simulation in discrete steps, evolve the Model and report results.
 */
template <class M = Model>
class Simulation
{

public:
    /**
     * @brief Create a simulation.
     * @param[in] t0 The starting time of the Simulation.
     * @param[in] model The Model to simulate.
     */
    Simulation(TimePoint t0, M&& model)
        : m_model(std::move(model))
        , m_t(t0)
        , m_dt(hours(1))
    {
    }

    /**
     * @brief Create a Simulation with an empty Model.
     * Model needs to be filled later.
     * @see Simulation::get_model
     * @param[in] t0 The starting time of the Simulation.
     */
    Simulation(TimePoint t0, size_t num_agegroups)
        : Simulation(t0, M(num_agegroups))
    {
    }

    /** 
     * @brief Run the Simulation from the current time to tmax.
     * @param[in] tmax Time to stop.
     * @param[in] history History object to log data of the Simulation.
     */
    template <typename... History>
    void advance(TimePoint tmax, History&... history)
    {
        //log initial system state
        (history.log(*this), ...);
        while (m_t < tmax) {
            evolve_model(tmax);
            (history.log(*this), ...);
        }
    }

    /**
     * @brief Get the current time of the Simulation.
     */
    TimePoint get_time() const
    {
        return m_t;
    }

    /**
     * @brief Get the Model that this Simulation evolves.
     */
    M& get_model()
    {
        return m_model;
    }
    const M& get_model() const
    {
        return m_model;
    }

private:
    void store_result_at(TimePoint t);
    void evolve_model(TimePoint tmax)
    {
        auto dt = std::min(m_dt, tmax - m_t);
        m_model.evolve(m_t, dt);
        m_t += m_dt;
    }

    M m_model; ///< The Model to simulate.
    TimePoint m_t; ///< The current TimePoint of the Simulation.
    TimeSpan m_dt; ///< The length of the time steps.
};

} // namespace abm
} // namespace mio

#endif
