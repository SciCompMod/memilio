/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#ifndef EPI_ABM_SIMULATOR_H
#define EPI_ABM_SIMULATOR_H

#include "abm/world.h"
#include "abm/time.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/history.h" // IWYU pragma: keep

namespace mio
{
namespace abm
{

/**
 * @brief Run the Simulation in discrete steps, evolve the World and report results.
 */
template<typename FP=double>
class Simulation
{
    using ResultVector = Eigen::Matrix<int, Eigen::Index(InfectionState::Count), 1>;

public:
    /**
     * @brief Create a simulation.
     * @param[in] t0 The starting time of the Simulation.
     * @param[in] world The World to simulate.
     */
    Simulation(TimePoint t, World<FP>&& world)
        : m_world(std::move(world))
        , m_result(Eigen::Index(InfectionState::Count))
        , m_t(t)
        , m_dt(hours(1))
    {
        initialize_locations(m_t);
    }

    /**
     * @brief Create a Simulation with an empty World.
     * World needs to be filled later.
     * @see Simulation::get_world
     * @param[in] t0 The starting time of the Simulation.
     */
    Simulation(TimePoint t0)
        : Simulation(t0, World<FP>())
    {
    }

    /** 
     * @brief Run the Simulation from the current time to tmax.
     * @param[in] tmax Time to stop.
     */
    void advance(TimePoint tmax)
    {
        //log initial system state
        initialize_locations(m_t);
        store_result_at(m_t);
        while (m_t < tmax) {
            evolve_world(tmax);
            store_result_at(m_t);
        }
    }

    /** 
     * @brief Run the Simulation from the current time to tmax.
     * @param[in] tmax Time to stop.
     * @param[in] history History object to log data of the Simulation.
     */
    template <typename History>
    void advance(TimePoint tmax, History& history)
    {
        //log initial system state
        initialize_locations(m_t);
        store_result_at(m_t);
        history.log(*this);
        while (m_t < tmax) {
            evolve_world(tmax);
            store_result_at(m_t);
            history.log(*this);
        }
    }

    /**
     * @brief Get the result of the Simulation.
     * Sum over all Location%s of the number of Person%s in an #InfectionState.
     */
    const TimeSeries<ScalarType>& get_result() const
    {
        return m_result;
    }

    /**
     * @brief Get the current time of the Simulation.
     */
    TimePoint get_time() const
    {
        return m_t;
    }

    /**
     * @brief Get the World that this Simulation evolves.
     */
    World<FP>& get_world()
    {
        return m_world;
    }
    const World<FP>& get_world() const
    {
        return m_world;
    }

private:
    void initialize_locations(TimePoint t)
    {
        for (auto& location : m_world.get_locations()) {
            location.initialize_subpopulations(t);
        }
    }

    void store_result_at(TimePoint t)
    {
        m_result.add_time_point(t.days());
        m_result.get_last_value().setZero();
        for (auto& location : m_world.get_locations()) {
            m_result.get_last_value() += location.get_subpopulations().get_last_value().template cast<ScalarType>();
        }
    }

    void evolve_world(TimePoint tmax)
    {
        auto dt = std::min(m_dt, tmax - m_t);
        m_world.evolve(m_t, dt);
        m_t += m_dt;
    }

    World<FP> m_world; ///< The World to simulate.
    TimeSeries<ScalarType> m_result; ///< The result of the Simulation.
    TimePoint m_t; ///< The current TimePoint of the Simulation.
    TimeSpan m_dt; ///< The length of the time steps.
};

} // namespace abm
} // namespace mio

#endif
