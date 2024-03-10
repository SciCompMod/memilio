/* 
* Copyright (C) 2020-2024 MEmilio
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
#ifndef EPI_ABM_SIMULATOR_H
#define EPI_ABM_SIMULATOR_H

#include "abm/world.h"
#include "abm/time.h"
#include "memilio/utils/time_series.h"
#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/io/history.h"

namespace mio
{
namespace abm
{

/**
 * @brief Run the Simulation in discrete steps, evolve the World and report results.
 */
class Simulation
{

public:
    /**
     * @brief Create a simulation.
     * @param[in] t0 The starting time of the Simulation.
     * @param[in] world The World to simulate.
     */
    Simulation(TimePoint t0, World&& world);

    /**
     * @brief Create a Simulation with an empty World.
     * World needs to be filled later.
     * @see Simulation::get_world
     * @param[in] t0 The starting time of the Simulation.
     */
    Simulation(TimePoint t0, size_t num_agegroups)
        : Simulation(t0, World(num_agegroups))
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
        double sum = 0;
        while (m_t < tmax) {
            evolve_world(tmax);
            double start, end;
            start = omp_get_wtime();
            (history.log(*this), ...);
            end = omp_get_wtime();
            sum = sum + (end - start);
        }
        std::cout << "Time spent on logging: " << sum << std::endl;
    }

    /**
     * @brief Get the current time of the Simulation.
     */
    TimePoint get_time() const
    {
        return m_t;
    }

    /**
     * @brief Get the previous time of the Simulation.
     */
    TimePoint get_prev_time() const
    {
        return m_prev_t;
    }

    /**
     * @brief Get the World that this Simulation evolves.
     */
    World& get_world()
    {
        return m_world;
    }
    const World& get_world() const
    {
        return m_world;
    }

private:
    void store_result_at(TimePoint t);
    void evolve_world(TimePoint tmax);

    World m_world; ///< The World to simulate.
    TimePoint m_t; ///< The current TimePoint of the Simulation.
    TimePoint m_prev_t; ///< The previous TimePoint of the Simulation.
    TimeSpan m_dt; ///< The length of the time steps.
};

} // namespace abm
} // namespace mio

#endif
