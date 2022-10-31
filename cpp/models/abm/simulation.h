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
#include "abm/output_results.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace abm
{

/**
 * run the simulation in discrete steps, evolve the world and report results.
 */
class Simulation
{
    using ResultVector = Eigen::Matrix<int, Eigen::Index(InfectionState::Count), 1>;

public:
    /**
     * Create a simulation.
     * @param t the starting time of the simulation
     * @param world the world to simulate
     */
    Simulation(TimePoint t, World&& world);

    /**
     * Create a simulation with an empty world.
     * World needs to be filled later.
     * @see Simulation::get_world
     * @param t the starting time of the simulation.
     */
    Simulation(TimePoint t)
        : Simulation(t, World())
    {
    }

    /** 
     * Run the simulation from the current time to tmax.
     * @param tmax time to stop
     */
    void advance(TimePoint tmax);

    /** 
     * Setup which data should get written to a txt file and which not.
     * @param print_data write aggregated data for each timestep
     * @param print_location_data write aggregated data for each location for each timestep
     */
    void set_print_data(bool print_data, bool print_location_data);

    /**
     * Get the Output Object that this simulation generates.
     */
    OutputResults& get_output()
    {
        return m_output;
    }

    /**
     * Get the World that this simulation evolves.
     * @{
     */
    World& get_world()
    {
        return m_world;
    }
    const World& get_world() const
    {
        return m_world;
    }
    /**@}*/

private:
    OutputResults m_output;
    World m_world;
    TimePoint m_t;
    TimeSpan m_dt;
};

} // namespace abm
} // namespace mio

#endif
