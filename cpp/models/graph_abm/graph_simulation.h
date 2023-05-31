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
#ifndef EPI_ABM_GRAPH_SIMULATOR_H
#define EPI_ABM_GRAPH_SIMULATOR_H

#include "models/graph_abm/graph_world.h"
#include "abm/simulation.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace graph_abm
{

/**
 * @brief Run the Simulation in discrete steps, evolve the World and report results.
 */
class GraphSimulation : public mio::abm::Simulation
{
    using Base = mio::abm::Simulation;

public:
    /**
     * @brief Create a graph simulation.
     * @param[in] t0 The starting time of the Simulation.
     * @param[in] graph_world The GraphWorld to simulate.
     */
    GraphSimulation(mio::abm::TimePoint t0, GraphWorld&& graph_world);

    /**
     * @brief Get the GraphWorld that this Simulation evolves.
     */
    GraphWorld& get_graph_world()
    {
        return m_graph_world;
    }
    const GraphWorld& get_graph_world() const
    {
        return m_graph_world;
    }

private:
    GraphWorld m_graph_world;
};

} // namespace graph_abm
} // namespace mio

#endif
