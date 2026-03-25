/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef MIO_MOBILITY_STAGE_ALIGNED_H
#define MIO_MOBILITY_STAGE_ALIGNED_H

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/math/stepper_wrapper.h"
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

namespace mio
{

/**
 * @brief Create a stage-aligned mobility simulation.
 *
 * This simulation is replacing the default adaptive integrator (Cash-Karp 5(4))
 * in each node with a fixed-step explicit RK4 whose step size matches the graph dt.
 * Combined with a model-specific calculate_mobility_returns overload that does
 * RK4 co-integration (e.g., from ode_seir/mobility.h). For more details, see:
 * 
 * Zunker et al. 2026, https://doi.org/10.48550/arXiv.2603.11275
 *
 * @param t0    Start time of the simulation.
 * @param dt    Time step between mobility events (also used as the RK4 step size).
 * @param graph The mobility graph with SimulationNodes and MobilityEdges.
 * @return A GraphSimulation configured for stage-aligned mobility updates.
 * @{
 */
template <typename FP, class Sim>
auto make_stage_aligned_mobility_sim(FP t0, FP dt, Graph<SimulationNode<FP, Sim>, MobilityEdge<FP>>& graph)
{
    // Switch all node integrators to fixed-step RK4
    for (auto& node : graph.nodes()) {
        node.property.get_simulation().set_integrator_core(
            std::make_unique<ExplicitStepperWrapper<FP, boost::numeric::odeint::runge_kutta4>>());
    }
    return make_mobility_sim<FP>(t0, dt, graph);
}

template <typename FP, class Sim>
auto make_stage_aligned_mobility_sim(FP t0, FP dt, Graph<SimulationNode<FP, Sim>, MobilityEdge<FP>>&& graph)
{
    // Switch all node integrators to fixed-step RK4
    for (auto& node : graph.nodes()) {
        node.property.get_simulation().set_integrator_core(
            std::make_unique<ExplicitStepperWrapper<FP, boost::numeric::odeint::runge_kutta4>>());
    }

    return make_mobility_sim<FP>(t0, dt, std::move(graph));
}
/** @} */

} // namespace mio

#endif // MIO_MOBILITY_STAGE_ALIGNED_H
