/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding, Henrik Zunker
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
#ifndef MIO_COMPARTMENTS_FLOW_SIMULATION_BASE_H
#define MIO_COMPARTMENTS_FLOW_SIMULATION_BASE_H

#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/simulation_base.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace details
{

/**
 * @brief Base class to define a FlowSimulation.
 * Provides a protected advance method that accepts the specified integrands, that must be exposed by the Derived
 * class. Also defines all relevant members and accessors for a FlowSimulation.
 * @tparam FP A floating point type, e.g. double.
 * @tparam M An implementation of FlowModel.
 * @tparam Integrands One or more function types used for defining the right hand side of a system of equations.
 */
template <typename FP, class M, class... Integrands>
class FlowSimulationBase : public SimulationBase<FP, M, Integrands...>
{
    static_assert(is_flow_model<FP, M>::value, "Template parameter must be a flow model.");

public:
    using Base  = SimulationBase<FP, M, Integrands...>;
    using Model = M;
    using Core  = IntegratorCore<FP, Integrands...>;

    /**
     * @brief Create a FlowSimulationBase.
     * @param[in] model An instance of a flow model.
     * @param[in] integrator_core A unique pointer to an object derived from IntegratorCore.
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration.
     */
    FlowSimulationBase(Model const& model, std::unique_ptr<Core>&& integrator_core, FP t0, FP dt)
        : Base(model, std::move(integrator_core), t0, dt)
        , m_flow_result(t0, model.get_initial_flows())
    {
    }

    /**
     * @brief Returns the simulation result describing the transitions between compartments for each time step.
     *
     * Which flows are used by the model is defined by the Flows template argument for the FlowModel.
     * To get the correct index for the flow between two compartments use FlowModel::get_flat_flow_index.
     *
     * @return A TimeSeries to represent a numerical solution for the flows in the model. 
     * For each simulated time step, the TimeSeries contains the value of each flow. 
     * @{
     */
    TimeSeries<FP>& get_flows()
    {
        return m_flow_result;
    }

    const TimeSeries<FP>& get_flows() const
    {
        return m_flow_result;
    }
    /** @} */

protected:
    /**
     * @brief Computes the distribution of the Population to the InfectionState%s based on the simulated flows.
     * Uses the same method as the DerivFunction used in advance to compute the population given the flows and initial
     * values. Adds time points to Base::m_result until it has the same number of time points as m_flow_result.
     * Does not recalculate older values.
     */
    void compute_population_results()
    {
        const auto& flows = get_flows();
        const auto& model = this->get_model();
        auto& result      = this->get_result();
        // take the last time point as base result (instead of the initial results), so that we use external changes
        const size_t last_tp = result.get_num_time_points() - 1;
        // calculate new time points
        for (Eigen::Index i = result.get_num_time_points(); i < flows.get_num_time_points(); i++) {
            result.add_time_point(flows.get_time(i));
            model.get_derivatives(flows.get_value(i) - flows.get_value(last_tp), result.get_value(i));
            result.get_value(i) += result.get_value(last_tp);
        }
    }

private:
    mio::TimeSeries<FP> m_flow_result; ///< Flow result of the simulation.
};

/// @brief Specialization of FlowSimulationBase that takes a SystemIntegrator instead of it's Integrands.
template <typename FP, class M, class... Integrands>
class FlowSimulationBase<FP, M, SystemIntegrator<FP, Integrands...>> : public FlowSimulationBase<FP, M, Integrands...>
{
    using FlowSimulationBase<FP, M, Integrands...>::FlowSimulationBase;
};

} // namespace details
} // namespace mio

#endif // MIO_COMPARTMENTS_FLOW_SIMULATION_BASE_H
