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
#ifndef MIO_COMPARTMENTS_FLOW_SIMULATION_H
#define MIO_COMPARTMENTS_FLOW_SIMULATION_H

#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/flow_simulation_base.h"
#include "memilio/compartments/simulation.h"
#include "memilio/math/integrator.h"

namespace mio
{

/**
 * @brief A class for simulating a FlowModel.
 * @tparam FP A floating point type, e.g. double.
 * @tparam M An implementation of a FlowModel.
 */
template <typename FP, IsFlowModel<FP> M>
class FlowSimulation : public details::FlowSimulationBase<FP, M, OdeIntegrator<FP>>
{
public:
    using Base  = details::FlowSimulationBase<FP, M, OdeIntegrator<FP>>;
    using Model = M;

    /**
     * @brief Set up the simulation with an ODE solver.
     * @param[in] model An instance of a flow model.
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration.
     */
    FlowSimulation(Model const& model, FP t0 = 0., FP dt = 0.1)
        : Base(model, std::make_unique<DefaultIntegratorCore<FP>>(), t0, dt)
        , m_pop(model.get_initial_values().size())
    {
    }

    /**
     * @brief Run the simulation up to a given time.
     * The time tmax must be greater than `get_result().get_last_time_point()`, which is used as the starting point. The
     * initial value is `get_result().get_last_value()`.
     * @param[in] tmax Next stopping point of the simulation.
     * @return The simulation result at tmax.
     */
    Eigen::Ref<Eigen::VectorX<FP>> advance(FP tmax)
    {
        // the derivfunktion (i.e. the lambda passed to m_integrator.advance below) requires that there are at least
        // as many entries in m_flow_result as in Base::m_result
        assert(Base::get_flows().get_num_time_points() == Base::get_result().get_num_time_points());
        const auto result = Base::advance(
            [this](auto&& flows, auto&& t, auto&& dflows_dt) {
                const auto& pop_result = this->get_result();
                const auto& model      = this->get_model();
                // compute current population
                //   flows contains the accumulated outflows of each compartment for each target compartment at time t.
                //   Using that the ODEs are linear expressions of the flows, get_derivatives can compute the total change
                //   in population from t0 to t.
                //   To incorporate external changes to the last values of pop_result (e.g. by applying mobility), we only
                //   calculate the change in population starting from the last available time point in m_result, instead
                //   of starting at t0. To do that, the following difference of flows is used.
                model.get_derivatives(flows - Base::get_flows().get_value(pop_result.get_num_time_points() - 1),
                                      m_pop); // note: overwrites values in pop
                //   add the "initial" value of the ODEs (using last available time point in pop_result)
                //     If no changes were made to the last value in m_result outside of FlowSimulation, the following
                //     line computes the same as `model.get_derivatives(flows, x); x += model.get_initial_values();`.
                m_pop += pop_result.get_last_value();
                // compute the current change in flows with respect to the current population
                dflows_dt.setZero();
                model.get_flows(m_pop, m_pop, t, dflows_dt); // this result is used by the integrator
            },
            tmax, Base::get_flows());
        Base::compute_population_results();
        return result;
    }

private:
    Eigen::VectorX<FP> m_pop; ///< pre-allocated temporary, used during computation of flow derivatives
};

/**
 * @brief Run a FlowSimulation of a FlowModel.
 * @param[in] t0 Start time.
 * @param[in] tmax End time.
 * @param[in] dt Initial step size of integration.
 * @param[in] model An instance of a FlowModel.
 * @param[in] integrator_core Optionally override the IntegratorCore used by the FlowSimulation.
 * @return The simulation result as two TimeSeries. The first describes the compartments at each time point,
 *         the second gives the corresponding flows that lead from t0 to each time point.
 * @tparam FP a floating point type, e.g., double
 * @tparam Model The particular Model derived from FlowModel to simulate.
 * @tparam Sim A FlowSimulation that can simulate the model.
 */
template <typename FP, class Model, class Sim = FlowSimulation<FP, Model>>
std::vector<TimeSeries<FP>> simulate_flows(FP t0, FP tmax, FP dt, Model const& model,
                                           std::unique_ptr<OdeIntegratorCore<FP>>&& integrator_core = nullptr)
{
    model.check_constraints();
    Sim sim(model, t0, dt);
    if (integrator_core) {
        sim.set_integrator_core(std::move(integrator_core));
    }
    sim.advance(tmax);
    return {sim.get_result(), sim.get_flows()};
}

} // namespace mio

#endif // MIO_COMPARTMENTS_FLOW_SIMULATION_H
