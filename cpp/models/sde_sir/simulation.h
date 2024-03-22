/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Nils Wassmuth, Rene Schmieding, Martin J. Kuehn
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
#ifndef MIO_SDE_SIR_SIMULATION_H
#define MIO_SDE_SIR_SIMULATION_H

#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "sde_sir/model.h"

namespace mio
{
namespace ssir
{

/// @brief A specialized Simulation for mio::ssir::Model.
class Simulation : public mio::Simulation<Model>
{
public:
    using mio::Simulation<Model>::Simulation;
    /**
     * @brief advance simulation to tmax
     * tmax must be greater than get_result().get_last_time_point()
     * @param tmax next stopping point of simulation
     */
    Eigen::Ref<Eigen::VectorXd> advance(double tmax)
    {
        return get_ode_integrator().advance(
            [this](auto&& y, auto&& t, auto&& dydt) {
                get_model().step_size = get_dt();
                get_model().get_derivatives(y, y, t, dydt);
            },
            tmax, get_dt(), get_result());
    }
};

/**
 * @brief Run a Simulation of a mio::ssir::Model.
 * @param[in] t0 Start time.
 * @param[in] tmax End time.
 * @param[in] dt Initial step size of integration.
 * @param[in] model An instance of mio::ssir::Model.
 * @param[in] integrator Optionally override the IntegratorCore used by the Simulation.
 * @return A TimeSeries to represent the final simulation result
 */
inline TimeSeries<ScalarType> simulate(double t0, double tmax, double dt, Model const& model,
                                       std::shared_ptr<IntegratorCore> integrator = nullptr)
{
    return mio::simulate<Model, Simulation>(t0, tmax, dt, model, integrator);
}

/// @brief A specialized FlowSimulation for mio::ssir::Model.
class FlowSimulation : public mio::FlowSimulation<Model>
{
public:
    using mio::FlowSimulation<Model>::FlowSimulation;
    /**
     * @brief Advance the simulation to tmax.
     * tmax must be greater than get_result().get_last_time_point().
     * @param[in] tmax Next stopping time of the simulation.
     */
    Eigen::Ref<Eigen::VectorXd> advance(double tmax)
    {
        assert(get_flows().get_num_time_points() == get_result().get_num_time_points());
        auto result = this->get_ode_integrator().advance(
            // see the general mio::FlowSimulation for more details on this DerivFunction
            [this](auto&& flows, auto&& t, auto&& dflows_dt) {
                const auto& pop_result = get_result();
                auto& model            = get_model();
                // compute current population
                model.get_derivatives(flows - get_flows().get_value(pop_result.get_num_time_points() - 1), m_pop);
                m_pop += pop_result.get_last_value();
                // compute the current change in flows with respect to the current population
                dflows_dt.setZero();
                model.step_size = get_dt(); // set the current step size
                model.get_flows(m_pop, m_pop, t, dflows_dt);
            },
            tmax, get_dt(), get_flows());
        compute_population_results();
        return result;
    }
};

/**
 * @brief Run a FlowSimulation of mio::ssir::Model.
 * @param[in] t0 Start time.
 * @param[in] tmax End time.
 * @param[in] dt Initial step size of integration.
 * @param[in] model An instance of mio::ssir::Model.
 * @param[in] integrator Optionally override the IntegratorCore used by the FlowSimulation.
 * @return The simulation result as two TimeSeries. The first describes the compartments at each time point,
 *         the second gives the corresponding flows that lead from t0 to each time point.
 */
inline std::vector<TimeSeries<ScalarType>> simulate_flows(double t0, double tmax, double dt, Model const& model,
                                                          std::shared_ptr<IntegratorCore> integrator = nullptr)
{
    return mio::simulate_flows<Model, FlowSimulation>(t0, tmax, dt, model, integrator);
}

} // namespace ssir
} // namespace mio

#endif // MIO_SDE_SIR_SIMULATION_H
