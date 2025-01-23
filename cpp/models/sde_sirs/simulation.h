/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef MIO_SDE_SIRS_SIMULATION_H
#define MIO_SDE_SIRS_SIMULATION_H

#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "sde_sirs/model.h"

namespace mio
{
namespace ssirs
{

/// @brief A specialized Simulation for mio::ssirs::Model.
class Simulation : public mio::Simulation<ScalarType, Model>
{
protected:
    using mio::Simulation<ScalarType, Model>::set_integrator;

public:
    /**
     * @brief Set up the simulation with an ODE solver.
     * @param[in] model An instance of mio::ssirs::Model.
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration.
     */
    Simulation(Model const& model, ScalarType t0 = 0., ScalarType dt = 0.1)
        : mio::Simulation<ScalarType, Model>(model, t0, dt)
    {
        auto integrator = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();
        set_integrator(integrator);
    }

    using mio::Simulation<ScalarType, Model>::Simulation;
    /**
     * @brief advance simulation to tmax
     * tmax must be greater than get_result().get_last_time_point()
     * @param tmax next stopping point of simulation
     */
    Eigen::Ref<Eigen::VectorX<ScalarType>> advance(ScalarType tmax)
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
 * @brief Run a Simulation of a mio::ssirs::Model.
 * @param[in] t0 Start time.
 * @param[in] tmax End time.
 * @param[in] dt Initial step size of integration.
 * @param[in] model An instance of mio::ssirs::Model.
 * @return A TimeSeries to represent the final simulation result
 */
inline TimeSeries<ScalarType> simulate(ScalarType t0, ScalarType tmax, ScalarType dt, Model const& model)
{
    model.check_constraints();
    Simulation sim(model, t0, dt);
    sim.advance(tmax);
    return sim.get_result();
}

/// @brief A specialized FlowSimulation for mio::ssirs::Model.
class FlowSimulation : public mio::FlowSimulation<ScalarType, Model>
{
protected:
    using mio::FlowSimulation<ScalarType, Model>::set_integrator;

public:
    /**
     * @brief Set up the simulation with an ODE solver.
     * @param[in] model An instance of mio::ssirs::Model.
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration.
     */
    FlowSimulation(Model const& model, ScalarType t0 = 0., ScalarType dt = 0.1)
        : mio::FlowSimulation<ScalarType, Model>(model, t0, dt)
    {
        auto integrator = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();
        set_integrator(integrator);
    }

    /**
     * @brief Advance the simulation to tmax.
     * tmax must be greater than get_result().get_last_time_point().
     * @param[in] tmax Next stopping time of the simulation.
     */
    Eigen::Ref<Eigen::VectorX<ScalarType>> advance(ScalarType tmax)
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
 * @brief Run a FlowSimulation of mio::ssirs::Model.
 * @param[in] t0 Start time.
 * @param[in] tmax End time.
 * @param[in] dt Initial step size of integration.
 * @param[in] model An instance of mio::ssirs::Model.
 * @return The simulation result as two TimeSeries. The first describes the compartments at each time point,
 *         the second gives the corresponding flows that lead from t0 to each time point.
 */
inline std::vector<TimeSeries<ScalarType>> simulate_flows(ScalarType t0, ScalarType tmax, ScalarType dt,
                                                          Model const& model)
{
    model.check_constraints();
    FlowSimulation sim(model, t0, dt);
    sim.advance(tmax);
    return {sim.get_result(), sim.get_flows()};
}

} // namespace ssirs
} // namespace mio

#endif // MIO_SDE_SIRS_SIMULATION_H
