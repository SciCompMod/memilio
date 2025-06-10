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
#ifndef MIO_SDE_SEIRVV_SIMULATION_H
#define MIO_SDE_SEIRVV_SIMULATION_H

#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "sde_seirvv/model.h"

namespace mio
{
namespace sseirvv
{

/// @brief A specialized Simulation for mio::sseirvv::Model.
template <typename FP>
class Simulation : public mio::Simulation<FP, Model<FP>>
{
protected:
    using mio::Simulation<FP, Model<FP>>::set_integrator;

public:
    /**
     * @brief Set up the simulation with an SDE solver.
     * @param[in] model An instance of mio::sseirvv::Model.
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration.
     */
    Simulation(Model<FP> const& model, FP t0 = 0., FP dt = 0.1)
        : mio::Simulation<FP, Model<FP>>(model, t0, dt)
    {
        auto integrator = std::make_shared<mio::EulerIntegratorCore<FP>>();
        set_integrator(integrator);
    }

    using mio::Simulation<FP, Model<FP>>::Simulation;
    /**
     * @brief Advance simulation to tmax.
     * tmax must be greater than get_result().get_last_time_point().
     * @param tmax Next stopping point of simulation.
     */
    Eigen::Ref<Eigen::VectorX<FP>> advance(FP tmax)
    {
        return this->get_ode_integrator().advance(
            [this](auto&& y, auto&& t, auto&& dydt) {
                this->get_model().step_size = this->get_dt();
                this->get_model().get_derivatives(y, y, t, dydt);
            },
            tmax, this->get_dt(), this->get_result());
    }
};

/**
 * @brief Run a Simulation of a mio::sseirvv::Model.
 * @param[in] t0 Start time.
 * @param[in] tmax End time.
 * @param[in] dt Initial step size of integration.
 * @param[in] model An instance of mio::sseirvv::Model.
 * @return A TimeSeries to represent the final simulation result
 */
template <typename FP>
inline TimeSeries<FP> simulate(FP t0, FP tmax, FP dt, Model<FP> const& model)
{
    model.check_constraints();
    Simulation sim(model, t0, dt);
    sim.advance(tmax);
    return sim.get_result();
}

/// @brief A specialized FlowSimulation for mio::sseirvv::Model.
template <typename FP>
class FlowSimulation : public mio::FlowSimulation<FP, Model<FP>>
{
protected:
    using mio::FlowSimulation<FP, Model<FP>>::set_integrator;

public:
    /**
     * @brief Set up the simulation with an ODE solver.
     * @param[in] model An instance of mio::sseirvv::Model.
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration.
     */
    FlowSimulation(Model<FP> const& model, FP t0 = 0., FP dt = 0.1)
        : mio::FlowSimulation<FP, Model<FP>>(model, t0, dt)
    {
        auto integrator = std::make_shared<mio::EulerIntegratorCore<FP>>();
        this->set_integrator(integrator);
    }

    /**
     * @brief Advance the simulation to tmax.
     * tmax must be greater than get_result().get_last_time_point().
     * @param[in] tmax Next stopping time of the simulation.
     */
    Eigen::Ref<Eigen::VectorX<FP>> advance(FP tmax)
    {
        assert(this->get_flows().get_num_time_points() == this->get_result().get_num_time_points());
        auto result = this->get_ode_integrator().advance(
            // See the general mio::FlowSimulation for more details on this DerivFunction.
            [this](auto&& flows, auto&& t, auto&& dflows_dt) {
                const auto& pop_result = this->get_result();
                auto& model            = this->get_model();
                // Compute current population.
                model.get_derivatives(flows - this->get_flows().get_value(pop_result.get_num_time_points() - 1),
                                      this->m_pop);
                this->m_pop += pop_result.get_last_value();
                // Compute the current change in flows with respect to the current population.
                dflows_dt.setZero();
                model.step_size = this->get_dt(); // set the current step size
                model.get_flows(this->m_pop, this->m_pop, t, dflows_dt);
            },
            tmax, this->get_dt(), this->get_flows());
        this->compute_population_results();
        return result;
    }
};

/**
 * @brief Run a FlowSimulation of mio::sseirvv::Model.
 * @param[in] t0 Start time.
 * @param[in] tmax End time.
 * @param[in] dt Initial step size of integration.
 * @param[in] model An instance of mio::sseirvv::Model.
 * @return The simulation result as two TimeSeries. The first describes the compartments at each time point,
 *         the second gives the corresponding flows that lead from t0 to each time point.
 */
template <typename FP>
inline std::vector<TimeSeries<FP>> simulate_flows(FP t0, FP tmax, FP dt, Model<FP> const& model)
{
    model.check_constraints();
    FlowSimulation sim(model, t0, dt);
    sim.advance(tmax);
    return {sim.get_result(), sim.get_flows()};
}

} // namespace sseirvv
} // namespace mio

#endif // MIO_SDE_SEIRVV_SIMULATION_H
