/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Jan Kleinert, Daniel Abele
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
#ifndef SIMULATION_H
#define SIMULATION_H

#include "memilio/config.h"
#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/utils/time_series.h"

namespace mio
{

template <typename FP = ScalarType>
using DefaultIntegratorCore = mio::ControlledStepperWrapper<FP, boost::numeric::odeint::runge_kutta_cash_karp54>;

/**
 * @brief A class for the simulation of a compartment model
 * @tparam FP floating point type, e.g., double
 * @tparam M a CompartmentModel type
 */
template <typename FP, class M>
class Simulation
{
    static_assert(is_compartment_model<FP, M>::value, "Template parameter must be a compartment model.");

public:
    using Model = M;

    /**
     * @brief Setup the simulation with an ODE solver.
     * @param[in] model An instance of a compartmental model
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration
     */
    Simulation(Model const& model, FP t0 = 0., FP dt = 0.1)
        : m_integratorCore(std::make_shared<DefaultIntegratorCore<FP>>())
        , m_model(std::make_unique<Model>(model))
        , m_integrator(m_integratorCore)
        , m_result(t0, m_model->get_initial_values())
        , m_dt(dt)
    {
    }

    Simulation(const Simulation& other)
        : m_integratorCore(other.m_integratorCore)
        , m_model(std::make_unique<Model>(*other.m_model))
        , m_integrator(m_integratorCore)
        , m_result(other.m_result)
        , m_dt(other.m_dt)
    {
    }

    /**
     * @brief Set the integrator core used in the simulation.
     * @param[in] integrator A shared pointer to an object derived from IntegratorCore.
     */
    void set_integrator(std::shared_ptr<IntegratorCore<FP>> integrator)
    {
        m_integratorCore = std::move(integrator);
        m_integrator.set_integrator(m_integratorCore);
    }

    /**
     * @brief Access the integrator core used in the simulation.
     * @return A reference to the integrator core used in the simulation
     * @{
     */
    IntegratorCore<FP>& get_integrator()
    {
        return *m_integratorCore;
    }

    /**
     * @brief get_integrator
     * @return reference to the core integrator used in the simulation
     */
    IntegratorCore<FP> const& get_integrator() const
    {
        return *m_integratorCore;
    }
    /** @} */

    /**
     * @brief advance simulation to tmax
     * tmax must be greater than get_result().get_last_time_point()
     * @param tmax next stopping point of simulation
     */
    Eigen::Ref<Vector<FP>> advance(FP tmax)
    {
        return m_integrator.advance(
            [this](auto&& y, auto&& t, auto&& dydt) {
                get_model().eval_right_hand_side(y, y, t, dydt);
            },
            tmax, m_dt, m_result);
    }

    /**
     * @brief Returns the simulation result describing the model population in each time step.
     *
     * Which compartments are used by the model is defined by the Comp template argument for the CompartmentalModel
     * (usually an enum named InfectionState).
     *
     * @return A TimeSeries to represent a numerical solution for the population of the model.
     * For each simulated time step, the TimeSeries contains the population size in each compartment.
     * @{
     */
    TimeSeries<FP>& get_result()
    {
        return m_result;
    }

    /**
     * @brief get_result returns the final simulation result
     * @return a TimeSeries to represent the final simulation result
     */
    const TimeSeries<FP>& get_result() const
    {
        return m_result;
    }
    /** @} */

    /**
     * @brief Get a reference to the model owned and used by the simulation.
     * @return The simulation model.
     * @{
     */
    const Model& get_model() const
    {
        return *m_model;
    }

    Model& get_model()
    {
        return *m_model;
    }
    /** @} */

    /**
     * @brief Returns the step size used by the integrator.
     * When using a integration scheme with adaptive time stepping, the integrator will store its estimate for the
     * next step size in this value.
     * @{
     */
    FP& get_dt()
    {
        return m_dt;
    }

    const FP& get_dt() const
    {
        return m_dt;
    }
    /** @} */

protected:
    /// @brief Get a reference to the integrater. Can be used to overwrite advance.
    OdeIntegrator<FP>& get_ode_integrator()
    {
        return m_integrator;
    }

private:
    std::shared_ptr<IntegratorCore<FP>> m_integratorCore; ///< Defines the integration scheme via its step function.
    std::unique_ptr<Model> m_model; ///< The model defining the ODE system and initial conditions.
    OdeIntegrator<FP> m_integrator; ///< Integrates the DerivFunction (see advance) and stores resutls in m_result.
    TimeSeries<FP> m_result; ///< The simulation results.
    FP m_dt; ///< The time step used (and possibly set) by m_integratorCore::step.
};

/**
 * Defines the return type of the `advance` member function of a type.
 * Template is invalid if this member function does not exist.
 *
 * @tparam FP floating point type, e.g., double
 * @tparam Sim a compartment model simulation type.
 */
template <typename FP, class Sim>
using advance_expr_t = decltype(std::declval<Sim>().advance(std::declval<FP>()));

/**
 * Template meta function to check if a type is a compartment model simulation. 
 * Defines a static constant of name `value`. 
 * The constant `value` will be equal to true if Sim is a valid compartment simulation type.
 * Otherwise, `value` will be equal to false.
 * @tparam FP floating point type, e.g., double
 * @tparam Sim a type that may or may not be a compartment model simulation.
 */
template <typename FP, class Sim>
using is_compartment_model_simulation =
    std::integral_constant<bool, (is_expression_valid<advance_expr_t, FP, Sim>::value &&
                                  is_compartment_model<FP, typename Sim::Model>::value)>;

/**
 * @brief Run a Simulation of a CompartmentalModel.
 * @param[in] t0 Start time.
 * @param[in] tmax End time.
 * @param[in] dt Initial step size of integration.
 * @param[in] model An instance of a CompartmentalModel.
 * @param[in] integrator Optionally override the IntegratorCore used by the Simulation.
 * @return A TimeSeries to represent the final Simulation result
 * @tparam FP floating point type, e.g., double
 * @tparam Model The particular Model derived from CompartmentModel to simulate.
 * @tparam Sim A Simulation that can simulate the model.
 */
template <typename FP, class Model, class Sim = Simulation<FP, Model>>
TimeSeries<FP> simulate(FP t0, FP tmax, FP dt, Model const& model,
                        std::shared_ptr<IntegratorCore<FP>> integrator = nullptr)
{
    model.check_constraints();
    Sim sim(model, t0, dt);
    if (integrator) {
        sim.set_integrator(integrator);
    }
    sim.advance(tmax);
    return sim.get_result();
}

} // namespace mio

#endif // SIMULATION_H
