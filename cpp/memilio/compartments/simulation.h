/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "memilio/math/euler.h"

namespace mio
{

/**
 * @brief A class for the simulation of a compartment model.
 * @tparam M a CompartmentModel type
 */
template <class M>
class Simulation
{
    static_assert(is_compartment_model<M>::value, "Template parameter must be a compartment model.");

public:
    using Model = M;

    /**
     * @brief setup the simulation with an ODE solver
     * @param[in] model: An instance of a compartmental model
     * @param[in] t0 start time
     * @param[in] dt initial step size of integration
     */
    Simulation(Model const& model, double t0 = 0., double dt = 0.1)
        : m_integratorCore(
              std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>())
        , m_model(std::make_unique<Model>(model))
        , m_integrator(
              [&model = *m_model](auto&& y, auto&& t, auto&& dydt) {
                  model.eval_right_hand_side(y, y, t, dydt);
              },
              t0, m_model->get_initial_values(), dt, m_integratorCore)
    {
    }

    /**
     * @brief set the core integrator used in the simulation
     */
    void set_integrator(std::shared_ptr<IntegratorCore> integrator)
    {
        m_integratorCore = std::move(integrator);
        m_integrator.set_integrator(m_integratorCore);
    }

    /**
     * @brief get_integrator
     * @return reference to the core integrator used in the simulation
     */
    IntegratorCore& get_integrator()
    {
        return *m_integratorCore;
    }

    /**
     * @brief get_integrator
     * @return reference to the core integrator used in the simulation
     */
    IntegratorCore const& get_integrator() const
    {
        return *m_integratorCore;
    }

    /**
     * @brief advance simulation to tmax
     * tmax must be greater than get_result().get_last_time_point()
     * @param tmax next stopping point of simulation
     */
    Eigen::Ref<Eigen::VectorXd> advance(double tmax)
    {
        return m_integrator.advance(tmax);
    }

    /**
     * @brief Get the result of the simulation.
     * Return the number of persons in all InfectionState%s.
     * @return The result of the simulation.
     */
    TimeSeries<ScalarType>& get_result()
    {
        return m_integrator.get_result();
    }

    /**
     * @brief get_result returns the final simulation result
     * @return a TimeSeries to represent the final simulation result
     */
    const TimeSeries<ScalarType>& get_result() const
    {
        return m_integrator.get_result();
    }

    /**
     * @brief returns the simulation model used in simulation
     */
    const Model& get_model() const
    {
        return *m_model;
    }

    /**
     * @brief returns the simulation model used in simulation
     */
    Model& get_model()
    {
        return *m_model;
    }

    /**
     * @brief returns the next time step chosen by integrator
    */
    double get_dt() const
    {
        return m_integrator.get_dt();
    }

private:
    std::shared_ptr<IntegratorCore> m_integratorCore;
    std::unique_ptr<Model> m_model;
    OdeIntegrator m_integrator;
}; // namespace mio

/**
 * Defines the return type of the `advance` member function of a type.
 * Template is invalid if this member function does not exist.
 * @tparam Sim a compartment model simulation type.
 */
template <class Sim>
using advance_expr_t = decltype(std::declval<Sim>().advance(std::declval<double>()));

/**
 * Template meta function to check if a type is a compartment model simulation. 
 * Defines a static constant of name `value`. 
 * The constant `value` will be equal to true if Sim is a valid compartment simulation type.
 * Otherwise, `value` will be equal to false.
 * @tparam Sim a type that may or may not be a compartment model simulation.
 */
template <class Sim>
using is_compartment_model_simulation =
    std::integral_constant<bool, (is_expression_valid<advance_expr_t, Sim>::value &&
                                  is_compartment_model<typename Sim::Model>::value)>;

/**
 * @brief simulate simulates a compartmental model
 * @param[in] t0 start time
 * @param[in] tmax end time
 * @param[in] dt initial step size of integration
 * @param[in] model: An instance of a compartmental model
 * @return a TimeSeries to represent the final simulation result
 * @tparam Model a compartment model type
 * @tparam Sim a simulation type that can simulate the model.
 */
template <class Model, class Sim = Simulation<Model>>
TimeSeries<ScalarType> simulate(double t0, double tmax, double dt, Model const& model,
                                std::shared_ptr<IntegratorCore> integrator = nullptr)
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
