/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef MIO_COMPARTMENTS_SIMULATION_H
#define MIO_COMPARTMENTS_SIMULATION_H

#include "memilio/compartments/compartmental_model.h"
#include "memilio/compartments/simulation_base.h"
#include "memilio/math/integrator.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/utils/time_series.h"

namespace mio
{

/// @brief The default integrator used by Simulation.
template <typename FP>
using DefaultIntegratorCore = mio::ControlledStepperWrapper<FP, boost::numeric::odeint::runge_kutta_cash_karp54>;

/**
 * @brief A class for simulating a CompartmentalModel.
 * @tparam FP A floating point type, e.g. double.
 * @tparam M An implementation of a CompartmentalModel.
 */
template <typename FP, class M>
class Simulation : public details::SimulationBase<FP, M, OdeIntegrator<FP>>
{
public:
    using Base  = details::SimulationBase<FP, M, OdeIntegrator<FP>>;
    using Model = M;

    /**
     * @brief Setup the simulation with an ODE solver.
     * @param[in] model An instance of a compartmental model
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration
     */
    Simulation(Model const& model, FP t0 = 0., FP dt = 0.1)
        : Base(model, std::make_shared<DefaultIntegratorCore<FP>>(), t0, dt)
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
        return Base::advance(
            [this](auto&& y, auto&& t, auto&& dydt) {
                Base::get_model().eval_right_hand_side(y, y, t, dydt);
            },
            tmax, Base::get_result());
    }
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
                        std::shared_ptr<OdeIntegratorCore<FP>> integrator = nullptr)
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

#endif // MIO_COMPARTMENTS_SIMULATION_H
