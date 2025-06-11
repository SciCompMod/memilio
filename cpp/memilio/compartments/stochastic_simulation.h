/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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
#ifndef MIO_COMPARTMENTS_STOCHASTIC_SIMULATION_H
#define MIO_COMPARTMENTS_STOCHASTIC_SIMULATION_H

#include "memilio/compartments/sde_model.h"
#include "memilio/compartments/simulation_base.h"
#include "memilio/math/euler_maruyama.h"
#include "memilio/utils/time_series.h"

namespace mio
{

/**
 * @brief A class for simulating a StochasticModel.
 * @tparam FP A floating point type, e.g. double.
 * @tparam M An implementation of a StochasticModel.
 */
template <typename FP, class M>
class StochasticSimulation : public SimulationBase<FP, M, 2>
{
    static_assert(is_stochsatic_model<FP, M>::value, "Template parameter must be a stochastic model.");

public:
    using Base  = SimulationBase<FP, M, 2>;
    using Model = M;

    /**
     * @brief Setup the simulation with an SDE solver.
     * @param[in] model An instance of a stochastic model.
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration.
     */
    StochasticSimulation(Model const& model, FP t0 = 0., FP dt = 0.1)
        : Base(model, std::make_shared<EulerMaruyamaIntegratorCore<FP>>(), t0, dt)
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
        return Base::advance({[this](auto&& y, auto&& t, auto&& dydt) {
                                  Base::get_model().eval_right_hand_side(y, y, t, dydt);
                              },
                              [this](auto&& y, auto&& t, auto&& dydt) {
                                  dydt.setZero();
                                  Base::get_model().get_noise(y, y, t, dydt);
                              }},
                             tmax, Base::get_result());
    }
};

template <typename FP, class Model, class Sim = StochasticSimulation<FP, Model>>
TimeSeries<FP> simulate_stochastic(FP t0, FP tmax, FP dt, Model const& model,
                                   std::shared_ptr<IntegratorCore<FP, 2>> integrator = nullptr)
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

#endif // MIO_COMPARTMENTS_STOCHASTIC_SIMULATION_H
