/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
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
#include "lct_secir/simulation.h"
#include "lct_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/stepper_wrapper.h"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"

#include <iostream>
#include <string>

namespace mio
{
namespace lsecir
{

Simulation::Simulation(Model const& model, ScalarType t0, ScalarType dt)
    : m_integratorCore(
          std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>())
    , m_model(std::make_unique<Model>(model))
    , m_integrator(
          [&model = *m_model](auto&& y, auto&& t, auto&& dydt) {
              model.eval_right_hand_side(y, t, dydt);
          },
          t0, m_model->get_initial_values(), dt, m_integratorCore)
{
}

void print_TimeSeries(const TimeSeries<ScalarType>& result, std::string heading)
{
    // print result after simulation
    std::cout << "# time | " + heading << std::endl;
    for (Eigen::Index i = 0; i < result.get_num_time_points(); ++i) {
        std::cout << result.get_time(i);
        for (Eigen::Index j = 0; j < result.get_num_elements(); ++j) {
            std::cout << " | " << std::fixed << std::setprecision(8) << result[i][j];
        }
        std::cout << "\n" << std::endl;
    }
}

TimeSeries<ScalarType> simulate(ScalarType t0, ScalarType tmax, ScalarType dt, Model const& model,
                                std::shared_ptr<IntegratorCore> integrator)
{
    model.check_constraints();
    Simulation sim(model, t0, dt);
    if (integrator) {
        sim.set_integrator(integrator);
    }
    sim.advance(tmax);
    return sim.get_result();
}

} // namespace lsecir
} // namespace mio