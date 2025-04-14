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
#include "benchmarks/integrator_step.h"
#include "benchmarks/secir_ageres_setups.h"

#include "memilio/math/adapt_rk.h"
#include "memilio/math/stepper_wrapper.h"

template <class Integrator>
void integrator_step(::benchmark::State& state)
{
    // suppress non-critical messages
    mio::set_log_level(mio::LogLevel::critical);
    // NOTE: make sure that yt has sensible values, e.g. by creating a simulation with the chosen model
    // with "num_agegroups" agegroups, and taking "yt" as the state of the simulation at "t_init"
    // NOTE: yt must have #agegroups * #compartments entries
    // benchmark setup
    const std::string config_path = mio::path_join(mio::memilio_dir(), "cpp/benchmarks/configs/integrator_step.config");
    auto cfg = mio::benchmark::IntegratorStepConfig::initialize(config_path);
    //auto cfg = mio::benchmark::IntegratorStepConfig::initialize();
    auto model = mio::benchmark::model::SecirAgeres(cfg.num_agegroups);
    // set deriv function and integrator
    mio::DerivFunction<ScalarType> f = [model](Eigen::Ref<const Eigen::VectorXd> x, double s,
                                               Eigen::Ref<Eigen::VectorXd> dxds) {
        model.eval_right_hand_side(x, x, s, dxds);
    };
    auto I = Integrator(cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);

    double t, dt;
    for (auto _ : state) {
        // This code gets timed
        t  = cfg.t_init;
        dt = cfg.dt_init;
        I.step(f, cfg.yt, t, dt, cfg.ytp1);
    }
}

// dummy runs to avoid large effects of cpu scaling on times of actual benchmarks
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore<ScalarType>)->Name("Dummy 1/3");
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore<ScalarType>)->Name("Dummy 2/3");
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore<ScalarType>)->Name("Dummy 3/3");
// register functions as a benchmarks and set a name
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore<ScalarType>)->Name("simulate SecirModel adapt_rk");
BENCHMARK_TEMPLATE(integrator_step,
                   mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>)
    ->Name("simulate SecirModel boost rk_ck54");
// BENCHMARK_TEMPLATE(integrator_step, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_dopri5>)
// ->Name("simulate SecirModel boost rk_dopri5"); // TODO: reenable once boost bug is fixed
BENCHMARK_TEMPLATE(integrator_step,
                   mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_fehlberg78>)
    ->Name("simulate SecirModel boost rkf78");
// run all benchmarks
BENCHMARK_MAIN();
