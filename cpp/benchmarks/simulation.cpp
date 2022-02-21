/* 
* Copyright (cfg) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "benchmarks/simulation.h"
#include "benchmarks/secir_ageres_setups.h"

#include "memilio/math/adapt_rk.h"
#include "memilio/math/adapt_rk_fast.h"
#include "memilio/math/vadapt_rk.h"
#include "memilio/math/vadapt_rk_opt.h"
#include "memilio/math/vadapt_rk_fast.h"
#include "memilio/math/stepper_wrapper.h"

template <class Integrator>
void simulation(::benchmark::State& state)
{
    // suppress non-critical messages
    mio::set_log_level(mio::LogLevel::critical);
    // setup benchmark parameters
    auto cfg = mio::benchmark::SimulationConfig::initialize("benchmarks/simulation.config");
    //auto cfg = mio::benchmark::SimulationConfig::initialize(10);
    auto model = mio::benchmark::model::SecirAgeres(cfg.num_agegroups);

    for (auto _ : state) {
        // This code gets timed
        std::shared_ptr<mio::IntegratorCore> I =
            std::make_shared<Integrator>(cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);
        simulate(cfg.t0, cfg.t_max, cfg.dt, model, I);
    }
}

// dummy runs to avoid large effects of cpu scaling on times of actual benchmarks
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore)->Name("Dummy 1/3");
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore)->Name("Dummy 2/3");
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore)->Name("Dummy 3/3");
// register functions as a benchmarks and set a name
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore)->Name("simulate SecirModel adapt_rk");
BENCHMARK_TEMPLATE(simulation, mio::VRKIntegratorCore)->Name("simulate SecirModel vadapt_rk");
BENCHMARK_TEMPLATE(simulation, mio::RKFastIntegratorCore)->Name("simulate SecirModel adapt_rk_fast");
BENCHMARK_TEMPLATE(simulation, mio::VRKFastIntegratorCore)->Name("simulate SecirModel vadapt_rk_fast");
BENCHMARK_TEMPLATE(simulation, mio::VRKOptIntegratorCore)->Name("simulate SecirModel vadapt_rk_opt");
BENCHMARK_TEMPLATE(simulation, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>)
    ->Name("simulate SecirModel boost rk_ck54");
BENCHMARK_TEMPLATE(simulation, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_dopri5>)
    ->Name("simulate SecirModel boost rk_dopri5");
BENCHMARK_TEMPLATE(simulation, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_fehlberg78>)
    ->Name("simulate SecirModel boost rkf78");
// run all benchmarks
BENCHMARK_MAIN();