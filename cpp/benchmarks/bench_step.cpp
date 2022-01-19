/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "bench_core.h"
#include "bench_model.h"

// google's benchmark makros require a function in the global scope, hence the following "using" statement
using mio::benchmark::integrator_step;

// dummy runs to avoid large effects of cpu scaling on times of actual benchmarks
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore)
    ->Name("Dummy 1/3");
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore)
    ->Name("Dummy 2/3");
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore)
    ->Name("Dummy 3/3");
// register functions as a benchmarks and set a name
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore)
    ->Name("simulate SecirModel adapt_rk");
BENCHMARK_TEMPLATE(integrator_step, mio::VRKIntegratorCore)
    ->Name("simulate SecirModel vadapt_rk");
BENCHMARK_TEMPLATE(integrator_step, mio::RKAllIntegratorCore)
    ->Name("simulate SecirModel adapt_rk_all");
BENCHMARK_TEMPLATE(integrator_step, mio::VRKAllIntegratorCore)
    ->Name("simulate SecirModel vadapt_rk_all");
BENCHMARK_TEMPLATE(integrator_step, mio::VRKOptIntegratorCore)
    ->Name("simulate SecirModel vadapt_rk_opt");
BENCHMARK_TEMPLATE(integrator_step, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>)
    ->Name("simulate SecirModel boost rk_ck54");
BENCHMARK_TEMPLATE(integrator_step, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_dopri5>)
    ->Name("simulate SecirModel boost rk_dopri5");
BENCHMARK_TEMPLATE(integrator_step, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_fehlberg78>)
    ->Name("simulate SecirModel boost rkf78");
// run all benchmarks
BENCHMARK_MAIN();