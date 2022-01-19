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

// google's benchmark templates require a function in the global scope, hence the following "using" statement
using mio::benchmark::simulation;
// set config path. the layout of the file depends on the initializer struct used to read it
const char configpath[] = "benchmarks/simulation_benchmark.config";
// specify a initilaizer which can read a config file for benchmark parameters.
// the 'filepath' template cannot be a string literal.
using init = mio::benchmark::simulation_file_init<configpath, mio::benchmark::model::SecirAgeresDampings>;

// dummy runs to avoid large effects of cpu scaling on times of actual benchmarks
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore, init)
    ->Name("Dummy 1/3");
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore, init)
    ->Name("Dummy 2/3");
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore, init)
    ->Name("Dummy 3/3");
// register functions as a benchmarks and set a name
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore, init)
    ->Name("simulate SecirModel adapt_rk");
BENCHMARK_TEMPLATE(simulation, mio::VRKIntegratorCore, init)
    ->Name("simulate SecirModel vadapt_rk");
BENCHMARK_TEMPLATE(simulation, mio::RKAllIntegratorCore, init)
    ->Name("simulate SecirModel adapt_rk_all");
BENCHMARK_TEMPLATE(simulation, mio::VRKAllIntegratorCore, init)
    ->Name("simulate SecirModel vadapt_rk_all");
BENCHMARK_TEMPLATE(simulation, mio::VRKOptIntegratorCore, init)
    ->Name("simulate SecirModel vadapt_rk_opt");
BENCHMARK_TEMPLATE(simulation, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>, init)
    ->Name("simulate SecirModel boost rk_ck54");
BENCHMARK_TEMPLATE(simulation, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_dopri5>, init)
    ->Name("simulate SecirModel boost rk_dopri5");
BENCHMARK_TEMPLATE(simulation, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_fehlberg78>, init)
    ->Name("simulate SecirModel boost rkf78");
// run all benchmarks
BENCHMARK_MAIN();