/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding, Daniel Abele, Martin J. Kuehn
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
#include "benchmarks/flow_simulation_ode_secirvvs.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "ode_secirvvs/model.h"
#include <string>

const std::string config_path = "../../benchmarks/simulation.config";

// simulation without flows (not in Model definition and not calculated by Simulation)
void flowless_sim(::benchmark::State& state)
{
    using Model = mio::benchmark::FlowlessModel;
    // suppress non-critical messages
    mio::set_log_level(mio::LogLevel::critical);
    // load config
    auto cfg = mio::benchmark::SimulationConfig::initialize(config_path);
    // create model
    Model model(cfg.num_agegroups);
    mio::benchmark::setup_model(model);
    // create simulation
    std::shared_ptr<mio::IntegratorCore<double>> I =
        std::make_shared<mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);
    // run benchmark
    for (auto _ : state) {
        // This code gets timed
        mio::benchmark::Simulation<mio::Simulation<double, Model>> sim(model, cfg.t0, cfg.dt);
        sim.set_integrator(I);
        // run sim
        sim.advance(cfg.t_max);
    }
}

// simulation with flows (in Model definition, but NOT calculated by Simulation)
void flow_sim_comp_only(::benchmark::State& state)
{
    using Model = mio::benchmark::FlowModel;
    // suppress non-critical messages
    mio::set_log_level(mio::LogLevel::critical);
    // load config
    auto cfg = mio::benchmark::SimulationConfig::initialize(config_path);
    // create model
    Model model(cfg.num_agegroups);
    mio::benchmark::setup_model(model);
    // create simulation
    std::shared_ptr<mio::IntegratorCore<double>> I =
        std::make_shared<mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);
    // run benchmark
    for (auto _ : state) {
        // This code gets timed
        mio::osecirvvs::Simulation<double, mio::Simulation<double, Model>> sim(model, cfg.t0, cfg.dt);
        sim.set_integrator(I);
        // run sim
        sim.advance(cfg.t_max);
    }
}

// simulation with flows (in Model definition and calculated by Simulation)
void flow_sim(::benchmark::State& state)
{
    using Model = mio::benchmark::FlowModel;
    // suppress non-critical messages
    mio::set_log_level(mio::LogLevel::critical);
    // load config
    auto cfg = mio::benchmark::SimulationConfig::initialize(config_path);
    // create model
    Model model(cfg.num_agegroups);
    mio::benchmark::setup_model(model);
    // create simulation
    std::shared_ptr<mio::IntegratorCore<double>> I =
        std::make_shared<mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);
    // run benchmark
    for (auto _ : state) {
        // This code gets timed
        mio::osecirvvs::Simulation<double, mio::FlowSimulation<double, Model>> sim(model, cfg.t0, cfg.dt);
        sim.set_integrator(I);
        // run sim
        sim.advance(cfg.t_max);
    }
}

// register functions as a benchmarks and set a name
// mitigate influence of cpu scaling
BENCHMARK(flowless_sim)->Name("Dummy 1/3");
BENCHMARK(flowless_sim)->Name("Dummy 2/3");
BENCHMARK(flowless_sim)->Name("Dummy 3/3");
// actual benchmarks
BENCHMARK(flowless_sim)
    ->Name(
        "osecirvvs::Simulation<mio::Simulation> on osecirvvs::Model (osecirvvs::* from pre 511 branch) without flows");
BENCHMARK(flow_sim_comp_only)->Name("osecirvvs::Simulation<mio::Simulation> on osecirvvs::Model with flows");
BENCHMARK(flow_sim)->Name("osecirvvs::Simulation<mio::FlowSimulation> on osecirvvs::Model with flows");
// run all benchmarks
BENCHMARK_MAIN();
