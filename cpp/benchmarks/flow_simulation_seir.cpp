/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "ode_secirvvs/model.h"
#include <string>

const std::string config_path = "../../benchmarks/simulation.config";

#include "memilio/compartments/simulation.h"
#include "models/ode_seir/model.h"

namespace mio
{
namespace benchmark
{

using FlowModel = ::mio::oseir::Model;

using namespace oseir;

// https://github.com/DLR-SC/memilio/blob/13555a6b23177d2d4633c393903461a27ce5762b/cpp/models/ode_seir/model.h
class FlowlessModel : public CompartmentalModel<InfectionState, Populations<InfectionState>, Parameters>
{
    using Base = CompartmentalModel<InfectionState, mio::Populations<InfectionState>, Parameters>;

public:
    FlowlessModel()
        : Base(Populations({InfectionState::Count}, 0.), ParameterSet())
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        auto& params     = this->parameters;
        double coeffStoE = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                           params.get<TransmissionProbabilityOnContact>() / populations.get_total();

        dydt[(size_t)InfectionState::Susceptible] =
            -coeffStoE * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected];
        dydt[(size_t)InfectionState::Exposed] =
            coeffStoE * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected] -
            (1.0 / params.get<TimeExposed>()) * y[(size_t)InfectionState::Exposed];
        dydt[(size_t)InfectionState::Infected] =
            (1.0 / params.get<TimeExposed>()) * y[(size_t)InfectionState::Exposed] -
            (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected];
        dydt[(size_t)InfectionState::Recovered] =
            (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected];
    }
};

template <class Model>
void setup_model(Model& model)
{
    double total_population                                                                            = 10000;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]   = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}]  = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}] = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}];
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    model.parameters.template set<mio::oseir::TimeExposed>(5.2);
    model.parameters.template set<mio::oseir::TimeInfected>(6);
    model.parameters.template set<mio::oseir::TransmissionProbabilityOnContact>(0.04);
    model.parameters.template get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 10;
}

} // namespace benchmark
} // namespace mio

// simulation without flows (not in Model definition and not calculated by Simulation)
void flowless_sim(::benchmark::State& state)
{
    using Model = mio::benchmark::FlowlessModel;
    // suppress non-critical messages
    mio::set_log_level(mio::LogLevel::critical);
    // load config
    auto cfg = mio::benchmark::SimulationConfig::initialize(config_path);
    // create model
    Model model;
    mio::benchmark::setup_model(model);
    // create simulation
    std::shared_ptr<mio::IntegratorCore> I =
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>(
            cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);
    mio::TimeSeries<ScalarType> results(static_cast<size_t>(Model::Compartments::Count));
    // run benchmark
    for (auto _ : state) {
        // This code gets timed
        results = mio::simulate(cfg.t0, cfg.t_max, cfg.dt, model, I);
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
    Model model;
    mio::benchmark::setup_model(model);
    // create simulation
    std::shared_ptr<mio::IntegratorCore> I =
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>(
            cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);
    mio::TimeSeries<ScalarType> results(static_cast<size_t>(Model::Compartments::Count));
    // run benchmark
    for (auto _ : state) {
        // This code gets timed
        results = mio::simulate(cfg.t0, cfg.t_max, cfg.dt, model, I);
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
    Model model;
    mio::benchmark::setup_model(model);
    // create simulation
    std::shared_ptr<mio::IntegratorCore> I =
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>(
            cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);
    mio::TimeSeries<ScalarType> results(static_cast<size_t>(Model::Compartments::Count));
    // run benchmark
    for (auto _ : state) {
        // This code gets timed
        results = mio::simulate_flows(cfg.t0, cfg.t_max, cfg.dt, model, I)[0];
    }
}

// simulation with flows (in Model definition and calculated by Simulation)
void flow_sim2(::benchmark::State& state)
{
    using Model = mio::benchmark::FlowModel;
    // suppress non-critical messages
    mio::set_log_level(mio::LogLevel::critical);
    // load config
    auto cfg = mio::benchmark::SimulationConfig::initialize(config_path);
    // create model
    Model model;
    mio::benchmark::setup_model(model);
    // create simulation
    std::shared_ptr<mio::IntegratorCore> I =
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>(
            cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);
    mio::TimeSeries<ScalarType> results(static_cast<size_t>(Model::Compartments::Count));
    // run benchmark
    for (auto _ : state) {
        // This code gets timed
        mio::SimulationFlows<Model> sim(model, cfg.t0, cfg.dt);
        sim.set_integrator(I);
        sim.advance(cfg.t_max);
        results = sim.get_result();
        // results = mio::simulate_flows(cfg.t0, cfg.t_max, cfg.dt, model, I)[0];
    }
}

// register functions as a benchmarks and set a name
// mitigate influence of cpu scaling
// BENCHMARK(flowless_sim)->Name("Dummy 1/3");
// BENCHMARK(flowless_sim)->Name("Dummy 2/3");
// BENCHMARK(flowless_sim)->Name("Dummy 3/3");
// actual benchmarks
BENCHMARK(flowless_sim)->Name("mio::Simulation on oseir::Model (pre 511 branch) without flows");
BENCHMARK(flow_sim_comp_only)->Name("mio::Simulation on oseir::Model with flows");
BENCHMARK(flow_sim)->Name("mio::FlowSimulation on oseir::Model with flows");
BENCHMARK(flow_sim2)->Name("mio::SimulationFlows on oseir::Model with flows");
// run all benchmarks
BENCHMARK_MAIN();
