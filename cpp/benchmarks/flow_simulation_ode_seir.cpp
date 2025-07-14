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
#include "memilio/compartments/flow_simulation.h"
#include "ode_seir/model.h"
#include <string>

const std::string config_path = "../../benchmarks/simulation.config";

#include "memilio/compartments/simulation.h"
#include "models/ode_seir/model.h"

namespace mio
{
namespace benchmark
{

using FlowModel = ::mio::oseir::Model<double>;

using namespace oseir;

// For comparison benchmarks, the ODE-SEIR model has been adapted into
// a compartmental model that does not rely on flows.
class FlowlessModel
    : public CompartmentalModel<double, oseir::InfectionState, Populations<double, AgeGroup, oseir::InfectionState>,
                                oseir::Parameters<double>>
{
    using Base =
        CompartmentalModel<double, oseir::InfectionState, mio::Populations<double, AgeGroup, oseir::InfectionState>,
                           oseir::Parameters<double>>;

public:
    FlowlessModel(int num_agegroups)
        : Base(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        auto& params                     = this->parameters;
        const Index<AgeGroup> age_groups = reduce_index<Index<AgeGroup>>(this->populations.size());

        for (auto i : mio::make_index_range(age_groups)) {
            size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            size_t Ei = this->populations.get_flat_index({i, InfectionState::Exposed});
            size_t Ii = this->populations.get_flat_index({i, InfectionState::Infected});
            size_t Ri = this->populations.get_flat_index({i, InfectionState::Recovered});

            for (auto j : make_index_range(age_groups)) {

                size_t Sj = this->populations.get_flat_index({j, InfectionState::Susceptible});
                size_t Ej = this->populations.get_flat_index({j, InfectionState::Exposed});
                size_t Ij = this->populations.get_flat_index({j, InfectionState::Infected});
                size_t Rj = this->populations.get_flat_index({j, InfectionState::Recovered});

                const double Nj_inv = 1.0 / (pop[Sj] + pop[Ej] + pop[Ij] + pop[Rj]);
                const double coeffStoE =
                    params.template get<ContactPatterns<double>>().get_cont_freq_mat().get_matrix_at(
                        SimulationTime<double>(t))(i.get(), j.get()) *
                    params.template get<TransmissionProbabilityOnContact<double>>()[i] * Nj_inv;

                dydt[Si] -= y[Si] * pop[Ij] * coeffStoE;
                dydt[Ei] += y[Si] * pop[Ij] * coeffStoE;
            }

            dydt[Ii] += (1.0 / params.get<TimeExposed<double>>()[i]) * y[Ei];
            dydt[Ii] -= (1.0 / params.get<TimeInfected<double>>()[i]) * y[Ii];
            dydt[Ri] = (1.0 / params.get<TimeInfected<double>>()[i]) * y[Ii];
        }
    }
};

template <class Model>
void setup_model(Model& model)
{
    const double total_population = 10000.0;
    const auto num_groups         = model.parameters.get_num_groups();
    for (AgeGroup i = 0; i < num_groups; i++) {
        model.populations[{i, oseir::InfectionState::Exposed}]   = 100.0 / static_cast<size_t>(num_groups);
        model.populations[{i, oseir::InfectionState::Infected}]  = 100.0 / static_cast<size_t>(num_groups);
        model.populations[{i, oseir::InfectionState::Recovered}] = 100.0 / static_cast<size_t>(num_groups);
        model.populations[{i, oseir::InfectionState::Susceptible}] =
            total_population / static_cast<size_t>(num_groups) -
            model.populations[{i, oseir::InfectionState::Exposed}] -
            model.populations[{i, oseir::InfectionState::Infected}] -
            model.populations[{i, oseir::InfectionState::Recovered}];
    }
    model.parameters.template set<mio::oseir::TimeExposed<double>>(5.2);
    model.parameters.template set<mio::oseir::TimeInfected<double>>(6);
    model.parameters.template set<mio::oseir::TransmissionProbabilityOnContact<double>>(0.04);
    mio::ContactMatrixGroup<double>& contact_matrix =
        model.parameters.template get<mio::oseir::ContactPatterns<double>>();
    contact_matrix[0].get_baseline().setConstant(10.0);
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
    Model model(cfg.num_agegroups);
    mio::benchmark::setup_model(model);
    // create simulation
    std::shared_ptr<mio::IntegratorCore<double>> I =
        std::make_shared<mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);
    mio::TimeSeries<double> results(static_cast<size_t>(Model::Compartments::Count));
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
    Model model(cfg.num_agegroups);
    mio::benchmark::setup_model(model);
    // create simulation
    std::shared_ptr<mio::IntegratorCore<double>> I =
        std::make_shared<mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);
    mio::TimeSeries<double> results(static_cast<size_t>(Model::Compartments::Count));
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
    Model model(cfg.num_agegroups);
    mio::benchmark::setup_model(model);
    // create simulation
    std::shared_ptr<mio::IntegratorCore<double>> I =
        std::make_shared<mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            cfg.abs_tol, cfg.rel_tol, cfg.dt_min, cfg.dt_max);
    mio::TimeSeries<double> results(static_cast<size_t>(Model::Compartments::Count));
    // run benchmark
    for (auto _ : state) {
        // This code gets timed
        results = mio::simulate_flows(cfg.t0, cfg.t_max, cfg.dt, model, I)[0];
    }
}

// register functions as a benchmarks and set a name
// mitigate influence of cpu scaling
BENCHMARK(flowless_sim)->Name("Dummy 1/3");
BENCHMARK(flowless_sim)->Name("Dummy 2/3");
BENCHMARK(flowless_sim)->Name("Dummy 3/3");
// actual benchmarks
BENCHMARK(flowless_sim)->Name("mio::Simulation on oseir::Model (pre 511 branch) without flows");
BENCHMARK(flow_sim_comp_only)->Name("mio::Simulation on oseir::Model with flows");
BENCHMARK(flow_sim)->Name("mio::FlowSimulation on oseir::Model with flows");
// run all benchmarks
BENCHMARK_MAIN();
