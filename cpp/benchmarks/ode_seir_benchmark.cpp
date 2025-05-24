/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik ZUnker
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

#include "benchmark/benchmark.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"
#include "ode_seir/model.h"
#include "examples/state_estimators.h"

namespace mio
{
namespace benchmark
{

using ModelType = mio::oseir::Model<ScalarType>;
using SimType   = mio::FlowSimulation<ScalarType, ModelType>;

void setup_model_benchmark(ModelType& model, size_t num_agegroups = 6)
{
    const double total_population = 10000.0;
    for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); i++) {
        model.populations[{i, mio::oseir::InfectionState::Exposed}]   = 100.0 / static_cast<double>(num_agegroups);
        model.populations[{i, mio::oseir::InfectionState::Infected}]  = 100.0 / static_cast<double>(num_agegroups);
        model.populations[{i, mio::oseir::InfectionState::Recovered}] = 100.0 / static_cast<double>(num_agegroups);
        model.populations[{i, mio::oseir::InfectionState::Susceptible}] =
            total_population / static_cast<double>(num_agegroups) -
            model.populations[{i, mio::oseir::InfectionState::Exposed}] -
            model.populations[{i, mio::oseir::InfectionState::Infected}] -
            model.populations[{i, mio::oseir::InfectionState::Recovered}];
    }
    model.parameters.template set<mio::oseir::TimeExposed<ScalarType>>(5.2);
    model.parameters.template set<mio::oseir::TimeInfected<ScalarType>>(6);
    model.parameters.template set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.04);
    mio::ContactMatrixGroup& contact_matrix = model.parameters.template get<mio::oseir::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(10.);
    model.apply_constraints();
}

const ScalarType t0               = 0.0;
const ScalarType t_max            = 50.0;
const ScalarType dt               = 0.5;
const size_t num_age_groups_bench = 1;

// Benchmark for Euler integration method
void benchmark_state_euler(::benchmark::State& state)
{
    mio::set_log_level(mio::LogLevel::critical);
    ModelType model(num_age_groups_bench);
    setup_model_benchmark(model);
    SimType sim(model, t0, dt / 10.0);
    sim.advance(t_max);
    const auto& seir_res = sim.get_result();

    Eigen::VectorXd initial_mobile_pop = seir_res.get_value(0) * 0.1;
    const auto step_size_ref           = seir_res.get_time(1) - seir_res.get_time(0);

    for (auto _ : state) {
        Eigen::VectorXd mobile_pop = initial_mobile_pop;
        for (ScalarType t = t0; t < t_max; t += dt) {
            if (t + dt > t_max + 1e-10) {
                break;
            }
            const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
            if (closest_idx_total >= static_cast<size_t>(seir_res.get_num_time_points())) {
                break;
            }
            const auto& total_pop = seir_res.get_value(closest_idx_total);
            mio::examples::integrate_mobile_population_euler(mobile_pop, sim, total_pop, t, dt);
        }
    }
}

// Benchmark for high order integration method
void benchmark_state_high_order(::benchmark::State& state)
{
    mio::set_log_level(mio::LogLevel::critical);
    ModelType model(num_age_groups_bench);
    setup_model_benchmark(model);

    SimType sim(model, t0, dt / 10.0);
    sim.advance(t_max);
    const auto& seir_res = sim.get_result();

    Eigen::VectorXd initial_mobile_pop = seir_res.get_value(0) * 0.1;
    const auto step_size_ref           = seir_res.get_time(1) - seir_res.get_time(0);

    for (auto _ : state) {
        Eigen::VectorXd mobile_pop = initial_mobile_pop;
        for (ScalarType t = t0; t < t_max; t += dt) {
            if (t + dt > t_max + 1e-10) {
                break;
            }
            const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
            if (closest_idx_total >= static_cast<size_t>(seir_res.get_num_time_points())) {
                break;
            }
            const auto& total_pop = seir_res.get_value(closest_idx_total);
            mio::examples::integrate_mobile_population_high_order(mobile_pop, sim, total_pop, t, dt);
        }
    }
}

// Benchmark for Flow-based mobility returns method
void benchmark_state_flow_based(::benchmark::State& state)
{
    mio::set_log_level(mio::LogLevel::critical);
    ModelType model(num_age_groups_bench);
    setup_model_benchmark(model);

    SimType sim(model, t0, dt);
    sim.advance(t_max);
    const auto& seir_res = sim.get_result();

    Eigen::VectorXd initial_mobile_pop = seir_res.get_value(0) * 0.1;
    const auto step_size_ref           = seir_res.get_time(1) - seir_res.get_time(0);

    for (auto _ : state) {
        Eigen::VectorXd mobile_pop = initial_mobile_pop;
        for (ScalarType t = t0; t < t_max; t += dt) {
            if (t + dt > t_max + 1e-10) {
                break;
            }
            const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
            if (closest_idx_total >= static_cast<size_t>(seir_res.get_num_time_points())) {
                break;
            }
            const auto& total_pop = seir_res.get_value(closest_idx_total);
            mio::examples::apply_flows_to_mobile_population_generic(mobile_pop, sim.get_model(), total_pop,
                                                                    sim.get_flows().get_value(closest_idx_total));
        }
    }
}

void benchmark_state_flow_based_model_specific(::benchmark::State& state)
{
    mio::set_log_level(mio::LogLevel::critical);
    ModelType model(num_age_groups_bench);
    setup_model_benchmark(model);

    SimType sim(model, t0, dt);
    sim.advance(t_max);
    const auto& seir_res = sim.get_result();

    Eigen::VectorXd initial_mobile_pop = seir_res.get_value(0) * 0.1;
    const auto step_size_ref           = seir_res.get_time(1) - seir_res.get_time(0);

    for (auto _ : state) {
        Eigen::VectorXd mobile_pop = initial_mobile_pop;
        for (ScalarType t = t0; t < t_max; t += dt) {
            if (t + dt > t_max + 1e-10) {
                break;
            }
            const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
            if (closest_idx_total >= static_cast<size_t>(seir_res.get_num_time_points())) {
                break;
            }
            const auto& total_pop = seir_res.get_value(closest_idx_total);
            mio::examples::apply_flows_to_mobile_population(mobile_pop, sim.get_model(), total_pop,
                                                            sim.get_flows().get_value(closest_idx_total));
        }
    }
}

// Benchmark for Probabilistic mobility returns method
void benchmark_state_probabilistic(::benchmark::State& state)
{
    mio::set_log_level(mio::LogLevel::critical);
    ModelType model(num_age_groups_bench);
    setup_model_benchmark(model);

    SimType sim(model, t0, dt / 10.0);
    sim.advance(t_max);
    const auto& seir_res = sim.get_result();

    Eigen::VectorXd initial_mobile_pop = seir_res.get_value(0) * 0.1;
    const auto step_size_ref           = seir_res.get_time(1) - seir_res.get_time(0);

    for (auto _ : state) {
        Eigen::VectorXd mobile_pop = initial_mobile_pop;
        for (ScalarType t = t0; t < t_max; t += dt) {
            if (t + dt > t_max + 1e-10) {
                break;
            }
            const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
            if (closest_idx_total >= static_cast<size_t>(seir_res.get_num_time_points())) {
                break;
            }
            const auto& total_pop = seir_res.get_value(closest_idx_total);
            mio::examples::probabilistic_mobility_returns(mobile_pop, sim, total_pop, t, dt);
        }
    }
}

} // namespace benchmark
} // namespace mio

BENCHMARK(mio::benchmark::benchmark_state_euler)->Name("Euler");
BENCHMARK(mio::benchmark::benchmark_state_high_order)->Name("High-order");
BENCHMARK(mio::benchmark::benchmark_state_flow_based)->Name("Flow-based");
BENCHMARK(mio::benchmark::benchmark_state_flow_based_model_specific)->Name("Flow-based Model Specific");
BENCHMARK(mio::benchmark::benchmark_state_probabilistic)->Name("Probabilistic");

// run all benchmarks
BENCHMARK_MAIN();
