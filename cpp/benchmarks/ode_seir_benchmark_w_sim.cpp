/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
namespace benchmark_mio
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
const size_t num_age_groups_bench = 6;

// Define the number of commuter groups for the benchmark
const std::vector<int> commuter_group_counts = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

static void bench_euler_by_groups(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
    mio::set_log_level(mio::LogLevel::critical);

    ModelType model(num_age_groups_bench);
    setup_model_benchmark(model);

    for (auto _ : state) {
        state.PauseTiming();
        SimType sim(model, t0, dt);
        auto integrator_rk =
            std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
        sim.set_integrator(integrator_rk);
        state.ResumeTiming();

        sim.advance(t_max);

        state.PauseTiming();
        const auto& seir_res               = sim.get_result();
        double mobile_population_fraction  = 0.1 * num_commuter_groups;
        Eigen::VectorXd initial_mobile_pop = seir_res.get_value(0) * (mobile_population_fraction / num_commuter_groups);
        const auto step_size_ref           = seir_res.get_time(1) - seir_res.get_time(0);

        std::vector<Eigen::VectorXd> mobile_pops(num_commuter_groups, initial_mobile_pop);
        state.ResumeTiming();

        for (ScalarType t = t0; t < t_max; t += dt) {
            if (t + dt > t_max + 1e-10)
                break;

            const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
            if (closest_idx_total >= static_cast<size_t>(seir_res.get_num_time_points()))
                break;

            const auto& total_pop = seir_res.get_value(closest_idx_total);

            // Für jede Pendlergruppe die Euler-Methode anwenden
            for (int i = 0; i < num_commuter_groups; ++i) {
                mio::examples::integrate_mobile_population_euler(mobile_pops[i], sim, total_pop, t, dt);
            }
        }
        benchmark::DoNotOptimize(mobile_pops); // Prevent mobile_pops calculation from being optimized away
    }
}

static void bench_high_order_by_groups(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
    mio::set_log_level(mio::LogLevel::critical);

    ModelType model(num_age_groups_bench);
    setup_model_benchmark(model);

    for (auto _ : state) {
        state.PauseTiming();
        SimType sim(model, t0, dt);
        auto integrator_rk =
            std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
        sim.set_integrator(integrator_rk);
        state.ResumeTiming();

        sim.advance(t_max);

        state.PauseTiming();
        const auto& seir_res               = sim.get_result();
        double mobile_population_fraction  = 0.1 * num_commuter_groups;
        Eigen::VectorXd initial_mobile_pop = seir_res.get_value(0) * (mobile_population_fraction / num_commuter_groups);
        const auto step_size_ref           = seir_res.get_time(1) - seir_res.get_time(0);

        std::vector<Eigen::VectorXd> mobile_pops(num_commuter_groups, initial_mobile_pop);
        state.ResumeTiming();

        for (ScalarType t = t0; t < t_max; t += dt) {
            if (t + dt > t_max + 1e-10)
                break;

            const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
            if (closest_idx_total >= static_cast<size_t>(seir_res.get_num_time_points()))
                break;

            const auto& total_pop = seir_res.get_value(closest_idx_total);

            for (int i = 0; i < num_commuter_groups; ++i) {
                mio::examples::integrate_mobile_population_high_order(mobile_pops[i], sim, total_pop, t, dt);
            }
        }
        benchmark::DoNotOptimize(mobile_pops); // Prevent mobile_pops calculation from being optimized away
    }
}

static void bench_flow_based_by_groups(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
    mio::set_log_level(mio::LogLevel::critical);

    ModelType model(num_age_groups_bench);
    setup_model_benchmark(model);

    for (auto _ : state) {
        state.PauseTiming();
        SimType sim(model, t0, dt);
        auto integrator_rk =
            std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
        sim.set_integrator(integrator_rk);
        state.ResumeTiming();

        sim.advance(t_max); // Measure advance

        state.PauseTiming();
        const auto& seir_res               = sim.get_result();
        double mobile_population_fraction  = 0.1 * num_commuter_groups;
        Eigen::VectorXd initial_mobile_pop = seir_res.get_value(0) * (mobile_population_fraction / num_commuter_groups);
        const auto step_size_ref           = seir_res.get_time(1) - seir_res.get_time(0);

        std::vector<Eigen::VectorXd> mobile_pops(num_commuter_groups, initial_mobile_pop);
        state.ResumeTiming();

        for (ScalarType t = t0; t < t_max; t += dt) {
            if (t + dt > t_max + 1e-10)
                break;

            const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
            if (closest_idx_total >= static_cast<size_t>(seir_res.get_num_time_points()))
                break;

            const auto& total_pop = seir_res.get_value(closest_idx_total);
            const auto& flows     = sim.get_flows().get_value(closest_idx_total);

            for (int i = 0; i < num_commuter_groups; ++i) {
                mio::examples::flow_based_mobility_returns(mobile_pops[i], sim, total_pop, t, dt);
            }
        }
        benchmark::DoNotOptimize(mobile_pops);
    }
}

static void bench_probabilistic_by_groups(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
    mio::set_log_level(mio::LogLevel::critical);

    ModelType model(num_age_groups_bench);
    setup_model_benchmark(model);

    for (auto _ : state) {
        state.PauseTiming();
        SimType sim(model, t0, dt);
        auto integrator_rk =
            std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
        sim.set_integrator(integrator_rk);
        state.ResumeTiming();

        sim.advance(t_max); // Measure advance

        state.PauseTiming();
        const auto& seir_res               = sim.get_result();
        double mobile_population_fraction  = 0.1 * num_commuter_groups;
        Eigen::VectorXd initial_mobile_pop = seir_res.get_value(0) * (mobile_population_fraction / num_commuter_groups);
        const auto step_size_ref           = seir_res.get_time(1) - seir_res.get_time(0);

        std::vector<Eigen::VectorXd> mobile_pops(num_commuter_groups, initial_mobile_pop);
        state.ResumeTiming();

        for (ScalarType t = t0; t < t_max; t += dt) {
            if (t + dt > t_max + 1e-10)
                break;

            const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
            if (closest_idx_total >= static_cast<size_t>(seir_res.get_num_time_points()))
                break;

            const auto& total_pop = seir_res.get_value(closest_idx_total);

            for (int i = 0; i < num_commuter_groups; ++i) {
                mio::examples::probabilistic_mobility_returns(mobile_pops[i], sim, total_pop, t, dt);
            }
        }
        benchmark::DoNotOptimize(mobile_pops);
    }
}

static void bench_explicit_by_groups(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
    mio::set_log_level(mio::LogLevel::critical);

    mio::examples::ExplicitModel model_explicit(num_age_groups_bench, num_commuter_groups);

    const double total_population = 10000.0;

    // Start with non commuting pop (50%)
    for (auto i = 0; i < num_age_groups_bench; i++) {
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::E}] =
            50.0 / static_cast<double>(num_age_groups_bench);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::I}] =
            50.0 / static_cast<double>(num_age_groups_bench);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::R}] =
            50.0 / static_cast<double>(num_age_groups_bench);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::S}] =
            (total_population * 0.5) / static_cast<double>(num_age_groups_bench) -
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::E}] -
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::I}] -
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::R}];
    }

    // Then the remaining commuting pop
    for (int c = 0; c < num_commuter_groups; c++) {
        auto commuter_type =
            static_cast<mio::examples::CommuterType>(static_cast<int>(mio::examples::CommuterType::CommuterBase) + c);
        double commuter_fraction = 0.5 / static_cast<double>(num_commuter_groups);

        for (auto i = 0; i < num_age_groups_bench; i++) {
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::E}] =
                50.0 / static_cast<double>(num_age_groups_bench * num_commuter_groups);
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::I}] =
                50.0 / static_cast<double>(num_age_groups_bench * num_commuter_groups);
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::R}] =
                50.0 / static_cast<double>(num_age_groups_bench * num_commuter_groups);

            double exposed   = model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                                   mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::E}];
            double infected  = model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                                    mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::I}];
            double recovered = model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                                     mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::R}];

            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::S}] =
                (total_population * commuter_fraction) / static_cast<double>(num_age_groups_bench) - exposed -
                infected - recovered;
        }
    }

    model_explicit.parameters.set<mio::oseir::TimeExposed<ScalarType>>(5.2);
    model_explicit.parameters.set<mio::oseir::TimeInfected<ScalarType>>(6);
    model_explicit.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.04);
    mio::ContactMatrixGroup& contact_matrix = model_explicit.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(10.);
    model_explicit.check_constraints();

    for (auto _ : state) {
        state.PauseTiming();
        mio::examples::ExplicitSim sim_explicit(model_explicit, t0, dt);
        auto integrator_rk =
            std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
        sim_explicit.set_integrator(integrator_rk);
        state.ResumeTiming();

        // Only measure the advance call
        sim_explicit.advance(t_max);
    }
}

} // namespace benchmark_mio
} // namespace mio

BENCHMARK(mio::benchmark_mio::bench_euler_by_groups)
    ->DenseRange(0, mio::benchmark_mio::commuter_group_counts.size() - 1)
    ->Name("Euler")
    ->Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_high_order_by_groups)
    ->DenseRange(0, mio::benchmark_mio::commuter_group_counts.size() - 1)
    ->Name("High-order")
    ->Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_flow_based_by_groups)
    ->DenseRange(0, mio::benchmark_mio::commuter_group_counts.size() - 1)
    ->Name("Flow-based")
    ->Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_probabilistic_by_groups)
    ->DenseRange(0, mio::benchmark_mio::commuter_group_counts.size() - 1)
    ->Name("Probabilistic")
    ->Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_explicit_by_groups)
    ->DenseRange(0, mio::benchmark_mio::commuter_group_counts.size() - 1)
    ->Name("Explicit")
    ->Unit(::benchmark::kMicrosecond);

// run all benchmarks
BENCHMARK_MAIN();
