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
#include <chrono>

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
    model.parameters.template set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.1);
    mio::ContactMatrixGroup& contact_matrix = model.parameters.template get<mio::oseir::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(10.);
    model.apply_constraints();
}

const ScalarType t0    = 0.0;
const ScalarType t_max = 50.0;
const ScalarType dt    = 0.5;

// Define the number of patches (fully connected) for the benchmark
// For N patches that are fully connected, we have N*(N-1) directed edges (e.g., for N=4 -> 12 edges).
const std::vector<int> number_patches = {4, 8, 16, 32, 64, 128, 256};

static void bench_euler_by_groups(::benchmark::State& state)
{
    const int num_patches       = number_patches[state.range(0)];
    const size_t num_age_groups = static_cast<size_t>(state.range(1));
    mio::set_log_level(mio::LogLevel::critical);

    ModelType model(num_age_groups);
    setup_model_benchmark(model);

    // accumulators for per-iteration averages
    double sum_advance_us  = 0.0;
    double sum_mobility_us = 0.0;

    for (auto _ : state) {
        state.PauseTiming();
        // Create one local simulation per patch
        std::vector<SimType> sims;
        sims.reserve(static_cast<size_t>(num_patches));
        for (int p = 0; p < num_patches; ++p) {
            sims.emplace_back(model, t0, dt);
            auto integrator_rk =
                std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
            sims.back().set_integrator(integrator_rk);
        }
        state.ResumeTiming();

        // Advance all local simulations
        const auto t_adv_begin = std::chrono::steady_clock::now();
        for (auto& s : sims) {
            s.advance(t_max);
        }
        const auto t_adv_end = std::chrono::steady_clock::now();
        sum_advance_us += std::chrono::duration<double, std::micro>(t_adv_end - t_adv_begin).count();

        state.PauseTiming();
        // Reference time series from the first patch to define the time grid
        const auto& seir_res0 = sims[0].get_result();
        double mobile_population_fraction =
            0.1 * static_cast<double>(num_patches) * (num_patches > 1 ? (num_patches - 1) : 1);
        const int num_edges = (num_patches > 1) ? (num_patches * (num_patches - 1)) : 0; // directed edges
        Eigen::VectorXd initial_mobile_pop =
            seir_res0.get_value(0) *
            (num_edges > 0 ? (mobile_population_fraction / static_cast<double>(num_edges)) : 0.0);
        const auto step_size_ref = seir_res0.get_time(1) - seir_res0.get_time(0);

        // Build edge list (u -> v) for a fully connected digraph without self-loops
        std::vector<std::pair<int, int>> edges;
        edges.reserve(static_cast<size_t>(num_edges));
        for (int u = 0; u < num_patches; ++u) {
            for (int v = 0; v < num_patches; ++v) {
                if (u != v) {
                    edges.emplace_back(u, v);
                }
            }
        }

        std::vector<Eigen::VectorXd> mobile_pops(static_cast<size_t>(edges.size()), initial_mobile_pop);
        state.ResumeTiming();

        const auto t_mob_begin = std::chrono::steady_clock::now();
        for (ScalarType t = t0; t < t_max; t += dt) {
            if (t + dt > t_max + 1e-10)
                break;

            const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
            if (closest_idx_total >= static_cast<size_t>(seir_res0.get_num_time_points()))
                break;

            // Update all directed edges (u -> v); use destination patch v for totals and model
            for (size_t e = 0; e < edges.size(); ++e) {
                const int v             = edges[e].second;
                const auto& total_pop_v = sims[static_cast<size_t>(v)].get_result().get_value(closest_idx_total);
                mio::examples::integrate_mobile_population_euler(mobile_pops[e], sims[static_cast<size_t>(v)],
                                                                 total_pop_v, t, dt);
            }
        }
        const auto t_mob_end = std::chrono::steady_clock::now();
        sum_mobility_us += std::chrono::duration<double, std::micro>(t_mob_end - t_mob_begin).count();
        benchmark::DoNotOptimize(mobile_pops); // Prevent mobile_pops calculation from being optimized away
    }

    // expose averages per iteration as GBM counters (microseconds)
    const double iters              = static_cast<double>(state.iterations());
    const double avg_adv_us         = (iters > 0.0) ? (sum_advance_us / iters) : 0.0;
    const double avg_mob_us         = (iters > 0.0) ? (sum_mobility_us / iters) : 0.0;
    const double avg_total_us       = avg_adv_us + avg_mob_us;
    const double adv_frac           = (avg_total_us > 0.0) ? (avg_adv_us / avg_total_us) : 0.0;
    const double mob_frac           = (avg_total_us > 0.0) ? (avg_mob_us / avg_total_us) : 0.0;
    state.counters["advance_us"]    = benchmark::Counter(avg_adv_us);
    state.counters["mobility_us"]   = benchmark::Counter(avg_mob_us);
    state.counters["advance_frac"]  = benchmark::Counter(adv_frac);
    state.counters["mobility_frac"] = benchmark::Counter(mob_frac);
}

static void bench_flow_based_by_groups(::benchmark::State& state)
{
    const int num_patches       = number_patches[state.range(0)];
    const size_t num_age_groups = static_cast<size_t>(state.range(1));
    mio::set_log_level(mio::LogLevel::critical);

    ModelType model(num_age_groups);
    setup_model_benchmark(model);

    // accumulators for per-iteration averages
    double sum_advance_us  = 0.0;
    double sum_mobility_us = 0.0;

    for (auto _ : state) {
        state.PauseTiming();
        // Create one local simulation per patch
        std::vector<SimType> sims;
        sims.reserve(static_cast<size_t>(num_patches));
        for (int p = 0; p < num_patches; ++p) {
            sims.emplace_back(model, t0, dt);
            auto integrator_rk =
                std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
            sims.back().set_integrator(integrator_rk);
        }
        state.ResumeTiming();

        // Advance all local simulations (measure)
        const auto t_adv_begin = std::chrono::steady_clock::now();
        for (auto& s : sims) {
            s.advance(t_max);
        }
        const auto t_adv_end = std::chrono::steady_clock::now();
        sum_advance_us += std::chrono::duration<double, std::micro>(t_adv_end - t_adv_begin).count();

        state.PauseTiming();
        const auto& seir_res0 = sims[0].get_result();
        const int num_edges   = (num_patches > 1) ? (num_patches * (num_patches - 1)) : 0; // directed edges
        double mobile_population_fraction =
            0.1 * static_cast<double>(num_patches) * (num_patches > 1 ? (num_patches - 1) : 1);
        Eigen::VectorXd initial_mobile_pop =
            seir_res0.get_value(0) *
            (num_edges > 0 ? (mobile_population_fraction / static_cast<double>(num_edges)) : 0.0);
        const auto step_size_ref = seir_res0.get_time(1) - seir_res0.get_time(0);

        // Build edge list (u -> v) for a fully connected digraph without self-loops
        std::vector<std::pair<int, int>> edges;
        edges.reserve(static_cast<size_t>(num_edges));
        for (int u = 0; u < num_patches; ++u) {
            for (int v = 0; v < num_patches; ++v) {
                if (u != v) {
                    edges.emplace_back(u, v);
                }
            }
        }

        std::vector<Eigen::VectorXd> mobile_pops(static_cast<size_t>(edges.size()), initial_mobile_pop);
        state.ResumeTiming();

        const auto t_mob_begin = std::chrono::steady_clock::now();
        for (ScalarType t = t0; t < t_max; t += dt) {
            if (t + dt > t_max + 1e-10)
                break;

            const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
            if (closest_idx_total >= static_cast<size_t>(seir_res0.get_num_time_points()))
                break;

            // Update all directed edges (u -> v); use destination patch v for totals and model/flows
            for (size_t e = 0; e < edges.size(); ++e) {
                const int v             = edges[e].second;
                const auto& total_pop_v = sims[static_cast<size_t>(v)].get_result().get_value(closest_idx_total);
                mio::examples::flow_based_mobility_returns(mobile_pops[e], sims[static_cast<size_t>(v)], total_pop_v, t,
                                                           dt);
            }
        }
        const auto t_mob_end = std::chrono::steady_clock::now();
        sum_mobility_us += std::chrono::duration<double, std::micro>(t_mob_end - t_mob_begin).count();
        benchmark::DoNotOptimize(mobile_pops);
    }

    // expose averages per iteration as GBM counters (microseconds)
    const double iters              = static_cast<double>(state.iterations());
    const double avg_adv_us         = (iters > 0.0) ? (sum_advance_us / iters) : 0.0;
    const double avg_mob_us         = (iters > 0.0) ? (sum_mobility_us / iters) : 0.0;
    const double avg_total_us       = avg_adv_us + avg_mob_us;
    const double adv_frac           = (avg_total_us > 0.0) ? (avg_adv_us / avg_total_us) : 0.0;
    const double mob_frac           = (avg_total_us > 0.0) ? (avg_mob_us / avg_total_us) : 0.0;
    state.counters["advance_us"]    = benchmark::Counter(avg_adv_us);
    state.counters["mobility_us"]   = benchmark::Counter(avg_mob_us);
    state.counters["advance_frac"]  = benchmark::Counter(adv_frac);
    state.counters["mobility_frac"] = benchmark::Counter(mob_frac);
}

static void bench_explicit_by_groups(::benchmark::State& state)
{
    const int num_commuter_groups = number_patches[state.range(0)];
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    mio::set_log_level(mio::LogLevel::critical);

    mio::examples::ExplicitModel model_explicit(num_age_groups, num_commuter_groups);

    const double total_population = 10000.0;

    // Start with non commuting pop (50%)
    for (auto i = 0; i < static_cast<int>(num_age_groups); i++) {
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::E}] =
            50.0 / static_cast<double>(num_age_groups);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::I}] =
            50.0 / static_cast<double>(num_age_groups);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::R}] =
            50.0 / static_cast<double>(num_age_groups);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::S}] =
            (total_population * 0.5) / static_cast<double>(num_age_groups) -
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

        for (auto i = 0; i < static_cast<int>(num_age_groups); i++) {
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::E}] =
                50.0 / static_cast<double>(num_age_groups * num_commuter_groups);
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::I}] =
                50.0 / static_cast<double>(num_age_groups * num_commuter_groups);
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::R}] =
                50.0 / static_cast<double>(num_age_groups * num_commuter_groups);

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
                (total_population * commuter_fraction) / static_cast<double>(num_age_groups) - exposed - infected -
                recovered;
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
    ->Apply([](auto* b) {
        for (int i = 0; i < static_cast<int>(mio::benchmark_mio::number_patches.size()); ++i) {
            for (int g = 1; g <= 6; ++g) {
                b->Args({i, g}); // range(0)=patch index, range(1)=num age groups
            }
        }
    }) -> Name("Euler") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_flow_based_by_groups)
    ->Apply([](auto* b) {
        for (int i = 0; i < static_cast<int>(mio::benchmark_mio::number_patches.size()); ++i) {
            for (int g = 1; g <= 6; ++g) {
                b->Args({i, g});
            }
        }
    }) -> Name("Flow-based") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_explicit_by_groups)
    ->Apply([](auto* b) {
        for (int i = 0; i < static_cast<int>(mio::benchmark_mio::number_patches.size()); ++i) {
            for (int g = 1; g <= 6; ++g) {
                b->Args({i, g});
            }
        }
    }) -> Name("Explicit") -> Unit(::benchmark::kMicrosecond);

// run all benchmarks
BENCHMARK_MAIN();
