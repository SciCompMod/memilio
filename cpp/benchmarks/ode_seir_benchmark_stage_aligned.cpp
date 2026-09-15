/*
* Copyright (C) 2020-2026 MEmilio
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
#include "ode_seir_benchmark_stage_aligned.h"
#include "ode_seir_runtime_explicit.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace mio::benchmark_mio
{

void set_counters(benchmark::State& state, const StageAlignedProblem& problem, int threads)
{
    state.counters["patches"]                = problem.patches;
    state.counters["edges"]                  = static_cast<double>(problem.edges());
    state.counters["age_groups"]             = problem.groups;
    state.counters["steps"]                  = integration_steps;
    state.counters["implementation_version"] = stage_aligned_implementation_version;
    state.counters["temporal_block_steps"]   = temporal_block_steps;
    if (threads > 0) {
        state.counters["threads"] = threads;
    }
    state.SetItemsProcessed(state.iterations() * integration_steps * static_cast<int64_t>(problem.edges()));
}

void benchmark_stage_aligned(benchmark::State& state, int threads)
{
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    StageAlignedProblem problem(static_cast<int>(state.range(0)), static_cast<int>(state.range(1)),
                                static_cast<int>(state.range(2)), threads);
    for (auto _ : state) {
        state.PauseTiming();
        problem.reset_state(threads);
        state.ResumeTiming();
        advance_stage_aligned(problem, step_size, threads, integration_steps);
        benchmark::DoNotOptimize(problem.totals.data());
        benchmark::DoNotOptimize(problem.travelers.data());
        benchmark::ClobberMemory();
    }
    set_counters(state, problem, threads);
}

void benchmark_stage_aligned_serial(benchmark::State& state)
{
    benchmark_stage_aligned(state, 0);
}

#ifdef _OPENMP
inline constexpr std::array<int, 5> scalability_thread_counts           = {1, 16, 32, 64, 128};
inline constexpr std::array<std::pair<int, int>, 5> weak_scaling_shapes = {
    std::pair{512, 1}, std::pair{2048, 16}, std::pair{2896, 32}, std::pair{4096, 64}, std::pair{5792, 128}};

void benchmark_stage_aligned_openmp(benchmark::State& state)
{
    benchmark_stage_aligned(state, static_cast<int>(state.range(3)));
}

#endif

void apply_problem_shapes(benchmark::internal::Benchmark* benchmark)
{
    for (const auto& [patches, travelers] : problem_shapes) {
        for (int groups : age_group_counts) {
            benchmark->Args({patches, travelers, groups});
        }
    }
    benchmark->Args({stage_aligned_strong_scaling_patches, stage_aligned_strong_scaling_patches - 1, 6});
}

#ifdef _OPENMP
void apply_openmp_shapes(benchmark::internal::Benchmark* benchmark)
{
    constexpr int maximum_threads = scalability_thread_counts.back();
    for (const auto& [patches, travelers] : problem_shapes) {
        for (int groups : age_group_counts) {
            benchmark->Args({patches, travelers, groups, maximum_threads});
        }
    }
    for (int threads : scalability_thread_counts) {
        benchmark->Args({stage_aligned_strong_scaling_patches, stage_aligned_strong_scaling_patches - 1, 6, threads});
    }
    for (const auto& [weak_patches, threads] : weak_scaling_shapes) {
        benchmark->Args({weak_patches, weak_patches - 1, 6, threads});
    }
}
#endif

template <int G>
bool validate_blocked(int threads, std::string& error)
{
    constexpr int validation_patches = 263;
    double maximum_absolute_error    = 0.0;
    double maximum_relative_error    = 0.0;
    bool values_match                = true;
    const auto compare               = [&](const auto& expected, const auto& actual) {
        if (expected.size() != actual.size()) {
            values_match = false;
            return;
        }
        for (size_t index = 0; index < expected.size(); ++index) {
            const double absolute_error = std::abs(expected[index] - actual[index]);
            const double relative_error = absolute_error / std::max(1.0, std::abs(expected[index]));
            maximum_absolute_error      = std::max(maximum_absolute_error, absolute_error);
            maximum_relative_error      = std::max(maximum_relative_error, relative_error);
            if (!std::isfinite(expected[index]) || !std::isfinite(actual[index]) ||
                absolute_error > 1e-10 * (1.0 + std::abs(expected[index]))) {
                values_match = false;
            }
        }
    };
    // Tiny and partial tiles, multiple patches and a partial second time block.
    for (const auto& [patches, travelers] : {std::pair{3, 2}, std::pair{validation_patches, 262}}) {
        for (int steps : {1, 3, temporal_block_steps, temporal_block_steps + 1}) {
            StageAlignedProblem reference(patches, travelers, G);
            StageAlignedProblem serial(patches, travelers, G);
            StageAlignedProblem parallel(patches, travelers, G, threads);
            // Check parallel initialization, then break the proportionality of
            // traveler compartments to exercise independent traveler dynamics.
            compare(reference.totals, parallel.totals);
            compare(reference.travelers, parallel.travelers);
            for (size_t index = 0; index < reference.travelers.size(); ++index) {
                reference.travelers[index] *= 0.9 + 0.01 * ((index * 7) % 13);
            }
            for (int group = 0; group < G; ++group) {
                reference.rate_exposed[group] *= 1.0 + 0.02 * group;
                reference.rate_infected[group] *= 1.0 + 0.03 * group;
            }
            serial   = reference;
            parallel = reference;
            advance_stage_aligned_reference(reference, step_size, steps);
            advance_stage_aligned(serial, step_size, 0, steps);
            advance_stage_aligned(parallel, step_size, threads, steps);
            compare(reference.totals, serial.totals);
            compare(reference.travelers, serial.travelers);
            compare(reference.totals, parallel.totals);
            compare(reference.travelers, parallel.travelers);
            compare(reference.stage_lambda, serial.stage_lambda);
            compare(reference.stage_lambda, parallel.stage_lambda);
            serial.reset_state();
            parallel.reset_state(threads);
            compare(serial.totals, parallel.totals);
            compare(serial.travelers, parallel.travelers);
        }
    }

    std::cout << "Stage-aligned stepwise/blocked validation N_G=" << G << ": max_abs=" << maximum_absolute_error
              << ", max_rel=" << maximum_relative_error << '\n';
    if (!values_match) {
        error = "Stepwise and blocked stage-aligned results differ for N_G=" + std::to_string(G) +
                " (max_abs=" + std::to_string(maximum_absolute_error) +
                ", max_rel=" + std::to_string(maximum_relative_error) + ").";
    }
    return values_match;
}

bool validate_blocked_all(std::string& error)
{
#ifdef _OPENMP
    const int threads = omp_get_max_threads();
#else
    const int threads = 0;
#endif
    return validate_blocked<1>(threads, error) && validate_blocked<3>(threads, error) &&
           validate_blocked<6>(threads, error) && validate_blocked<8>(threads, error);
}
namespace scenario = mio::runtime_scenario;
template <int G>
void runtime_benchmark_impl(benchmark::State& state, int threads)
{
    try {
        scenario::Inputs inputs(static_cast<int>(state.range(0)), G);
        scenario::ExplicitProblem problem(inputs, threads);
        const auto time = scenario::schedule();
        for (auto _ : state) {
            state.PauseTiming();
            problem.reset();
            state.ResumeTiming();
            scenario::advance_explicit<G>(problem, time, threads);
            benchmark::DoNotOptimize(problem.core.totals.data());
            benchmark::DoNotOptimize(problem.core.travelers.data());
            benchmark::ClobberMemory();
        }
        scenario::check_population(inputs, problem.core.totals);
        scenario::counters(state, inputs, true, threads > 0 ? threads : 1);
    }
    catch (const std::exception& error) {
        state.SkipWithError(error.what());
    }
}
void runtime_dispatch(benchmark::State& state, int threads)
{
    switch (static_cast<int>(state.range(1))) {
    case 1:
        runtime_benchmark_impl<1>(state, threads);
        break;
    case 3:
        runtime_benchmark_impl<3>(state, threads);
        break;
    case 6:
        runtime_benchmark_impl<6>(state, threads);
        break;
    case 8:
        runtime_benchmark_impl<8>(state, threads);
        break;
    default:
        state.SkipWithError("Unsupported runtime age groups.");
    }
}
void runtime_serial(benchmark::State& state)
{
    runtime_dispatch(state, 0);
}
void runtime_openmp(benchmark::State& state)
{
    runtime_dispatch(state, static_cast<int>(state.range(2)));
}

void runtime_phase_openmp(benchmark::State& state)
{
    try {
        if (state.range(1) != 6) {
            state.SkipWithError("Explicit phase diagnostics require six age groups.");
            return;
        }
        const int threads = static_cast<int>(state.range(2));
        scenario::Inputs inputs(static_cast<int>(state.range(0)), 6);
        scenario::ExplicitProblem problem(inputs, threads);
        const auto time = scenario::schedule();
        std::vector<double> control_residents(inputs.initial.size()), profiled_residents(inputs.initial.size());
        // Google Benchmark invokes this callback separately for each raw
        // repetition. Alternate ordering per shape/thread pair to expose, not
        // systematically favor, effects of the immediately preceding run.
        static std::map<std::pair<int, int>, size_t> repetitions;
        const bool control_first = repetitions[{inputs.p, threads}]++ % 2 == 0;
        scenario::ExplicitPhaseTimes phases;
        double control_seconds = 0.0;
        const auto run_control = [&] {
            problem.reset();
            const auto start = scenario::ExplicitPhaseClock::now();
            scenario::advance_explicit<6>(problem, time, threads);
            control_seconds =
                std::chrono::duration<double>(scenario::ExplicitPhaseClock::now() - start).count();
            benchmark::DoNotOptimize(problem.core.totals.data());
            benchmark::DoNotOptimize(problem.core.travelers.data());
            benchmark::ClobberMemory();
            control_residents = problem.core.totals;
        };
        for (auto _ : state) {
            state.PauseTiming();
            if (control_first)
                run_control();
            problem.reset();
            state.ResumeTiming();
            phases = scenario::advance_explicit_profiled<6>(problem, time, threads);
            benchmark::DoNotOptimize(problem.core.totals.data());
            benchmark::DoNotOptimize(problem.core.travelers.data());
            benchmark::ClobberMemory();
            state.PauseTiming();
            profiled_residents = problem.core.totals;
            if (!control_first)
                run_control();
            // Both are independent trajectories from the same initial state.
            // compare() also checks finite/nonnegative states and population
            // conservation. Copying and validation are outside both timers.
            scenario::compare(inputs, control_residents, profiled_residents);
            state.ResumeTiming();
        }
        if (state.iterations() != 1) {
            state.SkipWithError("Phase diagnostics require exactly one paired trajectory per repetition.");
            return;
        }
        scenario::counters(state, inputs, true, threads);
        const double integration_seconds = phases.home_seconds + phases.away_seconds;
        state.counters["home_seconds"] = phases.home_seconds;
        state.counters["away_seconds"] = phases.away_seconds;
        state.counters["integration_seconds"] = integration_seconds;
        state.counters["departure_seconds"] = phases.departure_seconds;
        state.counters["return_seconds"] = phases.return_seconds;
        state.counters["profiled_seconds"] = phases.profiled_seconds;
        state.counters["unattributed_seconds"] =
            phases.profiled_seconds - integration_seconds - phases.departure_seconds - phases.return_seconds;
        state.counters["control_seconds"] = control_seconds;
        state.counters["control_first"] = control_first ? 1 : 0;
        state.counters["phase_validation_passed"] = 1;
    }
    catch (const std::exception& error) {
        state.SkipWithError(error.what());
    }
}

} // namespace mio::benchmark_mio

BENCHMARK(mio::benchmark_mio::benchmark_stage_aligned_serial)
    ->Apply(mio::benchmark_mio::apply_problem_shapes)
    ->ArgNames({"patches", "travelers_per_patch", "age_groups"})
    ->Name("stage_aligned/serial")
    ->UseRealTime();

#ifdef _OPENMP
BENCHMARK(mio::benchmark_mio::benchmark_stage_aligned_openmp)
    ->Apply(mio::benchmark_mio::apply_openmp_shapes)
    ->ArgNames({"patches", "travelers_per_patch", "age_groups", "threads"})
    ->Name("stage_aligned/openmp")
    ->UseRealTime();
#endif

int main(int argc, char** argv)
{
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
#ifdef _OPENMP
#pragma omp parallel
    {
        LIKWID_MARKER_THREADINIT;
    }
#endif
#endif

    const bool runtime = mio::runtime_scenario::enabled();
    if (runtime) {
        mio::runtime_scenario::register_shapes("runtime/explicit/serial", mio::benchmark_mio::runtime_serial);
#ifdef _OPENMP
        const auto callback = mio::runtime_scenario::experiment() == "phases" ?
                                  mio::benchmark_mio::runtime_phase_openmp : mio::benchmark_mio::runtime_openmp;
        mio::runtime_scenario::register_shapes("runtime/explicit/openmp", callback, omp_get_max_threads());
#endif
    }
    bool needs_validation = true;
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument(argv[index]);
        if (mio::runtime_scenario::informational_argument(argument)) {
            needs_validation = false;
        }
    }
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
#ifdef LIKWID_PERFMON
        LIKWID_MARKER_CLOSE;
#endif
        return 1;
    }
    std::string error;
    try {
        if (runtime && needs_validation) {
            mio::runtime_scenario::validate_accuracy();
#ifdef _OPENMP
            const int threads = omp_get_max_threads();
#else
            const int threads = 0;
#endif
            const bool validate_phases = mio::runtime_scenario::experiment() == "phases";
            mio::runtime_scenario::validate_explicit_cpu<1>(threads, validate_phases);
            mio::runtime_scenario::validate_explicit_cpu<3>(threads, validate_phases);
            mio::runtime_scenario::validate_explicit_cpu<6>(threads, validate_phases);
            mio::runtime_scenario::validate_explicit_cpu<8>(threads, validate_phases);
        }
    }
    catch (const std::exception& exception) {
        std::cerr << "Runtime validation failed: " << exception.what() << '\n';
        ::benchmark::Shutdown();
        return 1;
    }
    if (!runtime && needs_validation && !mio::benchmark_mio::validate_blocked_all(error)) {
        std::cerr << error << '\n';
        ::benchmark::Shutdown();
#ifdef LIKWID_PERFMON
        LIKWID_MARKER_CLOSE;
#endif
        return 1;
    }
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_CLOSE;
#endif
    return 0;
}
