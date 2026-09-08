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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
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
    state.counters["patches"]    = problem.patches;
    state.counters["edges"]      = static_cast<double>(problem.edges());
    state.counters["age_groups"] = problem.groups;
    state.counters["steps"]      = integration_steps;
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
                                static_cast<int>(state.range(2)));
    for (auto _ : state) {
        state.PauseTiming();
        problem.reset_state();
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
inline constexpr std::array<int, 5> scalability_thread_counts = {1, 16, 32, 64, 128};
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
    benchmark->Args(
        {stage_aligned_strong_scaling_patches, stage_aligned_strong_scaling_patches - 1, 6});
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
        benchmark->Args(
            {stage_aligned_strong_scaling_patches, stage_aligned_strong_scaling_patches - 1, 6, threads});
    }
    for (const auto& [weak_patches, threads] : weak_scaling_shapes) {
        benchmark->Args({weak_patches, weak_patches - 1, 6, threads});
    }
}

template <int G>
bool validate_openmp(int threads, std::string& error)
{
    constexpr int validation_patches = 263;
    StageAlignedProblem serial(validation_patches, validation_patches - 1, G);
    StageAlignedProblem parallel(validation_patches, validation_patches - 1, G);
    advance_stage_aligned(serial, step_size, 0, integration_steps);
    advance_stage_aligned(parallel, step_size, threads, integration_steps);

    double maximum_absolute_error = 0.0;
    double maximum_relative_error = 0.0;
    bool values_match             = true;
    const auto compare = [&](const std::vector<double>& expected, const std::vector<double>& actual) {
        for (size_t index = 0; index < expected.size(); ++index) {
            const double absolute_error = std::abs(expected[index] - actual[index]);
            const double relative_error = absolute_error / std::max(1.0, std::abs(expected[index]));
            maximum_absolute_error      = std::max(maximum_absolute_error, absolute_error);
            maximum_relative_error      = std::max(maximum_relative_error, relative_error);
            if (!std::isfinite(actual[index]) || absolute_error > 1e-10 * (1.0 + std::abs(expected[index]))) {
                values_match = false;
            }
        }
    };
    compare(serial.totals, parallel.totals);
    compare(serial.travelers, parallel.travelers);

    std::cout << "Stage-aligned OpenMP validation N_G=" << G << ": max_abs=" << maximum_absolute_error
              << ", max_rel=" << maximum_relative_error << '\n';
    if (!values_match) {
        error = "Serial and OpenMP stage-aligned results differ for N_G=" + std::to_string(G) +
                " (max_abs=" + std::to_string(maximum_absolute_error) +
                ", max_rel=" + std::to_string(maximum_relative_error) + ").";
    }
    return values_match;
}

bool validate_openmp_all(std::string& error)
{
    const int threads = omp_get_max_threads();
    return validate_openmp<1>(threads, error) && validate_openmp<3>(threads, error) &&
           validate_openmp<6>(threads, error) && validate_openmp<8>(threads, error);
}
#endif

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

    bool needs_validation = true;
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument(argv[index]);
        if (argument.starts_with("--benchmark_list_tests") || argument == "--help" || argument == "-h") {
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
#ifdef _OPENMP
    std::string error;
    if (needs_validation && !mio::benchmark_mio::validate_openmp_all(error)) {
        std::cerr << error << '\n';
        ::benchmark::Shutdown();
#ifdef LIKWID_PERFMON
        LIKWID_MARKER_CLOSE;
#endif
        return 1;
    }
#endif
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_CLOSE;
#endif
    return 0;
}
