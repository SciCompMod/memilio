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
#include <cstdint>
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
void benchmark_stage_aligned_openmp(benchmark::State& state)
{
    benchmark_stage_aligned(state, static_cast<int>(state.range(3)));
}

std::vector<int> benchmark_thread_counts()
{
    const int maximum       = omp_get_max_threads();
    std::vector<int> counts = {1, std::max(1, maximum / 8), std::max(1, maximum / 4), std::max(1, maximum / 2),
                               maximum};
    std::sort(counts.begin(), counts.end());
    counts.erase(std::unique(counts.begin(), counts.end()), counts.end());
    return counts;
}
#endif

void apply_problem_shapes(benchmark::internal::Benchmark* benchmark)
{
    for (const auto& [patches, travelers] : problem_shapes) {
        for (int groups : age_group_counts) {
            benchmark->Args({patches, travelers, groups});
        }
    }
}

#ifdef _OPENMP
void apply_openmp_shapes(benchmark::internal::Benchmark* benchmark)
{
    const auto thread_counts  = benchmark_thread_counts();
    const int maximum_threads = thread_counts.back();
    for (const auto& [patches, travelers] : problem_shapes) {
        for (int groups : age_group_counts) {
            benchmark->Args({patches, travelers, groups, maximum_threads});
        }
    }
    const auto [patches, travelers] = problem_shapes.back();
    for (int threads : thread_counts) {
        if (threads != maximum_threads) {
            benchmark->Args({patches, travelers, 6, threads});
        }
    }
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

#ifdef LIKWID_PERFMON
int main(int argc, char** argv)
{
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
#ifdef _OPENMP
#pragma omp parallel
    {
        LIKWID_MARKER_THREADINIT;
    }
#endif
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
        LIKWID_MARKER_CLOSE;
        return 1;
    }
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
    LIKWID_MARKER_CLOSE;
    return 0;
}
#else
BENCHMARK_MAIN();
#endif
