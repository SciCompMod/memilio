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
#include "ode_seir_metapop_benchmark.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace mio::benchmark_mio
{

#ifdef _OPENMP
template <int G, int Stage>
void advance_openmp_stage(ImplicitProblem& problem, double dt)
{
    const int patches                  = problem.patches;
    const std::vector<double>& current = Stage == 0 ? problem.state : problem.stage_state;

#pragma omp for schedule(static)
    for (int patch = 0; patch < patches; ++patch) {
        std::array<double, G> aggregate{};
        const double* commuting =
            problem.commuting_transpose_row_major.data() + static_cast<size_t>(patch) * patches;
        if constexpr (G == 1) {
            double value = 0.0;
#pragma omp simd reduction(+ : value)
            for (int origin = 0; origin < patches; ++origin) {
                value += commuting[origin] * current[static_cast<size_t>(2) * patches + origin];
            }
            aggregate[0] = value;
        }
        else {
            for (int origin = 0; origin < patches; ++origin) {
                const double mobility = commuting[origin];
                for (int group = 0; group < G; ++group) {
                    aggregate[group] = std::fma(
                        mobility, current[static_cast<size_t>(3 * group + 2) * patches + origin], aggregate[group]);
                }
            }
        }
        for (int group = 0; group < G; ++group) {
            const size_t index = static_cast<size_t>(group) * patches + patch;
            problem.present_share[index] = aggregate[group] * problem.inverse_present_population[index];
        }
    }

#pragma omp for schedule(static)
    for (int patch = 0; patch < patches; ++patch) {
        std::array<double, G> aggregate{};
        const double* commuting = problem.commuting_row_major.data() + static_cast<size_t>(patch) * patches;
        if constexpr (G == 1) {
            double value = 0.0;
#pragma omp simd reduction(+ : value)
            for (int destination = 0; destination < patches; ++destination) {
                value += commuting[destination] * problem.present_share[destination];
            }
            aggregate[0] = value;
        }
        else {
            for (int destination = 0; destination < patches; ++destination) {
                const double mobility = commuting[destination];
                for (int group = 0; group < G; ++group) {
                    aggregate[group] = std::fma(
                        mobility, problem.present_share[static_cast<size_t>(group) * patches + destination],
                        aggregate[group]);
                }
            }
        }
        advance_cpu_patch<G, Stage>(problem, current, patch, aggregate, dt);
    }
}

template <int G>
void advance_openmp(ImplicitProblem& problem, int threads, double dt = step_size, int steps = integration_steps)
{
    if (threads <= 0) {
        throw std::invalid_argument("The OpenMP thread count must be positive.");
    }
#pragma omp parallel num_threads(threads)
    {
        for (int step = 0; step < steps; ++step) {
            advance_openmp_stage<G, 0>(problem, dt);
            advance_openmp_stage<G, 1>(problem, dt);
            advance_openmp_stage<G, 2>(problem, dt);
            advance_openmp_stage<G, 3>(problem, dt);
        }
    }
}
#endif

void set_counters(benchmark::State& state, const ImplicitProblem& problem, int threads = 0)
{
    state.counters["patches"]    = problem.patches;
    state.counters["age_groups"] = problem.groups;
    state.counters["steps"]      = integration_steps;
    if (threads > 0) {
        state.counters["threads"] = threads;
    }
}

template <int G>
void benchmark_serial_impl(benchmark::State& state)
{
    ImplicitProblem problem(static_cast<int>(state.range(0)), G);
    for (auto _ : state) {
        state.PauseTiming();
        problem.reset_state();
        state.ResumeTiming();
        advance_cpu<G>(problem);
        benchmark::DoNotOptimize(problem.state.data());
        benchmark::ClobberMemory();
    }
    set_counters(state, problem);
}

void benchmark_serial(benchmark::State& state)
{
    switch (static_cast<int>(state.range(1))) {
    case 1:
        benchmark_serial_impl<1>(state);
        break;
    case 3:
        benchmark_serial_impl<3>(state);
        break;
    case 6:
        benchmark_serial_impl<6>(state);
        break;
    case 8:
        benchmark_serial_impl<8>(state);
        break;
    default:
        state.SkipWithError("Unsupported age-group count.");
    }
}

#ifdef _OPENMP
template <int G>
void benchmark_openmp_impl(benchmark::State& state)
{
    const int threads = static_cast<int>(state.range(2));
    ImplicitProblem problem(static_cast<int>(state.range(0)), G);
    for (auto _ : state) {
        state.PauseTiming();
        problem.reset_state();
        state.ResumeTiming();
        advance_openmp<G>(problem, threads);
        benchmark::DoNotOptimize(problem.state.data());
        benchmark::ClobberMemory();
    }
    set_counters(state, problem, threads);
}

void benchmark_openmp(benchmark::State& state)
{
    switch (static_cast<int>(state.range(1))) {
    case 1:
        benchmark_openmp_impl<1>(state);
        break;
    case 3:
        benchmark_openmp_impl<3>(state);
        break;
    case 6:
        benchmark_openmp_impl<6>(state);
        break;
    case 8:
        benchmark_openmp_impl<8>(state);
        break;
    default:
        state.SkipWithError("Unsupported age-group count.");
    }
}
#endif

void apply_shapes(benchmark::internal::Benchmark* benchmark)
{
    for (int patches : patch_counts) {
        for (int groups : age_group_counts) {
            benchmark->Args({patches, groups});
        }
    }
}

#ifdef _OPENMP
std::vector<int> benchmark_thread_counts()
{
    const int maximum       = omp_get_max_threads();
    std::vector<int> counts = {1, std::max(1, maximum / 8), std::max(1, maximum / 4), std::max(1, maximum / 2),
                               maximum};
    std::sort(counts.begin(), counts.end());
    counts.erase(std::unique(counts.begin(), counts.end()), counts.end());
    return counts;
}

void apply_openmp_shapes(benchmark::internal::Benchmark* benchmark)
{
    const auto thread_counts  = benchmark_thread_counts();
    const int maximum_threads = thread_counts.back();
    for (int patches : patch_counts) {
        for (int groups : age_group_counts) {
            benchmark->Args({patches, groups, maximum_threads});
        }
    }
    for (int threads : thread_counts) {
        if (threads != maximum_threads) {
            benchmark->Args({patch_counts.back(), 6, threads});
        }
    }
}

template <int G>
bool validate_openmp(int threads, std::string& error)
{
    constexpr int validation_patches = 263;
    ImplicitProblem serial(validation_patches, G);
    ImplicitProblem parallel(validation_patches, G);
    advance_cpu<G>(serial);
    advance_openmp<G>(parallel, threads);

    double maximum_absolute_error = 0.0;
    double maximum_relative_error = 0.0;
    bool values_match             = true;
    for (size_t index = 0; index < serial.state.size(); ++index) {
        const double absolute_error = std::abs(serial.state[index] - parallel.state[index]);
        const double relative_error =
            absolute_error / std::max(validation_absolute_denominator, std::abs(serial.state[index]));
        maximum_absolute_error = std::max(maximum_absolute_error, absolute_error);
        maximum_relative_error = std::max(maximum_relative_error, relative_error);
        if (!std::isfinite(parallel.state[index]) ||
            absolute_error > validation_relative_tolerance * (1.0 + std::abs(serial.state[index]))) {
            values_match = false;
        }
    }
    std::cout << "OpenMP validation N_G=" << G << ": max_abs=" << maximum_absolute_error
              << ", max_rel=" << maximum_relative_error << '\n';
    if (!values_match) {
        error = "Serial and OpenMP results differ for N_G=" + std::to_string(G) +
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

BENCHMARK(mio::benchmark_mio::benchmark_serial)
    ->Apply(mio::benchmark_mio::apply_shapes)
    ->ArgNames({"patches", "age_groups"})
    ->Name("implicit/serial")
    ->UseRealTime();

#ifdef _OPENMP
BENCHMARK(mio::benchmark_mio::benchmark_openmp)
    ->Apply(mio::benchmark_mio::apply_openmp_shapes)
    ->ArgNames({"patches", "age_groups", "threads"})
    ->Name("implicit/openmp")
    ->UseRealTime();
#endif

int main(int argc, char** argv)
{
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
        return 1;
    }
#ifdef _OPENMP
    std::string error;
    if (!mio::benchmark_mio::validate_openmp_all(error)) {
        std::cerr << error << '\n';
        ::benchmark::Shutdown();
        return 1;
    }
#endif
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
    return 0;
}
