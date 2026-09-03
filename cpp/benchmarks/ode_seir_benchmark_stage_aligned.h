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
#ifndef MIO_ODE_SEIR_BENCHMARK_STAGE_ALIGNED_H
#define MIO_ODE_SEIR_BENCHMARK_STAGE_ALIGNED_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#endif

namespace mio::benchmark_mio
{

constexpr double step_size                                  = 0.5;
constexpr double total_traveler_fraction                    = 0.1;
constexpr int traveler_chunk_size                           = 256;
constexpr int integration_steps                             = 64;
constexpr std::array<int, 4> age_group_counts               = {1, 3, 6, 8};
constexpr std::array<std::pair<int, int>, 7> problem_shapes = {
    std::pair{16, 15},   std::pair{32, 31},   std::pair{64, 63},    std::pair{128, 127},
    std::pair{256, 255}, std::pair{512, 511}, std::pair{1024, 1023}};

/**
 * Synthetic, destination-oriented SEIR state used by all three benchmarks.
 * It represents the continuous integration between two mobility events; the
 * discrete event update itself is deliberately not part of this kernel benchmark.
 *
 * Aggregate states use the layout totals[compartment * patches + patch].
 * Traveler states use travelers[(patch * compartments + compartment) *
 * travelers_per_patch + traveler]. Thus, traveler states belonging to the
 * same compartment are contiguous, which permits SIMD loads on CPUs and
 * coalesced loads on GPUs.
 */
struct StageAlignedProblem {
    int patches;
    int travelers_per_patch;
    int groups;
    std::vector<double> totals;
    std::vector<double> travelers;
    std::vector<double> contact_beta;
    std::vector<double> rate_exposed;
    std::vector<double> rate_infected;
    std::vector<double> stage_lambda;

    StageAlignedProblem(int num_patches, int num_travelers_per_patch, int num_groups)
        : patches(num_patches)
        , travelers_per_patch(num_travelers_per_patch)
        , groups(num_groups)
        , totals(static_cast<size_t>(4 * groups) * patches)
        , travelers(static_cast<size_t>(patches) * 4 * groups * travelers_per_patch)
        , contact_beta(static_cast<size_t>(groups) * groups)
        , rate_exposed(static_cast<size_t>(groups), 1.0 / 5.2)
        , rate_infected(static_cast<size_t>(groups), 1.0 / 6.0)
        , stage_lambda(static_cast<size_t>(4 * groups) * patches)
    {
        if (patches <= 0 || travelers_per_patch <= 0 ||
            std::find(age_group_counts.begin(), age_group_counts.end(), groups) == age_group_counts.end()) {
            throw std::invalid_argument("Unsupported stage-aligned benchmark shape.");
        }
        for (int i = 0; i < groups; ++i) {
            for (int j = 0; j < groups; ++j) {
                contact_beta[static_cast<size_t>(i) * groups + j] = 0.27 / (1.0 + std::abs(i - j));
            }
        }
        reset_state();
    }

    int compartments() const
    {
        return 4 * groups;
    }

    size_t edges() const
    {
        return static_cast<size_t>(patches) * travelers_per_patch;
    }

    void reset_state()
    {
        const double exposed    = 100.0 / groups;
        const double infected   = 100.0 / groups;
        const double recovered  = 100.0 / groups;
        const double population = 10000.0 / groups;

        for (int p = 0; p < patches; ++p) {
            const double patch_position = static_cast<double>(p + 1) / static_cast<double>(patches + 1);
            const double patch_scale    = 0.95 + 0.1 * patch_position;
            double traveler_weight_sum  = 0.0;
            for (int traveler = 0; traveler < travelers_per_patch; ++traveler) {
                traveler_weight_sum += 1.0 + 0.05 * ((17 * traveler + 13 * p) % 11);
            }
            for (int g = 0; g < groups; ++g) {
                const double group_scale    = 0.9 + 0.2 * static_cast<double>(g + 1) / static_cast<double>(groups + 1);
                const double patch_exposed  = exposed * (0.8 + 0.4 * patch_position) / group_scale;
                const double patch_infected = infected * (1.2 - 0.4 * patch_position) * group_scale;
                const std::array<double, 4> initial = {
                    (population - patch_exposed - patch_infected - recovered) * patch_scale,
                    patch_exposed * patch_scale, patch_infected * patch_scale, recovered * patch_scale};
                for (int c = 0; c < 4; ++c) {
                    const int compartment                                  = 4 * g + c;
                    totals[static_cast<size_t>(compartment) * patches + p] = initial[c];
                    auto first                                             = travelers.begin() +
                                 static_cast<std::ptrdiff_t>((static_cast<size_t>(p) * compartments() + compartment) *
                                                             travelers_per_patch);
                    for (int traveler = 0; traveler < travelers_per_patch; ++traveler) {
                        const double weight = 1.0 + 0.05 * ((17 * traveler + 13 * p) % 11);
                        first[traveler]     = initial[c] * total_traveler_fraction * weight / traveler_weight_sum;
                    }
                }
            }
        }
    }
};

template <int G>
inline void integrate_totals(StageAlignedProblem& problem, int patch, double dt)
{
    std::array<double, G> s0, e0, i0, inverse_population, population, s, e, i, result_s, result_e, result_i;
    const int patches = problem.patches;
    for (int g = 0; g < G; ++g) {
        s0[g] = s[g] = result_s[g] = problem.totals[static_cast<size_t>(4 * g) * patches + patch];
        e0[g] = e[g] = result_e[g] = problem.totals[static_cast<size_t>(4 * g + 1) * patches + patch];
        i0[g] = i[g] = result_i[g] = problem.totals[static_cast<size_t>(4 * g + 2) * patches + patch];
        population[g] = s0[g] + e0[g] + i0[g] + problem.totals[static_cast<size_t>(4 * g + 3) * patches + patch];
        inverse_population[g] = population[g] > 1e-12 ? 1.0 / population[g] : 0.0;
    }

    const std::array<double, 4> weights     = {dt / 6.0, dt / 3.0, dt / 3.0, dt / 6.0};
    const std::array<double, 3> stage_steps = {dt * 0.5, dt * 0.5, dt};
    for (int stage = 0; stage < 4; ++stage) {
        for (int target = 0; target < G; ++target) {
            double lambda = 0.0;
            for (int source = 0; source < G; ++source) {
                lambda += problem.contact_beta[static_cast<size_t>(target) * G + source] * i[source] *
                          inverse_population[source];
            }
            problem.stage_lambda[static_cast<size_t>(stage * G + target) * patches + patch] = lambda;
        }

        for (int g = 0; g < G; ++g) {
            const double lambda  = problem.stage_lambda[static_cast<size_t>(stage * G + g) * patches + patch];
            const double flow_se = lambda * s[g];
            const double flow_ei = problem.rate_exposed[g] * e[g];
            const double flow_ir = problem.rate_infected[g] * i[g];
            const double ds      = -flow_se;
            const double de      = flow_se - flow_ei;
            const double di      = flow_ei - flow_ir;
            result_s[g] += weights[stage] * ds;
            result_e[g] += weights[stage] * de;
            result_i[g] += weights[stage] * di;
            if (stage < 3) {
                s[g] = s0[g] + stage_steps[stage] * ds;
                e[g] = e0[g] + stage_steps[stage] * de;
                i[g] = i0[g] + stage_steps[stage] * di;
            }
        }
    }

    for (int g = 0; g < G; ++g) {
        problem.totals[static_cast<size_t>(4 * g) * patches + patch]     = result_s[g];
        problem.totals[static_cast<size_t>(4 * g + 1) * patches + patch] = result_e[g];
        problem.totals[static_cast<size_t>(4 * g + 2) * patches + patch] = result_i[g];
        problem.totals[static_cast<size_t>(4 * g + 3) * patches + patch] =
            population[g] - result_s[g] - result_e[g] - result_i[g];
    }
}

template <int G>
inline void integrate_traveler_chunk(StageAlignedProblem& problem, int patch, int group, int first, int last, double dt)
{
    const int stride          = problem.travelers_per_patch;
    const size_t patch_offset = static_cast<size_t>(patch) * 4 * G * stride;
    double* susceptible       = problem.travelers.data() + patch_offset + static_cast<size_t>(4 * group) * stride;
    double* exposed           = susceptible + stride;
    double* infected          = exposed + stride;
    double* recovered         = infected + stride;
    const double rate_e       = problem.rate_exposed[group];
    const double rate_i       = problem.rate_infected[group];
    const double lambda1      = problem.stage_lambda[static_cast<size_t>(group) * problem.patches + patch];
    const double lambda2      = problem.stage_lambda[static_cast<size_t>(G + group) * problem.patches + patch];
    const double lambda3      = problem.stage_lambda[static_cast<size_t>(2 * G + group) * problem.patches + patch];
    const double lambda4      = problem.stage_lambda[static_cast<size_t>(3 * G + group) * problem.patches + patch];
    const double half_dt      = 0.5 * dt;
    const double sixth_dt     = dt / 6.0;
    const double third_dt     = dt / 3.0;

#ifdef _OPENMP
#pragma omp simd
#endif
    for (int traveler = first; traveler < last; ++traveler) {
        const double s0         = susceptible[traveler];
        const double e0         = exposed[traveler];
        const double i0         = infected[traveler];
        const double population = s0 + e0 + i0 + recovered[traveler];

        double flow_se  = lambda1 * s0;
        double flow_ei  = rate_e * e0;
        double flow_ir  = rate_i * i0;
        double result_s = s0 - sixth_dt * flow_se;
        double result_e = e0 + sixth_dt * (flow_se - flow_ei);
        double result_i = i0 + sixth_dt * (flow_ei - flow_ir);
        double s        = s0 - half_dt * flow_se;
        double e        = e0 + half_dt * (flow_se - flow_ei);
        double i        = i0 + half_dt * (flow_ei - flow_ir);

        flow_se = lambda2 * s;
        flow_ei = rate_e * e;
        flow_ir = rate_i * i;
        result_s -= third_dt * flow_se;
        result_e += third_dt * (flow_se - flow_ei);
        result_i += third_dt * (flow_ei - flow_ir);
        s = s0 - half_dt * flow_se;
        e = e0 + half_dt * (flow_se - flow_ei);
        i = i0 + half_dt * (flow_ei - flow_ir);

        flow_se = lambda3 * s;
        flow_ei = rate_e * e;
        flow_ir = rate_i * i;
        result_s -= third_dt * flow_se;
        result_e += third_dt * (flow_se - flow_ei);
        result_i += third_dt * (flow_ei - flow_ir);
        s = s0 - dt * flow_se;
        e = e0 + dt * (flow_se - flow_ei);
        i = i0 + dt * (flow_ei - flow_ir);

        flow_se = lambda4 * s;
        flow_ei = rate_e * e;
        flow_ir = rate_i * i;
        result_s -= sixth_dt * flow_se;
        result_e += sixth_dt * (flow_se - flow_ei);
        result_i += sixth_dt * (flow_ei - flow_ir);
        susceptible[traveler] = result_s;
        exposed[traveler]     = result_e;
        infected[traveler]    = result_i;
        recovered[traveler]   = population - result_s - result_e - result_i;
    }
}

template <int G>
inline void advance_stage_aligned_impl(StageAlignedProblem& problem, double dt, int threads, int steps)
{
    const int chunks      = (problem.travelers_per_patch + traveler_chunk_size - 1) / traveler_chunk_size;
    const long long tasks = static_cast<long long>(problem.patches) * G * chunks;

#ifdef _OPENMP
    if (threads > 0) {
#pragma omp parallel num_threads(threads)
        {
#ifdef LIKWID_PERFMON
            LIKWID_MARKER_START("stage_aligned_openmp");
#endif
            for (int step = 0; step < steps; ++step) {
#pragma omp for schedule(static)
                for (int patch = 0; patch < problem.patches; ++patch) {
                    integrate_totals<G>(problem, patch, dt);
                }
#pragma omp for schedule(static)
                for (long long task = 0; task < tasks; ++task) {
                    const int chunk             = static_cast<int>(task % chunks);
                    const long long patch_group = task / chunks;
                    const int group             = static_cast<int>(patch_group % G);
                    const int patch             = static_cast<int>(patch_group / G);
                    const int first             = chunk * traveler_chunk_size;
                    const int last              = std::min(first + traveler_chunk_size, problem.travelers_per_patch);
                    integrate_traveler_chunk<G>(problem, patch, group, first, last, dt);
                }
            }
#ifdef LIKWID_PERFMON
            LIKWID_MARKER_STOP("stage_aligned_openmp");
#endif
        }
        return;
    }
#else
    if (threads > 0) {
        throw std::invalid_argument("OpenMP benchmark requested without OpenMP support.");
    }
#endif

#ifdef LIKWID_PERFMON
    LIKWID_MARKER_START("stage_aligned_serial");
#endif
    for (int step = 0; step < steps; ++step) {
        for (int patch = 0; patch < problem.patches; ++patch) {
            integrate_totals<G>(problem, patch, dt);
        }
        for (long long task = 0; task < tasks; ++task) {
            const int chunk             = static_cast<int>(task % chunks);
            const long long patch_group = task / chunks;
            const int group             = static_cast<int>(patch_group % G);
            const int patch             = static_cast<int>(patch_group / G);
            const int first             = chunk * traveler_chunk_size;
            const int last              = std::min(first + traveler_chunk_size, problem.travelers_per_patch);
            integrate_traveler_chunk<G>(problem, patch, group, first, last, dt);
        }
    }
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_STOP("stage_aligned_serial");
#endif
}

inline void advance_stage_aligned(StageAlignedProblem& problem, double dt = step_size, int threads = 0, int steps = 1)
{
    if (steps <= 0) {
        throw std::invalid_argument("The number of integration steps must be positive.");
    }
    switch (problem.groups) {
    case 1:
        advance_stage_aligned_impl<1>(problem, dt, threads, steps);
        break;
    case 3:
        advance_stage_aligned_impl<3>(problem, dt, threads, steps);
        break;
    case 6:
        advance_stage_aligned_impl<6>(problem, dt, threads, steps);
        break;
    case 8:
        advance_stage_aligned_impl<8>(problem, dt, threads, steps);
        break;
    default:
        throw std::invalid_argument("Unsupported number of age groups.");
    }
}

} // namespace mio::benchmark_mio

#endif // MIO_ODE_SEIR_BENCHMARK_STAGE_ALIGNED_H
