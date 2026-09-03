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

#ifndef MIO_ODE_SEIR_METAPOP_BENCHMARK_H
#define MIO_ODE_SEIR_METAPOP_BENCHMARK_H

#include "Eigen/Dense"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace mio::benchmark_mio
{

inline constexpr double step_size                       = 0.5;
inline constexpr int integration_steps                  = 20;
inline constexpr int maximum_groups                     = 8;
inline constexpr std::array<int, 4> age_group_counts    = {1, 3, 6, 8};
inline constexpr std::array<int, 7> patch_counts        = {17, 33, 65, 129, 257, 513, 1025};
inline constexpr double validation_relative_tolerance   = 1e-10;
inline constexpr double validation_absolute_denominator = 1.0;

/**
 * Synthetic fixed-H problem for the dense implicit-mobility benchmarks.
 * Dynamic S/E/I values use a compartment-grouped layout; R is reconstructed
 * from the conserved age-group population. Setup and transfers are not timed.
 */
struct ImplicitProblem {
    int patches;
    int groups;
    std::vector<double> initial_state;
    std::vector<double> state;
    std::vector<double> stage_state;
    std::vector<double> result;
    std::vector<double> population;
    std::vector<double> inverse_population;
    std::vector<double> inverse_present_population;
    std::vector<double> commuting_row_major;
    std::vector<double> commuting_transpose_row_major;
    std::vector<double> infection_coefficients;
    std::vector<double> rate_exposed;
    std::vector<double> rate_infected;
    std::vector<double> present_share;
    std::vector<double> mobile_share;

    ImplicitProblem(int num_patches, int num_groups)
        : patches(num_patches)
        , groups(num_groups)
        , initial_state(static_cast<size_t>(3 * groups) * patches)
        , state(initial_state.size())
        , stage_state(initial_state.size())
        , result(initial_state.size())
        , population(static_cast<size_t>(groups) * patches)
        , inverse_population(population.size())
        , inverse_present_population(population.size())
        , commuting_row_major(static_cast<size_t>(patches) * patches)
        , commuting_transpose_row_major(commuting_row_major.size())
        , infection_coefficients(static_cast<size_t>(groups) * groups)
        , rate_exposed(groups)
        , rate_infected(groups)
        , present_share(population.size())
        , mobile_share(population.size())
    {
        if (patches < 2 ||
            std::find(age_group_counts.begin(), age_group_counts.end(), groups) == age_group_counts.end()) {
            throw std::invalid_argument("Unsupported implicit benchmark shape.");
        }
        initialize_state();
        initialize_parameters();
        initialize_commuting_matrix();
        initialize_present_population();
        reset_state();
    }

    size_t state_index(int compartment, int group, int patch) const
    {
        return static_cast<size_t>(3 * group + compartment) * patches + patch;
    }

    size_t group_patch_index(int group, int patch) const
    {
        return static_cast<size_t>(group) * patches + patch;
    }

    void initialize_state()
    {
        for (int patch = 0; patch < patches; ++patch) {
            const double patch_position = static_cast<double>(patch + 1) / static_cast<double>(patches + 1);
            for (int group = 0; group < groups; ++group) {
                const double group_position = static_cast<double>(group + 1) / static_cast<double>(groups + 1);
                const double total     = 10000.0 * (0.9 + 0.2 * patch_position) * (0.9 + 0.2 * group_position) / groups;
                const double exposed   = total * (0.007 + 0.003 * patch_position) / (1.0 + 0.1 * group);
                const double infected  = total * (0.012 - 0.004 * patch_position) * (1.0 + 0.08 * group);
                const double recovered = total * (0.02 + 0.002 * group);
                const size_t gp        = group_patch_index(group, patch);
                population[gp]         = total;
                inverse_population[gp] = 1.0 / total;
                initial_state[state_index(0, group, patch)] = total - exposed - infected - recovered;
                initial_state[state_index(1, group, patch)] = exposed;
                initial_state[state_index(2, group, patch)] = infected;
            }
        }
    }

    void initialize_parameters()
    {
        for (int target = 0; target < groups; ++target) {
            const double transmission_probability = 0.035 + 0.002 * target;
            for (int source = 0; source < groups; ++source) {
                const double contacts = 9.0 * (1.0 + 0.02 * target + 0.01 * source) / (1.0 + std::abs(target - source));
                // Includes the model's factor 1/2 and rho_i * phi_{i,j}.
                infection_coefficients[static_cast<size_t>(target) * groups + source] =
                    0.5 * transmission_probability * contacts;
            }
            rate_exposed[target]  = 1.0 / (5.2 + 0.1 * target);
            rate_infected[target] = 1.0 / (6.0 + 0.1 * target);
        }
    }

    void initialize_commuting_matrix()
    {
        for (int origin = 0; origin < patches; ++origin) {
            double weight_sum = 0.0;
            for (int patch = 0; patch < patches; ++patch) {
                if (patch != origin) {
                    weight_sum += 1.0 + 0.02 * ((17 * origin + 13 * patch) % 23);
                }
            }
            const double mobile_fraction = 0.12 + 0.04 * static_cast<double>(origin % 7) / 6.0;
            for (int patch = 0; patch < patches; ++patch) {
                const double value =
                    patch == origin ? 1.0 - mobile_fraction
                                    : mobile_fraction * (1.0 + 0.02 * ((17 * origin + 13 * patch) % 23)) / weight_sum;
                commuting_row_major[static_cast<size_t>(origin) * patches + patch]           = value;
                commuting_transpose_row_major[static_cast<size_t>(patch) * patches + origin] = value;
            }
        }
    }

    void initialize_present_population()
    {
        for (int patch = 0; patch < patches; ++patch) {
            const double* commuting = commuting_transpose_row_major.data() + static_cast<size_t>(patch) * patches;
            for (int group = 0; group < groups; ++group) {
                double present_population = 0.0;
                for (int origin = 0; origin < patches; ++origin) {
                    present_population = std::fma(commuting[origin], population[group_patch_index(group, origin)],
                                                  present_population);
                }
                inverse_present_population[group_patch_index(group, patch)] = 1.0 / present_population;
            }
        }
    }

    void reset_state()
    {
        state = initial_state;
    }
};

template <int Stage>
inline void update_rk_value(double base, double derivative, double dt, double& accumulated, double& stage_value,
                            double& final_value)
{
    if constexpr (Stage == 0) {
        accumulated = std::fma(dt / 6.0, derivative, base);
        stage_value = std::fma(dt / 2.0, derivative, base);
    }
    else if constexpr (Stage == 1) {
        accumulated = std::fma(dt / 3.0, derivative, accumulated);
        stage_value = std::fma(dt / 2.0, derivative, base);
    }
    else if constexpr (Stage == 2) {
        accumulated = std::fma(dt / 3.0, derivative, accumulated);
        stage_value = std::fma(dt, derivative, base);
    }
    else {
        final_value = std::fma(dt / 6.0, derivative, accumulated);
    }
}

template <int G, int Stage>
inline void advance_cpu_patch(ImplicitProblem& problem, const std::vector<double>& current, int patch,
                              const std::array<double, G>& mobile_infectious_share, double dt)
{
    const int patches = problem.patches;
    std::array<double, G> home_infectious_share{};
    for (int source = 0; source < G; ++source) {
        home_infectious_share[source] = current[static_cast<size_t>(3 * source + 2) * patches + patch] *
                                        problem.inverse_population[static_cast<size_t>(source) * patches + patch];
    }
    for (int target = 0; target < G; ++target) {
        double pressure = 0.0;
        for (int source = 0; source < G; ++source) {
            pressure = std::fma(problem.infection_coefficients[static_cast<size_t>(target) * G + source],
                                home_infectious_share[source] + mobile_infectious_share[source],
                                pressure);
        }
        const size_t susceptible = static_cast<size_t>(3 * target) * patches + patch;
        const size_t exposed     = susceptible + patches;
        const size_t infected    = exposed + patches;
        const double flow_se     = pressure * current[susceptible];
        const double flow_ei     = problem.rate_exposed[target] * current[exposed];
        const double flow_ir     = problem.rate_infected[target] * current[infected];
        update_rk_value<Stage>(problem.state[susceptible], -flow_se, dt, problem.result[susceptible],
                               problem.stage_state[susceptible], problem.state[susceptible]);
        update_rk_value<Stage>(problem.state[exposed], flow_se - flow_ei, dt, problem.result[exposed],
                               problem.stage_state[exposed], problem.state[exposed]);
        update_rk_value<Stage>(problem.state[infected], flow_ei - flow_ir, dt, problem.result[infected],
                               problem.stage_state[infected], problem.state[infected]);
    }
}

template <int G, int Stage>
inline void advance_cpu_patch(ImplicitProblem& problem, const std::vector<double>& current, int patch, double dt)
{
    std::array<double, G> mobile_infectious_share{};
    for (int source = 0; source < G; ++source) {
        mobile_infectious_share[source] =
            problem.mobile_share[static_cast<size_t>(source) * problem.patches + patch];
    }
    advance_cpu_patch<G, Stage>(problem, current, patch, mobile_infectious_share, dt);
}

template <int G, int Stage>
void advance_cpu_stage(ImplicitProblem& problem, double dt)
{
    const int patches                  = problem.patches;
    const std::vector<double>& current = Stage == 0 ? problem.state : problem.stage_state;

    using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using ColMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    const Eigen::Map<const RowMatrix> commuting(problem.commuting_row_major.data(), patches, patches);
    const Eigen::Map<const ColMatrix, 0, Eigen::OuterStride<>> infected(
        current.data() + 2 * patches, patches, G, Eigen::OuterStride<>(3 * patches));
    Eigen::Map<ColMatrix> present_share(problem.present_share.data(), patches, G);
    const Eigen::Map<const ColMatrix> inverse_present(problem.inverse_present_population.data(), patches, G);
    Eigen::Map<ColMatrix> mobile_share(problem.mobile_share.data(), patches, G);
    present_share.noalias() = commuting.transpose() * infected;
    present_share.array() *= inverse_present.array();
    mobile_share.noalias() = commuting * present_share;

    for (int patch = 0; patch < patches; ++patch) {
        advance_cpu_patch<G, Stage>(problem, current, patch, dt);
    }
}

template <int G>
void advance_cpu(ImplicitProblem& problem, double dt = step_size, int steps = integration_steps)
{
    for (int step = 0; step < steps; ++step) {
        advance_cpu_stage<G, 0>(problem, dt);
        advance_cpu_stage<G, 1>(problem, dt);
        advance_cpu_stage<G, 2>(problem, dt);
        advance_cpu_stage<G, 3>(problem, dt);
    }
}

} // namespace mio::benchmark_mio

#endif // MIO_ODE_SEIR_METAPOP_BENCHMARK_H
