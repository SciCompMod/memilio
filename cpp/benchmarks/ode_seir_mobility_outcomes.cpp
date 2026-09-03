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

#include "models/ode_seir_metapop/model.h"
#include "models/ode_seir_metapop/parameters.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr int num_compartments = 4;
constexpr int susceptible      = 0;
constexpr int exposed          = 1;
constexpr int infected         = 2;
constexpr int recovered        = 3;
constexpr double time_tolerance = 1e-11;

struct ExperimentConfig {
    int rows                 = 1;
    int columns              = 4;
    int patches              = rows * columns;
    int groups               = 1;
    double theta             = 0.5;
    double step_size         = 1.0 / 64.0;
    double output_interval   = 0.25;
    double end_time          = 160.0;
    std::vector<double> mobility;
    std::vector<double> contact_rates;
    std::vector<double> infection_risk;
    std::vector<double> time_exposed;
    std::vector<double> time_infected;
    std::vector<double> initial_state;
    std::vector<double> resident_population;
};

struct Trajectory {
    std::string model;
    double cycle_days = 0.0;
    std::vector<double> times;
    std::vector<std::vector<double>> states;
    double max_population_error       = 0.0;
    double max_stage_consistency_error = 0.0;
};

struct ComparisonMetrics {
    double peak_infected_fraction    = 0.0;
    double peak_time_days            = 0.0;
    double attack_rate               = 0.0;
    double max_normalized_patch_l1   = 0.0;
    double time_of_max_difference    = 0.0;
    double relative_infectious_l2    = 0.0;
};

struct ExplicitWorkspace {
    std::vector<double> totals;
    std::vector<double> k1;
    std::vector<double> k2;
    std::vector<double> k3;
    std::vector<double> k4;
    std::vector<double> temporary;
    std::vector<double> updated_totals;
    std::vector<double> checked_totals;
    std::vector<double> stage_lambda;

    explicit ExplicitWorkspace(const ExperimentConfig& config)
    {
        const auto aggregate_size = static_cast<size_t>(config.patches * config.groups * num_compartments);
        totals.resize(aggregate_size);
        k1.resize(aggregate_size);
        k2.resize(aggregate_size);
        k3.resize(aggregate_size);
        k4.resize(aggregate_size);
        temporary.resize(aggregate_size);
        updated_totals.resize(aggregate_size);
        checked_totals.resize(aggregate_size);
        stage_lambda.resize(static_cast<size_t>(4 * config.patches * config.groups));
    }
};

size_t resident_index(const ExperimentConfig& config, int patch, int group, int compartment)
{
    return static_cast<size_t>(((patch * config.groups + group) * num_compartments) + compartment);
}

size_t population_index(const ExperimentConfig& config, int patch, int group)
{
    return static_cast<size_t>(patch * config.groups + group);
}

size_t mobility_index(const ExperimentConfig& config, int origin, int location)
{
    return static_cast<size_t>(origin * config.patches + location);
}

size_t located_index(const ExperimentConfig& config, int origin, int location, int group, int compartment)
{
    return static_cast<size_t>((((origin * config.patches + location) * config.groups + group) * num_compartments) +
                               compartment);
}

size_t lambda_index(const ExperimentConfig& config, int stage, int patch, int group)
{
    return static_cast<size_t>((stage * config.patches + patch) * config.groups + group);
}

bool equal_time(double left, double right)
{
    return std::abs(left - right) <= time_tolerance * std::max({1.0, std::abs(left), std::abs(right)});
}

ExperimentConfig make_config()
{
    ExperimentConfig config;
    config.mobility.resize(static_cast<size_t>(config.patches * config.patches), 0.0);
    config.contact_rates = {6.0};
    config.infection_risk = {0.05};
    config.time_exposed   = {4.0};
    config.time_infected  = {6.0};

    config.initial_state.resize(
        static_cast<size_t>(config.patches * config.groups * num_compartments), 0.0);
    config.resident_population.resize(static_cast<size_t>(config.patches * config.groups), 0.0);

    for (int patch = 0; patch < config.patches; ++patch) {
        constexpr double population = 10000.0;
        const double initial_exposed  = patch == 0 ? 80.0 : 0.0;
        const double initial_infected = patch == 0 ? 20.0 : 0.0;
        config.resident_population[population_index(config, patch, 0)] = population;
        config.initial_state[resident_index(config, patch, 0, susceptible)] =
            population - initial_exposed - initial_infected;
        config.initial_state[resident_index(config, patch, 0, exposed)]  = initial_exposed;
        config.initial_state[resident_index(config, patch, 0, infected)] = initial_infected;
    }

    for (int patch = 0; patch < config.patches; ++patch) {
        config.mobility[mobility_index(config, patch, patch)] =
            patch == 0 || patch == config.patches - 1 ? 0.9 : 0.8;
        if (patch > 0) {
            config.mobility[mobility_index(config, patch, patch - 1)] = 0.1;
        }
        if (patch + 1 < config.patches) {
            config.mobility[mobility_index(config, patch, patch + 1)] = 0.1;
        }
    }

    for (int origin = 0; origin < config.patches; ++origin) {
        double row_sum = 0.0;
        for (int location = 0; location < config.patches; ++location) {
            const double value = config.mobility[mobility_index(config, origin, location)];
            if (value < 0.0 || value > 1.0) {
                throw std::runtime_error("Mobility fractions must lie in [0, 1].");
            }
            row_sum += value;
        }
        if (std::abs(row_sum - 1.0) > 1e-12) {
            throw std::runtime_error("Each mobility-matrix row must sum to one.");
        }
    }
    return config;
}

using ImplicitModel = mio::oseirmetapop::Model<double>;

void configure_implicit_model(ImplicitModel& model, const ExperimentConfig& config)
{
    using mio::AgeGroup;
    using mio::regions::Region;
    using namespace mio::oseirmetapop;

    for (int patch = 0; patch < config.patches; ++patch) {
        for (int group = 0; group < config.groups; ++group) {
            model.populations[{Region(patch), AgeGroup(group), InfectionState::Susceptible}] =
                config.initial_state[resident_index(config, patch, group, susceptible)];
            model.populations[{Region(patch), AgeGroup(group), InfectionState::Exposed}] =
                config.initial_state[resident_index(config, patch, group, exposed)];
            model.populations[{Region(patch), AgeGroup(group), InfectionState::Infected}] =
                config.initial_state[resident_index(config, patch, group, infected)];
            model.populations[{Region(patch), AgeGroup(group), InfectionState::Recovered}] =
                config.initial_state[resident_index(config, patch, group, recovered)];
            model.parameters.get<TransmissionProbabilityOnContact<>>()[{Region(patch), AgeGroup(group)}] =
                config.infection_risk[static_cast<size_t>(group)];
            model.parameters.get<TimeExposed<>>()[{Region(patch), AgeGroup(group)}] =
                config.time_exposed[static_cast<size_t>(group)];
            model.parameters.get<TimeInfected<>>()[{Region(patch), AgeGroup(group)}] =
                config.time_infected[static_cast<size_t>(group)];
        }
    }

    auto& contact_matrix = model.parameters.get<ContactPatterns<>>().get_cont_freq_mat()[0].get_baseline();
    for (int target = 0; target < config.groups; ++target) {
        for (int source = 0; source < config.groups; ++source) {
            contact_matrix(target, source) =
                config.contact_rates[static_cast<size_t>(target * config.groups + source)];
        }
    }

    Eigen::MatrixXd mobility(config.patches, config.patches);
    for (int origin = 0; origin < config.patches; ++origin) {
        for (int location = 0; location < config.patches; ++location) {
            mobility(origin, location) = config.mobility[mobility_index(config, origin, location)];
        }
    }
    model.set_commuting_strengths(mobility);
    if (model.check_constraints()) {
        throw std::runtime_error("The implicit model configuration violates a model constraint.");
    }
}

void rk4_implicit_step(const ImplicitModel& model, Eigen::VectorXd& state, double time, double step)
{
    Eigen::VectorXd k1(state.size());
    Eigen::VectorXd k2(state.size());
    Eigen::VectorXd k3(state.size());
    Eigen::VectorXd k4(state.size());
    Eigen::VectorXd temporary(state.size());

    model.eval_right_hand_side(state, state, time, k1);
    temporary = state + 0.5 * step * k1;
    model.eval_right_hand_side(temporary, temporary, time + 0.5 * step, k2);
    temporary = state + 0.5 * step * k2;
    model.eval_right_hand_side(temporary, temporary, time + 0.5 * step, k3);
    temporary = state + step * k3;
    model.eval_right_hand_side(temporary, temporary, time + step, k4);
    state += (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

std::vector<double> extract_implicit_state(const ImplicitModel& model, const ExperimentConfig& config,
                                           const Eigen::VectorXd& state)
{
    using mio::AgeGroup;
    using mio::regions::Region;
    using mio::oseirmetapop::InfectionState;
    const std::array<InfectionState, num_compartments> infection_states = {
        InfectionState::Susceptible, InfectionState::Exposed, InfectionState::Infected, InfectionState::Recovered};
    std::vector<double> result(config.initial_state.size());
    for (int patch = 0; patch < config.patches; ++patch) {
        for (int group = 0; group < config.groups; ++group) {
            for (int compartment = 0; compartment < num_compartments; ++compartment) {
                const auto model_index = model.populations.get_flat_index(
                    {Region(patch), AgeGroup(group), infection_states[static_cast<size_t>(compartment)]});
                result[resident_index(config, patch, group, compartment)] = state[model_index];
            }
        }
    }
    return result;
}

double population_error(const ExperimentConfig& config, const std::vector<double>& resident_state)
{
    double maximum = 0.0;
    for (int patch = 0; patch < config.patches; ++patch) {
        for (int group = 0; group < config.groups; ++group) {
            double current = 0.0;
            for (int compartment = 0; compartment < num_compartments; ++compartment) {
                const double value = resident_state[resident_index(config, patch, group, compartment)];
                if (!std::isfinite(value) || value < -1e-8) {
                    throw std::runtime_error("A simulated compartment value is invalid or negative.");
                }
                current += value;
            }
            const double reference = config.resident_population[population_index(config, patch, group)];
            maximum = std::max(maximum, std::abs(current - reference) / std::max(1.0, reference));
        }
    }
    return maximum;
}

Trajectory simulate_implicit(const ExperimentConfig& config)
{
    ImplicitModel model(config.patches, config.groups);
    configure_implicit_model(model, config);
    Eigen::VectorXd state = model.get_initial_values();
    Trajectory trajectory{"implicit", 0.0, {}, {}, 0.0, 0.0};
    trajectory.times.push_back(0.0);
    trajectory.states.push_back(extract_implicit_state(model, config, state));

    double time        = 0.0;
    double next_output = config.output_interval;
    while (time < config.end_time - time_tolerance) {
        const double target = std::min(next_output, config.end_time);
        while (time < target - time_tolerance) {
            const double step = std::min(config.step_size, target - time);
            rk4_implicit_step(model, state, time, step);
            time += step;
        }
        time = target;
        auto snapshot = extract_implicit_state(model, config, state);
        trajectory.max_population_error =
            std::max(trajectory.max_population_error, population_error(config, snapshot));
        trajectory.times.push_back(time);
        trajectory.states.push_back(std::move(snapshot));
        next_output += config.output_interval;
    }
    return trajectory;
}

std::vector<double> residents_from_locations(const ExperimentConfig& config, const std::vector<double>& located)
{
    std::vector<double> residents(config.initial_state.size(), 0.0);
    for (int origin = 0; origin < config.patches; ++origin) {
        for (int location = 0; location < config.patches; ++location) {
            for (int group = 0; group < config.groups; ++group) {
                for (int compartment = 0; compartment < num_compartments; ++compartment) {
                    residents[resident_index(config, origin, group, compartment)] +=
                        located[located_index(config, origin, location, group, compartment)];
                }
            }
        }
    }
    return residents;
}

void set_home_locations(const ExperimentConfig& config, const std::vector<double>& residents,
                        std::vector<double>& located)
{
    std::fill(located.begin(), located.end(), 0.0);
    for (int origin = 0; origin < config.patches; ++origin) {
        for (int group = 0; group < config.groups; ++group) {
            for (int compartment = 0; compartment < num_compartments; ++compartment) {
                located[located_index(config, origin, origin, group, compartment)] =
                    residents[resident_index(config, origin, group, compartment)];
            }
        }
    }
}

void depart(const ExperimentConfig& config, std::vector<double>& located)
{
    const auto residents = residents_from_locations(config, located);
    std::fill(located.begin(), located.end(), 0.0);
    for (int origin = 0; origin < config.patches; ++origin) {
        for (int location = 0; location < config.patches; ++location) {
            const double fraction = config.mobility[mobility_index(config, origin, location)];
            for (int group = 0; group < config.groups; ++group) {
                for (int compartment = 0; compartment < num_compartments; ++compartment) {
                    located[located_index(config, origin, location, group, compartment)] =
                        fraction * residents[resident_index(config, origin, group, compartment)];
                }
            }
        }
    }
}

void return_home(const ExperimentConfig& config, std::vector<double>& located)
{
    const auto residents = residents_from_locations(config, located);
    set_home_locations(config, residents, located);
    for (int origin = 0; origin < config.patches; ++origin) {
        for (int location = 0; location < config.patches; ++location) {
            if (origin == location) {
                continue;
            }
            for (int group = 0; group < config.groups; ++group) {
                for (int compartment = 0; compartment < num_compartments; ++compartment) {
                    if (std::abs(located[located_index(config, origin, location, group, compartment)]) > 1e-12) {
                        throw std::runtime_error("An off-diagonal traveler state remained after return.");
                    }
                }
            }
        }
    }
}

void totals_from_locations(const ExperimentConfig& config, const std::vector<double>& located,
                           std::vector<double>& totals)
{
    std::fill(totals.begin(), totals.end(), 0.0);
    for (int origin = 0; origin < config.patches; ++origin) {
        for (int location = 0; location < config.patches; ++location) {
            for (int group = 0; group < config.groups; ++group) {
                for (int compartment = 0; compartment < num_compartments; ++compartment) {
                    totals[resident_index(config, location, group, compartment)] +=
                        located[located_index(config, origin, location, group, compartment)];
                }
            }
        }
    }
}

void aggregate_rhs(const ExperimentConfig& config, const std::vector<double>& state, std::vector<double>& derivative,
                   double* stage_lambda)
{
    std::fill(derivative.begin(), derivative.end(), 0.0);
    for (int patch = 0; patch < config.patches; ++patch) {
        for (int target = 0; target < config.groups; ++target) {
            double force_of_infection = 0.0;
            for (int source = 0; source < config.groups; ++source) {
                double population = 0.0;
                for (int compartment = 0; compartment < num_compartments; ++compartment) {
                    population += state[resident_index(config, patch, source, compartment)];
                }
                const double infectious_share =
                    population > 1e-14 ? state[resident_index(config, patch, source, infected)] / population : 0.0;
                const double beta = config.infection_risk[static_cast<size_t>(target)] *
                                    config.contact_rates[static_cast<size_t>(target * config.groups + source)];
                force_of_infection += beta * infectious_share;
            }
            stage_lambda[static_cast<size_t>(patch * config.groups + target)] = force_of_infection;
            const double s = state[resident_index(config, patch, target, susceptible)];
            const double e = state[resident_index(config, patch, target, exposed)];
            const double i = state[resident_index(config, patch, target, infected)];
            const double flow_se = force_of_infection * s;
            const double flow_ei = e / config.time_exposed[static_cast<size_t>(target)];
            const double flow_ir = i / config.time_infected[static_cast<size_t>(target)];
            derivative[resident_index(config, patch, target, susceptible)] = -flow_se;
            derivative[resident_index(config, patch, target, exposed)]     = flow_se - flow_ei;
            derivative[resident_index(config, patch, target, infected)]    = flow_ei - flow_ir;
            derivative[resident_index(config, patch, target, recovered)]   = flow_ir;
        }
    }
}

std::array<double, num_compartments> traveler_rhs(const ExperimentConfig& config, int group, double lambda,
                                                   const std::array<double, num_compartments>& state)
{
    const double flow_se = lambda * state[susceptible];
    const double flow_ei = state[exposed] / config.time_exposed[static_cast<size_t>(group)];
    const double flow_ir = state[infected] / config.time_infected[static_cast<size_t>(group)];
    return {-flow_se, flow_se - flow_ei, flow_ei - flow_ir, flow_ir};
}

double stage_aligned_step(const ExperimentConfig& config, std::vector<double>& located, double step,
                          ExplicitWorkspace& workspace)
{
    totals_from_locations(config, located, workspace.totals);
    const size_t aggregate_size = workspace.totals.size();

    aggregate_rhs(config, workspace.totals, workspace.k1, workspace.stage_lambda.data());
    for (size_t index = 0; index < aggregate_size; ++index) {
        workspace.temporary[index] = workspace.totals[index] + 0.5 * step * workspace.k1[index];
    }
    aggregate_rhs(config, workspace.temporary, workspace.k2,
                  workspace.stage_lambda.data() + static_cast<size_t>(config.patches * config.groups));
    for (size_t index = 0; index < aggregate_size; ++index) {
        workspace.temporary[index] = workspace.totals[index] + 0.5 * step * workspace.k2[index];
    }
    aggregate_rhs(config, workspace.temporary, workspace.k3,
                  workspace.stage_lambda.data() + static_cast<size_t>(2 * config.patches * config.groups));
    for (size_t index = 0; index < aggregate_size; ++index) {
        workspace.temporary[index] = workspace.totals[index] + step * workspace.k3[index];
    }
    aggregate_rhs(config, workspace.temporary, workspace.k4,
                  workspace.stage_lambda.data() + static_cast<size_t>(3 * config.patches * config.groups));
    for (size_t index = 0; index < aggregate_size; ++index) {
        workspace.updated_totals[index] = workspace.totals[index] +
                                          (step / 6.0) * (workspace.k1[index] + 2.0 * workspace.k2[index] +
                                                          2.0 * workspace.k3[index] + workspace.k4[index]);
    }

    for (int origin = 0; origin < config.patches; ++origin) {
        for (int location = 0; location < config.patches; ++location) {
            for (int group = 0; group < config.groups; ++group) {
                std::array<double, num_compartments> initial{};
                for (int compartment = 0; compartment < num_compartments; ++compartment) {
                    initial[static_cast<size_t>(compartment)] =
                        located[located_index(config, origin, location, group, compartment)];
                }
                const double lambda1 = workspace.stage_lambda[lambda_index(config, 0, location, group)];
                const double lambda2 = workspace.stage_lambda[lambda_index(config, 1, location, group)];
                const double lambda3 = workspace.stage_lambda[lambda_index(config, 2, location, group)];
                const double lambda4 = workspace.stage_lambda[lambda_index(config, 3, location, group)];
                const auto first     = traveler_rhs(config, group, lambda1, initial);
                std::array<double, num_compartments> temporary{};
                for (int compartment = 0; compartment < num_compartments; ++compartment) {
                    temporary[static_cast<size_t>(compartment)] =
                        initial[static_cast<size_t>(compartment)] +
                        0.5 * step * first[static_cast<size_t>(compartment)];
                }
                const auto second = traveler_rhs(config, group, lambda2, temporary);
                for (int compartment = 0; compartment < num_compartments; ++compartment) {
                    temporary[static_cast<size_t>(compartment)] =
                        initial[static_cast<size_t>(compartment)] +
                        0.5 * step * second[static_cast<size_t>(compartment)];
                }
                const auto third = traveler_rhs(config, group, lambda3, temporary);
                for (int compartment = 0; compartment < num_compartments; ++compartment) {
                    temporary[static_cast<size_t>(compartment)] =
                        initial[static_cast<size_t>(compartment)] + step * third[static_cast<size_t>(compartment)];
                }
                const auto fourth = traveler_rhs(config, group, lambda4, temporary);
                for (int compartment = 0; compartment < num_compartments; ++compartment) {
                    located[located_index(config, origin, location, group, compartment)] =
                        initial[static_cast<size_t>(compartment)] +
                        (step / 6.0) * (first[static_cast<size_t>(compartment)] +
                                        2.0 * second[static_cast<size_t>(compartment)] +
                                        2.0 * third[static_cast<size_t>(compartment)] +
                                        fourth[static_cast<size_t>(compartment)]);
                }
            }
        }
    }

    totals_from_locations(config, located, workspace.checked_totals);
    double consistency_error = 0.0;
    for (size_t index = 0; index < aggregate_size; ++index) {
        consistency_error =
            std::max(consistency_error, std::abs(workspace.checked_totals[index] - workspace.updated_totals[index]) /
                                            std::max(1.0, std::abs(workspace.updated_totals[index])));
    }
    return consistency_error;
}

Trajectory simulate_explicit(const ExperimentConfig& config, double cycle_days)
{
    if (cycle_days <= 0.0) {
        throw std::invalid_argument("The explicit mobility-cycle length must be positive.");
    }
    std::vector<double> located(
        static_cast<size_t>(config.patches * config.patches * config.groups * num_compartments), 0.0);
    set_home_locations(config, config.initial_state, located);
    ExplicitWorkspace workspace(config);
    Trajectory trajectory{"explicit", cycle_days, {}, {}, 0.0, 0.0};
    trajectory.times.push_back(0.0);
    trajectory.states.push_back(config.initial_state);

    double time             = 0.0;
    double next_output      = config.output_interval;
    double next_departure   = config.theta * cycle_days;
    double next_return      = cycle_days;
    bool away               = false;

    while (time < config.end_time - time_tolerance) {
        const double target = std::min({next_output, next_departure, next_return, config.end_time});
        while (time < target - time_tolerance) {
            const double step = std::min(config.step_size, target - time);
            trajectory.max_stage_consistency_error =
                std::max(trajectory.max_stage_consistency_error,
                         stage_aligned_step(config, located, step, workspace));
            time += step;
        }
        time = target;

        if (equal_time(time, next_departure)) {
            if (away) {
                throw std::runtime_error("A departure event occurred during an away phase.");
            }
            depart(config, located);
            away = true;
            next_departure += cycle_days;
        }
        if (equal_time(time, next_return)) {
            if (!away) {
                throw std::runtime_error("A return event occurred during a home phase.");
            }
            return_home(config, located);
            away = false;
            next_return += cycle_days;
        }
        if (equal_time(time, next_output) || equal_time(time, config.end_time)) {
            auto snapshot = residents_from_locations(config, located);
            trajectory.max_population_error =
                std::max(trajectory.max_population_error, population_error(config, snapshot));
            trajectory.times.push_back(time);
            trajectory.states.push_back(std::move(snapshot));
            next_output += config.output_interval;
        }
    }
    return trajectory;
}

double total_population(const ExperimentConfig& config)
{
    return std::accumulate(config.resident_population.begin(), config.resident_population.end(), 0.0);
}

double compartment_total(const ExperimentConfig& config, const std::vector<double>& state, int compartment)
{
    double result = 0.0;
    for (int patch = 0; patch < config.patches; ++patch) {
        for (int group = 0; group < config.groups; ++group) {
            result += state[resident_index(config, patch, group, compartment)];
        }
    }
    return result;
}

ComparisonMetrics compare_with_implicit(const ExperimentConfig& config, const Trajectory& implicit,
                                        const Trajectory& trajectory)
{
    if (implicit.times.size() != trajectory.times.size()) {
        throw std::runtime_error("Compared trajectories use different output grids.");
    }
    const double population = total_population(config);
    double implicit_peak_infected = 0.0;
    for (const auto& state : implicit.states) {
        implicit_peak_infected = std::max(implicit_peak_infected, compartment_total(config, state, infected));
    }

    ComparisonMetrics metrics;
    double l2_numerator   = 0.0;
    double l2_denominator = 0.0;
    double peak_infected  = -1.0;
    for (size_t time_index = 0; time_index < trajectory.times.size(); ++time_index) {
        if (!equal_time(implicit.times[time_index], trajectory.times[time_index])) {
            throw std::runtime_error("Compared trajectories use different output times.");
        }
        const auto& reference = implicit.states[time_index];
        const auto& state     = trajectory.states[time_index];
        const double current_infected = compartment_total(config, state, infected);
        if (current_infected > peak_infected) {
            peak_infected         = current_infected;
            metrics.peak_time_days = trajectory.times[time_index];
        }

        double patch_l1 = 0.0;
        for (int patch = 0; patch < config.patches; ++patch) {
            double reference_patch = 0.0;
            double current_patch   = 0.0;
            for (int group = 0; group < config.groups; ++group) {
                const double reference_value = reference[resident_index(config, patch, group, infected)];
                const double current_value   = state[resident_index(config, patch, group, infected)];
                reference_patch += reference_value;
                current_patch += current_value;
                const double difference = current_value - reference_value;
                l2_numerator += difference * difference;
                l2_denominator += reference_value * reference_value;
            }
            patch_l1 += std::abs(current_patch - reference_patch);
        }
        const double normalized_patch_l1 = patch_l1 / std::max(1.0, implicit_peak_infected);
        if (normalized_patch_l1 > metrics.max_normalized_patch_l1) {
            metrics.max_normalized_patch_l1 = normalized_patch_l1;
            metrics.time_of_max_difference  = trajectory.times[time_index];
        }
    }

    metrics.peak_infected_fraction = peak_infected / population;
    metrics.attack_rate = (compartment_total(config, trajectory.states.front(), susceptible) -
                           compartment_total(config, trajectory.states.back(), susceptible)) /
                          population;
    metrics.relative_infectious_l2 = std::sqrt(l2_numerator / std::max(1.0, l2_denominator));
    return metrics;
}

double max_l1_fraction(const ExperimentConfig& config, const Trajectory& first, const Trajectory& second)
{
    if (first.times.size() != second.times.size()) {
        throw std::runtime_error("Validation trajectories use different output grids.");
    }
    double maximum = 0.0;
    const double population = total_population(config);
    for (size_t time_index = 0; time_index < first.times.size(); ++time_index) {
        double difference = 0.0;
        for (size_t state_index = 0; state_index < first.states[time_index].size(); ++state_index) {
            difference += std::abs(first.states[time_index][state_index] - second.states[time_index][state_index]);
        }
        maximum = std::max(maximum, difference / population);
    }
    return maximum;
}

void write_timeseries(const std::filesystem::path& filename, const ExperimentConfig& config,
                      const std::vector<Trajectory>& trajectories)
{
    std::ofstream output(filename);
    if (!output) {
        throw std::runtime_error("Could not open trajectory output file: " + filename.string());
    }
    output << "model,cycle_days,time_days,patch,row,column,age_group,population,S,E,I,R\n";
    output << std::setprecision(17);
    for (const auto& trajectory : trajectories) {
        for (size_t time_index = 0; time_index < trajectory.times.size(); ++time_index) {
            for (int patch = 0; patch < config.patches; ++patch) {
                for (int group = 0; group < config.groups; ++group) {
                    output << trajectory.model << ',' << trajectory.cycle_days << ',' << trajectory.times[time_index]
                           << ',' << patch << ',' << patch / config.columns << ',' << patch % config.columns << ','
                           << group << ',' << config.resident_population[population_index(config, patch, group)];
                    for (int compartment = 0; compartment < num_compartments; ++compartment) {
                        output << ',' << trajectory.states[time_index]
                                                  [resident_index(config, patch, group, compartment)];
                    }
                    output << '\n';
                }
            }
        }
    }
}

void write_summary(const std::filesystem::path& filename, const ExperimentConfig& config,
                   const Trajectory& implicit, const std::vector<Trajectory>& explicit_trajectories,
                   double implicit_refinement, double explicit_refinement)
{
    std::ofstream output(filename);
    if (!output) {
        throw std::runtime_error("Could not open summary output file: " + filename.string());
    }
    output << "model,cycle_days,peak_infected_fraction,peak_time_days,attack_rate,max_normalized_patch_l1,"
              "time_of_max_difference,relative_infectious_l2,max_population_error,max_stage_consistency_error,"
              "step_refinement_max_l1_fraction\n";
    output << std::setprecision(17);
    const auto implicit_metrics = compare_with_implicit(config, implicit, implicit);
    output << "implicit,0," << implicit_metrics.peak_infected_fraction << ',' << implicit_metrics.peak_time_days << ','
           << implicit_metrics.attack_rate << ",0,0,0," << implicit.max_population_error << ",0,"
           << implicit_refinement << '\n';
    for (const auto& trajectory : explicit_trajectories) {
        const auto metrics = compare_with_implicit(config, implicit, trajectory);
        const double refinement = equal_time(trajectory.cycle_days, 1.0)
                                      ? explicit_refinement
                                      : std::numeric_limits<double>::quiet_NaN();
        output << "explicit," << trajectory.cycle_days << ',' << metrics.peak_infected_fraction << ','
               << metrics.peak_time_days << ',' << metrics.attack_rate << ',' << metrics.max_normalized_patch_l1 << ','
               << metrics.time_of_max_difference << ',' << metrics.relative_infectious_l2 << ','
               << trajectory.max_population_error << ',' << trajectory.max_stage_consistency_error << ','
               << refinement << '\n';
    }
}

void print_summary(const ExperimentConfig& config, const Trajectory& implicit,
                   const std::vector<Trajectory>& explicit_trajectories)
{
    const auto implicit_metrics = compare_with_implicit(config, implicit, implicit);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Implicit: peak=" << 100.0 * implicit_metrics.peak_infected_fraction << "% at day "
              << implicit_metrics.peak_time_days << ", attack rate=" << 100.0 * implicit_metrics.attack_rate << "%\n";
    for (const auto& trajectory : explicit_trajectories) {
        const auto metrics = compare_with_implicit(config, implicit, trajectory);
        std::cout << "Explicit delta=" << trajectory.cycle_days << " d: peak="
                  << 100.0 * metrics.peak_infected_fraction << "% at day " << metrics.peak_time_days
                  << ", attack rate=" << 100.0 * metrics.attack_rate
                  << "%, max spatial L1=" << 100.0 * metrics.max_normalized_patch_l1 << "% of implicit peak\n";
    }
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::filesystem::path output_directory = argc > 1 ? argv[1] : "hpc_paper/results";
        std::filesystem::create_directories(output_directory);

        const ExperimentConfig config = make_config();
        const Trajectory implicit      = simulate_implicit(config);
        const std::array<double, 5> cycle_lengths = {0.25, 0.5, 1.0, 2.0, 4.0};
        std::vector<Trajectory> explicit_trajectories;
        explicit_trajectories.reserve(cycle_lengths.size());
        for (double cycle_days : cycle_lengths) {
            explicit_trajectories.push_back(simulate_explicit(config, cycle_days));
        }

        ExperimentConfig no_mobility = config;
        no_mobility.end_time          = 4.0;
        std::fill(no_mobility.mobility.begin(), no_mobility.mobility.end(), 0.0);
        for (int patch = 0; patch < no_mobility.patches; ++patch) {
            no_mobility.mobility[mobility_index(no_mobility, patch, patch)] = 1.0;
        }
        const auto no_mobility_implicit = simulate_implicit(no_mobility);
        const auto no_mobility_explicit = simulate_explicit(no_mobility, 1.0);
        const double no_mobility_difference =
            max_l1_fraction(no_mobility, no_mobility_implicit, no_mobility_explicit);
        if (no_mobility_difference > 1e-10) {
            throw std::runtime_error("Implicit and explicit models disagree in the no-mobility limit.");
        }

        ExperimentConfig refined = config;
        refined.step_size *= 0.5;
        const auto refined_implicit = simulate_implicit(refined);
        const auto refined_explicit = simulate_explicit(refined, 1.0);
        const double implicit_refinement = max_l1_fraction(config, implicit, refined_implicit);
        const auto daily = std::find_if(explicit_trajectories.begin(), explicit_trajectories.end(),
                                        [](const auto& trajectory) { return equal_time(trajectory.cycle_days, 1.0); });
        if (daily == explicit_trajectories.end()) {
            throw std::runtime_error("The daily explicit trajectory is missing.");
        }
        const double explicit_refinement = max_l1_fraction(config, *daily, refined_explicit);

        std::vector<Trajectory> all_trajectories;
        all_trajectories.reserve(1 + explicit_trajectories.size());
        all_trajectories.push_back(implicit);
        all_trajectories.insert(all_trajectories.end(), explicit_trajectories.begin(), explicit_trajectories.end());
        write_timeseries(output_directory / "mobility_outcomes_timeseries.csv", config, all_trajectories);
        write_summary(output_directory / "mobility_outcomes_summary.csv", config, implicit, explicit_trajectories,
                      implicit_refinement, explicit_refinement);

        print_summary(config, implicit, explicit_trajectories);
        std::cout << std::scientific << std::setprecision(3)
                  << "No-mobility validation difference: " << no_mobility_difference << '\n'
                  << "Implicit h/2 refinement difference: " << implicit_refinement << '\n'
                  << "Explicit h/2 refinement difference: " << explicit_refinement << '\n'
                  << "Wrote results to " << std::filesystem::absolute(output_directory) << '\n';
    }
    catch (const std::exception& exception) {
        std::cerr << "mobility outcome experiment failed: " << exception.what() << '\n';
        return 1;
    }
    return 0;
}
