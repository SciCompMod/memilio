/* Copyright (C) 2026 MEmilio. Licensed under the Apache License, Version 2.0. */
#ifndef MIO_ODE_SEIR_RUNTIME_SCENARIO_H
#define MIO_ODE_SEIR_RUNTIME_SCENARIO_H

#include "benchmark/benchmark.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mio::runtime_scenario
{
inline constexpr int version                = 1;
// Separate the registered grid from the unchanged physical scenario. Version 1
// contained only the original 28 shapes; version 2 adds four large G=6 shapes.
inline constexpr int shape_set_version      = 2;
inline constexpr double accuracy_tolerance  = 1e-6; // fraction of initial resident population
inline constexpr std::array<int, 7> patches = {16, 32, 64, 128, 256, 512, 1024};
inline constexpr std::array<int, 4> large_patches = {2048, 4096, 8192, 16384};
inline constexpr std::array<int, 4> groups  = {1, 3, 6, 8};
inline constexpr int scaling_shape_set_version = 3;
inline constexpr int strong_patches = 8192;
inline constexpr std::array<std::pair<int, int>, 5> weak_scaling_shapes = {
    std::pair{1, 512}, std::pair{16, 2048}, std::pair{32, 2896}, std::pair{64, 4096}, std::pair{128, 5792}};
inline double maximum_refinement_error      = 0.0;
inline bool accuracy_checked                = false;

inline std::string experiment()
{
    const char* value = std::getenv("RUNTIME_EXPERIMENT");
    const std::string selected = value ? value : "runtime";
    if (selected != "runtime" && selected != "strong" && selected != "weak" && selected != "scaling")
        throw std::invalid_argument("RUNTIME_EXPERIMENT must be runtime, strong, weak or scaling.");
    return selected;
}

inline void check_openmp_team()
{
#ifdef _OPENMP
    const int requested = omp_get_max_threads();
    int actual = 0;
#pragma omp parallel num_threads(requested)
    {
#pragma omp single
        actual = omp_get_num_threads();
    }
    if (actual != requested)
        throw std::runtime_error("OpenMP created " + std::to_string(actual) + " threads, expected " +
                                 std::to_string(requested) + "; check OMP_THREAD_LIMIT and allocation.");
#endif
}

inline bool informational_argument(std::string_view argument)
{
    return argument == "--benchmark_list_tests" || argument == "--benchmark_list_tests=true" ||
           argument == "--benchmark_list_tests=1" || argument == "--help" || argument == "-h";
}

inline bool enabled()
{
    const char* value = std::getenv("MEMILIO_RUNTIME_SCENARIO");
    return value && std::string(value) == "1";
}

inline int positive_environment(const char* key, int fallback, int maximum)
{
    const char* value = std::getenv(key);
    if (!value)
        return fallback;
    std::string text(value);
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos) {
        throw std::invalid_argument(std::string(key) + " must be a positive integer.");
    }
    const auto number = std::stoll(text);
    if (number < 1 || number > maximum)
        throw std::invalid_argument(std::string(key) + " is out of range.");
    return static_cast<int>(number);
}

struct Schedule {
    int days           = 32;
    int half_day_steps = 32;
    double dt() const
    {
        return 0.5 / half_day_steps;
    }
    int total_steps() const
    {
        return 2 * days * half_day_steps;
    }
};

inline Schedule schedule()
{
    return {positive_environment("RUNTIME_DAYS", 32, 366), positive_environment("RUNTIME_HALF_DAY_STEPS", 32, 1024)};
}

// One source of physical inputs for implicit/explicit and CPU/CUDA. States use
// [(4 * group + compartment) * patches + resident_patch]; H[origin, destination].
struct Inputs {
    int p;
    int g;
    std::vector<double> initial, population, h, ht, beta, rate_e, rate_i;
    Inputs(int num_patches, int num_groups, bool identity_mobility = false)
        : p(num_patches)
        , g(num_groups)
        , initial(static_cast<size_t>(4 * g) * p)
        , population(static_cast<size_t>(g) * p)
        , h(static_cast<size_t>(p) * p)
        , ht(h.size())
        , beta(static_cast<size_t>(g) * g)
        , rate_e(g)
        , rate_i(g)
    {
        if (p < 2 || std::find(groups.begin(), groups.end(), g) == groups.end()) {
            throw std::invalid_argument("Unsupported runtime scenario shape.");
        }
        for (int origin = 0; origin < p; ++origin) {
            double weight = 0.0;
            for (int dest = 0; dest < p; ++dest) {
                if (dest != origin)
                    weight += 1.0 + 0.02 * ((17 * origin + 13 * dest) % 23);
            }
            const double mobile = identity_mobility ? 0.0 : 0.10 + 0.10 * (origin % 7) / 6.0;
            for (int dest = 0; dest < p; ++dest) {
                const double value =
                    dest == origin ? 1.0 - mobile : mobile * (1.0 + 0.02 * ((17 * origin + 13 * dest) % 23)) / weight;
                h[static_cast<size_t>(origin) * p + dest] = ht[static_cast<size_t>(dest) * p + origin] = value;
            }
            const double position = static_cast<double>(origin + 1) / (p + 1);
            for (int group = 0; group < g; ++group) {
                const double group_position = static_cast<double>(group + 1) / (g + 1);
                const double total          = 10000.0 * (0.9 + 0.2 * position) * (0.9 + 0.2 * group_position) / g;
                const double e              = total * (0.007 + 0.003 * position) / (1.0 + 0.1 * group);
                const double i              = total * (0.012 - 0.004 * position) * (1.0 + 0.08 * group);
                const double r              = total * (0.02 + 0.002 * group);
                population[static_cast<size_t>(group) * p + origin] = total;
                initial[index(origin, group, 0)]                    = total - e - i - r;
                initial[index(origin, group, 1)]                    = e;
                initial[index(origin, group, 2)]                    = i;
                initial[index(origin, group, 3)]                    = r;
            }
        }
        for (int target = 0; target < g; ++target) {
            rate_e[target] = 1.0 / (5.2 + 0.1 * target);
            rate_i[target] = 1.0 / (6.0 + 0.1 * target);
            for (int source = 0; source < g; ++source) {
                beta[static_cast<size_t>(target) * g + source] = (0.035 + 0.002 * target) * 9.0 *
                                                                 (1.0 + 0.02 * target + 0.01 * source) /
                                                                 (1.0 + std::abs(target - source));
            }
        }
    }
    size_t index(int patch, int group, int compartment) const
    {
        return static_cast<size_t>(4 * group + compartment) * p + patch;
    }
};

template <class Problem>
void configure_implicit(Problem& problem, const Inputs& inputs)
{
    problem.population                    = inputs.population;
    problem.commuting_row_major           = inputs.h;
    problem.commuting_transpose_row_major = inputs.ht;
    problem.rate_exposed                  = inputs.rate_e;
    problem.rate_infected                 = inputs.rate_i;
    for (size_t i = 0; i < inputs.beta.size(); ++i)
        problem.infection_coefficients[i] = 0.5 * inputs.beta[i];
    for (size_t i = 0; i < inputs.population.size(); ++i)
        problem.inverse_population[i] = 1.0 / inputs.population[i];
    for (int group = 0; group < inputs.g; ++group) {
        for (int c = 0; c < 3; ++c) {
            for (int patch = 0; patch < inputs.p; ++patch) {
                problem.initial_state[static_cast<size_t>(3 * group + c) * inputs.p + patch] =
                    inputs.initial[inputs.index(patch, group, c)];
            }
        }
    }
    problem.initialize_present_population();
    problem.reset_state();
}

inline std::vector<double> implicit_residents(const Inputs& inputs, const std::vector<double>& state)
{
    auto result = inputs.initial;
    for (int group = 0; group < inputs.g; ++group) {
        for (int patch = 0; patch < inputs.p; ++patch) {
            double recovered = inputs.population[static_cast<size_t>(group) * inputs.p + patch];
            for (int c = 0; c < 3; ++c) {
                const double value                    = state[static_cast<size_t>(3 * group + c) * inputs.p + patch];
                result[inputs.index(patch, group, c)] = value;
                recovered -= value;
            }
            result[inputs.index(patch, group, 3)] = recovered;
        }
    }
    return result;
}

inline void check_population(const Inputs& inputs, const std::vector<double>& residents)
{
    if (residents.size() != inputs.initial.size())
        throw std::runtime_error("Resident-state size mismatch.");
    for (int group = 0; group < inputs.g; ++group) {
        for (int patch = 0; patch < inputs.p; ++patch) {
            const double initial = inputs.population[static_cast<size_t>(group) * inputs.p + patch];
            double total         = 0.0;
            for (int c = 0; c < 4; ++c) {
                const double value = residents[inputs.index(patch, group, c)];
                if (!std::isfinite(value) || value < -1e-10 * initial) {
                    throw std::runtime_error("Non-finite/negative scenario population.");
                }
                total += value;
            }
            if (std::abs(total - initial) > 1e-9 * initial) {
                throw std::runtime_error("Resident population was not conserved through mobility events.");
            }
        }
    }
}

inline double compare(const Inputs& inputs, const std::vector<double>& expected, const std::vector<double>& actual,
                      double tolerance = 1e-10)
{
    check_population(inputs, expected);
    check_population(inputs, actual);
    double maximum = 0.0;
    for (int g = 0; g < inputs.g; ++g) {
        for (int p = 0; p < inputs.p; ++p) {
            for (int c = 0; c < 4; ++c) {
                maximum = std::max(maximum, std::abs(expected[inputs.index(p, g, c)] - actual[inputs.index(p, g, c)]) /
                                                inputs.population[static_cast<size_t>(g) * inputs.p + p]);
            }
        }
    }
    if (maximum > tolerance)
        throw std::runtime_error("Scenario validation error: " + std::to_string(maximum));
    return maximum;
}

// Independent, deliberately simple full-state RK4 reference for validation only.
// Explicit states include the diagonal and reconstruct all aggregates at EVERY
// RK stage. No optimized traveler kernel or event implementation is used here.
class Reference
{
public:
    Reference(const Inputs& inputs, bool explicit_model)
        : in(inputs)
        , explicit_mobility(explicit_model)
        , state(inputs.initial.size() * (explicit_model ? inputs.p : 1), 0.0)
        , k1(state.size())
        , k2(state.size())
        , k3(state.size())
        , k4(state.size())
        , temporary(state.size())
        , totals(inputs.initial.size())
        , pressure(static_cast<size_t>(inputs.p) * inputs.g)
        , present(pressure.size())
    {
        if (explicit_mobility) {
            for (int p = 0; p < in.p; ++p) {
                for (int g = 0; g < in.g; ++g) {
                    for (int c = 0; c < 4; ++c)
                        state[located(p, p, g, c)] = in.initial[in.index(p, g, c)];
                }
            }
        }
        else
            state = in.initial;
    }
    std::vector<double> run(Schedule times)
    {
        for (int day = 0; day < times.days; ++day) {
            for (int n = 0; n < times.half_day_steps; ++n)
                step(times.dt());
            if (explicit_mobility)
                depart();
            for (int n = 0; n < times.half_day_steps; ++n)
                step(times.dt());
            if (explicit_mobility)
                return_home(); // final return is included
        }
        return residents();
    }

private:
    size_t located(int origin, int destination, int group, int c) const
    {
        return (static_cast<size_t>(destination) * 4 * in.g + 4 * group + c) * in.p + origin;
    }
    std::vector<double> residents() const
    {
        if (!explicit_mobility)
            return state;
        std::vector<double> result(in.initial.size(), 0.0);
        for (int origin = 0; origin < in.p; ++origin) {
            for (int destination = 0; destination < in.p; ++destination) {
                for (int g = 0; g < in.g; ++g) {
                    for (int c = 0; c < 4; ++c)
                        result[in.index(origin, g, c)] += state[located(origin, destination, g, c)];
                }
            }
        }
        return result;
    }
    void depart()
    {
        const auto home = residents();
        for (int origin = 0; origin < in.p; ++origin) {
            for (int dest = 0; dest < in.p; ++dest) {
                for (int g = 0; g < in.g; ++g) {
                    for (int c = 0; c < 4; ++c) {
                        state[located(origin, dest, g, c)] =
                            in.h[static_cast<size_t>(origin) * in.p + dest] * home[in.index(origin, g, c)];
                    }
                }
            }
        }
    }
    void return_home()
    {
        const auto home = residents();
        std::fill(state.begin(), state.end(), 0.0);
        for (int p = 0; p < in.p; ++p) {
            for (int g = 0; g < in.g; ++g) {
                for (int c = 0; c < 4; ++c)
                    state[located(p, p, g, c)] = home[in.index(p, g, c)];
            }
        }
    }
    void rhs(const std::vector<double>& y, std::vector<double>& derivative)
    {
        std::fill(totals.begin(), totals.end(), 0.0);
        if (explicit_mobility) {
            for (int dest = 0; dest < in.p; ++dest) {
                for (int origin = 0; origin < in.p; ++origin) {
                    for (int g = 0; g < in.g; ++g) {
                        for (int c = 0; c < 4; ++c)
                            totals[in.index(dest, g, c)] += y[located(origin, dest, g, c)];
                    }
                }
            }
        }
        else {
            for (int dest = 0; dest < in.p; ++dest) {
                for (int g = 0; g < in.g; ++g) {
                    double infectious = 0.0, population = 0.0;
                    for (int origin = 0; origin < in.p; ++origin) {
                        const double fraction = in.h[static_cast<size_t>(origin) * in.p + dest];
                        infectious += fraction * y[in.index(origin, g, 2)];
                        population += fraction * in.population[static_cast<size_t>(g) * in.p + origin];
                    }
                    present[static_cast<size_t>(g) * in.p + dest] = infectious / population;
                }
            }
        }
        for (int p = 0; p < in.p; ++p) {
            for (int target = 0; target < in.g; ++target) {
                double lambda = 0.0;
                for (int source = 0; source < in.g; ++source) {
                    double fraction = 0.0;
                    if (explicit_mobility) {
                        double population = 0.0;
                        for (int c = 0; c < 4; ++c)
                            population += totals[in.index(p, source, c)];
                        if (population > 1e-12)
                            fraction = totals[in.index(p, source, 2)] / population;
                    }
                    else {
                        fraction =
                            0.5 * y[in.index(p, source, 2)] / in.population[static_cast<size_t>(source) * in.p + p];
                        for (int dest = 0; dest < in.p; ++dest) {
                            fraction += 0.5 * in.h[static_cast<size_t>(p) * in.p + dest] *
                                        present[static_cast<size_t>(source) * in.p + dest];
                        }
                    }
                    lambda += in.beta[static_cast<size_t>(target) * in.g + source] * fraction;
                }
                pressure[static_cast<size_t>(target) * in.p + p] = lambda;
                for (int origin = 0; origin < (explicit_mobility ? in.p : 1); ++origin) {
                    const auto at = [&](int c) {
                        return explicit_mobility ? located(origin, p, target, c) : in.index(p, target, c);
                    };
                    const double se = lambda * y[at(0)], ei = in.rate_e[target] * y[at(1)],
                                 ir   = in.rate_i[target] * y[at(2)];
                    derivative[at(0)] = -se;
                    derivative[at(1)] = se - ei;
                    derivative[at(2)] = ei - ir;
                    derivative[at(3)] = ir;
                }
            }
        }
    }
    void step(double dt)
    {
        rhs(state, k1);
        for (size_t i = 0; i < state.size(); ++i)
            temporary[i] = state[i] + 0.5 * dt * k1[i];
        rhs(temporary, k2);
        for (size_t i = 0; i < state.size(); ++i)
            temporary[i] = state[i] + 0.5 * dt * k2[i];
        rhs(temporary, k3);
        for (size_t i = 0; i < state.size(); ++i)
            temporary[i] = state[i] + dt * k3[i];
        rhs(temporary, k4);
        for (size_t i = 0; i < state.size(); ++i)
            state[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
    const Inputs& in;
    bool explicit_mobility;
    std::vector<double> state, k1, k2, k3, k4, temporary, totals, pressure, present;
};

inline void validate_accuracy()
{
    // Before any measured work: requested resource counts must be real teams.
    check_openmp_team();
    accuracy_checked  = false;
    const auto coarse = schedule();
    const Schedule fine{coarse.days, 2 * coarse.half_day_steps};
    for (int g : groups) {
        Inputs in(8, g);
        for (bool explicit_model : {false, true}) {
            auto expected            = Reference(in, explicit_model).run(fine);
            auto actual              = Reference(in, explicit_model).run(coarse);
            const double error       = compare(in, expected, actual, accuracy_tolerance);
            maximum_refinement_error = std::max(maximum_refinement_error, error);
            std::cout << "Runtime step refinement " << (explicit_model ? "explicit" : "implicit") << " N_G=" << g
                      << " max_population_fraction_error=" << error << '\n';
        }
        Inputs no_mobility(3, g, true);
        compare(no_mobility, Reference(no_mobility, false).run({2, coarse.half_day_steps}),
                Reference(no_mobility, true).run({2, coarse.half_day_steps}));
    }
    accuracy_checked = true;
}

inline void counters(benchmark::State& state, const Inputs& in, bool explicit_model, int threads)
{
    if (!accuracy_checked)
        throw std::runtime_error("The runtime accuracy gate has not been executed.");
    const auto time                                 = schedule();
    state.counters["scenario_version"]              = version;
    state.counters["shape_set_version"]             =
        experiment() == "runtime" ? shape_set_version : scaling_shape_set_version;
    state.counters["patches"]                       = in.p;
    state.counters["age_groups"]                    = in.g;
    state.counters["cpu_threads"]                   = threads;
    state.counters["simulation_days"]               = time.days;
    state.counters["step_days"]                     = time.dt();
    state.counters["half_day_steps"]                = time.half_day_steps;
    state.counters["steps"]                         = time.total_steps();
    state.counters["event_interval_days"]           = 0.5;
    state.counters["mobility_events"]               = explicit_model ? 2 * time.days : 0;
    state.counters["events_timed"]                  = explicit_model ? 1 : 0;
    state.counters["home_fraction"]                 = 0.5;
    state.counters["validation_passed"]             = 1;
    state.counters["refinement_error"]              = maximum_refinement_error;
    state.counters["accuracy_tolerance"]            = accuracy_tolerance;
    state.counters["refinement_validation_patches"] = 8;
}

inline void register_shapes(const char* name, void (*function)(benchmark::State&), int threads = 0)
{
    const auto selected = experiment();
    if (selected != "runtime") {
        // Scaling uses the exact same OpenMP functions, including OMP1. Never
        // register serial/GPU baselines as scaling cases.
        if (threads <= 0)
            return;
        const auto shape = std::find_if(weak_scaling_shapes.begin(), weak_scaling_shapes.end(),
                                       [threads](const auto& item) { return item.first == threads; });
        if (shape == weak_scaling_shapes.end())
            throw std::invalid_argument("Daily scaling supports 1, 16, 32, 64 or 128 OpenMP threads.");
        for (const char* kind : {"strong", "weak"}) {
            if (selected != "scaling" && selected != kind)
                continue;
            const auto scaling_name = std::string(kind) + std::string(name).substr(7);
            benchmark::RegisterBenchmark(scaling_name.c_str(), function)
                ->Args({std::string_view(kind) == "strong" ? strong_patches : shape->second, 6, threads})
                ->ArgNames({"patches", "age_groups", "threads"})
                ->UseRealTime();
        }
        return;
    }
    auto* b = benchmark::RegisterBenchmark(name, function);
    for (int p : patches) {
        for (int g : groups) {
            if (threads > 0)
                b->Args({p, g, threads});
            else
                b->Args({p, g});
        }
    }
    // Extend the size range at G=6 only, not the full patch/age-group product.
    for (int p : large_patches) {
        if (threads > 0)
            b->Args({p, 6, threads});
        else
            b->Args({p, 6});
    }
    if (threads > 0)
        b->ArgNames({"patches", "age_groups", "threads"});
    else
        b->ArgNames({"patches", "age_groups"});
    b->UseRealTime();
}
} // namespace mio::runtime_scenario
#endif
