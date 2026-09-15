/* Copyright (C) 2026 MEmilio. Licensed under the Apache License, Version 2.0. */
#ifndef MIO_ODE_SEIR_RUNTIME_EXPLICIT_H
#define MIO_ODE_SEIR_RUNTIME_EXPLICIT_H
#include "ode_seir_benchmark_stage_aligned.h"
#include "ode_seir_runtime_scenario.h"
#include <chrono>

namespace mio::runtime_scenario
{
struct ExplicitProblem {
    benchmark_mio::StageAlignedProblem core;
    std::vector<double> residents;
    const Inputs& inputs;
    ExplicitProblem(const Inputs& in, int threads = 0)
        : core(in.p, in.p - 1, in.g, threads)
        , residents(in.initial.size())
        , inputs(in)
    {
        core.contact_beta  = in.beta;
        core.rate_exposed  = in.rate_e;
        core.rate_infected = in.rate_i;
        reset();
    }
    void reset()
    {
        core.totals = inputs.initial;
        // At home the traveler buffer is inactive. Every entry is overwritten
        // at departure before it can be read; no stale traveler state is used.
    }
};

template <int G>
inline void explicit_home_workshare(ExplicitProblem& p, Schedule time)
{
    auto& core = p.core;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int patch = 0; patch < core.patches; ++patch) {
        for (int step = 0; step < time.half_day_steps; ++step)
            benchmark_mio::integrate_totals<G>(core, patch, time.dt());
    }
}

template <int G>
inline void explicit_departure_workshare(ExplicitProblem& p)
{
    auto& core        = p.core;
    const int patches = core.patches, travelers = patches - 1;
    // Departure must read an immutable resident snapshot while destination
    // totals and off-diagonal travelers are populated in parallel.
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (size_t i = 0; i < core.totals.size(); ++i)
        p.residents[i] = core.totals[i];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int dest = 0; dest < patches; ++dest) {
        for (int c = 0; c < 4 * G; ++c) {
            double total = p.inputs.ht[static_cast<size_t>(dest) * patches + dest] *
                           p.residents[static_cast<size_t>(c) * patches + dest];
            auto* row = core.travelers.data() + (static_cast<size_t>(dest) * 4 * G + c) * travelers;
#ifdef _OPENMP
#pragma omp simd reduction(+ : total)
#endif
            for (int t = 0; t < travelers; ++t) {
                const int origin = t + (t >= dest);
                row[t]           = p.inputs.ht[static_cast<size_t>(dest) * patches + origin] *
                         p.residents[static_cast<size_t>(c) * patches + origin];
                total += row[t];
            }
            core.totals[static_cast<size_t>(c) * patches + dest] = total;
        }
    }
}

template <int G>
inline void explicit_away_workshare(ExplicitProblem& p, Schedule time)
{
    auto& core = p.core;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int patch = 0; patch < core.patches; ++patch) {
        // Never cross the return event, even if the kernel's time tile is larger.
        benchmark_mio::advance_stage_aligned_patch<G>(core, patch, time.dt(), time.half_day_steps);
    }
}

template <int G>
inline void explicit_return_workshare(ExplicitProblem& p)
{
    auto& core        = p.core;
    const int patches = core.patches, travelers = patches - 1;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int origin = 0; origin < patches; ++origin) {
        for (int c = 0; c < 4 * G; ++c) {
            const auto* incoming = core.travelers.data() + (static_cast<size_t>(origin) * 4 * G + c) * travelers;
            double home          = core.totals[static_cast<size_t>(c) * patches + origin];
            for (int t = 0; t < travelers; ++t)
                home -= incoming[t];
            for (int dest = 0; dest < patches; ++dest) {
                if (dest != origin) {
                    const int t = origin - (origin > dest);
                    home += core.travelers[(static_cast<size_t>(dest) * 4 * G + c) * travelers + t];
                }
            }
            p.residents[static_cast<size_t>(c) * patches + origin] = home;
        }
    }
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (size_t i = 0; i < core.totals.size(); ++i)
        core.totals[i] = p.residents[i];
}

template <int G>
inline void explicit_day_workshare(ExplicitProblem& p, Schedule time)
{
    explicit_home_workshare<G>(p, time);
    explicit_departure_workshare<G>(p);
    explicit_away_workshare<G>(p, time);
    explicit_return_workshare<G>(p);
}

template <int G>
void advance_explicit(ExplicitProblem& problem, Schedule time, int threads = 0)
{
#ifndef _OPENMP
    if (threads > 0)
        throw std::invalid_argument("Runtime OpenMP requires OpenMP support.");
#else
#pragma omp parallel if (threads > 0) num_threads(threads > 0 ? threads : 1)
#endif
    {
        for (int day = 0; day < time.days; ++day)
            explicit_day_workshare<G>(problem, time);
    }
}

struct ExplicitPhaseTimes {
    double home_seconds      = 0.0;
    double departure_seconds = 0.0;
    double away_seconds      = 0.0;
    double return_seconds    = 0.0;
    double profiled_seconds  = 0.0;
};

using ExplicitPhaseClock = std::chrono::steady_clock;

template <class Work>
inline void profile_explicit_workshare(double& seconds, Work&& work)
{
    // Each thread owns its start variable, but only the same master thread
    // writes and reads its timestamp. Workers cannot start before this clock.
    ExplicitPhaseClock::time_point start;
#ifdef _OPENMP
#pragma omp master
#endif
    {
        start = ExplicitPhaseClock::now();
    }
#ifdef _OPENMP
#pragma omp barrier
#endif
    work(); // Its existing final workshare barrier waits for every worker.
#ifdef _OPENMP
#pragma omp master
#endif
    {
        seconds += std::chrono::duration<double>(ExplicitPhaseClock::now() - start).count();
    }
}

template <int G>
ExplicitPhaseTimes advance_explicit_profiled(ExplicitProblem& problem, Schedule time, int threads = 0)
{
#ifndef _OPENMP
    if (threads > 0)
        throw std::invalid_argument("Runtime OpenMP requires OpenMP support.");
#endif
    ExplicitPhaseTimes result;
    const auto start = ExplicitPhaseClock::now();
#ifdef _OPENMP
#pragma omp parallel if (threads > 0) num_threads(threads > 0 ? threads : 1)
#endif
    {
        for (int day = 0; day < time.days; ++day) {
            profile_explicit_workshare(result.home_seconds, [&] { explicit_home_workshare<G>(problem, time); });
            profile_explicit_workshare(result.departure_seconds, [&] { explicit_departure_workshare<G>(problem); });
            profile_explicit_workshare(result.away_seconds, [&] { explicit_away_workshare<G>(problem, time); });
            profile_explicit_workshare(result.return_seconds, [&] { explicit_return_workshare<G>(problem); });
        }
    }
    result.profiled_seconds = std::chrono::duration<double>(ExplicitPhaseClock::now() - start).count();
    return result;
}

template <int G>
void validate_explicit_cpu(int threads, bool validate_phases = false)
{
    for (bool no_mobility : {false, true}) {
        Inputs in(3, G, no_mobility);
        for (int phase_steps : {1, 3, 32, 65}) {
            Schedule time{2, phase_steps};
            const auto expected = Reference(in, true).run(time);
            ExplicitProblem actual(in, threads);
            for (int replay = 0; replay < 2; ++replay) {
                actual.reset();
                advance_explicit<G>(actual, time, threads);
                compare(in, expected, actual.core.totals);
                if (validate_phases) {
                    const auto control = actual.core.totals;
                    actual.reset();
                    const auto phases = advance_explicit_profiled<G>(actual, time, threads);
                    compare(in, expected, actual.core.totals);
                    compare(in, control, actual.core.totals);
                    for (double elapsed : {phases.home_seconds, phases.departure_seconds, phases.away_seconds,
                                           phases.return_seconds, phases.profiled_seconds}) {
                        if (!std::isfinite(elapsed) || elapsed < 0.0)
                            throw std::runtime_error("Invalid explicit phase duration.");
                    }
                    const double accounted = phases.home_seconds + phases.departure_seconds + phases.away_seconds +
                                             phases.return_seconds;
                    if (phases.profiled_seconds <= 0.0 || accounted > phases.profiled_seconds * (1.0 + 1e-10))
                        throw std::runtime_error("Explicit phase durations exceed the profiled trajectory.");
                }
            }
        }
    }
    std::cout << "Runtime explicit full-state reference N_G=" << G << ": passed\n";
    if (validate_phases)
        std::cout << "Explicit phase/control reference and accounting N_G=" << G << ": passed\n";
}
} // namespace mio::runtime_scenario
#endif
