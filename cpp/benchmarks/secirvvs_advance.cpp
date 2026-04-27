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

/**
 * Benchmark comparing the overhead of the model-specific osecirvvs::Simulation::advance()
 * (which uses substeps and calls apply_vaccination / apply_variant / dynamic NPI checks)
 * versus the generic mio::Simulation::advance() (a single integrator call for the entire range).
 */

#include "ode_secirvvs/model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "benchmark/benchmark.h"

using FP    = double;
using Model = mio::osecirvvs::Model<FP>;

static Model make_model(bool with_npis = false)
{
    constexpr int tmax = 10;
    Model model(1);
    model.populations[{mio::AgeGroup(0), mio::osecirvvs::InfectionState::InfectedSymptomsNaive}] = 100.0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecirvvs::InfectionState::SusceptibleNaive},
                                                10000.0);
    model.parameters.get<mio::osecirvvs::DailyPartialVaccinations<FP>>().resize(mio::SimulationDay(size_t(tmax + 1)));
    model.parameters.get<mio::osecirvvs::DailyFullVaccinations<FP>>().resize(mio::SimulationDay(size_t(tmax + 1)));
    if (with_npis) {
        auto& npis = model.parameters.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms<FP>>();
        npis.set_threshold(0.01 * 100'000, {mio::DampingSampling<FP>{1.0,
                                                                     mio::DampingLevel(0),
                                                                     mio::DampingType(0),
                                                                     mio::SimulationTime<FP>(0),
                                                                     {0},
                                                                     Eigen::VectorXd::Ones(1)}});
        npis.set_duration(mio::SimulationTime<FP>(14.0));
        npis.set_base_value(100'000);
    }
    return model;
}

// Generic advance: single integrator call for the full [t0, tmax] range.
// Uses the fixed-step Euler integrator so that the number of integration steps
// is deterministic and independent of contact-rate changes.
static void BM_generic(benchmark::State& state)
{
    mio::set_log_level(mio::LogLevel::off);
    auto model = make_model();
    for (auto _ : state) {
        mio::simulate<FP, Model>(0., 10., 0.1, model, std::make_unique<mio::EulerIntegratorCore<FP>>());
    }
}

// Model-specific advance without dynamic NPIs: 1-day loop with apply_vaccination + apply_variant,
// dynamic NPI threshold check is skipped (thresholds empty).
static void BM_secirvvs_no_npis(benchmark::State& state)
{
    mio::set_log_level(mio::LogLevel::off);
    auto model = make_model(/*with_npis=*/false);
    for (auto _ : state) {
        mio::osecirvvs::simulate<FP>(0., 10., 0.1, model, std::make_unique<mio::EulerIntegratorCore<FP>>());
    }
}

// Model-specific advance with dynamic NPIs: same as above plus get_infections_relative +
// threshold comparison on every day step.
static void BM_secirvvs_with_npis(benchmark::State& state)
{
    mio::set_log_level(mio::LogLevel::off);
    auto model = make_model(/*with_npis=*/true);
    for (auto _ : state) {
        mio::osecirvvs::simulate<FP>(0., 10., 0.1, model, std::make_unique<mio::EulerIntegratorCore<FP>>());
    }
}

BENCHMARK(BM_generic)->Name("SECIRVVS generic advance");
BENCHMARK(BM_secirvvs_no_npis)->Name("SECIRVVS model-specific advance (no dynamic NPIs)");
BENCHMARK(BM_secirvvs_with_npis)->Name("SECIRVVS model-specific advance (with dynamic NPIs)");
BENCHMARK_MAIN();
