/*
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/compartments/flow_simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"
#include "ode_seir/model.h"
#include "ode_seir/state_estimators.h"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/euler.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/euler.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#endif

namespace mio
{
namespace oseir
{

struct SEIRIndexCache {
    std::vector<int> S, E, I, R;
    size_t G = 0;
};

inline SEIRIndexCache make_seir_index_cache(const mio::oseir::Model<ScalarType>& model)
{
    SEIRIndexCache ic;
    ic.G = static_cast<size_t>(model.parameters.get_num_groups());
    ic.S.resize(ic.G);
    ic.E.resize(ic.G);
    ic.I.resize(ic.G);
    ic.R.resize(ic.G);
    using IS = mio::oseir::InfectionState;
    for (size_t g = 0; g < ic.G; ++g) {
        auto Gg = mio::AgeGroup(static_cast<int>(g));
        ic.S[g] = model.populations.get_flat_index({Gg, IS::Susceptible});
        ic.E[g] = model.populations.get_flat_index({Gg, IS::Exposed});
        ic.I[g] = model.populations.get_flat_index({Gg, IS::Infected});
        ic.R[g] = model.populations.get_flat_index({Gg, IS::Recovered});
    }
    return ic;
}

} // namespace oseir

namespace benchmark_mio
{

using ModelType = mio::oseir::Model<ScalarType>;
using SimType   = mio::FlowSimulation<ScalarType, ModelType>;

void setup_model_benchmark(ModelType& model, size_t num_agegroups = 6)
{
    const double total_population = 10000.0;
    for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); i++) {
        model.populations[{i, mio::oseir::InfectionState::Exposed}]   = 100.0 / static_cast<double>(num_agegroups);
        model.populations[{i, mio::oseir::InfectionState::Infected}]  = 100.0 / static_cast<double>(num_agegroups);
        model.populations[{i, mio::oseir::InfectionState::Recovered}] = 100.0 / static_cast<double>(num_agegroups);
        model.populations[{i, mio::oseir::InfectionState::Susceptible}] =
            total_population / static_cast<double>(num_agegroups) -
            model.populations[{i, mio::oseir::InfectionState::Exposed}] -
            model.populations[{i, mio::oseir::InfectionState::Infected}] -
            model.populations[{i, mio::oseir::InfectionState::Recovered}];
    }
    model.parameters.template set<mio::oseir::TimeExposed<ScalarType>>(5.2);
    model.parameters.template set<mio::oseir::TimeInfected<ScalarType>>(6);
    model.parameters.template set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.1);
    mio::ContactMatrixGroup& contact_matrix = model.parameters.template get<mio::oseir::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(2.7);
    model.apply_constraints();
}

void setup_explicit_model_benchmark(mio::oseir::StandardModelLagrangian& model_explicit, size_t num_age_groups,
                                    int num_commuter_groups, const double total_population = 10000.0)
{

    // Start with non commuting pop
    for (auto i = 0; i < static_cast<int>(num_age_groups); i++) {
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType, mio::oseir::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::oseir::CommuterType::NonCommuter, mio::oseir::InfectionStateExplicit::E}] =
            50.0 / static_cast<double>(num_age_groups);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType, mio::oseir::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::oseir::CommuterType::NonCommuter, mio::oseir::InfectionStateExplicit::I}] =
            50.0 / static_cast<double>(num_age_groups);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType, mio::oseir::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::oseir::CommuterType::NonCommuter, mio::oseir::InfectionStateExplicit::R}] =
            50.0 / static_cast<double>(num_age_groups);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType, mio::oseir::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::oseir::CommuterType::NonCommuter, mio::oseir::InfectionStateExplicit::S}] =
            (total_population * 0.5) / static_cast<double>(num_age_groups) -
            model_explicit
                .populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType, mio::oseir::InfectionStateExplicit>{
                    mio::AgeGroup(i), mio::oseir::CommuterType::NonCommuter, mio::oseir::InfectionStateExplicit::E}] -
            model_explicit
                .populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType, mio::oseir::InfectionStateExplicit>{
                    mio::AgeGroup(i), mio::oseir::CommuterType::NonCommuter, mio::oseir::InfectionStateExplicit::I}] -
            model_explicit
                .populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType, mio::oseir::InfectionStateExplicit>{
                    mio::AgeGroup(i), mio::oseir::CommuterType::NonCommuter, mio::oseir::InfectionStateExplicit::R}];
    }

    // Then the remaining commuting pop
    for (int c = 0; c < num_commuter_groups; c++) {
        auto commuter_type =
            static_cast<mio::oseir::CommuterType>(static_cast<int>(mio::oseir::CommuterType::CommuterBase) + c);
        double commuter_fraction = 0.5 / static_cast<double>(num_commuter_groups);

        for (auto i = 0; i < static_cast<int>(num_age_groups); i++) {
            model_explicit
                .populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType, mio::oseir::InfectionStateExplicit>{
                    mio::AgeGroup(i), commuter_type, mio::oseir::InfectionStateExplicit::E}] =
                50.0 / static_cast<double>(num_age_groups * num_commuter_groups);
            model_explicit
                .populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType, mio::oseir::InfectionStateExplicit>{
                    mio::AgeGroup(i), commuter_type, mio::oseir::InfectionStateExplicit::I}] =
                50.0 / static_cast<double>(num_age_groups * num_commuter_groups);
            model_explicit
                .populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType, mio::oseir::InfectionStateExplicit>{
                    mio::AgeGroup(i), commuter_type, mio::oseir::InfectionStateExplicit::R}] =
                50.0 / static_cast<double>(num_age_groups * num_commuter_groups);

            double exposed   = model_explicit.populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType,
                                                                   mio::oseir::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::oseir::InfectionStateExplicit::E}];
            double infected  = model_explicit.populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType,
                                                                    mio::oseir::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::oseir::InfectionStateExplicit::I}];
            double recovered = model_explicit.populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType,
                                                                     mio::oseir::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::oseir::InfectionStateExplicit::R}];

            model_explicit
                .populations[mio::Index<mio::AgeGroup, mio::oseir::CommuterType, mio::oseir::InfectionStateExplicit>{
                    mio::AgeGroup(i), commuter_type, mio::oseir::InfectionStateExplicit::S}] =
                (total_population * commuter_fraction) / static_cast<double>(num_age_groups) - exposed - infected -
                recovered;
        }
    }

    model_explicit.parameters.set<mio::oseir::TimeExposed<ScalarType>>(5.2);
    model_explicit.parameters.set<mio::oseir::TimeInfected<ScalarType>>(6);
    model_explicit.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.04);
    mio::ContactMatrixGroup& contact_matrix = model_explicit.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(1.0);
    model_explicit.check_constraints();
}

const ScalarType t0 = 0.0;
// const ScalarType t_max = 0.2;
// const ScalarType dt    = 0.1;

const ScalarType t_max = 50.0;
const ScalarType dt    = 0.5;

const ScalarType t_max_phi = 1;
const ScalarType dt_phi    = 0.5;

namespace
{
struct BenchSetupPrinter {
    BenchSetupPrinter()
    {
        std::cout << "Benchmark setup:\n";
        std::cout << "  t_max (normal)    = " << t_max << "\n";
        std::cout << "  dt (normal)       = " << dt << "\n";
        std::cout << "  t_max_phi (phi)   = " << t_max_phi << "\n";
        std::cout << "  dt_phi (phi)      = " << dt_phi << "\n";
        std::cout << std::flush;
    }
};
static BenchSetupPrinter bench_setup_printer;
} // namespace
// Define the number of commuter groups for the benchmark
// const std::vector<int> commuter_group_counts = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
// const std::vector<int> age_group_counts      = {1, 2, 3, 4, 5, 6};

const std::vector<int> commuter_group_counts = {256, 512, 1024, 2048, 4096, 8192, 16384};
const std::vector<int> age_group_counts      = {1, 3, 6}; //, 8, 12, 16};
const std::vector<int> core_counts           = {1, 7, 14, 28, 56};

static void bench_auxiliary_euler(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();
            ModelType model(num_age_groups);
            setup_model_benchmark(model);
            SimType sim(model, t0, dt);
            auto integrator_rk =
                std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
            sim.set_integrator(integrator_rk);
            state.ResumeTiming();

            sim.advance(t_max);

            state.PauseTiming();
            const auto& seir_res              = sim.get_result();
            double mobile_population_fraction = 0.1 * num_commuter_groups;
            Eigen::VectorXd initial_mobile_pop =
                seir_res.get_value(0) * (mobile_population_fraction / num_commuter_groups);
            const auto step_size_ref = seir_res.get_time(1) - seir_res.get_time(0);

            std::vector<Eigen::VectorXd> mobile_pops(num_commuter_groups, initial_mobile_pop);
            state.ResumeTiming();

            for (ScalarType t = t0; t < t_max; t += dt) {
                if (t + dt > t_max + 1e-10)
                    break;

                const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
                if (closest_idx_total >= static_cast<size_t>(seir_res.get_num_time_points()))
                    break;

                const auto& total_pop = seir_res.get_value(closest_idx_total);

                for (int i = 0; i < num_commuter_groups; ++i) {
                    mio::oseir::integrate_mobile_population_euler(mobile_pops[i], sim, total_pop, t, dt);
                }
            }
            benchmark::DoNotOptimize(mobile_pops);
        }
    }
}

static void bench_stage_aligned_rk4(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    struct RK4StageCache {
        std::vector<Eigen::VectorXd> y;
        std::vector<Eigen::VectorXd> k_tot;
        std::vector<std::vector<double>> lambda;

        void resize(size_t nc, size_t num_age_groups)
        {
            y.resize(4, Eigen::VectorXd(nc));
            k_tot.resize(4, Eigen::VectorXd(nc));
            lambda.resize(4, std::vector<double>(num_age_groups, 0.0));
        }
    };

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();
            ModelType model(num_age_groups);
            setup_model_benchmark(model);

            Eigen::VectorXd current_totals = model.populations.get_compartments();

            double mobile_fraction = 0.1 * num_commuter_groups;
            Eigen::VectorXd initial_mobile =
                current_totals * (num_commuter_groups > 0 ? (mobile_fraction / num_commuter_groups) : 0.0);
            std::vector<Eigen::VectorXd> mobile_pops(num_commuter_groups, initial_mobile);

            const size_t NC = current_totals.size();
            RK4StageCache cache;
            cache.resize(NC, num_age_groups);

            auto ic = mio::oseir::make_seir_index_cache(model);

            std::vector<double> rate_E(num_age_groups, 0.0), rate_I(num_age_groups, 0.0);
            for (size_t g = 0; g < num_age_groups; ++g) {
                double t_E =
                    model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(static_cast<int>(g))];
                double t_I =
                    model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(static_cast<int>(g))];
                rate_E[g] = (t_E > 1e-10) ? (1.0 / t_E) : 0.0;
                rate_I[g] = (t_I > 1e-10) ? (1.0 / t_I) : 0.0;
            }

            auto compute_lambda = [&](const Eigen::VectorXd& y, double tt, std::vector<double>& lambda_out) {
                lambda_out.assign(num_age_groups, 0.0);
                for (size_t j = 0; j < num_age_groups; ++j) {
                    const double Nj    = y[ic.S[j]] + y[ic.E[j]] + y[ic.I[j]] + y[ic.R[j]];
                    const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;
                    for (size_t i = 0; i < num_age_groups; ++i) {
                        const double coeff =
                            model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>()
                                .get_cont_freq_mat()
                                .get_matrix_at(tt)(i, j) *
                            model.parameters
                                .get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(
                                    static_cast<int>(i))] *
                            divNj;
                        lambda_out[i] += coeff * y[ic.I[j]];
                    }
                }
            };

            Eigen::VectorXd k1_com(NC), k2_com(NC), k3_com(NC), k4_com(NC);
            Eigen::VectorXd Xc2(NC), Xc3(NC), Xc4(NC);

            state.ResumeTiming();

            for (ScalarType t = t0; t < t_max; t += dt) {

                // Integrate totals
                // Stage 1
                cache.y[0] = current_totals;
                model.get_derivatives(cache.y[0], cache.y[0], t, cache.k_tot[0]);
                compute_lambda(cache.y[0], t, cache.lambda[0]);

                // Stage 2
                cache.y[1] = cache.y[0] + (dt * 0.5) * cache.k_tot[0];
                model.get_derivatives(cache.y[1], cache.y[1], t + 0.5 * dt, cache.k_tot[1]);
                compute_lambda(cache.y[1], t + 0.5 * dt, cache.lambda[1]);

                // Stage 3
                cache.y[2] = cache.y[0] + (dt * 0.5) * cache.k_tot[1];
                model.get_derivatives(cache.y[2], cache.y[2], t + 0.5 * dt, cache.k_tot[2]);
                compute_lambda(cache.y[2], t + 0.5 * dt, cache.lambda[2]);

                // Stage 4
                cache.y[3] = cache.y[0] + dt * cache.k_tot[2];
                model.get_derivatives(cache.y[3], cache.y[3], t + dt, cache.k_tot[3]);
                compute_lambda(cache.y[3], t + dt, cache.lambda[3]);

                Eigen::VectorXd next_totals = current_totals + (dt / 6.0) * (cache.k_tot[0] + 2 * cache.k_tot[1] +
                                                                             2 * cache.k_tot[2] + cache.k_tot[3]);

                // Reconstruction
                for (int cg = 0; cg < num_commuter_groups; ++cg) {
                    Eigen::Ref<Eigen::VectorXd> Xc = mobile_pops[cg];

                    // Stage 1
                    for (size_t g = 0; g < num_age_groups; ++g) {
                        int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        double fSE = cache.lambda[0][g] * Xc[iS];
                        double fEI = rate_E[g] * Xc[iE];
                        double fIR = rate_I[g] * Xc[iI];
                        k1_com[iS] = -fSE;
                        k1_com[iE] = fSE - fEI;
                        k1_com[iI] = fEI - fIR;
                        k1_com[iR] = fIR;
                    }

                    // Stage 2
                    Xc2 = Xc + (dt * 0.5) * k1_com;
                    for (size_t g = 0; g < num_age_groups; ++g) {
                        int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        double fSE = cache.lambda[1][g] * Xc2[iS];
                        double fEI = rate_E[g] * Xc2[iE];
                        double fIR = rate_I[g] * Xc2[iI];
                        k2_com[iS] = -fSE;
                        k2_com[iE] = fSE - fEI;
                        k2_com[iI] = fEI - fIR;
                        k2_com[iR] = fIR;
                    }

                    // Stage 3
                    Xc3 = Xc + (dt * 0.5) * k2_com;
                    for (size_t g = 0; g < num_age_groups; ++g) {
                        int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        double fSE = cache.lambda[2][g] * Xc3[iS];
                        double fEI = rate_E[g] * Xc3[iE];
                        double fIR = rate_I[g] * Xc3[iI];
                        k3_com[iS] = -fSE;
                        k3_com[iE] = fSE - fEI;
                        k3_com[iI] = fEI - fIR;
                        k3_com[iR] = fIR;
                    }

                    // Stage 4
                    Xc4 = Xc + dt * k3_com;
                    for (size_t g = 0; g < num_age_groups; ++g) {
                        int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        double fSE = cache.lambda[3][g] * Xc4[iS];
                        double fEI = rate_E[g] * Xc4[iE];
                        double fIR = rate_I[g] * Xc4[iI];
                        k4_com[iS] = -fSE;
                        k4_com[iE] = fSE - fEI;
                        k4_com[iI] = fEI - fIR;
                        k4_com[iR] = fIR;
                    }

                    Xc += (dt / 6.0) * (k1_com + 2.0 * k2_com + 2.0 * k3_com + k4_com);
                }
                current_totals = next_totals;
            }
            benchmark::DoNotOptimize(mobile_pops);
        }
    }
}

static void bench_standard_lagrangian_rk4(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();
            mio::oseir::StandardModelLagrangian model_explicit(num_age_groups, num_commuter_groups);
            setup_explicit_model_benchmark(model_explicit, num_age_groups, num_commuter_groups);
            mio::oseir::StandardLagrangianSim sim_explicit(model_explicit, t0, dt);
            auto integrator_rk =
                std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
            sim_explicit.set_integrator(integrator_rk);
            state.ResumeTiming();

            sim_explicit.advance(t_max);
        }
    }
}

static void bench_standard_lagrangian_euler(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();
            mio::oseir::StandardModelLagrangian model_explicit(num_age_groups, num_commuter_groups);
            setup_explicit_model_benchmark(model_explicit, num_age_groups, num_commuter_groups);
            mio::oseir::StandardLagrangianSim sim_explicit(model_explicit, t0, dt);
            auto integrator_euler = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();
            sim_explicit.set_integrator(integrator_euler);
            state.ResumeTiming();

            sim_explicit.advance(t_max);
        }
    }
}

// ============================================================================
// OpenMP-parallel stage-aligned RK4: the outer patch loop runs in parallel.
//
// All NP = num_commuter_groups + 1 patches share the same model type
// and are independent of each other, so they can be computed simultaneously
// across all available CPU threads.
// ============================================================================
#ifdef _OPENMP
static void bench_stage_aligned_rk4_omp(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const int num_threads         = static_cast<int>(state.range(2));
    const int num_patches         = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);
    omp_set_num_threads(num_threads);

    // Per-patch runtime state
    struct PatchState {
        Eigen::VectorXd current_totals;
        std::vector<Eigen::VectorXd> mobile_pops;
        // RK4 stage vectors for the totals ODE
        Eigen::VectorXd y0, y1, y2, y3; // stage inputs
        Eigen::VectorXd k0, k1, k2, k3; // stage derivatives
        // lambda[stage][age_group]
        std::vector<std::vector<double>> lambda;
        // RK4 workspace for commuter ODE
        Eigen::VectorXd k1c, k2c, k3c, k4c;
        Eigen::VectorXd Xc2, Xc3, Xc4;
    };
    std::vector<PatchState> patches(static_cast<size_t>(num_patches));

    for (auto _ : state) {

        // ---- NOT measured: model setup and data initialisation ----
        state.PauseTiming();

        ModelType model(num_age_groups);
        setup_model_benchmark(model);
        const auto ic               = make_seir_index_cache(model);
        const Eigen::VectorXd proto = model.populations.get_compartments();
        const size_t NC             = static_cast<size_t>(proto.size());

        std::vector<double> rate_E(num_age_groups, 0.0);
        std::vector<double> rate_I(num_age_groups, 0.0);
        for (size_t g = 0; g < num_age_groups; ++g) {
            const double tE = model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(g)];
            const double tI = model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(g)];
            rate_E[g]       = (tE > 1e-10) ? 1.0 / tE : 0.0;
            rate_I[g]       = (tI > 1e-10) ? 1.0 / tI : 0.0;
        }

        const double mobile_fraction = 0.1 * num_commuter_groups;
        const Eigen::VectorXd init_mob =
            proto * (num_commuter_groups > 0 ? mobile_fraction / num_commuter_groups : 0.0);

#pragma omp parallel for schedule(static)
        for (int p = 0; p < num_patches; ++p) {
            PatchState& ps    = patches[static_cast<size_t>(p)];
            ps.current_totals = proto;
            ps.mobile_pops.assign(num_commuter_groups, init_mob);
            ps.lambda.assign(4, std::vector<double>(num_age_groups, 0.0));

            ps.y0.resize(NC);
            ps.y0.setZero();
            ps.y1.resize(NC);
            ps.y1.setZero();
            ps.y2.resize(NC);
            ps.y2.setZero();
            ps.y3.resize(NC);
            ps.y3.setZero();
            ps.k0.resize(NC);
            ps.k0.setZero();
            ps.k1.resize(NC);
            ps.k1.setZero();
            ps.k2.resize(NC);
            ps.k2.setZero();
            ps.k3.resize(NC);
            ps.k3.setZero();
            ps.k1c.resize(NC);
            ps.k1c.setZero();
            ps.k2c.resize(NC);
            ps.k2c.setZero();
            ps.k3c.resize(NC);
            ps.k3c.setZero();
            ps.k4c.resize(NC);
            ps.k4c.setZero();
            ps.Xc2.resize(NC);
            ps.Xc2.setZero();
            ps.Xc3.resize(NC);
            ps.Xc3.setZero();
            ps.Xc4.resize(NC);
            ps.Xc4.setZero();
        }

        // ---- Pre-compute contact*beta matrix ----
        // Doing this once before the parallel region avoids calling
        // get_matrix_at() inside the loop, where all 36 OMP threads would
        // simultaneously call it.
        std::vector<double> contact_beta(num_age_groups * num_age_groups, 0.0);
        for (size_t i = 0; i < num_age_groups; ++i) {
            const double beta_i =
                model.parameters.get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(i)];
            for (size_t j = 0; j < num_age_groups; ++j) {
                contact_beta[i * num_age_groups + j] =
                    model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>().get_cont_freq_mat().get_matrix_at(
                        0.0)(i, j) *
                    beta_i;
            }
        }

        auto compute_lambda = [&](const Eigen::VectorXd& y, std::vector<double>& lam_out) {
            for (size_t i = 0; i < num_age_groups; ++i)
                lam_out[i] = 0.0;
            for (size_t j = 0; j < num_age_groups; ++j) {
                const double Nj    = y[ic.S[j]] + y[ic.E[j]] + y[ic.I[j]] + y[ic.R[j]];
                const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;
                const double wj    = divNj * y[ic.I[j]];
                for (size_t i = 0; i < num_age_groups; ++i)
                    lam_out[i] += contact_beta[i * num_age_groups + j] * wj;
            }
        };

        // Computes lambda AND the totals RHS and
        // avoids calling model.get_derivatives() inside the parallel part.
        auto compute_rhs_and_lambda = [&](const Eigen::VectorXd& y, Eigen::VectorXd& k_out,
                                          std::vector<double>& lam_out) {
            compute_lambda(y, lam_out);
            for (size_t g = 0; g < num_age_groups; ++g) {
                const int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                const double fSE = lam_out[g] * y[iS];
                const double fEI = rate_E[g] * y[iE];
                const double fIR = rate_I[g] * y[iI];
                k_out[iS]        = -fSE;
                k_out[iE]        = fSE - fEI;
                k_out[iI]        = fEI - fIR;
                k_out[iR]        = fIR;
            }
        };

        // ---- MEASURED: parallel time integration over all patches ----
        state.ResumeTiming();
#ifdef LIKWID_PERFMON
        LIKWID_MARKER_START("omp_rk4_hot");
#endif

#pragma omp parallel for schedule(static)
        for (int p = 0; p < num_patches; ++p) {
            PatchState& ps = patches[static_cast<size_t>(p)];

            for (ScalarType t = t0; t < t_max; t += dt) {

                // --- Totals RK4 + lambda per stage (no heap allocations) ---
                ps.y0 = ps.current_totals;
                compute_rhs_and_lambda(ps.y0, ps.k0, ps.lambda[0]);

                ps.y1 = ps.y0 + (dt * 0.5) * ps.k0;
                compute_rhs_and_lambda(ps.y1, ps.k1, ps.lambda[1]);

                ps.y2 = ps.y0 + (dt * 0.5) * ps.k1;
                compute_rhs_and_lambda(ps.y2, ps.k2, ps.lambda[2]);

                ps.y3 = ps.y0 + dt * ps.k2;
                compute_rhs_and_lambda(ps.y3, ps.k3, ps.lambda[3]);

                const Eigen::VectorXd next_totals =
                    ps.current_totals + (dt / 6.0) * (ps.k0 + 2.0 * ps.k1 + 2.0 * ps.k2 + ps.k3);

                // --- Commuter RK4 ---
                for (int cg = 0; cg < num_commuter_groups; ++cg) {
                    Eigen::Ref<Eigen::VectorXd> Xc = ps.mobile_pops[static_cast<size_t>(cg)];

                    for (size_t g = 0; g < num_age_groups; ++g) {
                        const int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        const double fSE = ps.lambda[0][g] * Xc[iS];
                        const double fEI = rate_E[g] * Xc[iE];
                        const double fIR = rate_I[g] * Xc[iI];
                        ps.k1c[iS]       = -fSE;
                        ps.k1c[iE]       = fSE - fEI;
                        ps.k1c[iI]       = fEI - fIR;
                        ps.k1c[iR]       = fIR;
                    }

                    ps.Xc2 = Xc + (dt * 0.5) * ps.k1c;
                    for (size_t g = 0; g < num_age_groups; ++g) {
                        const int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        const double fSE = ps.lambda[1][g] * ps.Xc2[iS];
                        const double fEI = rate_E[g] * ps.Xc2[iE];
                        const double fIR = rate_I[g] * ps.Xc2[iI];
                        ps.k2c[iS]       = -fSE;
                        ps.k2c[iE]       = fSE - fEI;
                        ps.k2c[iI]       = fEI - fIR;
                        ps.k2c[iR]       = fIR;
                    }

                    ps.Xc3 = Xc + (dt * 0.5) * ps.k2c;
                    for (size_t g = 0; g < num_age_groups; ++g) {
                        const int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        const double fSE = ps.lambda[2][g] * ps.Xc3[iS];
                        const double fEI = rate_E[g] * ps.Xc3[iE];
                        const double fIR = rate_I[g] * ps.Xc3[iI];
                        ps.k3c[iS]       = -fSE;
                        ps.k3c[iE]       = fSE - fEI;
                        ps.k3c[iI]       = fEI - fIR;
                        ps.k3c[iR]       = fIR;
                    }

                    ps.Xc4 = Xc + dt * ps.k3c;
                    for (size_t g = 0; g < num_age_groups; ++g) {
                        const int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        const double fSE = ps.lambda[3][g] * ps.Xc4[iS];
                        const double fEI = rate_E[g] * ps.Xc4[iE];
                        const double fIR = rate_I[g] * ps.Xc4[iI];
                        ps.k4c[iS]       = -fSE;
                        ps.k4c[iE]       = fSE - fEI;
                        ps.k4c[iI]       = fEI - fIR;
                        ps.k4c[iR]       = fIR;
                    }

                    Xc += (dt / 6.0) * (ps.k1c + 2.0 * ps.k2c + 2.0 * ps.k3c + ps.k4c);
                }
                ps.current_totals = next_totals;
            }
        } // end omp parallel for
#ifdef LIKWID_PERFMON
        LIKWID_MARKER_STOP("omp_rk4_hot");
#endif

        benchmark::DoNotOptimize(patches);
    }
}
#endif // _OPENMP

static void bench_stage_aligned_euler(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();
            ModelType model(num_age_groups);
            setup_model_benchmark(model);

            Eigen::VectorXd current_totals = model.populations.get_compartments();

            double mobile_fraction = 0.1 * num_commuter_groups;
            Eigen::VectorXd initial_mobile =
                current_totals * (num_commuter_groups > 0 ? (mobile_fraction / num_commuter_groups) : 0.0);
            std::vector<Eigen::VectorXd> mobile_pops(num_commuter_groups, initial_mobile);

            const size_t NC = current_totals.size();
            auto ic         = mio::oseir::make_seir_index_cache(model);

            std::vector<double> rate_E(num_age_groups, 0.0), rate_I(num_age_groups, 0.0);
            for (size_t g = 0; g < num_age_groups; ++g) {
                double t_E =
                    model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(static_cast<int>(g))];
                double t_I =
                    model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(static_cast<int>(g))];
                rate_E[g] = (t_E > 1e-10) ? (1.0 / t_E) : 0.0;
                rate_I[g] = (t_I > 1e-10) ? (1.0 / t_I) : 0.0;
            }

            auto compute_lambda = [&](const Eigen::VectorXd& y, double tt, std::vector<double>& lambda_out) {
                lambda_out.assign(num_age_groups, 0.0);
                for (size_t j = 0; j < num_age_groups; ++j) {
                    const double Nj    = y[ic.S[j]] + y[ic.E[j]] + y[ic.I[j]] + y[ic.R[j]];
                    const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;
                    for (size_t i = 0; i < num_age_groups; ++i) {
                        const double coeff =
                            model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>()
                                .get_cont_freq_mat()
                                .get_matrix_at(tt)(i, j) *
                            model.parameters
                                .get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(
                                    static_cast<int>(i))] *
                            divNj;
                        lambda_out[i] += coeff * y[ic.I[j]];
                    }
                }
            };

            Eigen::VectorXd k_tot(NC);
            std::vector<double> current_lambda(num_age_groups, 0.0);

            state.ResumeTiming();

            for (ScalarType t = t0; t < t_max; t += dt) {

                // Solve totals with Euler
                model.get_derivatives(current_totals, current_totals, t, k_tot);
                compute_lambda(current_totals, t, current_lambda);
                Eigen::VectorXd next_totals = current_totals + dt * k_tot;

                // reconstruction commuter groups
                for (int cg = 0; cg < num_commuter_groups; ++cg) {
                    Eigen::Ref<Eigen::VectorXd> Xc = mobile_pops[cg];

                    for (size_t g = 0; g < num_age_groups; ++g) {
                        int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];

                        double fSE = current_lambda[g] * Xc[iS];
                        double fEI = rate_E[g] * Xc[iE];
                        double fIR = rate_I[g] * Xc[iI];

                        Xc[iS] += dt * (-fSE);
                        Xc[iE] += dt * (fSE - fEI);
                        Xc[iI] += dt * (fEI - fIR);
                        Xc[iR] += dt * (fIR);
                    }
                }
                current_totals = next_totals;
            }
            benchmark::DoNotOptimize(mobile_pops);
        }
    }
}

static void bench_stage_aligned_hybrid(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {

            state.PauseTiming();
            ModelType model(num_age_groups);
            setup_model_benchmark(model);
            SimType sim(model, t0, dt);

            auto integrator_rk =
                std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
            sim.set_integrator(integrator_rk);
            state.ResumeTiming();

            // 1. Totals mit RK-4 integrieren (MEmilio sim.advance)
            sim.advance(t_max);

            state.PauseTiming();

            const auto& seir_res = sim.get_result();
            const size_t n_steps = static_cast<size_t>(seir_res.get_num_time_points()) > 0
                                       ? static_cast<size_t>(seir_res.get_num_time_points()) - 1
                                       : 0;

            double mobile_fraction = 0.1 * num_commuter_groups;
            Eigen::VectorXd initial_mobile =
                seir_res.get_value(0) * (num_commuter_groups > 0 ? (mobile_fraction / num_commuter_groups) : 0.0);
            std::vector<Eigen::VectorXd> mobile_pops(static_cast<size_t>(num_commuter_groups), initial_mobile);

            auto ic = mio::oseir::make_seir_index_cache(model);

            std::vector<double> rate_E(num_age_groups, 0.0), rate_I(num_age_groups, 0.0);
            for (size_t g = 0; g < num_age_groups; ++g) {
                double t_E =
                    model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(static_cast<int>(g))];
                double t_I =
                    model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(static_cast<int>(g))];
                rate_E[g] = (t_E > 1e-10) ? (1.0 / t_E) : 0.0;
                rate_I[g] = (t_I > 1e-10) ? (1.0 / t_I) : 0.0;
            }

            auto compute_lambda = [&](const Eigen::VectorXd& y, double tt, std::vector<double>& lambda_out) {
                lambda_out.assign(num_age_groups, 0.0);
                for (size_t j = 0; j < num_age_groups; ++j) {
                    const double Nj    = y[ic.S[j]] + y[ic.E[j]] + y[ic.I[j]] + y[ic.R[j]];
                    const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;
                    for (size_t i = 0; i < num_age_groups; ++i) {
                        const double coeff =
                            model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>()
                                .get_cont_freq_mat()
                                .get_matrix_at(tt)(i, j) *
                            model.parameters
                                .get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(
                                    static_cast<int>(i))] *
                            divNj;
                        lambda_out[i] += coeff * y[ic.I[j]];
                    }
                }
            };

            std::vector<double> current_lambda(num_age_groups, 0.0);

            state.ResumeTiming();

            // 2. Reconstruction (Euler Update)
            for (size_t s = 0; s < n_steps; ++s) {
                const Eigen::VectorXd& tot_t = seir_res.get_value(s);

                double current_t  = seir_res.get_time(s);
                double current_dt = seir_res.get_time(s + 1) - current_t;

                compute_lambda(tot_t, current_t, current_lambda);

                for (int cg = 0; cg < num_commuter_groups; ++cg) {
                    Eigen::VectorXd& x = mobile_pops[static_cast<size_t>(cg)];

                    for (size_t g = 0; g < num_age_groups; ++g) {
                        const size_t iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];

                        const double fSE = current_lambda[g] * x[static_cast<Eigen::Index>(iS)];
                        const double fEI = rate_E[g] * x[static_cast<Eigen::Index>(iE)];
                        const double fIR = rate_I[g] * x[static_cast<Eigen::Index>(iI)];

                        x[static_cast<Eigen::Index>(iS)] += current_dt * (-fSE);
                        x[static_cast<Eigen::Index>(iE)] += current_dt * (fSE - fEI);
                        x[static_cast<Eigen::Index>(iI)] += current_dt * (fEI - fIR);
                        x[static_cast<Eigen::Index>(iR)] += current_dt * (fIR);
                    }
                }
            }

            benchmark::DoNotOptimize(mobile_pops);
        }
    }
}

template <class Integrator>
static void bench_matrix_phi_reconstruction(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {

            state.PauseTiming();
            ModelType model(num_age_groups);
            setup_model_benchmark(model);

            // Dimensions
            size_t NC          = (size_t)model.populations.get_num_compartments();
            size_t system_size = NC + NC * NC;

            // Initial state: y0 = [z0, Identity]
            Eigen::VectorXd y0(system_size);
            y0.head(NC)        = model.populations.get_compartments();
            Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(NC, NC);
            y0.tail(NC * NC)   = Eigen::Map<Eigen::VectorXd>(Id.data(), NC * NC);

            // Integrator Setup
            Integrator stepper;
            mio::oseir::AugmentedPhiSystem sys(model);

            // Commuter initialization (as a matrix)
            double mobile_fraction = 0.1 * num_commuter_groups;
            Eigen::VectorXd initial_mobile_pop =
                y0.head(NC) * (num_commuter_groups > 0 ? (mobile_fraction / num_commuter_groups) : 0.0);

            // Matrix X0: size (NC x num_commuter_groups)
            Eigen::MatrixXd X0(NC, num_commuter_groups);
            for (int cg = 0; cg < num_commuter_groups; ++cg) {
                X0.col(cg) = initial_mobile_pop;
            }

            // Result matrix (overwritten in each step)
            Eigen::MatrixXd Xt(NC, num_commuter_groups);

            state.ResumeTiming();

            // 1. Integrate (z and Phi) over the entire period t_max

            double t          = t0;
            Eigen::VectorXd y = y0;
            while (t < t_max_phi - 1e-10) {
                double dt_eff = std::min(dt_phi, t_max_phi - t);
                stepper.do_step(sys, y, t, dt_eff);
                t += dt_eff;
            }

            // 2. Reconstruction at tmax
            // Update all commuter groups: X(t) = Phi(t, t0) * X(t0)
            const auto Phi_final = Eigen::Map<const Eigen::MatrixXd>(y.tail(NC * NC).data(), NC, NC);
            Xt.noalias()         = Phi_final * X0;

            benchmark::DoNotOptimize(Xt);
        }
    }
}

} // namespace benchmark_mio
} // namespace mio

namespace mio
{
namespace benchmark_mio
{

/**
 * @brief Block-diagonal Augmented ODE system for Phi.
 * * Exploits the fact that the fundamental matrix Phi is completely decoupled between age groups.
 * State vector: y = [z (NC), vec(Phi_0) (16), vec(Phi_1) (16), ..., vec(Phi_G-1) (16)].
 * ODE size drops from NC + NC^2 to NC + 16 * G.
 */
struct AugmentedPhiSystemBlockDiag {
    const mio::oseir::Model<ScalarType>& model;
    size_t G;
    size_t NC;

    explicit AugmentedPhiSystemBlockDiag(const mio::oseir::Model<ScalarType>& m)
        : model(m)
        , G(static_cast<size_t>(m.parameters.get_num_groups()))
        , NC(static_cast<size_t>(m.populations.get_num_compartments()))
    {
    }

    void operator()(const Eigen::VectorXd& y, Eigen::VectorXd& dydt, double t)
    {
        const auto z = y.head(NC);

        // 1. dz/dt (standard SEIR RHS)
        Eigen::VectorXd dz(NC);
        model.get_derivatives(z, z, t, dz);
        dydt.head(NC) = dz;

        // 2. Compute force of infection for each age group
        using IS            = mio::oseir::InfectionState;
        Eigen::VectorXd foi = Eigen::VectorXd::Zero(G);
        for (size_t j = 0; j < G; ++j) {
            const size_t Sj = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Susceptible});
            const size_t Ej = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Exposed});
            const size_t Ij = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Infected});
            const size_t Rj = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Recovered});

            const double Nj    = z[Sj] + z[Ej] + z[Ij] + z[Rj];
            const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;

            for (size_t i = 0; i < G; ++i) {
                const double coeff =
                    model.parameters.template get<mio::oseir::ContactPatterns<ScalarType>>()
                        .get_cont_freq_mat()
                        .get_matrix_at(t)(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) *
                    model.parameters
                        .template get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(
                            static_cast<int>(i))] *
                    divNj;
                foi[i] += coeff * z[Ij];
            }
        }

        // 3. Compute dPhi_g/dt for each isolated age group block
        for (size_t g = 0; g < G; ++g) {
            const double time_exposed =
                model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(static_cast<int>(g))];
            const double time_infected =
                model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(static_cast<int>(g))];

            const double rate_E = (time_exposed > 1e-10) ? (1.0 / time_exposed) : 0.0;
            const double rate_I = (time_infected > 1e-10) ? (1.0 / time_infected) : 0.0;
            const double lambda = foi[g];

            // Build 4x4 transition matrix A_g for this age group
            Eigen::Matrix4d Ag = Eigen::Matrix4d::Zero();
            Ag(0, 0)           = -lambda;
            Ag(1, 0)           = lambda;
            Ag(1, 1)           = -rate_E;
            Ag(2, 1)           = rate_E;
            Ag(2, 2)           = -rate_I;
            Ag(3, 2)           = rate_I;

            // Extract current 4x4 Phi block from the state vector
            const auto Phi_g = Eigen::Map<const Eigen::Matrix4d>(y.data() + NC + g * 16);

            // Compute derivative and map directly to dydt
            Eigen::Map<Eigen::Matrix4d>(dydt.data() + NC + g * 16) = Ag * Phi_g;
        }
    }
};

template <class Integrator>
static void bench_matrix_phi_reconstruction_blockdiag(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();
            ModelType model(num_age_groups);
            setup_model_benchmark(model);

            size_t G  = num_age_groups;
            size_t NC = 4 * G;

            // NC + G * 16
            size_t system_size = NC + G * 16;

            Eigen::VectorXd y0(system_size);
            y0.head(NC) = model.populations.get_compartments();

            // Initialize the diagonal blocks with Identity
            for (size_t g = 0; g < G; ++g) {
                Eigen::Map<Eigen::Matrix4d>(y0.data() + NC + g * 16) = Eigen::Matrix4d::Identity();
            }

            Integrator stepper;
            AugmentedPhiSystemBlockDiag sys(model);

            double mobile_fraction = 0.1 * num_commuter_groups;
            Eigen::VectorXd initial_mobile_pop =
                y0.head(NC) * (num_commuter_groups > 0 ? (mobile_fraction / num_commuter_groups) : 0.0);

            Eigen::MatrixXd X0(NC, num_commuter_groups);
            for (int cg = 0; cg < num_commuter_groups; ++cg) {
                X0.col(cg) = initial_mobile_pop;
            }
            Eigen::MatrixXd Xt(NC, num_commuter_groups);

            state.ResumeTiming();

            // 1. Integrate over the entire period t_max
            double t          = t0;
            Eigen::VectorXd y = y0;
            while (t < t_max_phi - 1e-10) {
                double dt_eff = std::min(dt_phi, t_max_phi - t);
                stepper.do_step(sys, y, t, dt_eff);
                t += dt_eff;
            }

            // 2. Block-Diagonal Reconstruction
            for (size_t g = 0; g < G; ++g) {
                const auto Phi_g = Eigen::Map<const Eigen::Matrix4d>(y.data() + NC + g * 16);
                // MEmilio stores age groups sequentially: S_g, E_g, I_g, R_g are contiguous (size 4)
                Xt.block(4 * g, 0, 4, num_commuter_groups).noalias() =
                    Phi_g * X0.block(4 * g, 0, 4, num_commuter_groups);
            }

            benchmark::DoNotOptimize(Xt);
        }
    }
}

} // namespace benchmark_mio
} // namespace mio

// static void bench_matrix_phi_reconstruction_blockdiag_rk4(benchmark::State& state)
// {
//     mio::benchmark_mio::bench_matrix_phi_reconstruction_blockdiag<boost::numeric::odeint::runge_kutta4<
//         Eigen::VectorXd, double, Eigen::VectorXd, double, boost::numeric::odeint::vector_space_algebra>>(state);
// }

// BENCHMARK(bench_matrix_phi_reconstruction_blockdiag_rk4)
//     ->Apply([](auto* b) {
//         for (int i : mio::benchmark_mio::commuter_group_counts) {
//             for (int g : mio::benchmark_mio::age_group_counts) {
//                 b->Args({i, g});
//             }
//         }
//     }) -> Name("matrix_phi_blockdiag(RK4)") -> Unit(::benchmark::kMicrosecond);

// static void bench_matrix_phi_reconstruction_blockdiag_euler(benchmark::State& state)
// {
//     mio::benchmark_mio::bench_matrix_phi_reconstruction_blockdiag<boost::numeric::odeint::euler<
//         Eigen::VectorXd, double, Eigen::VectorXd, double, boost::numeric::odeint::vector_space_algebra>>(state);
// }

// BENCHMARK(bench_matrix_phi_reconstruction_blockdiag_euler)
//     ->Apply([](auto* b) {
//         for (int i : mio::benchmark_mio::commuter_group_counts) {
//             for (int g : mio::benchmark_mio::age_group_counts) {
//                 b->Args({i, g});
//             }
//         }
//     }) -> Name("matrix_phi_blockdiag(Euler)") -> Unit(::benchmark::kMicrosecond);

// static void bench_matrix_phi_reconstruction_rk4(benchmark::State& state)
// {
//     mio::benchmark_mio::bench_matrix_phi_reconstruction<boost::numeric::odeint::runge_kutta4<
//         Eigen::VectorXd, double, Eigen::VectorXd, double, boost::numeric::odeint::vector_space_algebra>>(state);
// }

// BENCHMARK(bench_matrix_phi_reconstruction_rk4)
//     ->Apply([](auto* b) {
//         for (int i : mio::benchmark_mio::commuter_group_counts) {
//             for (int g : mio::benchmark_mio::age_group_counts) {
//                 b->Args({i, g});
//             }
//         }
//     }) -> Name("matrix_phi_reconstruction(RK4)") -> Unit(::benchmark::kMicrosecond);

// static void bench_matrix_phi_reconstruction_euler(benchmark::State& state)
// {
//     mio::benchmark_mio::bench_matrix_phi_reconstruction<boost::numeric::odeint::euler<
//         Eigen::VectorXd, double, Eigen::VectorXd, double, boost::numeric::odeint::vector_space_algebra>>(state);
// }

// BENCHMARK(bench_matrix_phi_reconstruction_euler)
//     ->Apply([](auto* b) {
//         for (int i : mio::benchmark_mio::commuter_group_counts) {
//             for (int g : mio::benchmark_mio::age_group_counts) {
//                 b->Args({i, g});
//             }
//         }
//     }) -> Name("matrix_phi_reconstruction(Euler)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_stage_aligned_rk4)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts) {
            for (int g : mio::benchmark_mio::age_group_counts) {
                b->Args({i, g});
            }
        }
    }) -> Name("stage-aligned(RK4)") -> Unit(::benchmark::kMicrosecond);

#ifdef _OPENMP
BENCHMARK(mio::benchmark_mio::bench_stage_aligned_rk4_omp)
    ->Apply([](auto* b) {
        for (int t : mio::benchmark_mio::core_counts) {
            for (int i : mio::benchmark_mio::commuter_group_counts) {
                for (int g : mio::benchmark_mio::age_group_counts) {
                    b->Args({i, g, t});
                }
            }
        }
    }) -> Name("stage-aligned-omp(RK4)") -> Unit(::benchmark::kMicrosecond);
#endif

// BENCHMARK(mio::benchmark_mio::bench_stage_aligned_euler)
//     ->Apply([](auto* b) {
//         for (int i : mio::benchmark_mio::commuter_group_counts) {
//             for (int g : mio::benchmark_mio::age_group_counts) {
//                 b->Args({i, g});
//             }
//         }
//     }) -> Name("stage-aligned(Euler)") -> Unit(::benchmark::kMicrosecond);

// BENCHMARK(mio::benchmark_mio::bench_auxiliary_euler)
//     ->Apply([](auto* b) {
//         for (int i : mio::benchmark_mio::commuter_group_counts) {
//             for (int g : mio::benchmark_mio::age_group_counts) {
//                 b->Args({i, g});
//             }
//         }
//     }) -> Name("auxiliary_Euler") -> Unit(::benchmark::kMicrosecond);

// BENCHMARK(mio::benchmark_mio::bench_stage_aligned_hybrid)
//     ->Apply([](auto* b) {
//         for (int i : mio::benchmark_mio::commuter_group_counts) {
//             for (int g : mio::benchmark_mio::age_group_counts) {
//                 b->Args({i, g});
//             }
//         }
//     }) -> Name("stage-aligned(hybrid)") -> Unit(::benchmark::kMicrosecond);

// BENCHMARK(mio::benchmark_mio::bench_standard_lagrangian_rk4)
//     ->Apply([](auto* b) {
//         for (int i : mio::benchmark_mio::commuter_group_counts) {

//             for (int g : mio::benchmark_mio::age_group_counts) {
//                 b->Args({i, g});
//             }
//         }
//     }) -> Name("lagrange_rk4") -> Unit(::benchmark::kMicrosecond);

// BENCHMARK(mio::benchmark_mio::bench_standard_lagrangian_euler)
//     ->Apply([](auto* b) {
//         for (int i : mio::benchmark_mio::commuter_group_counts) {
//             for (int g : mio::benchmark_mio::age_group_counts) {
//                 b->Args({i, g});
//             }
//         }
//     }) -> Name("lagrange_euler") -> Unit(::benchmark::kMicrosecond);

// run all benchmarks
// When compiled with LIKWID Marker API (-DLIKWID_PERFMON), replace the
// standard BENCHMARK_MAIN() macro with a custom main that calls
// LIKWID_MARKER_INIT / LIKWID_MARKER_CLOSE around the benchmark runner.
// Without LIKWID the behaviour is identical to BENCHMARK_MAIN().
#ifdef LIKWID_PERFMON
int main(int argc, char** argv)
{
    LIKWID_MARKER_INIT;
#pragma omp parallel
    {
        LIKWID_MARKER_THREADINIT;
    }
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
