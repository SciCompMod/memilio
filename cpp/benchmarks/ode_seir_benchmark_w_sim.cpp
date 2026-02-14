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
#include "examples/state_estimators.h"
#include "examples/standard_lagrangian.h"

namespace mio
{
namespace examples
{

struct SEIRIndexCache {
    std::vector<int> S, E, I, R; // flat indices in totals/commuter vectors
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

struct StageCacheRK4 {
    double dt = 0.0;
    // totals at RK4 stages y1 (=t), y2 (=t+dt/2), y3 (=t+dt/2), y4 (=t+dt)
    Eigen::VectorXd y1, y2, y3, y4; // length = 4*G
    // totals flows at RK4 stages: per group [f_SE, f_EI, f_IR], length = 3*G
    Eigen::VectorXd f1, f2, f3, f4; // length = 3*G
};

template <class SimType>
inline StageCacheRK4 make_stage_cache_rk4(const SimType& sim, size_t step_index)
{
    // Totals an Schrittgrenzen
    const auto& res = sim.get_result();
    const double t  = res.get_time(step_index);
    const double tp = res.get_time(step_index + 1);
    const double dt = tp - t;

    const Eigen::VectorXd y1 = res.get_value(step_index);
    const size_t NC          = static_cast<size_t>(y1.size());

    const auto& model = sim.get_model();

    // Hilfs-Lambdas für RHS/Flows der Totals
    auto totals_rhs = [&](const Eigen::VectorXd& y, double tt) {
        Eigen::VectorXd dy(NC);
        model.get_derivatives(y, y, tt, dy);
        return dy;
    };
    auto totals_flows = [&](const Eigen::VectorXd& y, double tt) {
        Eigen::VectorXd f(3 * (size_t)model.parameters.get_num_groups());
        f.setZero();
        model.get_flows(y, y, tt, f);
        return f;
    };

    // RK4-Stufen für Totals aufbauen (nur einmal pro Schritt)
    const Eigen::VectorXd k1 = totals_rhs(y1, t);
    const Eigen::VectorXd y2 = y1 + (dt * 0.5) * k1;
    const Eigen::VectorXd k2 = totals_rhs(y2, t + 0.5 * dt);
    const Eigen::VectorXd y3 = y1 + (dt * 0.5) * k2;
    const Eigen::VectorXd k3 = totals_rhs(y3, t + 0.5 * dt);
    const Eigen::VectorXd y4 = y1 + dt * k3;

    StageCacheRK4 cache;
    cache.dt = dt;
    cache.y1 = y1;
    cache.y2 = y2;
    cache.y3 = y3;
    cache.y4 = y4;

    cache.f1 = totals_flows(cache.y1, t);
    cache.f2 = totals_flows(cache.y2, t + 0.5 * dt);
    cache.f3 = totals_flows(cache.y3, t + 0.5 * dt);
    cache.f4 = totals_flows(cache.y4, t + dt);

    return cache;
}

/**
 * @brief Stagewise Left-Share Update (klassisches RK4, 4 Stufen) aus StageCache.
 *
 * Commuter-Ableitungen werden JE STUFE ausschließlich durch Left-Share der Totals-Flüsse
 * f_SE, f_EI, f_IR berechnet. Die Stufen-Anteile xi berechnen sich elementweise als
 * xi = Xc_stage / Xt_stage (sicherheits-halber in [0,1] geklemmt).
 *
 * Layout: commuter & totals sind beide 4*G lang in Reihenfolge [S,E,I,R] je Altersgruppe.
 */
inline void
flow_based_mobility_returns_leftshare_rk4_from_cache(Eigen::Ref<Eigen::VectorXd> commuter, // Xc(t) -> Xc(t+dt)
                                                     const SEIRIndexCache& ic, // Indizes S,E,I,R je Gruppe
                                                     const StageCacheRK4& st) // Totals-Stufenzustände & -Flüsse
{
    const size_t G  = ic.G;
    const size_t NC = 4 * G;

    // schnelle Zugriffe
    const auto& y1  = st.y1;
    const auto& y2  = st.y2;
    const auto& y3  = st.y3;
    const auto& y4  = st.y4;
    const auto& f1  = st.f1;
    const auto& f2  = st.f2;
    const auto& f3  = st.f3;
    const auto& f4  = st.f4;
    const double dt = st.dt;

    auto clamp01 = [](double v) {
        return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
    };

    // Stufen-Ableitungen für Commuter
    Eigen::VectorXd k1 = Eigen::VectorXd::Zero(NC);
    Eigen::VectorXd k2 = Eigen::VectorXd::Zero(NC);
    Eigen::VectorXd k3 = Eigen::VectorXd::Zero(NC);
    Eigen::VectorXd k4 = Eigen::VectorXd::Zero(NC);

    // ---- Stage 1 ----
    {
        // xi1 = Xc0 / y1
        // und k1 nur aus f1 via Left-Share
        for (size_t g = 0; g < G; ++g) {
            const int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
            const int bf = 3 * static_cast<int>(g);

            const double xiS = (y1[iS] > 0.0) ? clamp01(commuter[iS] / y1[iS]) : 0.0;
            const double xiE = (y1[iE] > 0.0) ? clamp01(commuter[iE] / y1[iE]) : 0.0;
            const double xiI = (y1[iI] > 0.0) ? clamp01(commuter[iI] / y1[iI]) : 0.0;

            const double dSE = f1[bf + 0];
            const double dEI = f1[bf + 1];
            const double dIR = f1[bf + 2];

            k1[iS] = -xiS * dSE;
            k1[iE] = xiS * dSE - xiE * dEI;
            k1[iI] = xiE * dEI - xiI * dIR;
            k1[iR] = xiI * dIR;
        }
    }

    // ---- Stage 2 ----
    Eigen::VectorXd Xc2 = commuter + (dt * 0.5) * k1;
    {
        for (size_t g = 0; g < G; ++g) {
            const int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
            const int bf = 3 * static_cast<int>(g);

            const double xiS = (y2[iS] > 0.0) ? clamp01(Xc2[iS] / y2[iS]) : 0.0;
            const double xiE = (y2[iE] > 0.0) ? clamp01(Xc2[iE] / y2[iE]) : 0.0;
            const double xiI = (y2[iI] > 0.0) ? clamp01(Xc2[iI] / y2[iI]) : 0.0;

            const double dSE = f2[bf + 0];
            const double dEI = f2[bf + 1];
            const double dIR = f2[bf + 2];

            k2[iS] = -xiS * dSE;
            k2[iE] = xiS * dSE - xiE * dEI;
            k2[iI] = xiE * dEI - xiI * dIR;
            k2[iR] = xiI * dIR;
        }
    }

    // ---- Stage 3 ----
    Eigen::VectorXd Xc3 = commuter + (dt * 0.5) * k2;
    {
        for (size_t g = 0; g < G; ++g) {
            const int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
            const int bf = 3 * static_cast<int>(g);

            const double xiS = (y3[iS] > 0.0) ? clamp01(Xc3[iS] / y3[iS]) : 0.0;
            const double xiE = (y3[iE] > 0.0) ? clamp01(Xc3[iE] / y3[iE]) : 0.0;
            const double xiI = (y3[iI] > 0.0) ? clamp01(Xc3[iI] / y3[iI]) : 0.0;

            const double dSE = f3[bf + 0];
            const double dEI = f3[bf + 1];
            const double dIR = f3[bf + 2];

            k3[iS] = -xiS * dSE;
            k3[iE] = xiS * dSE - xiE * dEI;
            k3[iI] = xiE * dEI - xiI * dIR;
            k3[iR] = xiI * dIR;
        }
    }

    // ---- Stage 4 ----
    Eigen::VectorXd Xc4 = commuter + dt * k3;
    {
        for (size_t g = 0; g < G; ++g) {
            const int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
            const int bf = 3 * static_cast<int>(g);

            const double xiS = (y4[iS] > 0.0) ? clamp01(Xc4[iS] / y4[iS]) : 0.0;
            const double xiE = (y4[iE] > 0.0) ? clamp01(Xc4[iE] / y4[iE]) : 0.0;
            const double xiI = (y4[iI] > 0.0) ? clamp01(Xc4[iI] / y4[iI]) : 0.0;

            const double dSE = f4[bf + 0];
            const double dEI = f4[bf + 1];
            const double dIR = f4[bf + 2];

            k4[iS] = -xiS * dSE;
            k4[iE] = xiS * dSE - xiE * dEI;
            k4[iI] = xiE * dEI - xiI * dIR;
            k4[iR] = xiI * dIR;
        }
    }

    // RK4-Kombination
    Eigen::VectorXd Xc_next = commuter + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

    // Numerische Sicherungen
    for (int i = 0; i < (int)NC; ++i) {
        if (Xc_next[i] < 0.0)
            Xc_next[i] = 0.0;
        // Optional: keine strikte Begrenzung auf totals(t+dt) hier, weil totals extern geführt werden.
        // Falls du streng Xc<=Totals(t+dt) willst, gib totals(t+dt) in den Cache und klemm hier.
    }

    commuter = Xc_next;
}

} // namespace examples

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

void setup_explicit_model_benchmark(mio::examples::ExplicitModel& model_explicit, size_t num_age_groups,
                                    int num_commuter_groups, const double total_population = 10000.0)
{

    // Start with non commuting pop (50%)
    for (auto i = 0; i < static_cast<int>(num_age_groups); i++) {
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::E}] =
            50.0 / static_cast<double>(num_age_groups);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::I}] =
            50.0 / static_cast<double>(num_age_groups);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::R}] =
            50.0 / static_cast<double>(num_age_groups);
        model_explicit
            .populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType, mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::S}] =
            (total_population * 0.5) / static_cast<double>(num_age_groups) -
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::E}] -
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::I}] -
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), mio::examples::CommuterType::NonCommuter, mio::examples::InfectionStateExplicit::R}];
    }

    // Then the remaining commuting pop
    for (int c = 0; c < num_commuter_groups; c++) {
        auto commuter_type =
            static_cast<mio::examples::CommuterType>(static_cast<int>(mio::examples::CommuterType::CommuterBase) + c);
        double commuter_fraction = 0.5 / static_cast<double>(num_commuter_groups);

        for (auto i = 0; i < static_cast<int>(num_age_groups); i++) {
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::E}] =
                50.0 / static_cast<double>(num_age_groups * num_commuter_groups);
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::I}] =
                50.0 / static_cast<double>(num_age_groups * num_commuter_groups);
            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::R}] =
                50.0 / static_cast<double>(num_age_groups * num_commuter_groups);

            double exposed   = model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                                   mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::E}];
            double infected  = model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                                    mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::I}];
            double recovered = model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                                     mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::R}];

            model_explicit.populations[mio::Index<mio::AgeGroup, mio::examples::CommuterType,
                                                  mio::examples::InfectionStateExplicit>{
                mio::AgeGroup(i), commuter_type, mio::examples::InfectionStateExplicit::S}] =
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

const ScalarType t0    = 0.0;
const ScalarType t_max = 50.0;
const ScalarType dt    = 0.5;

// Define the number of commuter groups for the benchmark
const std::vector<int> commuter_group_counts = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

static void bench_auxiliary_euler(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
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

                // Für jede Pendlergruppe die Euler-Methode anwenden
                for (int i = 0; i < num_commuter_groups; ++i) {
                    mio::examples::integrate_mobile_population_euler(mobile_pops[i], sim, total_pop, t, dt);
                }
            }
            benchmark::DoNotOptimize(mobile_pops); // Prevent mobile_pops calculation from being optimized away
        }
    }
}

static void bench_stage_aligned_rk4_neu(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    // Hilfs-Struct für den Cache (um Allokationen in der Loop zu sparen)
    struct RK4StageCache {
        // Totals an den Stufen: y1(t), y2(t+h/2), y3(t+h/2), y4(t+h)
        std::vector<Eigen::VectorXd> y;
        // Flows an den Stufen: f1, f2, f3, f4
        std::vector<Eigen::VectorXd> flows;
        // Totals-Ableitungen (k_tot): k1, k2, k3, k4
        std::vector<Eigen::VectorXd> k_tot;

        void resize(size_t nc, size_t flow_dim)
        {
            y.resize(4, Eigen::VectorXd(nc));
            flows.resize(4, Eigen::VectorXd(flow_dim));
            k_tot.resize(4, Eigen::VectorXd(nc));
        }
    };

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();
            ModelType model(num_age_groups);
            setup_model_benchmark(model);

            // Wir nutzen hier keinen Standard-Integrator, sondern steppen manuell,
            // um Zugriff auf die Zwischenstufen zu haben.
            // Startwerte
            Eigen::VectorXd current_totals = model.populations.get_compartments();

            // Pendler Init
            double mobile_fraction = 0.1 * num_commuter_groups;
            Eigen::VectorXd initial_mobile =
                current_totals * (num_commuter_groups > 0 ? (mobile_fraction / num_commuter_groups) : 0.0);
            std::vector<Eigen::VectorXd> mobile_pops(num_commuter_groups, initial_mobile);

            // Cache einmalig allokieren
            const size_t NC       = current_totals.size();
            const size_t flow_dim = 3 * num_age_groups;
            RK4StageCache cache;
            cache.resize(NC, flow_dim);

            // Pre-Compute Indizes für Algebra-Loop
            auto ic = mio::examples::make_seir_index_cache(model);

            // Temporäre Vektoren für Commuter Updates (um new/delete zu vermeiden)
            Eigen::VectorXd k1_com(NC), k2_com(NC), k3_com(NC), k4_com(NC);
            Eigen::VectorXd Xc2(NC), Xc3(NC), Xc4(NC);

            state.ResumeTiming();

            // Zeitschleife (Manuelles RK4 Stepping)
            for (ScalarType t = t0; t < t_max; t += dt) {

                // ---------------------------------------------------------
                // 1. PHASE: PHYSIK (Teuer, aber nur 1x pro Zeitschritt)
                // ---------------------------------------------------------

                // Stage 1
                cache.y[0] = current_totals;
                model.get_derivatives(cache.y[0], cache.y[0], t, cache.k_tot[0]);
                model.get_flows(cache.y[0], cache.y[0], t, cache.flows[0]);

                // Stage 2
                cache.y[1] = cache.y[0] + (dt * 0.5) * cache.k_tot[0];
                model.get_derivatives(cache.y[1], cache.y[1], t + 0.5 * dt, cache.k_tot[1]);
                model.get_flows(cache.y[1], cache.y[1], t + 0.5 * dt, cache.flows[1]);

                // Stage 3
                cache.y[2] = cache.y[0] + (dt * 0.5) * cache.k_tot[1];
                model.get_derivatives(cache.y[2], cache.y[2], t + 0.5 * dt, cache.k_tot[2]);
                model.get_flows(cache.y[2], cache.y[2], t + 0.5 * dt, cache.flows[2]);

                // Stage 4
                cache.y[3] = cache.y[0] + dt * cache.k_tot[2];
                model.get_derivatives(cache.y[3], cache.y[3], t + dt, cache.k_tot[3]);
                model.get_flows(cache.y[3], cache.y[3], t + dt, cache.flows[3]);

                // Totals Update für nächsten Schritt vorbereiten
                Eigen::VectorXd next_totals = current_totals + (dt / 6.0) * (cache.k_tot[0] + 2 * cache.k_tot[1] +
                                                                             2 * cache.k_tot[2] + cache.k_tot[3]);

                // ---------------------------------------------------------
                // 2. PHASE: ALGEBRA (Billig, N-mal pro Zeitschritt)
                // ---------------------------------------------------------
                // Hier passiert der Speedup: Wir nutzen die Cache-Werte wieder

                for (int cg = 0; cg < num_commuter_groups; ++cg) {
                    Eigen::Ref<Eigen::VectorXd> Xc = mobile_pops[cg];

                    // -- Stage 1 Update --
                    // Berechne k1_com basierend auf cache.flows[0] und Shares (Xc / cache.y[0])
                    // (Logik analog zu commuter_rhs_from_shares, aber Inline für Speed)
                    for (size_t g = 0; g < num_age_groups; ++g) {
                        int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        int off = 3 * (int)g;
                        // Share berechnen
                        double S_tot = cache.y[0][iS];
                        double xiS   = (S_tot > 1e-10) ? (Xc[iS] / S_tot) : 0.0;
                        double E_tot = cache.y[0][iE];
                        double xiE   = (E_tot > 1e-10) ? (Xc[iE] / E_tot) : 0.0;
                        double I_tot = cache.y[0][iI];
                        double xiI   = (I_tot > 1e-10) ? (Xc[iI] / I_tot) : 0.0;

                        // Ableitungen
                        k1_com[iS] = -xiS * cache.flows[0][off];
                        k1_com[iE] = xiS * cache.flows[0][off] - xiE * cache.flows[0][off + 1];
                        k1_com[iI] = xiE * cache.flows[0][off + 1] - xiI * cache.flows[0][off + 2];
                        k1_com[iR] = xiI * cache.flows[0][off + 2];
                    }

                    // -- Stage 2 Update --
                    Xc2 = Xc + (dt * 0.5) * k1_com; // Predictor
                    for (size_t g = 0; g < num_age_groups; ++g) {
                        int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        int off    = 3 * (int)g;
                        double xiS = (cache.y[1][iS] > 1e-10) ? (Xc2[iS] / cache.y[1][iS]) : 0.0;
                        double xiE = (cache.y[1][iE] > 1e-10) ? (Xc2[iE] / cache.y[1][iE]) : 0.0;
                        double xiI = (cache.y[1][iI] > 1e-10) ? (Xc2[iI] / cache.y[1][iI]) : 0.0;

                        k2_com[iS] = -xiS * cache.flows[1][off];
                        k2_com[iE] = xiS * cache.flows[1][off] - xiE * cache.flows[1][off + 1];
                        k2_com[iI] = xiE * cache.flows[1][off + 1] - xiI * cache.flows[1][off + 2];
                        k2_com[iR] = xiI * cache.flows[1][off + 2];
                    }

                    // -- Stage 3 Update --
                    Xc3 = Xc + (dt * 0.5) * k2_com;
                    for (size_t g = 0; g < num_age_groups; ++g) {
                        int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        int off    = 3 * (int)g;
                        double xiS = (cache.y[2][iS] > 1e-10) ? (Xc3[iS] / cache.y[2][iS]) : 0.0;
                        double xiE = (cache.y[2][iE] > 1e-10) ? (Xc3[iE] / cache.y[2][iE]) : 0.0;
                        double xiI = (cache.y[2][iI] > 1e-10) ? (Xc3[iI] / cache.y[2][iI]) : 0.0;

                        k3_com[iS] = -xiS * cache.flows[2][off];
                        k3_com[iE] = xiS * cache.flows[2][off] - xiE * cache.flows[2][off + 1];
                        k3_com[iI] = xiE * cache.flows[2][off + 1] - xiI * cache.flows[2][off + 2];
                        k3_com[iR] = xiI * cache.flows[2][off + 2];
                    }

                    // -- Stage 4 Update --
                    Xc4 = Xc + dt * k3_com;
                    for (size_t g = 0; g < num_age_groups; ++g) {
                        int iS = ic.S[g], iE = ic.E[g], iI = ic.I[g], iR = ic.R[g];
                        int off    = 3 * (int)g;
                        double xiS = (cache.y[3][iS] > 1e-10) ? (Xc4[iS] / cache.y[3][iS]) : 0.0;
                        double xiE = (cache.y[3][iE] > 1e-10) ? (Xc4[iE] / cache.y[3][iE]) : 0.0;
                        double xiI = (cache.y[3][iI] > 1e-10) ? (Xc4[iI] / cache.y[3][iI]) : 0.0;

                        k4_com[iS] = -xiS * cache.flows[3][off];
                        k4_com[iE] = xiS * cache.flows[3][off] - xiE * cache.flows[3][off + 1];
                        k4_com[iI] = xiE * cache.flows[3][off + 1] - xiI * cache.flows[3][off + 2];
                        k4_com[iR] = xiI * cache.flows[3][off + 2];
                    }

                    // Final Update Commuter
                    Xc = Xc + (dt / 6.0) * (k1_com + 2.0 * k2_com + 2.0 * k3_com + k4_com);

                    // Simple numeric guard (optional for benchmark but fair)
                    Xc = Xc.cwiseMax(0.0).cwiseMin(next_totals);
                }

                // Advance Totals
                current_totals = next_totals;
            }
            benchmark::DoNotOptimize(mobile_pops);
        }
    }
}

static void bench_stage_aligned_rk4(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
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
            // state.ResumeTiming();

            // sim.advance(t_max); // Measure advance

            // state.PauseTiming();
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

                const auto& total_pop_const = seir_res.get_value(closest_idx_total);
                const auto& flows           = sim.get_flows().get_value(closest_idx_total);

                Eigen::VectorXd totals = total_pop_const;

                for (int i = 0; i < num_commuter_groups; ++i) {
                    // pass (commuter, totals, model, t, dt)
                    mio::examples::flow_based_mobility_returns_rk4(mobile_pops[i], totals, sim.get_model(), t, dt);
                }
            }
            benchmark::DoNotOptimize(mobile_pops);
        }
    }
}

static void bench_standard_lagrangian_rk4(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();
            mio::examples::StandardModelLagrangian model_explicit(num_age_groups, num_commuter_groups);
            setup_explicit_model_benchmark(model_explicit, num_age_groups, num_commuter_groups);
            mio::examples::StandardLagrangianSim sim_explicit(model_explicit, t0, dt);
            auto integrator_rk =
                std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
            sim_explicit.set_integrator(integrator_rk);
            state.ResumeTiming();

            // Only measure the advance call
            sim_explicit.advance(t_max);
        }
    }
}

static void bench_standard_lagrangian_euler(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();
            mio::examples::StandardModelLagrangian model_explicit(num_age_groups, num_commuter_groups);
            setup_explicit_model_benchmark(model_explicit, num_age_groups, num_commuter_groups);
            mio::examples::StandardLagrangianSim sim_explicit(model_explicit, t0, dt);
            auto integrator_euler = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();
            sim_explicit.set_integrator(integrator_euler);
            state.ResumeTiming();

            // Only measure the advance call
            sim_explicit.advance(t_max);
        }
    }
}

static void bench_stage_aligned_euler(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {

            state.PauseTiming();
            ModelType model(num_age_groups);
            setup_model_benchmark(model);
            SimType sim(model, t0, dt);
            auto integrator_euler = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();
            sim.set_integrator(integrator_euler);
            state.ResumeTiming();

            sim.advance(t_max);

            state.PauseTiming();

            const auto& seir_res = sim.get_result();
            const auto& flows_ts = sim.get_flows();
            const size_t n_steps = static_cast<size_t>(seir_res.get_num_time_points()) > 0
                                       ? static_cast<size_t>(seir_res.get_num_time_points()) - 1
                                       : 0;

            double mobile_fraction = 0.1 * num_commuter_groups;
            Eigen::VectorXd initial_mobile =
                seir_res.get_value(0) * (num_commuter_groups > 0 ? (mobile_fraction / num_commuter_groups) : 0.0);
            std::vector<Eigen::VectorXd> mobile_pops(static_cast<size_t>(num_commuter_groups), initial_mobile);

            using IS = mio::oseir::InfectionState;
            std::vector<size_t> idxS(num_age_groups), idxE(num_age_groups), idxI(num_age_groups), idxR(num_age_groups);
            for (size_t g = 0; g < num_age_groups; ++g) {
                auto G  = mio::AgeGroup(static_cast<int>(g));
                idxS[g] = model.populations.get_flat_index({G, IS::Susceptible});
                idxE[g] = model.populations.get_flat_index({G, IS::Exposed});
                idxI[g] = model.populations.get_flat_index({G, IS::Infected});
                idxR[g] = model.populations.get_flat_index({G, IS::Recovered});
            }

            state.ResumeTiming();

            // 3) Pendler-Update je Schritt: rein algebraisch, left-share
            for (size_t s = 0; s < n_steps; ++s) {
                const Eigen::VectorXd& tot_t = seir_res.get_value(s);

                const Eigen::VectorXd& F0 = flows_ts.get_value(s);
                const Eigen::VectorXd& F1 = flows_ts.get_value(s + 1);

                // Vorfaktoren 1/X_tot(t) für S,E,I je Altersgruppe (Divisionen sparen)
                std::vector<double> invS(num_age_groups), invE(num_age_groups), invI(num_age_groups);
                for (size_t g = 0; g < num_age_groups; ++g) {
                    const double sTot = tot_t[idxS[g]], eTot = tot_t[idxE[g]], iTot = tot_t[idxI[g]];
                    invS[g] = (sTot > 0.0) ? 1.0 / sTot : 0.0;
                    invE[g] = (eTot > 0.0) ? 1.0 / eTot : 0.0;
                    invI[g] = (iTot > 0.0) ? 1.0 / iTot : 0.0;
                }

                // Für jede Pendlergruppe: dF aufteilen und anwenden
                for (int cg = 0; cg < num_commuter_groups; ++cg) {
                    Eigen::VectorXd& x = mobile_pops[static_cast<size_t>(cg)];

                    for (size_t g = 0; g < num_age_groups; ++g) {
                        const size_t base = g * 3; // Reihenfolge der SEIR-Flüsse im memilio-SEIR: S->E, E->I, I->R

                        const double dSE =
                            F1[static_cast<Eigen::Index>(base + 0)] - F0[static_cast<Eigen::Index>(base + 0)];
                        const double dEI =
                            F1[static_cast<Eigen::Index>(base + 1)] - F0[static_cast<Eigen::Index>(base + 1)];
                        const double dIR =
                            F1[static_cast<Eigen::Index>(base + 2)] - F0[static_cast<Eigen::Index>(base + 2)];

                        const size_t iS = idxS[g], iE = idxE[g], iI = idxI[g], iR = idxR[g];

                        // Anteile am linken Rand (ta): xi = X_commuter / X_total
                        const double xiS = (invS[g] > 0.0) ? x[static_cast<Eigen::Index>(iS)] * invS[g] : 0.0;
                        const double xiE = (invE[g] > 0.0) ? x[static_cast<Eigen::Index>(iE)] * invE[g] : 0.0;
                        const double xiI = (invI[g] > 0.0) ? x[static_cast<Eigen::Index>(iI)] * invI[g] : 0.0;

                        // Aufteilen der Total-Flüsse auf die Gruppe
                        const double fSE = dSE * xiS;
                        const double fEI = dEI * xiE;
                        const double fIR = dIR * xiI;

                        // Massenhaltend updaten
                        x[static_cast<Eigen::Index>(iS)] -= fSE;
                        x[static_cast<Eigen::Index>(iE)] += fSE - fEI;
                        x[static_cast<Eigen::Index>(iI)] += fEI - fIR;
                        x[static_cast<Eigen::Index>(iR)] += fIR;
                    }
                }
            }

            benchmark::DoNotOptimize(mobile_pops);
        }
    }
}

static void bench_stage_aligned_hybrid(::benchmark::State& state)
{
    const int num_commuter_groups = commuter_group_counts[state.range(0)];
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

            const auto& seir_res = sim.get_result();
            const auto& flows_ts = sim.get_flows();
            const size_t n_steps = static_cast<size_t>(seir_res.get_num_time_points()) > 0
                                       ? static_cast<size_t>(seir_res.get_num_time_points()) - 1
                                       : 0;

            double mobile_fraction = 0.1 * num_commuter_groups;
            Eigen::VectorXd initial_mobile =
                seir_res.get_value(0) * (num_commuter_groups > 0 ? (mobile_fraction / num_commuter_groups) : 0.0);
            std::vector<Eigen::VectorXd> mobile_pops(static_cast<size_t>(num_commuter_groups), initial_mobile);

            using IS = mio::oseir::InfectionState;
            std::vector<size_t> idxS(num_age_groups), idxE(num_age_groups), idxI(num_age_groups), idxR(num_age_groups);
            for (size_t g = 0; g < num_age_groups; ++g) {
                auto G  = mio::AgeGroup(static_cast<int>(g));
                idxS[g] = model.populations.get_flat_index({G, IS::Susceptible});
                idxE[g] = model.populations.get_flat_index({G, IS::Exposed});
                idxI[g] = model.populations.get_flat_index({G, IS::Infected});
                idxR[g] = model.populations.get_flat_index({G, IS::Recovered});
            }

            state.ResumeTiming();

            // 3) Pendler-Update je Schritt: rein algebraisch, left-share
            for (size_t s = 0; s < n_steps; ++s) {
                const Eigen::VectorXd& tot_t = seir_res.get_value(s);

                const Eigen::VectorXd& F0 = flows_ts.get_value(s);
                const Eigen::VectorXd& F1 = flows_ts.get_value(s + 1);

                // Vorfaktoren 1/X_tot(t) für S,E,I je Altersgruppe (Divisionen sparen)
                std::vector<double> invS(num_age_groups), invE(num_age_groups), invI(num_age_groups);
                for (size_t g = 0; g < num_age_groups; ++g) {
                    const double sTot = tot_t[idxS[g]], eTot = tot_t[idxE[g]], iTot = tot_t[idxI[g]];
                    invS[g] = (sTot > 0.0) ? 1.0 / sTot : 0.0;
                    invE[g] = (eTot > 0.0) ? 1.0 / eTot : 0.0;
                    invI[g] = (iTot > 0.0) ? 1.0 / iTot : 0.0;
                }

                // Für jede Pendlergruppe: dF aufteilen und anwenden
                for (int cg = 0; cg < num_commuter_groups; ++cg) {
                    Eigen::VectorXd& x = mobile_pops[static_cast<size_t>(cg)];

                    for (size_t g = 0; g < num_age_groups; ++g) {
                        const size_t base = g * 3; // Reihenfolge der SEIR-Flüsse im memilio-SEIR: S->E, E->I, I->R

                        const double dSE =
                            F1[static_cast<Eigen::Index>(base + 0)] - F0[static_cast<Eigen::Index>(base + 0)];
                        const double dEI =
                            F1[static_cast<Eigen::Index>(base + 1)] - F0[static_cast<Eigen::Index>(base + 1)];
                        const double dIR =
                            F1[static_cast<Eigen::Index>(base + 2)] - F0[static_cast<Eigen::Index>(base + 2)];

                        const size_t iS = idxS[g], iE = idxE[g], iI = idxI[g], iR = idxR[g];

                        // Anteile am linken Rand (ta): xi = X_commuter / X_total
                        const double xiS = (invS[g] > 0.0) ? x[static_cast<Eigen::Index>(iS)] * invS[g] : 0.0;
                        const double xiE = (invE[g] > 0.0) ? x[static_cast<Eigen::Index>(iE)] * invE[g] : 0.0;
                        const double xiI = (invI[g] > 0.0) ? x[static_cast<Eigen::Index>(iI)] * invI[g] : 0.0;

                        // Aufteilen der Total-Flüsse auf die Gruppe
                        const double fSE = dSE * xiS;
                        const double fEI = dEI * xiE;
                        const double fIR = dIR * xiI;

                        // Massenhaltend updaten
                        x[static_cast<Eigen::Index>(iS)] -= fSE;
                        x[static_cast<Eigen::Index>(iE)] += fSE - fEI;
                        x[static_cast<Eigen::Index>(iI)] += fEI - fIR;
                        x[static_cast<Eigen::Index>(iR)] += fIR;
                    }
                }
            }

            benchmark::DoNotOptimize(mobile_pops);
        }
    }
}

static void bench_flow_based_exact_stagecache_by_groups(::benchmark::State& state)
{
    const int num_commuter_groups = mio::benchmark_mio::commuter_group_counts[state.range(0)];
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    using ModelType = mio::benchmark_mio::ModelType;
    using SimType   = mio::benchmark_mio::SimType;

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {

            state.PauseTiming();

            ModelType model(num_age_groups);
            mio::benchmark_mio::setup_model_benchmark(model);

            // 1) Totals einmal integrieren (wird gemessen – konsistent zu deinen anderen Benches)
            SimType sim(model, mio::benchmark_mio::t0, mio::benchmark_mio::dt);
            auto integrator_rk =
                std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
            sim.set_integrator(integrator_rk);

            state.ResumeTiming();
            sim.advance(mio::benchmark_mio::t_max);
            state.PauseTiming();

            const auto& res = sim.get_result();
            const size_t nsteps =
                (res.get_num_time_points() > 1) ? static_cast<size_t>(res.get_num_time_points() - 1) : 0;

            // Indizes für schnellen Zugriff
            const auto ic = mio::examples::make_seir_index_cache(sim.get_model());

            std::vector<mio::examples::StageCacheRK4> caches;
            caches.reserve(nsteps);
            for (size_t s = 0; s < nsteps; ++s) {
                caches.emplace_back(mio::examples::make_stage_cache_rk4(sim, s));
            }

            // Initiale Pendler je Gruppe (gleichmäßig verteilte Gesamtfraktion, identisch zu deinen anderen Benches)
            const double mobile_fraction = 0.1 * num_commuter_groups;
            const Eigen::VectorXd initial_mobile =
                res.get_value(0) *
                (num_commuter_groups > 0 ? (mobile_fraction / static_cast<double>(num_commuter_groups)) : 0.0);
            std::vector<Eigen::VectorXd> mobile_pops(static_cast<size_t>(num_commuter_groups), initial_mobile);

            state.ResumeTiming();

            // 3) Pro Schritt alle Gruppen mit dem Left-Share-je-Stufe-Update updaten (nur Vektor-Operationen)
            for (size_t s = 0; s < nsteps; ++s) {
                const auto& cache = caches[s];
                for (int cg = 0; cg < num_commuter_groups; ++cg) {
                    mio::examples::flow_based_mobility_returns_leftshare_rk4_from_cache(
                        mobile_pops[static_cast<size_t>(cg)], ic, cache);
                }
            }

            benchmark::DoNotOptimize(mobile_pops);
        }
    }
}

} // namespace benchmark_mio
} // namespace mio

BENCHMARK(mio::benchmark_mio::bench_stage_aligned_rk4_neu)
    ->Apply([](auto* b) {
        for (int i = 0; i < static_cast<int>(mio::benchmark_mio::commuter_group_counts.size()); ++i) {
            for (int g = 1; g <= 6; ++g) {
                b->Args({i, g});
            }
        }
    }) -> Name("stage-aligned(RK4)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_auxiliary_euler)
    ->Apply([](auto* b) {
        for (int i = 0; i < static_cast<int>(mio::benchmark_mio::commuter_group_counts.size()); ++i) {
            for (int g = 1; g <= 6; ++g) {
                b->Args({i, g});
            }
        }
    }) -> Name("auxiliary_Euler") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_stage_aligned_euler)
    ->Apply([](auto* b) {
        for (int i = 0; i < static_cast<int>(mio::benchmark_mio::commuter_group_counts.size()); ++i) {
            for (int g = 1; g <= 6; ++g) {
                b->Args({i, g});
            }
        }
    }) -> Name("stage-aligned(Euler)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_stage_aligned_hybrid)
    ->Apply([](auto* b) {
        for (int i = 0; i < static_cast<int>(mio::benchmark_mio::commuter_group_counts.size()); ++i) {
            for (int g = 1; g <= 6; ++g) {
                b->Args({i, g});
            }
        }
    }) -> Name("stage-aligned(hybrid)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_standard_lagrangian_rk4)
    ->Apply([](auto* b) {
        for (int i = 0; i < static_cast<int>(mio::benchmark_mio::commuter_group_counts.size()); ++i) {
            for (int g = 1; g <= 6; ++g) {
                b->Args({i, g});
            }
        }
    }) -> Name("lagrange_rk4") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_standard_lagrangian_euler)
    ->Apply([](auto* b) {
        for (int i = 0; i < static_cast<int>(mio::benchmark_mio::commuter_group_counts.size()); ++i) {
            for (int g = 1; g <= 6; ++g) {
                b->Args({i, g});
            }
        }
    }) -> Name("lagrange_euler") -> Unit(::benchmark::kMicrosecond);

// run all benchmarks
BENCHMARK_MAIN();
