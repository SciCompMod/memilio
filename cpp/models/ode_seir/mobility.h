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
#ifndef OSEIR_MOBILITY_H
#define OSEIR_MOBILITY_H

#include "ode_seir/model.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "memilio/math/integrator.h"
#include "memilio/math/eigen.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/time_series.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <memory>
#include <vector>

namespace mio
{
namespace oseir
{

/**
 * @brief Pre-computed SEIR rate evaluator for efficient stage-aligned mobility.
 *
 * Exploits the factorization  f(x) = B · D(z,t) · x  (cf. Zunker et al. 2026, 
 * https://doi.org/10.48550/arXiv.2603.11275) where the rate matrix D depends
 * only on the total (aggregated) population z, not on the individual commuter
 * subpopulation.
 *
 * @tparam FP Floating point type.
 */
template <typename FP>
struct SeirCoefficientEvaluator {
    const Model<FP>& model;
    const size_t NG; ///< Number of age groups.
    const size_t NC; ///< Total number of compartments (= 4 · NG for SEIR).

    // Pre-cached constant transition rates (independent of state z)
    std::vector<FP> rate_E; ///< 1 / TimeExposed per age group.
    std::vector<FP> rate_I; ///< 1 / TimeInfected per age group.
    std::vector<FP> beta; ///< TransmissionProbabilityOnContact per age group.

    // Pre-cached flat population indices
    std::vector<size_t> idxS, idxE, idxI, idxR;

    explicit SeirCoefficientEvaluator(const Model<FP>& m)
        : model(m)
        , NG(static_cast<size_t>(m.parameters.get_num_groups()))
        , NC(static_cast<size_t>(m.populations.get_num_compartments()))
        , rate_E(NG)
        , rate_I(NG)
        , beta(NG)
        , idxS(NG)
        , idxE(NG)
        , idxI(NG)
        , idxR(NG)
    {
        for (size_t g = 0; g < NG; ++g) {
            auto ag   = AgeGroup(static_cast<int>(g));
            idxS[g]   = m.populations.get_flat_index({ag, InfectionState::Susceptible});
            idxE[g]   = m.populations.get_flat_index({ag, InfectionState::Exposed});
            idxI[g]   = m.populations.get_flat_index({ag, InfectionState::Infected});
            idxR[g]   = m.populations.get_flat_index({ag, InfectionState::Recovered});
            rate_E[g] = FP(1) / m.parameters.template get<TimeExposed<FP>>()[ag];
            rate_I[g] = FP(1) / m.parameters.template get<TimeInfected<FP>>()[ag];
            beta[g]   = m.parameters.template get<TransmissionProbabilityOnContact<FP>>()[ag];
        }
    }

    /**
     * @brief Compute force of infection lambda_i from total population z.
     *
     * @param[in]  z      Total population vector (length NC).
     * @param[in]  t      Current time.
     * @param[out] lambda Force of infection per age group (length NG).
     */
    void compute_foi(const Eigen::VectorX<FP>& z, FP t, std::vector<FP>& lambda) const
    {
        std::fill(lambda.begin(), lambda.end(), FP(0));
        auto cm = model.parameters.template get<ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(
            SimulationTime<FP>(t));
        for (size_t j = 0; j < NG; ++j) {
            const FP Nj       = z[idxS[j]] + z[idxE[j]] + z[idxI[j]] + z[idxR[j]];
            const FP Ij_divNj = (Nj < Limits<FP>::zero_tolerance()) ? FP(0) : z[idxI[j]] / Nj;
            for (size_t i = 0; i < NG; ++i) {
                lambda[i] += cm(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) * beta[i] * Ij_divNj;
            }
        }
    }

    /**
     * @brief Apply SEIR dynamics to any subpopulation x using pre-computed lambda.
     *
     * @param[in]  lambda Pre-computed force of infection (from compute_foi).
     * @param[in]  x      Population vector (totals or commuters, length NC).
     * @param[out] dxdt   Derivative output vector (length NC).
     */
    void apply_rhs(const std::vector<FP>& lambda, const Eigen::VectorX<FP>& x, Eigen::VectorX<FP>& dxdt) const
    {
        dxdt.setZero();
        for (size_t g = 0; g < NG; ++g) {
            const FP flow_SE = lambda[g] * x[idxS[g]];
            const FP flow_EI = rate_E[g] * x[idxE[g]];
            const FP flow_IR = rate_I[g] * x[idxI[g]];
            dxdt[idxS[g]]    = -flow_SE;
            dxdt[idxE[g]]    = flow_SE - flow_EI;
            dxdt[idxI[g]]    = flow_EI - flow_IR;
            dxdt[idxR[g]]    = flow_IR;
        }
    }
};

/**
 * @brief SEIR-specific fixed-step RK4 integrator that caches lambda at every stage.
 *
 * During each `step()`, the classical 4-stage RK4 is performed using
 * SeirCoefficientEvaluator instead of the generic DerivFunction.
 * The force-of-infection lambda_s computed at each stage is stored
 * in an internal cache.
 *
 * After the node integration, `integrate_commuters()` reuses these cached
 * values to advance commuter subpopulations without any additional compute_foi
 * calls.
 *
 * The DerivFunction parameter in step() is ignored. The integrator uses its
 *       own SeirCoefficientEvaluator which must match the model attached to the
 *       simulation. Construct with the same model reference.
 *
 * @tparam FP Floating point type.
 */
template <typename FP>
class SeirRK4IntegratorCore : public OdeIntegratorCore<FP>
{
public:
    /// Cached lambda-vectors at the 4 RK stages.
    struct StageCache {
        std::array<std::vector<FP>, 4> lambda; ///< lambda_s for stages s=0..3.
        FP t0      = FP(0); ///< Start time of the cached step.
        FP dt      = FP(0); ///< Step size of the cached step.
        bool valid = false; ///< True after a successful step().
    };

    /**
     * @brief Construct from an oseir::Model.
     * @param model The SEIR model whose parameters define the RHS.
     *              Must outlive this integrator.
     */
    explicit SeirRK4IntegratorCore(const Model<FP>& model)
        : OdeIntegratorCore<FP>(FP{}, FP{}) // fixed-step: no dt bounds
        , m_eval(model)
    {
        const auto NC = static_cast<Eigen::Index>(m_eval.NC);
        m_k1.resize(NC);
        m_k2.resize(NC);
        m_k3.resize(NC);
        m_k4.resize(NC);
        m_tmp.resize(NC);
        for (auto& l : m_cache.lambda) {
            l.resize(m_eval.NG);
        }
    }

    std::unique_ptr<OdeIntegratorCore<FP>> clone() const override
    {
        return std::make_unique<SeirRK4IntegratorCore>(*this);
    }

    /**
     * @brief One fixed-step explicit RK4 integration step with lambda-caching.
     */
    bool step(const DerivFunction<FP>& /*f*/, Eigen::Ref<const Eigen::VectorX<FP>> yt, FP& t, FP& dt,
              Eigen::Ref<Eigen::VectorX<FP>> ytp1) const override
    {
        const FP half_dt = FP(0.5) * dt;

        // Stage 1 (t)
        m_eval.compute_foi(yt, t, m_cache.lambda[0]);
        m_eval.apply_rhs(m_cache.lambda[0], yt, m_k1);

        // Stage 2 (t + dt/2)
        m_tmp.noalias() = yt + half_dt * m_k1;
        m_eval.compute_foi(m_tmp, t + half_dt, m_cache.lambda[1]);
        m_eval.apply_rhs(m_cache.lambda[1], m_tmp, m_k2);

        // Stage 3 (t + dt/2)
        m_tmp.noalias() = yt + half_dt * m_k2;
        m_eval.compute_foi(m_tmp, t + half_dt, m_cache.lambda[2]);
        m_eval.apply_rhs(m_cache.lambda[2], m_tmp, m_k3);

        // Stage 4 (t + dt)
        m_tmp.noalias() = yt + dt * m_k3;
        m_eval.compute_foi(m_tmp, t + dt, m_cache.lambda[3]);
        m_eval.apply_rhs(m_cache.lambda[3], m_tmp, m_k4);

        // Final update
        ytp1.noalias() = yt + (dt / FP(6)) * (m_k1 + FP(2) * m_k2 + FP(2) * m_k3 + m_k4);

        // Fill cache
        m_cache.t0    = t;
        m_cache.dt    = dt;
        m_cache.valid = true;

        t += dt;
        return true;
    }

    /**
     * @brief Integrate a commuter subpopulation using cached lambda.
     *
     * @param[in,out] c  Commuter state vector; updated in-place to t + dt.
     * @param[in]     dt Step size (must match the cached step size).
     */
    void integrate_commuters(Eigen::Ref<typename TimeSeries<FP>::Vector> c, FP dt) const
    {
        assert(m_cache.valid && "integrate_commuters called without a preceding step()");
        const FP half_dt = FP(0.5) * dt;

        // Stage 1: apply cached lambda[0] to commuters
        m_eval.apply_rhs(m_cache.lambda[0], c, m_k1);

        // Stage 2: apply cached lambda[1]
        m_tmp.noalias() = c + half_dt * m_k1;
        m_eval.apply_rhs(m_cache.lambda[1], m_tmp, m_k2);

        // Stage 3: apply cached lambda[2]
        m_tmp.noalias() = c + half_dt * m_k2;
        m_eval.apply_rhs(m_cache.lambda[2], m_tmp, m_k3);

        // Stage 4: apply cached lambda[3]
        m_tmp.noalias() = c + dt * m_k3;
        m_eval.apply_rhs(m_cache.lambda[3], m_tmp, m_k4);

        // Final update
        c += (dt / FP(6)) * (m_k1 + FP(2) * m_k2 + FP(2) * m_k3 + m_k4);
    }

    /// @brief Whether the lambda-cache is valid (filled by the last step()).
    bool has_valid_cache() const
    {
        return m_cache.valid;
    }

    /// @brief Access the stage cache (e.g., for testing).
    const StageCache& get_cache() const
    {
        return m_cache;
    }

    /// @brief Access the evaluator (e.g., for testing).
    const SeirCoefficientEvaluator<FP>& get_evaluator() const
    {
        return m_eval;
    }

private:
    SeirCoefficientEvaluator<FP> m_eval;
    mutable StageCache m_cache;
    mutable Eigen::VectorX<FP> m_k1, m_k2, m_k3, m_k4, m_tmp;
};

namespace detail
{

/**
 * @brief Stage-aligned RK4 co-integration of commuter and total population.
 *
 * Fallback path when the node does NOT use SeirRK4IntegratorCore.
 */
template <typename FP>
void calculate_mobility_returns_seir(const SeirCoefficientEvaluator<FP>& eval,
                                     Eigen::Ref<typename TimeSeries<FP>::Vector> mobile_population,
                                     Eigen::Ref<const typename TimeSeries<FP>::Vector> total, FP t, FP dt)
{
    const size_t NC = eval.NC;
    const size_t NG = eval.NG;

    Eigen::VectorX<FP> z = total;
    Eigen::VectorX<FP> c = mobile_population;

    Eigen::VectorX<FP> k1_z(NC), k2_z(NC), k3_z(NC), k4_z(NC);
    Eigen::VectorX<FP> k1_c(NC), k2_c(NC), k3_c(NC), k4_c(NC);
    Eigen::VectorX<FP> z_s(NC), c_s(NC);
    std::vector<FP> lambda(NG);

    const FP half_dt = FP(0.5) * dt;

    // ---- Stage 1 (t) ----
    eval.compute_foi(z, t, lambda);
    eval.apply_rhs(lambda, z, k1_z);
    eval.apply_rhs(lambda, c, k1_c);

    // ---- Stage 2 (t + dt/2) ----
    z_s.noalias() = z + half_dt * k1_z;
    c_s.noalias() = c + half_dt * k1_c;
    eval.compute_foi(z_s, t + half_dt, lambda);
    eval.apply_rhs(lambda, z_s, k2_z);
    eval.apply_rhs(lambda, c_s, k2_c);

    // ---- Stage 3 (t + dt/2) ----
    z_s.noalias() = z + half_dt * k2_z;
    c_s.noalias() = c + half_dt * k2_c;
    eval.compute_foi(z_s, t + half_dt, lambda);
    eval.apply_rhs(lambda, z_s, k3_z);
    eval.apply_rhs(lambda, c_s, k3_c);

    // ---- Stage 4 (t + dt) ----
    z_s.noalias() = z + dt * k3_z;
    c_s.noalias() = c + dt * k3_c;
    eval.compute_foi(z_s, t + dt, lambda);
    eval.apply_rhs(lambda, z_s, k4_z);
    eval.apply_rhs(lambda, c_s, k4_c);

    // ---- Final update (only commuters) ----
    mobile_population = c + (dt / FP(6)) * (k1_c + FP(2) * k2_c + FP(2) * k3_c + k4_c);
}

} // namespace detail

/**
 * @brief Stage-aligned mobility returns for SEIR Simulation nodes.
 *
 * Ideal (SeirRK4IntegratorCore on node):
 *   Uses cached lambda from the node's integration -> only 4× apply_rhs,
 *   zero compute_foi calls, regardless of how many edges return to this node.
 *
 * Fallback (any other integrator):
 *   Full co-integration with 4× compute_foi + 8× apply_rhs.
 */
template <typename FP>
void calculate_mobility_returns(Eigen::Ref<typename TimeSeries<FP>::Vector> mobile_population,
                                const Simulation<FP, Model<FP>>& sim,
                                Eigen::Ref<const typename TimeSeries<FP>::Vector> total, FP t, FP dt)
{
    // Try ideal path: reuse cached lambda from node's SeirRK4IntegratorCore
    const auto* seir_core = dynamic_cast<const SeirRK4IntegratorCore<FP>*>(&sim.get_integrator_core());
    if (seir_core != nullptr && seir_core->has_valid_cache()) {
        seir_core->integrate_commuters(mobile_population, dt);
        return;
    }

    // Fallback: full co-integration (e.g. when using adaptive integrator)
    SeirCoefficientEvaluator<FP> eval(sim.get_model());
    detail::calculate_mobility_returns_seir<FP>(eval, mobile_population, total, t, dt);
}

/**
 * @brief Stage-aligned mobility returns for SEIR FlowSimulation nodes.
 * @see Simulation overload above.
 */
template <typename FP>
void calculate_mobility_returns(Eigen::Ref<typename TimeSeries<FP>::Vector> mobile_population,
                                const FlowSimulation<FP, Model<FP>>& sim,
                                Eigen::Ref<const typename TimeSeries<FP>::Vector> total, FP t, FP dt)
{
    // FlowSimulation: currently no caching (integrator is different), always fallback
    SeirCoefficientEvaluator<FP> eval(sim.get_model());
    detail::calculate_mobility_returns_seir<FP>(eval, mobile_population, total, t, dt);
}

/**
 * @brief Create a stage-aligned SEIR mobility simulation with lambda-caching.
 *
 * Sets SeirRK4IntegratorCore on every node and synchronizes the node's
 * internal dt with the graph dt.
 *
 * @param t0    Start time.
 * @param dt    Graph time step (= RK4 step size = mobility interval).
 * @param graph Mobility graph (will be moved).
 * @return A GraphSimulation ready to advance.
 */
template <typename FP>
auto make_seir_mobility_sim(FP t0, FP dt,
                            Graph<SimulationNode<FP, Simulation<FP, Model<FP>>>, MobilityEdge<FP>>&& graph)
{
    for (auto& node : graph.nodes()) {
        auto& sim = node.property.get_simulation();
        // Set the custom SEIR RK4 integrator with lambda-caching
        sim.set_integrator_core(std::make_unique<SeirRK4IntegratorCore<FP>>(sim.get_model()));
        // Synchronize node dt with graph dt so exactly one RK4 step is taken per graph step
        sim.get_dt() = dt;
    }
    return make_mobility_sim<FP>(t0, dt, std::move(graph));
}

/// @overload For lvalue graphs.
template <typename FP>
auto make_seir_mobility_sim(FP t0, FP dt, Graph<SimulationNode<FP, Simulation<FP, Model<FP>>>, MobilityEdge<FP>>& graph)
{
    for (auto& node : graph.nodes()) {
        auto& sim = node.property.get_simulation();
        sim.set_integrator_core(std::make_unique<SeirRK4IntegratorCore<FP>>(sim.get_model()));
        sim.get_dt() = dt;
    }
    return make_mobility_sim<FP>(t0, dt, graph);
}

} // namespace oseir
} // namespace mio

#endif // OSEIR_MOBILITY_H
