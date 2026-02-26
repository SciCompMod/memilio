#ifndef MIO_OSEIR_STATE_ESTIMATORS_H
#define MIO_OSEIR_STATE_ESTIMATORS_H

#include "memilio/compartments/flow_simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/type_list.h"
#include "memilio/math/eigen.h"
#include "memilio/math/stepper_wrapper.h"
#include <vector>
#include "ode_seir/model.h"
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

namespace mio
{
namespace oseir
{

using ScalarType = double;

// Enum class for dynamic commuter groups
enum class CommuterType
{
    NonCommuter = 0,
    CommuterBase,
};

enum class InfectionStateExplicit
{
    S = 0,
    E,
    I,
    R,
    Count
};

using ParametersExplicit = mio::oseir::Parameters<ScalarType>;

// Helper struct to dynamically create flows for multiple commuter groups
template <typename FP = ScalarType>
struct ExplicitFlowsCreator {
    static auto create_flows()
    {
        // Base flows for SEIR model (same for all commuter types)
        using BaseFlows = mio::TypeList<mio::Flow<InfectionStateExplicit::S, InfectionStateExplicit::E>,
                                        mio::Flow<InfectionStateExplicit::E, InfectionStateExplicit::I>,
                                        mio::Flow<InfectionStateExplicit::I, InfectionStateExplicit::R>>;

        return BaseFlows{};
    }
};

// Define flows using the helper struct
using FlowsExplicit = decltype(ExplicitFlowsCreator<>::create_flows());

template <typename FP = ScalarType>
class ModelExplicit : public mio::FlowModel<FP, InfectionStateExplicit,
                                            mio::Populations<FP, mio::AgeGroup, CommuterType, InfectionStateExplicit>,
                                            ParametersExplicit, FlowsExplicit>
{
    using Base = mio::FlowModel<FP, InfectionStateExplicit,
                                mio::Populations<FP, mio::AgeGroup, CommuterType, InfectionStateExplicit>,
                                ParametersExplicit, FlowsExplicit>;

public:
    using typename Base::ParameterSet;
    using typename Base::Populations;
    using CommuterIndex = CommuterType;

protected:
    int m_num_commuter_groups;

public:
    ModelExplicit(int num_agegroups, int num_commuter_groups)
        : Base(Populations({mio::AgeGroup(num_agegroups), CommuterIndex(num_commuter_groups + 1),
                            InfectionStateExplicit::Count}),
               ParametersExplicit(num_agegroups))
        , m_num_commuter_groups(num_commuter_groups)
    {
    }

    int get_num_commuter_groups() const
    {
        return m_num_commuter_groups;
    }

    // overwrite initial version
    Eigen::VectorX<FP> get_initial_flows() const
    {
        size_t num_flows = 3 * static_cast<size_t>(this->parameters.get_num_groups()) * (m_num_commuter_groups + 1);
        return Eigen::VectorX<FP>::Zero(num_flows);
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const auto& params = this->parameters;

        // E->I and I->R flows
        for (mio::AgeGroup g(0); g < mio::AgeGroup(params.get_num_groups()); ++g) {
            {
                CommuterType commuter_type = CommuterType::NonCommuter;
                const size_t E_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::E});
                const size_t I_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::I});

                flows[Base::template get_flat_flow_index<InfectionStateExplicit::E, InfectionStateExplicit::I>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeExposed<FP>>()[g]) * y[E_idx];

                flows[Base::template get_flat_flow_index<InfectionStateExplicit::I, InfectionStateExplicit::R>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeInfected<FP>>()[g]) * y[I_idx];
            }

            for (int c = 0; c < this->m_num_commuter_groups; ++c) {
                CommuterType commuter_type =
                    static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c);

                const size_t E_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::E});
                const size_t I_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::I});

                flows[Base::template get_flat_flow_index<InfectionStateExplicit::E, InfectionStateExplicit::I>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeExposed<FP>>()[g]) * y[E_idx];

                flows[Base::template get_flat_flow_index<InfectionStateExplicit::I, InfectionStateExplicit::R>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeInfected<FP>>()[g]) * y[I_idx];
            }
        }

        for (mio::AgeGroup i(0); i < mio::AgeGroup(params.get_num_groups()); ++i) {

            std::vector<FP> N_per_age_local(static_cast<size_t>(params.get_num_groups()));
            std::vector<FP> I_per_age_local(static_cast<size_t>(params.get_num_groups()));

            for (mio::AgeGroup j(0); j < mio::AgeGroup(params.get_num_groups()); ++j) {
                FP N_j = 0;
                FP I_j = 0;

                // NonCommuter
                for (size_t state = 0; state < static_cast<size_t>(InfectionStateExplicit::Count); ++state) {
                    N_j += pop[this->populations.get_flat_index(
                        {j, CommuterType::NonCommuter, static_cast<InfectionStateExplicit>(state)})];
                }
                I_j = pop[this->populations.get_flat_index({j, CommuterType::NonCommuter, InfectionStateExplicit::I})];

                // Commuters
                for (int c2 = 0; c2 < this->m_num_commuter_groups; ++c2) {
                    CommuterType ct_j = static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c2);
                    for (size_t state = 0; state < static_cast<size_t>(InfectionStateExplicit::Count); ++state) {
                        N_j += pop[this->populations.get_flat_index(
                            {j, ct_j, static_cast<InfectionStateExplicit>(state)})];
                    }
                    I_j += pop[this->populations.get_flat_index({j, ct_j, InfectionStateExplicit::I})];
                }

                N_per_age_local[j.get()] = N_j;
                I_per_age_local[j.get()] = I_j;
            }

            // Non-commuter flow
            {
                CommuterType commuter_type_i = CommuterType::NonCommuter;
                const size_t Si = this->populations.get_flat_index({i, commuter_type_i, InfectionStateExplicit::S});
                auto S_to_E_flow_idx =
                    Base::template get_flat_flow_index<InfectionStateExplicit::S, InfectionStateExplicit::E>(
                        {i, commuter_type_i});

                flows[S_to_E_flow_idx] = 0;

                for (mio::AgeGroup j(0); j < mio::AgeGroup(params.get_num_groups()); ++j) {
                    std::vector<FP> temp_div(1);
                    temp_div[0] = (N_per_age_local[j.get()] < mio::Limits<FP>::zero_tolerance())
                                      ? 0.0
                                      : FP(1.0) / N_per_age_local[j.get()];
                    const FP divN_j = temp_div[0];

                    const FP coeffStoE =
                        params.template get<mio::oseir::ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(t)(
                            i.get(), j.get()) *
                        params.template get<mio::oseir::TransmissionProbabilityOnContact<FP>>()[i] * divN_j;

                    flows[S_to_E_flow_idx] += coeffStoE * y[Si] * I_per_age_local[j.get()];
                }
            }

            // Commuter flows
            for (int c = 0; c < this->m_num_commuter_groups; ++c) {
                CommuterType commuter_type_i =
                    static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c);
                const size_t Si = this->populations.get_flat_index({i, commuter_type_i, InfectionStateExplicit::S});
                auto S_to_E_flow_idx =
                    Base::template get_flat_flow_index<InfectionStateExplicit::S, InfectionStateExplicit::E>(
                        {i, commuter_type_i});

                flows[S_to_E_flow_idx] = 0;

                for (mio::AgeGroup j(0); j < mio::AgeGroup(params.get_num_groups()); ++j) {
                    std::vector<FP> temp_div(1);
                    temp_div[0] = (N_per_age_local[j.get()] < mio::Limits<FP>::zero_tolerance())
                                      ? 0.0
                                      : FP(1.0) / N_per_age_local[j.get()];
                    const FP divN_j = temp_div[0];

                    const FP coeffStoE =
                        params.template get<mio::oseir::ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(t)(
                            i.get(), j.get()) *
                        params.template get<mio::oseir::TransmissionProbabilityOnContact<FP>>()[i] * divN_j;

                    flows[S_to_E_flow_idx] += coeffStoE * y[Si] * I_per_age_local[j.get()];
                }
            }
        }
    }
};

using StandardModelLagrangian = ModelExplicit<ScalarType>;
using StandardLagrangianSim   = mio::FlowSimulation<ScalarType, StandardModelLagrangian>;

using StandardModel = mio::oseir::Model<ScalarType>;
using StandardSim   = mio::FlowSimulation<ScalarType, StandardModel>;

template <class ModelType, size_t FlowIdx = 0>
void apply_flows_recursive(Eigen::Ref<Eigen::VectorXd> mobile_population,
                           const Eigen::VectorXd& initial_mobile_population, const ModelType& model,
                           const typename ModelType::FlowIndex& flow_indices, Eigen::Ref<const Eigen::VectorXd> total,
                           Eigen::Ref<const Eigen::VectorXd> diff);

template <class ModelType>
void apply_flows_to_mobile_population_generic(Eigen::Ref<Eigen::VectorXd> mobile_population, const ModelType& model,
                                              Eigen::Ref<const Eigen::VectorXd> total, const Eigen::VectorXd& diff)
{
    using FlowIndex = typename ModelType::FlowIndex;

    // Copy of mobile population to calculate the scaling factors
    Eigen::VectorXd initial_mobile_population = mobile_population;

    // Iterate through all flow indices
    for (FlowIndex flow_indices : mio::make_index_range(mio::reduce_index<FlowIndex>(model.populations.size()))) {
        apply_flows_recursive<ModelType>(mobile_population, initial_mobile_population, model, flow_indices, total,
                                         diff);
    }
}

// Recursive helper function to iterate through flows
template <class ModelType, size_t FlowIdx>
void apply_flows_recursive(Eigen::Ref<Eigen::VectorXd> mobile_population,
                           const Eigen::VectorXd& initial_mobile_population, const ModelType& model,
                           const typename ModelType::FlowIndex& flow_indices, Eigen::Ref<const Eigen::VectorXd> total,
                           Eigen::Ref<const Eigen::VectorXd> diff)
{
    using FlowsList   = typename ModelType::Flow_List;
    using CurrentFlow = mio::type_at_index_t<FlowIdx, FlowsList>;
    using PopIndex    = typename ModelType::Populations::Index;

    constexpr auto source_comp = CurrentFlow::source;
    constexpr auto target_comp = CurrentFlow::target;

    // Get flat indices for population compartments
    size_t source_pop_idx =
        model.populations.get_flat_index(mio::extend_index<PopIndex>(flow_indices, (size_t)source_comp));
    size_t target_pop_idx =
        model.populations.get_flat_index(mio::extend_index<PopIndex>(flow_indices, (size_t)target_comp));

    // Calculate compartment-specific scaling factor for the source compartment using initial population.
    ScalarType scale_source = (total[source_pop_idx] > mio::Limits<ScalarType>::zero_tolerance())
                                  ? (initial_mobile_population[source_pop_idx] / total[source_pop_idx])
                                  : 0.0;
    scale_source = std::max((ScalarType)0.0, std::min((ScalarType)1.0, scale_source));

    size_t flow_flat_idx = model.template get_flat_flow_index<source_comp, target_comp>(flow_indices);

    if (flow_flat_idx < static_cast<size_t>(diff.size())) {
        ScalarType flow_mobile = diff[flow_flat_idx] * scale_source;
        mobile_population[source_pop_idx] -= flow_mobile;
        mobile_population[target_pop_idx] += flow_mobile;
    }
    else {
        mio::log_warning("Flow index {} out of bounds for diff vector size {}.", flow_flat_idx, diff.size());
    }

    // Recursively call for the next flow
    if constexpr (FlowIdx + 1 < FlowsList::size()) {
        apply_flows_recursive<ModelType, FlowIdx + 1>(mobile_population, initial_mobile_population, model, flow_indices,
                                                      total, diff);
    }
}

using SimType = mio::FlowSimulation<ScalarType, mio::oseir::Model<ScalarType>>;

void integrate_mobile_population_euler(Eigen::Ref<Eigen::VectorXd> mobile_population, const SimType& sim,
                                       Eigen::Ref<const Eigen::VectorXd> total, ScalarType t, ScalarType dt)
{
    auto y0 = mobile_population.eval();
    auto y1 = mobile_population;
    mio::EulerIntegratorCore<ScalarType>().step(
        [&](auto&& y, auto&& t_, auto&& dydt) {
            sim.get_model().get_derivatives(total, y, t_, dydt);
        },
        y0, t, dt, y1);
}

// Helper function to apply flows to mobile population (Specific to SEIR model structure).
void apply_flows_to_mobile_population(Eigen::Ref<Eigen::VectorXd> mobile_population, const SimType& sim,
                                      Eigen::Ref<const Eigen::VectorXd> total, const Eigen::VectorXd& diff)
{
    const auto& model       = sim.get_model();
    const size_t num_groups = static_cast<size_t>(model.parameters.get_num_groups());
    using InfState          = mio::oseir::InfectionState;

    std::vector<size_t> S_idx(num_groups), E_idx(num_groups), I_idx(num_groups), R_idx(num_groups);
    for (size_t g = 0; g < num_groups; ++g) {
        mio::AgeGroup group(g);
        S_idx[g] = model.populations.get_flat_index({group, InfState::Susceptible});
        E_idx[g] = model.populations.get_flat_index({group, InfState::Exposed});
        I_idx[g] = model.populations.get_flat_index({group, InfState::Infected});
        R_idx[g] = model.populations.get_flat_index({group, InfState::Recovered});
    }

    for (size_t g = 0; g < num_groups; ++g) {
        size_t flow_base_idx = g * 3;

        double scale_S = (total[S_idx[g]] > 1e-10) ? (mobile_population[S_idx[g]] / total[S_idx[g]]) : 0.0;
        double scale_E = (total[E_idx[g]] > 1e-10) ? (mobile_population[E_idx[g]] / total[E_idx[g]]) : 0.0;
        double scale_I = (total[I_idx[g]] > 1e-10) ? (mobile_population[I_idx[g]] / total[I_idx[g]]) : 0.0;

        scale_S = std::max(0.0, std::min(1.0, scale_S));
        scale_E = std::max(0.0, std::min(1.0, scale_E));
        scale_I = std::max(0.0, std::min(1.0, scale_I));

        double flow_SE_mobile = diff[flow_base_idx] * scale_S;
        double flow_EI_mobile = diff[flow_base_idx + 1] * scale_E;
        double flow_IR_mobile = diff[flow_base_idx + 2] * scale_I;

        mobile_population[S_idx[g]] -= flow_SE_mobile;
        mobile_population[E_idx[g]] += flow_SE_mobile - flow_EI_mobile;
        mobile_population[I_idx[g]] += flow_EI_mobile - flow_IR_mobile;
        mobile_population[R_idx[g]] += flow_IR_mobile;
    }
}

void flow_based_mobility_returns(Eigen::Ref<Eigen::VectorXd> mobile_population, const SimType& sim,
                                 Eigen::Ref<const Eigen::VectorXd> total, ScalarType t, ScalarType dt)
{
    const auto& flows = sim.get_flows();

    const auto num_time_points = flows.get_num_time_points();

    // Check if there are at least two flow data points; if not, return
    if (num_time_points < 2) {
        mio::log_warning("Flow data not available for mobility return calculation.");
        return;
    }

    // assuming equidistant time points in flows
    auto dt_flows = flows.get_time(1) - flows.get_time(0);
    if (dt_flows <= 0) {
        mio::log_warning("Flow time step is non-positive.");
        return;
    }
    auto to_idx = [&](double tau) {
        return static_cast<size_t>(std::floor((tau + 1e-12) / dt_flows));
    };
    size_t current_idx = to_idx(t);
    size_t next_idx    = to_idx(t + dt);

    if (current_idx >= static_cast<size_t>(num_time_points) || next_idx >= static_cast<size_t>(num_time_points)) {
        mio::log_warning("Calculated flow indices out of bounds.");
        return;
    }

    if (current_idx != next_idx) {
        auto current_flow    = flows.get_value(current_idx);
        auto next_flow       = flows.get_value(next_idx);
        Eigen::VectorXd diff = next_flow - current_flow;
        apply_flows_to_mobile_population_generic(mobile_population, sim.get_model(), total, diff);
    }
    else {
        mio::log_warning("Could not find appropriate flow data points dt={} apart at time t={}. dt_flows={}", dt, t,
                         dt_flows);
    }
}

mio::TimeSeries<ScalarType> calculate_mobile_population(
    ScalarType t0, ScalarType tmax, ScalarType step_size,
    const std::function<void(Eigen::Ref<Eigen::VectorXd>, const SimType&, Eigen::Ref<const Eigen::VectorXd>, ScalarType,
                             ScalarType)>& method,
    const SimType& sim, const mio::TimeSeries<ScalarType>& seir_res, const Eigen::VectorXd& initial_mobile_pop)
{
    Eigen::VectorXd mobile_pop = initial_mobile_pop;
    mio::TimeSeries<ScalarType> solution(mobile_pop.size());
    solution.add_time_point(t0, mobile_pop);
    const auto step_size_ref = seir_res.get_time(1) - seir_res.get_time(0);
    for (ScalarType t = t0; t < tmax; t += step_size) {
        if (t + step_size > tmax + 1e-10) {
            break;
        }
        const auto closest_idx_total = static_cast<size_t>(std::round(t / step_size_ref));
        const auto& pop              = seir_res.get_value(closest_idx_total);
        method(mobile_pop, sim, pop, t, step_size);
        auto ti = seir_res.get_time(closest_idx_total) + step_size;
        mio::unused(ti);
        solution.add_time_point(seir_res.get_time(closest_idx_total) + step_size, mobile_pop);
    }

    return solution;
}

namespace
{

struct ModelEvaluator {
    const mio::oseir::Model<double>& model;
    std::size_t NG;
    std::size_t NC;

    // Precomputed rates
    std::vector<double> rate_E;
    std::vector<double> rate_I;

    // Index maps
    std::vector<int> idxS, idxE, idxI, idxR;

    ModelEvaluator(const mio::oseir::Model<double>& m)
        : model(m)
    {
        using IS = mio::oseir::InfectionState;
        NG       = static_cast<std::size_t>(model.parameters.get_num_groups());
        NC       = 4 * NG;

        rate_E.resize(NG);
        rate_I.resize(NG);
        idxS.resize(NG);
        idxE.resize(NG);
        idxI.resize(NG);
        idxR.resize(NG);

        for (std::size_t g = 0; g < NG; ++g) {
            auto ag = mio::AgeGroup(static_cast<int>(g));
            // Cache indices
            idxS[g] = (int)model.populations.get_flat_index({ag, IS::Susceptible});
            idxE[g] = (int)model.populations.get_flat_index({ag, IS::Exposed});
            idxI[g] = (int)model.populations.get_flat_index({ag, IS::Infected});
            idxR[g] = (int)model.populations.get_flat_index({ag, IS::Recovered});

            // Cache rates
            double t_E = model.parameters.get<mio::oseir::TimeExposed<double>>()[ag];
            double t_I = model.parameters.get<mio::oseir::TimeInfected<double>>()[ag];
            rate_E[g]  = (t_E > 1e-10) ? (1.0 / t_E) : 0.0;
            rate_I[g]  = (t_I > 1e-10) ? (1.0 / t_I) : 0.0;
        }
    }

    // 1. Totals RHS (Standard SEIR)
    void get_totals_rhs(const Eigen::VectorXd& y, double t, Eigen::VectorXd& dy) const
    {
        dy.resize(NC);
        model.get_derivatives(y, y, t, dy);
    }

    // 2. Compute Intensity
    void get_lambda(const Eigen::VectorXd& y, double t, std::vector<double>& lambda_out) const
    {
        lambda_out.assign(NG, 0.0);
        for (std::size_t j = 0; j < NG; ++j) {
            const double Nj    = y[idxS[j]] + y[idxE[j]] + y[idxI[j]] + y[idxR[j]];
            const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;

            for (std::size_t i = 0; i < NG; ++i) {
                const double coeff =
                    model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat().get_matrix_at(t)(
                        i, j) *
                    model.parameters
                        .get<mio::oseir::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup((int)i)] *
                    divNj;
                lambda_out[i] += coeff * y[idxI[j]];
            }
        }
    }

    // 3. Commuter RHS
    void get_commuter_rhs(const std::vector<double>& lambda, const Eigen::VectorXd& Xc, Eigen::VectorXd& dxc) const
    {
        dxc.setZero(NC);
        for (std::size_t g = 0; g < NG; ++g) {
            const int iS = idxS[g], iE = idxE[g], iI = idxI[g], iR = idxR[g];

            const double fSE = lambda[g] * Xc[iS];
            const double fEI = rate_E[g] * Xc[iE];
            const double fIR = rate_I[g] * Xc[iI];

            dxc[iS] = -fSE;
            dxc[iE] = fSE - fEI;
            dxc[iI] = fEI - fIR;
            dxc[iR] = fIR;
        }
    }
};

} // namespace

inline void flow_based_mobility_returns_rk2(Eigen::Ref<Eigen::VectorXd> commuter, Eigen::Ref<Eigen::VectorXd> totals,
                                            const mio::oseir::Model<double>& model, double t, double dt)
{
    // Setup Helper
    ModelEvaluator eval(model);
    const std::size_t NC = eval.NC;
    const std::size_t NG = eval.NG;

    // Cache
    Eigen::VectorXd k1_tot(NC), k2_tot(NC);
    Eigen::VectorXd k1_com(NC), k2_com(NC);
    std::vector<double> lambda(NG);

    // --- Stage 1 (t) ---
    eval.get_totals_rhs(totals, t, k1_tot);
    eval.get_lambda(totals, t, lambda);
    eval.get_commuter_rhs(lambda, commuter, k1_com);

    // --- Stage 2 (t + dt/2) ---
    Eigen::VectorXd totals_mid = totals + (dt * 0.5) * k1_tot;
    Eigen::VectorXd comm_mid   = commuter + (dt * 0.5) * k1_com;

    eval.get_totals_rhs(totals_mid, t + 0.5 * dt, k2_tot);
    eval.get_lambda(totals_mid, t + 0.5 * dt, lambda);
    eval.get_commuter_rhs(lambda, comm_mid, k2_com);

    // --- Final Update (Standard RK2) ---
    totals   = totals + dt * k2_tot;
    commuter = commuter + dt * k2_com;
}

inline void flow_based_mobility_returns_rk3(Eigen::Ref<Eigen::VectorXd> commuter, Eigen::Ref<Eigen::VectorXd> totals,
                                            const mio::oseir::Model<double>& model, double t, double dt)
{
    // Setup Helper
    ModelEvaluator eval(model);
    const std::size_t NC = eval.NC;
    const std::size_t NG = eval.NG;

    // Cache
    Eigen::VectorXd k1_tot(NC), k2_tot(NC), k3_tot(NC);
    Eigen::VectorXd k1_com(NC), k2_com(NC), k3_com(NC);
    std::vector<double> lambda(NG);

    // --- Stage 1 (t) ---
    eval.get_totals_rhs(totals, t, k1_tot);
    eval.get_lambda(totals, t, lambda);
    eval.get_commuter_rhs(lambda, commuter, k1_com);

    // --- Stage 2 (t + dt/2) ---
    Eigen::VectorXd totals_2 = totals + (dt * 0.5) * k1_tot;
    Eigen::VectorXd comm_2   = commuter + (dt * 0.5) * k1_com;

    eval.get_totals_rhs(totals_2, t + 0.5 * dt, k2_tot);
    eval.get_lambda(totals_2, t + 0.5 * dt, lambda);
    eval.get_commuter_rhs(lambda, comm_2, k2_com);

    // --- Stage 3 (t + dt) using Kutta's formula: y - dt*k1 + 2*dt*k2 ---
    Eigen::VectorXd totals_3 = totals + dt * (-k1_tot + 2.0 * k2_tot);
    Eigen::VectorXd comm_3   = commuter + dt * (-k1_com + 2.0 * k2_com);

    eval.get_totals_rhs(totals_3, t + dt, k3_tot);
    eval.get_lambda(totals_3, t + dt, lambda);
    eval.get_commuter_rhs(lambda, comm_3, k3_com);

    // --- Final Update (Standard RK3: 1/6 * (k1 + 4k2 + k3)) ---
    totals   = totals + (dt / 6.0) * (k1_tot + 4.0 * k2_tot + k3_tot);
    commuter = commuter + (dt / 6.0) * (k1_com + 4.0 * k2_com + k3_com);
}

inline void flow_based_mobility_returns_rk4(Eigen::Ref<Eigen::VectorXd> commuter, Eigen::Ref<Eigen::VectorXd> totals,
                                            const mio::oseir::Model<double>& model, double t, double dt)
{
    ModelEvaluator eval(model);
    const std::size_t NC = eval.NC;
    const std::size_t NG = eval.NG;

    Eigen::VectorXd k1_tot(NC), k2_tot(NC), k3_tot(NC), k4_tot(NC);
    Eigen::VectorXd k1_com(NC), k2_com(NC), k3_com(NC), k4_com(NC);
    std::vector<double> lambda(NG);

    // Stage 1
    eval.get_totals_rhs(totals, t, k1_tot);
    eval.get_lambda(totals, t, lambda);
    eval.get_commuter_rhs(lambda, commuter, k1_com);

    // Stage 2
    Eigen::VectorXd y2 = totals + (dt * 0.5) * k1_tot;
    Eigen::VectorXd c2 = commuter + (dt * 0.5) * k1_com;
    eval.get_totals_rhs(y2, t + 0.5 * dt, k2_tot);
    eval.get_lambda(y2, t + 0.5 * dt, lambda);
    eval.get_commuter_rhs(lambda, c2, k2_com);

    // Stage 3
    Eigen::VectorXd y3 = totals + (dt * 0.5) * k2_tot;
    Eigen::VectorXd c3 = commuter + (dt * 0.5) * k2_com;
    eval.get_totals_rhs(y3, t + 0.5 * dt, k3_tot);
    eval.get_lambda(y3, t + 0.5 * dt, lambda);
    eval.get_commuter_rhs(lambda, c3, k3_com);

    // Stage 4
    Eigen::VectorXd y4 = totals + dt * k3_tot;
    Eigen::VectorXd c4 = commuter + dt * k3_com;
    eval.get_totals_rhs(y4, t + dt, k4_tot);
    eval.get_lambda(y4, t + dt, lambda);
    eval.get_commuter_rhs(lambda, c4, k4_com);

    // Final Update
    totals += (dt / 6.0) * (k1_tot + 2.0 * k2_tot + 2.0 * k3_tot + k4_tot);
    commuter += (dt / 6.0) * (k1_com + 2.0 * k2_com + 2.0 * k3_com + k4_com);
}

// ============================================================================
// Fundamental Matrix (Phi) Methods
// ============================================================================

/**
 * @brief Builds the Jacobian matrix A(z,t) for the SEIR model.
 * Block-diagonal form per age group g:
 * [ -lambda_g     0         0     0 ]
 * [  lambda_g  -1/T_E       0     0 ]
 * [     0       1/T_E   -1/T_I    0 ]
 * [     0         0      1/T_I    0 ]
 *
 * @param A     NC x NC output matrix (written by this function)
 * @param model SEIR model (provides parameters and population indices)
 * @param z     current total population vector (length NC = 4*num_groups)
 * @param t     current time
 */
inline void build_seir_matrix_A(Eigen::MatrixXd& A, const mio::oseir::Model<double>& model, const Eigen::VectorXd& z,
                                double t)
{
    using IS                = mio::oseir::InfectionState;
    const size_t num_groups = static_cast<size_t>(model.parameters.get_num_groups());

    // 1. Compute force of infection lambda_i for each age group i
    Eigen::VectorXd force_of_infection = Eigen::VectorXd::Zero(num_groups);
    for (size_t j = 0; j < num_groups; ++j) {
        const size_t Sj = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Susceptible});
        const size_t Ej = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Exposed});
        const size_t Ij = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Infected});
        const size_t Rj = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Recovered});

        const double Nj    = z[Sj] + z[Ej] + z[Ij] + z[Rj];
        const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;

        for (size_t i = 0; i < num_groups; ++i) {
            const double coeffStoE =
                model.parameters.template get<mio::oseir::ContactPatterns<ScalarType>>()
                    .get_cont_freq_mat()
                    .get_matrix_at(t)(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) *
                model.parameters.template get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(
                    static_cast<int>(i))] *
                divNj;
            force_of_infection[i] += coeffStoE * z[Ij];
        }
    }

    A.setZero();

    // 2. Fill matrix A blockwise (one 4x4 block per age group)
    for (size_t g = 0; g < num_groups; ++g) {
        const double time_exposed  = model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(g)];
        const double time_infected = model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(g)];

        const double rate_E = (time_exposed > 1e-10) ? (1.0 / time_exposed) : 0.0;
        const double rate_I = (time_infected > 1e-10) ? (1.0 / time_infected) : 0.0;
        const double lambda = force_of_infection[g];

        const size_t iS = model.populations.get_flat_index({mio::AgeGroup(g), IS::Susceptible});
        const size_t iE = model.populations.get_flat_index({mio::AgeGroup(g), IS::Exposed});
        const size_t iI = model.populations.get_flat_index({mio::AgeGroup(g), IS::Infected});
        const size_t iR = model.populations.get_flat_index({mio::AgeGroup(g), IS::Recovered});

        // S -> E
        A(iS, iS) = -lambda;
        A(iE, iS) = lambda;
        // E -> I
        A(iE, iE) = -rate_E;
        A(iI, iE) = rate_E;
        // I -> R
        A(iI, iI) = -rate_I;
        A(iR, iI) = rate_I;
    }
}

/**
 * @brief Augmented ODE system coupling total population dynamics z with the fundamental matrix Phi.
 *
 * State vector: y = [z (NC), vec(Phi) (NC*NC)].
 * Phi satisfies dPhi/dt = A(z,t)*Phi with Phi(t0,t0) = I.
 * Once integrated to time t: X_c(t) = Phi(t, t0) * X_c(t0).
 */
struct AugmentedPhiSystem {
    const mio::oseir::Model<ScalarType>& model;
    size_t NC;
    Eigen::MatrixXd A_mem;

    explicit AugmentedPhiSystem(const mio::oseir::Model<ScalarType>& m)
        : model(m)
        , NC(static_cast<size_t>(m.populations.get_num_compartments()))
    {
        A_mem = Eigen::MatrixXd::Zero(NC, NC);
    }

    void operator()(const Eigen::VectorXd& y, Eigen::VectorXd& dydt, double t)
    {
        const auto z = y.head(NC);

        // dz/dt (standard SEIR RHS)
        Eigen::VectorXd dz(NC);
        model.get_derivatives(z, z, t, dz);
        dydt.head(NC) = dz;

        // Build A(z,t) and compute dPhi/dt = A * Phi
        build_seir_matrix_A(A_mem, model, z, t);
        const auto Phi       = Eigen::Map<const Eigen::MatrixXd>(y.tail(NC * NC).data(), NC, NC);
        Eigen::MatrixXd dPhi = A_mem * Phi;
        dydt.tail(NC * NC)   = Eigen::Map<const Eigen::VectorXd>(dPhi.data(), NC * NC);
    }
};

/**
 * @brief Commuter update via fundamental matrix Phi (single RK4 step on augmented system).
 *
 * Integrates the coupled system (z, Phi) from t to t+dt with a single classical RK4 step,
 * then applies the linear map: X_c(t+dt) = Phi(t+dt, t) * X_c(t).
 *
 * @param commuter  in/out: commuter population X_c(t) -> X_c(t+dt)
 * @param totals    in/out: total population z(t) -> z(t+dt)
 * @param model     SEIR model (provides RHS and parameters)
 * @param t         current time
 * @param dt        step size
 */
inline void flow_based_mobility_returns_phi_rk4(Eigen::Ref<Eigen::VectorXd> commuter,
                                                Eigen::Ref<Eigen::VectorXd> totals,
                                                const mio::oseir::Model<double>& model, double t, double dt)
{
    const std::size_t NC      = static_cast<std::size_t>(model.populations.get_num_compartments());
    const std::size_t sys_dim = NC + NC * NC;

    // Augmented initial state: y = [z(t), vec(I_NC)]
    Eigen::VectorXd y(sys_dim);
    y.head(NC)               = totals;
    const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(NC, NC);
    y.tail(NC * NC)          = Eigen::Map<const Eigen::VectorXd>(Id.data(), NC * NC);

    // Single classical RK4 step on the augmented system
    AugmentedPhiSystem sys(model);
    boost::numeric::odeint::runge_kutta4<Eigen::VectorXd, double, Eigen::VectorXd, double,
                                         boost::numeric::odeint::vector_space_algebra>
        stepper;
    stepper.do_step(sys, y, t, dt);

    // Extract updated totals z(t+dt)
    Eigen::VectorXd z_next = y.head(NC);

    // Apply fundamental matrix: X_c(t+dt) = Phi * X_c(t)
    const auto Phi            = Eigen::Map<const Eigen::MatrixXd>(y.tail(NC * NC).data(), NC, NC);
    const Eigen::VectorXd Xc0 = commuter;
    Eigen::VectorXd Xc_next   = Phi * Xc0;

    // Numerical guards: non-negativity and commuter <= totals
    for (int i = 0; i < (int)NC; ++i) {
        if (z_next[i] < 0.0)
            z_next[i] = 0.0;
        if (Xc_next[i] < 0.0)
            Xc_next[i] = 0.0;
        if (Xc_next[i] > z_next[i])
            Xc_next[i] = z_next[i];
    }

    totals   = z_next;
    commuter = Xc_next;
}

/**
 * @brief Helper: extract totals and apply Phi-map with numerical guards.
 * Internal use by phi_euler/rk2/rk3.
 */
inline void phi_apply_and_guard(Eigen::Ref<Eigen::VectorXd> commuter, Eigen::Ref<Eigen::VectorXd> totals,
                                const Eigen::VectorXd& y_next, std::size_t NC)
{
    Eigen::VectorXd z_next    = y_next.head(NC);
    const auto Phi            = Eigen::Map<const Eigen::MatrixXd>(y_next.tail(NC * NC).data(), NC, NC);
    const Eigen::VectorXd Xc0 = commuter;
    Eigen::VectorXd Xc_next   = Phi * Xc0;

    for (int i = 0; i < (int)NC; ++i) {
        if (z_next[i] < 0.0)
            z_next[i] = 0.0;
        if (Xc_next[i] < 0.0)
            Xc_next[i] = 0.0;
        if (Xc_next[i] > z_next[i])
            Xc_next[i] = z_next[i];
    }
    totals   = z_next;
    commuter = Xc_next;
}

/**
 * @brief Commuter update via fundamental matrix Phi – explicit Euler step.
 *
 * Integrates the augmented system (z, Phi) with a single explicit Euler step,
 * then applies X_c(t+dt) = Phi(t+dt, t) * X_c(t).
 */
inline void flow_based_mobility_returns_phi_euler(Eigen::Ref<Eigen::VectorXd> commuter,
                                                  Eigen::Ref<Eigen::VectorXd> totals,
                                                  const mio::oseir::Model<double>& model, double t, double dt)
{
    const std::size_t NC      = static_cast<std::size_t>(model.populations.get_num_compartments());
    const std::size_t sys_dim = NC + NC * NC;

    Eigen::VectorXd y(sys_dim);
    y.head(NC)         = totals;
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(NC, NC);
    y.tail(NC * NC)    = Eigen::Map<const Eigen::VectorXd>(Id.data(), NC * NC);

    AugmentedPhiSystem sys(model);
    Eigen::VectorXd k1(sys_dim);
    sys(y, k1, t);
    Eigen::VectorXd y_next = y + dt * k1;

    phi_apply_and_guard(commuter, totals, y_next, NC);
}

/**
 * @brief Commuter update via fundamental matrix Phi – RK2 (midpoint) step.
 *
 * Integrates the augmented system (z, Phi) with a single classical midpoint step,
 * then applies X_c(t+dt) = Phi(t+dt, t) * X_c(t).
 */
inline void flow_based_mobility_returns_phi_rk2(Eigen::Ref<Eigen::VectorXd> commuter,
                                                Eigen::Ref<Eigen::VectorXd> totals,
                                                const mio::oseir::Model<double>& model, double t, double dt)
{
    const std::size_t NC      = static_cast<std::size_t>(model.populations.get_num_compartments());
    const std::size_t sys_dim = NC + NC * NC;

    Eigen::VectorXd y(sys_dim);
    y.head(NC)         = totals;
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(NC, NC);
    y.tail(NC * NC)    = Eigen::Map<const Eigen::VectorXd>(Id.data(), NC * NC);

    AugmentedPhiSystem sys(model);
    Eigen::VectorXd k1(sys_dim), k2(sys_dim);

    sys(y, k1, t);
    sys(y + (dt * 0.5) * k1, k2, t + 0.5 * dt);

    Eigen::VectorXd y_next = y + dt * k2;

    phi_apply_and_guard(commuter, totals, y_next, NC);
}

/**
 * @brief Commuter update via fundamental matrix Phi – RK3 (Kutta) step.
 *
 * Integrates the augmented system (z, Phi) with Kutta's classical 3rd-order method:
 *   k1 = f(t, y)
 *   k2 = f(t + dt/2, y + dt/2 * k1)
 *   k3 = f(t + dt,   y - dt*k1 + 2*dt*k2)
 *   y_{n+1} = y + dt/6*(k1 + 4*k2 + k3)
 * then applies X_c(t+dt) = Phi(t+dt, t) * X_c(t).
 */
inline void flow_based_mobility_returns_phi_rk3(Eigen::Ref<Eigen::VectorXd> commuter,
                                                Eigen::Ref<Eigen::VectorXd> totals,
                                                const mio::oseir::Model<double>& model, double t, double dt)
{
    const std::size_t NC      = static_cast<std::size_t>(model.populations.get_num_compartments());
    const std::size_t sys_dim = NC + NC * NC;

    Eigen::VectorXd y(sys_dim);
    y.head(NC)         = totals;
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(NC, NC);
    y.tail(NC * NC)    = Eigen::Map<const Eigen::VectorXd>(Id.data(), NC * NC);

    AugmentedPhiSystem sys(model);
    Eigen::VectorXd k1(sys_dim), k2(sys_dim), k3(sys_dim);

    sys(y, k1, t);
    sys(y + (dt * 0.5) * k1, k2, t + 0.5 * dt);
    sys(y + dt * (-k1 + 2.0 * k2), k3, t + dt);

    Eigen::VectorXd y_next = y + (dt / 6.0) * (k1 + 4.0 * k2 + k3);

    phi_apply_and_guard(commuter, totals, y_next, NC);
}

} // namespace oseir
} // namespace mio

#endif // MIO_OSEIR_STATE_ESTIMATORS_H
