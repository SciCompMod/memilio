#ifndef MIO_EXAMPLES_STATE_ESTIMATORS_H
#define MIO_EXAMPLES_STATE_ESTIMATORS_H

#include "memilio/compartments/flow_simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/type_list.h"
#include "memilio/math/eigen.h"
#include "memilio/math/stepper_wrapper.h"
#include "ode_seir/model.h"

namespace mio
{
namespace examples
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

private:
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

        // Calculate and apply flows for each age group and commuter type
        for (mio::AgeGroup g(0); g < mio::AgeGroup(params.get_num_groups()); ++g) {
            // Process E->I and I->R flows for non-commuter group first
            {
                CommuterType commuter_type = CommuterType::NonCommuter;
                const size_t E_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::E});
                const size_t I_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::I});

                // E->I flows
                flows[Base::template get_flat_flow_index<InfectionStateExplicit::E, InfectionStateExplicit::I>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeExposed<FP>>()[g]) * y[E_idx];

                // I->R flows
                flows[Base::template get_flat_flow_index<InfectionStateExplicit::I, InfectionStateExplicit::R>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeInfected<FP>>()[g]) * y[I_idx];
            }

            // Process E->I and I->R flows for commuter groups
            for (int c = 0; c < m_num_commuter_groups; ++c) {
                CommuterType commuter_type =
                    static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c);

                const size_t E_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::E});
                const size_t I_idx = this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::I});

                // E->I flows
                flows[Base::template get_flat_flow_index<InfectionStateExplicit::E, InfectionStateExplicit::I>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeExposed<FP>>()[g]) * y[E_idx];

                // I->R flows
                flows[Base::template get_flat_flow_index<InfectionStateExplicit::I, InfectionStateExplicit::R>(
                    {g, commuter_type})] = (FP(1) / params.template get<mio::oseir::TimeInfected<FP>>()[g]) * y[I_idx];
            }
        }

        for (mio::AgeGroup i(0); i < mio::AgeGroup(params.get_num_groups()); ++i) {
            // For each commuter type of age group i
            // First non-commuter
            CommuterType commuter_type_i = CommuterType::NonCommuter;
            const size_t Si = this->populations.get_flat_index({i, commuter_type_i, InfectionStateExplicit::S});
            auto S_to_E_flow_idx =
                Base::template get_flat_flow_index<InfectionStateExplicit::S, InfectionStateExplicit::E>(
                    {i, commuter_type_i});

            // Initialize flow to zero
            flows[S_to_E_flow_idx] = 0;

            // Sum up contributions from all age groups j
            for (mio::AgeGroup j(0); j < mio::AgeGroup(params.get_num_groups()); ++j) {
                // Calculate total N_j for group j (across all commuter types)
                FP N_j = 0;

                // Add non-commuter population for group j
                for (size_t state = 0; state < static_cast<size_t>(InfectionStateExplicit::Count); ++state) {
                    N_j += pop[this->populations.get_flat_index(
                        {j, CommuterType::NonCommuter, static_cast<InfectionStateExplicit>(state)})];
                }

                // Add commuter populations for group j
                for (int c = 0; c < m_num_commuter_groups; ++c) {
                    CommuterType commuter_type_j =
                        static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c);
                    for (size_t state = 0; state < static_cast<size_t>(InfectionStateExplicit::Count); ++state) {
                        N_j += pop[this->populations.get_flat_index(
                            {j, commuter_type_j, static_cast<InfectionStateExplicit>(state)})];
                    }
                }

                // Calculate divN_j (1/N_j) with safety check
                const FP divN_j = (N_j < mio::Limits<FP>::zero_tolerance()) ? 0.0 : FP(1.0) / N_j;

                // Get transmission coefficient
                const FP coeffStoE =
                    params.template get<mio::oseir::ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(t)(
                        i.get(), j.get()) *
                    params.template get<mio::oseir::TransmissionProbabilityOnContact<FP>>()[i] * divN_j;

                // Sum up infected from all commuter types in group j
                FP total_infected_j =
                    pop[this->populations.get_flat_index({j, CommuterType::NonCommuter, InfectionStateExplicit::I})];
                for (int c = 0; c < m_num_commuter_groups; ++c) {
                    CommuterType commuter_type_j =
                        static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c);
                    total_infected_j +=
                        pop[this->populations.get_flat_index({j, commuter_type_j, InfectionStateExplicit::I})];
                }

                // Add contribution to S->E flow
                flows[S_to_E_flow_idx] += coeffStoE * y[Si] * total_infected_j;
            }

            // Now do the same for each commuter type
            for (int c = 0; c < m_num_commuter_groups; ++c) {
                CommuterType commuter_type_i_inner =
                    static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c);
                const size_t Si_inner =
                    this->populations.get_flat_index({i, commuter_type_i_inner, InfectionStateExplicit::S});
                auto S_to_E_flow_idx_inner =
                    Base::template get_flat_flow_index<InfectionStateExplicit::S, InfectionStateExplicit::E>(
                        {i, commuter_type_i_inner});

                // Initialize flow to zero
                flows[S_to_E_flow_idx_inner] = 0;

                // Sum up contributions from all age groups j
                for (mio::AgeGroup j(0); j < mio::AgeGroup(params.get_num_groups()); ++j) {
                    // Calculate total N_j for group j (across all commuter types)
                    FP N_j = 0;

                    // Add non-commuter population for group j
                    for (size_t state = 0; state < static_cast<size_t>(InfectionStateExplicit::Count); ++state) {
                        N_j += pop[this->populations.get_flat_index(
                            {j, CommuterType::NonCommuter, static_cast<InfectionStateExplicit>(state)})];
                    }

                    // Add commuter populations for group j
                    for (int c2 = 0; c2 < m_num_commuter_groups; ++c2) {
                        CommuterType commuter_type_j =
                            static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c2);
                        for (size_t state = 0; state < static_cast<size_t>(InfectionStateExplicit::Count); ++state) {
                            N_j += pop[this->populations.get_flat_index(
                                {j, commuter_type_j, static_cast<InfectionStateExplicit>(state)})];
                        }
                    }

                    const FP divN_j = (N_j < mio::Limits<FP>::zero_tolerance()) ? 0.0 : FP(1.0) / N_j;

                    const FP coeffStoE =
                        params.template get<mio::oseir::ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(t)(
                            i.get(), j.get()) *
                        params.template get<mio::oseir::TransmissionProbabilityOnContact<FP>>()[i] * divN_j;

                    // Sum up infected from all commuter types in group j
                    FP total_infected_j = pop[this->populations.get_flat_index(
                        {j, CommuterType::NonCommuter, InfectionStateExplicit::I})];
                    for (int c2 = 0; c2 < m_num_commuter_groups; ++c2) {
                        CommuterType commuter_type_j =
                            static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c2);
                        total_infected_j +=
                            pop[this->populations.get_flat_index({j, commuter_type_j, InfectionStateExplicit::I})];
                    }

                    // Add contribution to S->E flow
                    flows[S_to_E_flow_idx_inner] += coeffStoE * y[Si_inner] * total_infected_j;
                }
            }
        }
    }
};

using ExplicitModel = ModelExplicit<ScalarType>;
using ExplicitSim   = mio::FlowSimulation<ScalarType, ExplicitModel>;

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

inline void
flow_based_mobility_returns_rk4(Eigen::Ref<Eigen::VectorXd> commuter, // in/out: X_c(t) -> X_c(t+dt)
                                Eigen::Ref<Eigen::VectorXd> totals, // in/out: \widetilde X(t) -> \widetilde X(t+dt)
                                const mio::oseir::Model<double>& model, // model supplies totals-RHS and totals-flows
                                double t, // current time
                                double dt) // step size
{
    using IS             = mio::oseir::InfectionState;
    const std::size_t G  = static_cast<std::size_t>(model.parameters.get_num_groups());
    const std::size_t NC = 4 * G; // S,E,I,R per age group

    assert((std::size_t)totals.size() == NC && "totals size must be 4*G");
    assert((std::size_t)commuter.size() == NC && "commuter size must be 4*G");

    // ------- Helpers: index maps S,E,I,R for each age group g -------
    auto idxS = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Susceptible});
    };
    auto idxE = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Exposed});
    };
    auto idxI = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Infected});
    };
    auto idxR = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Recovered});
    };

    // ------- lambdas for totals derivatives & flows -------
    auto totals_rhs = [&](const Eigen::VectorXd& y, double tt, Eigen::VectorXd& dy) {
        dy.resize(NC);
        // In the SEIR integrators in MEmilio: get_derivatives(pop, y, t, dy)
        // Here, locals+commuters are already aggregated in y.
        model.get_derivatives(y, y, tt, dy);
    };

    // flows vector stores (per group): [ f_SE, f_EI, f_IR ], length = 3*G
    auto totals_flows = [&](const Eigen::VectorXd& y, double tt, Eigen::VectorXd& flows) {
        flows.resize(3 * (int)G);
        flows.setZero();
        model.get_flows(y, y, tt, flows);
    };

    // ------- Commuter derivative from Left-Share at a given stage -------
    // Given totals flows (f_SE, f_EI, f_IR per group) and shares ξ = Xc / Xt at that stage,
    // compute dXc/dt per compartment using only scaled totals flows:
    //
    //   S: dXc_S = - ξ_S * f_SE
    //   E: dXc_E = + ξ_S * f_SE - ξ_E * f_EI
    //   I: dXc_I = + ξ_E * f_EI - ξ_I * f_IR
    //   R: dXc_R = + ξ_I * f_IR
    //
    auto commuter_rhs_from_shares = [&](const Eigen::VectorXd& flows, const Eigen::VectorXd& xi, Eigen::VectorXd& dxc) {
        dxc.setZero(NC);
        for (std::size_t g = 0; g < G; ++g) {
            const int oS = idxS(g), oE = idxE(g), oI = idxI(g), oR = idxR(g);
            const int of = 3 * (int)g;

            const double fSE = flows[of + 0];
            const double fEI = flows[of + 1];
            const double fIR = flows[of + 2];

            const double xiS = xi[oS];
            const double xiE = xi[oE];
            const double xiI = xi[oI];
            // xiR not needed in RHS (no outflow from R in SEIR)

            dxc[oS] = -xiS * fSE;
            dxc[oE] = xiS * fSE - xiE * fEI;
            dxc[oI] = xiE * fEI - xiI * fIR;
            dxc[oR] = xiI * fIR;
        }
    };

    // ------- share vector ξ = Xc / Xt  (safe division, clamped to [0,1]) -------
    auto compute_shares = [&](const Eigen::VectorXd& Xc, const Eigen::VectorXd& Xt, Eigen::VectorXd& xi) {
        xi.resize(NC);
        for (int i = 0; i < (int)NC; ++i) {
            const double denom = Xt[i];
            double v           = (denom > 0.0) ? (Xc[i] / denom) : 0.0;
            // numerical safety
            if (v < 0.0)
                v = 0.0;
            if (v > 1.0)
                v = 1.0;
            xi[i] = v;
        }
    };

    // ------- Storage for RK4 stages -------
    Eigen::VectorXd k1_tot(NC), k2_tot(NC), k3_tot(NC), k4_tot(NC);
    Eigen::VectorXd f1(3 * (int)G), f2(3 * (int)G), f3(3 * (int)G), f4(3 * (int)G);

    Eigen::VectorXd k1_com(NC), k2_com(NC), k3_com(NC), k4_com(NC);
    Eigen::VectorXd xi1(NC), xi2(NC), xi3(NC), xi4(NC);

    // stage states
    Eigen::VectorXd y1  = totals;
    Eigen::VectorXd Xc0 = commuter;

    // --- Stage 1 (t) ---
    totals_rhs(y1, t, k1_tot);
    totals_flows(y1, t, f1);
    compute_shares(Xc0, y1, xi1);
    commuter_rhs_from_shares(f1, xi1, k1_com);

    // --- Stage 2 (t + dt/2) ---
    Eigen::VectorXd y2 = y1 + (dt * 0.5) * k1_tot;
    // commuter stage-2 state via RK4 relation (NO rhs eval for Xc needed here)
    Eigen::VectorXd Xc2 = Xc0 + (dt * 0.5) * k1_com;
    totals_rhs(y2, t + 0.5 * dt, k2_tot);
    totals_flows(y2, t + 0.5 * dt, f2);
    compute_shares(Xc2, y2, xi2);
    commuter_rhs_from_shares(f2, xi2, k2_com);

    // --- Stage 3 (t + dt/2) ---
    Eigen::VectorXd y3  = y1 + (dt * 0.5) * k2_tot;
    Eigen::VectorXd Xc3 = Xc0 + (dt * 0.5) * k2_com;
    totals_rhs(y3, t + 0.5 * dt, k3_tot);
    totals_flows(y3, t + 0.5 * dt, f3);
    compute_shares(Xc3, y3, xi3);
    commuter_rhs_from_shares(f3, xi3, k3_com);

    // --- Stage 4 (t + dt) ---
    Eigen::VectorXd y4  = y1 + dt * k3_tot;
    Eigen::VectorXd Xc4 = Xc0 + dt * k3_com;
    totals_rhs(y4, t + dt, k4_tot);
    totals_flows(y4, t + dt, f4);
    compute_shares(Xc4, y4, xi4);
    commuter_rhs_from_shares(f4, xi4, k4_com);

    // --- Combine (classical RK4 weights) ---
    totals = y1 + (dt / 6.0) * (k1_tot + 2.0 * k2_tot + 2.0 * k3_tot + k4_tot);

    Eigen::VectorXd Xc_next = Xc0 + (dt / 6.0) * (k1_com + 2.0 * k2_com + 2.0 * k3_com + k4_com);

    // Numerical guards: non-negativity and commuter <= totals (component-wise)
    for (int i = 0; i < (int)NC; ++i) {
        if (Xc_next[i] < 0.0)
            Xc_next[i] = 0.0;
        // enforce Xc <= totals (avoid tiny drifts)
        if (totals[i] < 0.0)
            totals[i] = 0.0; // paranoid clamp for totals as well
        if (Xc_next[i] > totals[i])
            Xc_next[i] = totals[i];
    }

    commuter = Xc_next;
}

// ============================================================================
// RK2 Implementation
// ============================================================================
// Classical Runge-Kutta 2nd order method (Midpoint):
//   k1 = f(t, y)
//   k2 = f(t + dt/2, y + dt/2 * k1)
//   y_{n+1} = y_n + dt * k2
//
// For the coupled totals+commuter system:
//   - Totals: standard RK2 on model.get_derivatives()
//   - Commuter: RK2 on share-scaled flows (ξ computed at each stage)
//
inline void
flow_based_mobility_returns_rk2(Eigen::Ref<Eigen::VectorXd> commuter, // in/out: X_c(t) -> X_c(t+dt)
                                Eigen::Ref<Eigen::VectorXd> totals, // in/out: \widetilde X(t) -> \widetilde X(t+dt)
                                const mio::oseir::Model<double>& model, // model supplies totals-RHS and totals-flows
                                double t, // current time
                                double dt) // step size
{
    using IS             = mio::oseir::InfectionState;
    const std::size_t G  = static_cast<std::size_t>(model.parameters.get_num_groups());
    const std::size_t NC = 4 * G; // S,E,I,R per age group

    assert((std::size_t)totals.size() == NC && "totals size must be 4*G");
    assert((std::size_t)commuter.size() == NC && "commuter size must be 4*G");

    // ------- Helpers: index maps S,E,I,R for each age group g -------
    auto idxS = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Susceptible});
    };
    auto idxE = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Exposed});
    };
    auto idxI = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Infected});
    };
    auto idxR = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Recovered});
    };

    // ------- lambdas for totals derivatives & flows -------
    auto totals_rhs = [&](const Eigen::VectorXd& y, double tt, Eigen::VectorXd& dy) {
        dy.resize(NC);
        model.get_derivatives(y, y, tt, dy);
    };

    // flows vector stores (per group): [ f_SE, f_EI, f_IR ], length = 3*G
    auto totals_flows = [&](const Eigen::VectorXd& y, double tt, Eigen::VectorXd& flows) {
        flows.resize(3 * (int)G);
        flows.setZero();
        model.get_flows(y, y, tt, flows);
    };

    // ------- Commuter derivative from Left-Share at a given stage -------
    // Given totals flows (f_SE, f_EI, f_IR per group) and shares ξ = Xc / Xt at that stage,
    // compute dXc/dt per compartment using only scaled totals flows:
    //
    //   S: dXc_S = - ξ_S * f_SE
    //   E: dXc_E = + ξ_S * f_SE - ξ_E * f_EI
    //   I: dXc_I = + ξ_E * f_EI - ξ_I * f_IR
    //   R: dXc_R = + ξ_I * f_IR
    //
    auto commuter_rhs_from_shares = [&](const Eigen::VectorXd& flows, const Eigen::VectorXd& xi, Eigen::VectorXd& dxc) {
        dxc.setZero(NC);
        for (std::size_t g = 0; g < G; ++g) {
            const int oS = idxS(g), oE = idxE(g), oI = idxI(g), oR = idxR(g);
            const int of = 3 * (int)g;

            const double fSE = flows[of + 0];
            const double fEI = flows[of + 1];
            const double fIR = flows[of + 2];

            const double xiS = xi[oS];
            const double xiE = xi[oE];
            const double xiI = xi[oI];

            dxc[oS] = -xiS * fSE;
            dxc[oE] = xiS * fSE - xiE * fEI;
            dxc[oI] = xiE * fEI - xiI * fIR;
            dxc[oR] = xiI * fIR;
        }
    };

    // ------- share vector ξ = Xc / Xt  (safe division, clamped to [0,1]) -------
    auto compute_shares = [&](const Eigen::VectorXd& Xc, const Eigen::VectorXd& Xt, Eigen::VectorXd& xi) {
        xi.resize(NC);
        for (int i = 0; i < (int)NC; ++i) {
            const double denom = Xt[i];
            double v           = (denom > 0.0) ? (Xc[i] / denom) : 0.0;
            // numerical safety: clamp to [0,1]
            if (v < 0.0)
                v = 0.0;
            if (v > 1.0)
                v = 1.0;
            xi[i] = v;
        }
    };

    // ------- Storage for RK2 stages -------
    Eigen::VectorXd k1_tot(NC), k2_tot(NC);
    Eigen::VectorXd f1(3 * (int)G), f2(3 * (int)G);
    Eigen::VectorXd k1_com(NC), k2_com(NC);
    Eigen::VectorXd xi1(NC), xi2(NC);

    // Initial states
    Eigen::VectorXd y0  = totals;
    Eigen::VectorXd Xc0 = commuter;

    // --- Stage 1: Evaluate at (t, y0, Xc0) ---
    totals_rhs(y0, t, k1_tot);
    totals_flows(y0, t, f1);
    compute_shares(Xc0, y0, xi1);
    commuter_rhs_from_shares(f1, xi1, k1_com);

    // --- Stage 2: Evaluate at midpoint (t + dt/2, y_mid, Xc_mid) ---
    // Midpoint prediction for totals
    Eigen::VectorXd y_mid = y0 + (dt * 0.5) * k1_tot;

    // Midpoint prediction for commuter
    Eigen::VectorXd Xc_mid = Xc0 + (dt * 0.5) * k1_com;

    // Evaluate derivatives at midpoint
    totals_rhs(y_mid, t + 0.5 * dt, k2_tot);
    totals_flows(y_mid, t + 0.5 * dt, f2);
    compute_shares(Xc_mid, y_mid, xi2);
    commuter_rhs_from_shares(f2, xi2, k2_com);

    // --- Final update: use k2 (derivative at midpoint) for full step ---
    totals                  = y0 + dt * k2_tot;
    Eigen::VectorXd Xc_next = Xc0 + dt * k2_com;

    // Numerical guards: non-negativity and commuter <= totals (component-wise)
    for (int i = 0; i < (int)NC; ++i) {
        if (Xc_next[i] < 0.0)
            Xc_next[i] = 0.0;
        if (totals[i] < 0.0)
            totals[i] = 0.0;
        if (Xc_next[i] > totals[i])
            Xc_next[i] = totals[i];
    }

    commuter = Xc_next;
}

// ============================================================================
// RK3 (Kutta's Third-Order Method) Implementation
// ============================================================================
// Classical Runge-Kutta 3rd order method:
//   k1 = f(t, y)
//   k2 = f(t + dt/2, y + dt/2 * k1)
//   k3 = f(t + dt, y - dt*k1 + 2*dt*k2)
//   y_{n+1} = y_n + dt/6 * (k1 + 4*k2 + k3)
//
// For the coupled totals+commuter system:
//   - Totals: standard RK3 on model.get_derivatives()
//   - Commuter: RK3 on share-scaled flows (ξ computed at each stage)
//
inline void
flow_based_mobility_returns_rk3(Eigen::Ref<Eigen::VectorXd> commuter, // in/out: X_c(t) -> X_c(t+dt)
                                Eigen::Ref<Eigen::VectorXd> totals, // in/out: \widetilde X(t) -> \widetilde X(t+dt)
                                const mio::oseir::Model<double>& model, // model supplies totals-RHS and totals-flows
                                double t, // current time
                                double dt) // step size
{
    using IS             = mio::oseir::InfectionState;
    const std::size_t G  = static_cast<std::size_t>(model.parameters.get_num_groups());
    const std::size_t NC = 4 * G; // S,E,I,R per age group

    assert((std::size_t)totals.size() == NC && "totals size must be 4*G");
    assert((std::size_t)commuter.size() == NC && "commuter size must be 4*G");

    // ------- Helpers: index maps S,E,I,R for each age group g -------
    auto idxS = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Susceptible});
    };
    auto idxE = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Exposed});
    };
    auto idxI = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Infected});
    };
    auto idxR = [&](std::size_t g) {
        return (int)model.populations.get_flat_index({mio::AgeGroup(g), IS::Recovered});
    };

    // ------- lambdas for totals derivatives & flows -------
    auto totals_rhs = [&](const Eigen::VectorXd& y, double tt, Eigen::VectorXd& dy) {
        dy.resize(NC);
        model.get_derivatives(y, y, tt, dy);
    };

    // flows vector stores (per group): [ f_SE, f_EI, f_IR ], length = 3*G
    auto totals_flows = [&](const Eigen::VectorXd& y, double tt, Eigen::VectorXd& flows) {
        flows.resize(3 * (int)G);
        flows.setZero();
        model.get_flows(y, y, tt, flows);
    };

    // ------- Commuter derivative from Left-Share at a given stage -------
    // Given totals flows (f_SE, f_EI, f_IR per group) and shares ξ = Xc / Xt at that stage,
    // compute dXc/dt per compartment using only scaled totals flows:
    //
    //   S: dXc_S = - ξ_S * f_SE
    //   E: dXc_E = + ξ_S * f_SE - ξ_E * f_EI
    //   I: dXc_I = + ξ_E * f_EI - ξ_I * f_IR
    //   R: dXc_R = + ξ_I * f_IR
    //
    auto commuter_rhs_from_shares = [&](const Eigen::VectorXd& flows, const Eigen::VectorXd& xi, Eigen::VectorXd& dxc) {
        dxc.setZero(NC);
        for (std::size_t g = 0; g < G; ++g) {
            const int oS = idxS(g), oE = idxE(g), oI = idxI(g), oR = idxR(g);
            const int of = 3 * (int)g;

            const double fSE = flows[of + 0];
            const double fEI = flows[of + 1];
            const double fIR = flows[of + 2];

            const double xiS = xi[oS];
            const double xiE = xi[oE];
            const double xiI = xi[oI];

            dxc[oS] = -xiS * fSE;
            dxc[oE] = xiS * fSE - xiE * fEI;
            dxc[oI] = xiE * fEI - xiI * fIR;
            dxc[oR] = xiI * fIR;
        }
    };

    // ------- share vector ξ = Xc / Xt  (safe division, clamped to [0,1]) -------
    auto compute_shares = [&](const Eigen::VectorXd& Xc, const Eigen::VectorXd& Xt, Eigen::VectorXd& xi) {
        xi.resize(NC);
        for (int i = 0; i < (int)NC; ++i) {
            const double denom = Xt[i];
            double v           = (denom > 0.0) ? (Xc[i] / denom) : 0.0;
            // numerical safety: clamp to [0,1]
            if (v < 0.0)
                v = 0.0;
            if (v > 1.0)
                v = 1.0;
            xi[i] = v;
        }
    };

    // ------- Storage for RK3 stages -------
    Eigen::VectorXd k1_tot(NC), k2_tot(NC), k3_tot(NC);
    Eigen::VectorXd f1(3 * (int)G), f2(3 * (int)G), f3(3 * (int)G);
    Eigen::VectorXd k1_com(NC), k2_com(NC), k3_com(NC);
    Eigen::VectorXd xi1(NC), xi2(NC), xi3(NC);

    // Initial states
    Eigen::VectorXd y0  = totals;
    Eigen::VectorXd Xc0 = commuter;

    // --- Stage 1: Evaluate at (t, y0, Xc0) ---
    totals_rhs(y0, t, k1_tot);
    totals_flows(y0, t, f1);
    compute_shares(Xc0, y0, xi1);
    commuter_rhs_from_shares(f1, xi1, k1_com);

    // --- Stage 2: Evaluate at midpoint (t + dt/2, y2, Xc2) ---
    // RK3 uses midpoint prediction like RK2 for stage 2
    Eigen::VectorXd y2  = y0 + (dt * 0.5) * k1_tot;
    Eigen::VectorXd Xc2 = Xc0 + (dt * 0.5) * k1_com;

    totals_rhs(y2, t + 0.5 * dt, k2_tot);
    totals_flows(y2, t + 0.5 * dt, f2);
    compute_shares(Xc2, y2, xi2);
    commuter_rhs_from_shares(f2, xi2, k2_com);

    // --- Stage 3: Evaluate at endpoint (t + dt, y3, Xc3) ---
    // RK3 formula: y3 = y0 - dt*k1 + 2*dt*k2 = y0 + dt*(-k1 + 2*k2)
    Eigen::VectorXd y3  = y0 + dt * (-k1_tot + 2.0 * k2_tot);
    Eigen::VectorXd Xc3 = Xc0 + dt * (-k1_com + 2.0 * k2_com);

    totals_rhs(y3, t + dt, k3_tot);
    totals_flows(y3, t + dt, f3);
    compute_shares(Xc3, y3, xi3);
    commuter_rhs_from_shares(f3, xi3, k3_com);

    // --- Final update: weighted combination of all three stages ---
    // RK3 weights: (k1 + 4*k2 + k3) / 6
    totals                  = y0 + (dt / 6.0) * (k1_tot + 4.0 * k2_tot + k3_tot);
    Eigen::VectorXd Xc_next = Xc0 + (dt / 6.0) * (k1_com + 4.0 * k2_com + k3_com);

    // Numerical guards: non-negativity and commuter <= totals (component-wise)
    for (int i = 0; i < (int)NC; ++i) {
        if (Xc_next[i] < 0.0)
            Xc_next[i] = 0.0;
        if (totals[i] < 0.0)
            totals[i] = 0.0;
        if (Xc_next[i] > totals[i])
            Xc_next[i] = totals[i];
    }

    commuter = Xc_next;
}

} // namespace examples
} // namespace mio

#endif // MIO_EXAMPLES_STATE_ESTIMATORS_H
