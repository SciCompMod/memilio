#ifndef MIO_EXAMPLES_STATE_ESTIMATORS_H
#define MIO_EXAMPLES_STATE_ESTIMATORS_H

#include "memilio/compartments/flow_simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/type_list.h"
#include "memilio/math/eigen.h"
#include "ode_seir/model.h"

namespace mio
{
namespace examples
{

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

        // Calculate total population across all commuter types
        FP N = 0;
        for (mio::AgeGroup g(0); g < mio::AgeGroup(params.get_num_groups()); ++g) {
            // Start with non-commuter group
            N += pop[this->populations.get_flat_index({g, CommuterType::NonCommuter, InfectionStateExplicit::S})];
            N += pop[this->populations.get_flat_index({g, CommuterType::NonCommuter, InfectionStateExplicit::E})];
            N += pop[this->populations.get_flat_index({g, CommuterType::NonCommuter, InfectionStateExplicit::I})];
            N += pop[this->populations.get_flat_index({g, CommuterType::NonCommuter, InfectionStateExplicit::R})];

            // Handle commuter groups
            for (int i = 0; i < m_num_commuter_groups; ++i) {
                CommuterType commuter_type =
                    static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + i);
                for (size_t state = 0; state < static_cast<size_t>(InfectionStateExplicit::Count); ++state) {
                    N += pop[this->populations.get_flat_index(
                        {g, commuter_type, static_cast<InfectionStateExplicit>(state)})];
                }
            }
        }

        FP invN = (N < mio::Limits<FP>::zero_tolerance()) ? FP(0) : FP(1) / N;

        // Calculate total infected population
        FP total_infected = 0;
        for (mio::AgeGroup g(0); g < mio::AgeGroup(params.get_num_groups()); ++g) {
            // Add infected from non-commuter group
            total_infected +=
                pop[this->populations.get_flat_index({g, CommuterType::NonCommuter, InfectionStateExplicit::I})];

            // Add infected from commuter groups
            for (int i = 0; i < m_num_commuter_groups; ++i) {
                CommuterType commuter_type =
                    static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + i);
                total_infected += pop[this->populations.get_flat_index({g, commuter_type, InfectionStateExplicit::I})];
            }
        }

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
                CommuterType commuter_type_i =
                    static_cast<CommuterType>(static_cast<int>(CommuterType::CommuterBase) + c);
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
                    flows[S_to_E_flow_idx] += coeffStoE * y[Si] * total_infected_j;
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

void integrate_mobile_population_high_order(Eigen::Ref<Eigen::VectorXd> mobile_population, const SimType& sim,
                                            Eigen::Ref<const Eigen::VectorXd> total, ScalarType t, ScalarType dt)
{
    auto y0         = mobile_population.eval();
    auto y1         = mobile_population;
    auto integrator = std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
    integrator->step(
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
    size_t current_idx = static_cast<size_t>(std::round(t / dt_flows));
    size_t next_idx    = static_cast<size_t>(std::round((t + dt) / dt_flows));

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

void probabilistic_mobility_returns(Eigen::Ref<Eigen::VectorXd> mobile_population, const SimType& sim,
                                    Eigen::Ref<const Eigen::VectorXd> total, ScalarType t, ScalarType dt)
{
    mio::unused(t);
    const auto& model       = sim.get_model();
    const size_t num_groups = static_cast<size_t>(model.parameters.get_num_groups());
    using InfState          = mio::oseir::InfectionState;

    Eigen::VectorXd mobile_pop_initial = mobile_population;

    for (mio::AgeGroup i(0); i < mio::AgeGroup(num_groups); ++i) {
        size_t S_i = model.populations.get_flat_index({i, InfState::Susceptible});
        size_t E_i = model.populations.get_flat_index({i, InfState::Exposed});
        size_t I_i = model.populations.get_flat_index({i, InfState::Infected});
        size_t R_i = model.populations.get_flat_index({i, InfState::Recovered});

        ScalarType time_exposed  = model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[i];
        ScalarType time_infected = model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[i];

        ScalarType force_of_infection = 0.0;
        for (mio::AgeGroup j(0); j < mio::AgeGroup(num_groups); ++j) {
            const size_t Ij = model.populations.get_flat_index({j, InfState::Infected});

            ScalarType Nj = 0;
            for (size_t comp = 0; comp < static_cast<size_t>(InfState::Count); ++comp) {
                Nj += total[model.populations.get_flat_index({j, static_cast<InfState>(comp)})];
            }

            const ScalarType divNj = (Nj < mio::Limits<ScalarType>::zero_tolerance()) ? 0.0 : 1.0 / Nj;
            const ScalarType coeffStoE =
                model.parameters.template get<mio::oseir::ContactPatterns<ScalarType>>()
                    .get_cont_freq_mat()
                    .get_matrix_at(t)(i.get(), j.get()) *
                model.parameters.template get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[i] * divNj;

            force_of_infection += coeffStoE * total[Ij];
        }

        double p_S_to_S = std::exp(-force_of_infection * dt);
        double p_E_to_E = (time_exposed > 1e-10) ? std::exp(-dt / time_exposed) : 1.0;
        double p_I_to_I = (time_infected > 1e-10) ? std::exp(-dt / time_infected) : 1.0;

        p_S_to_S = std::max(0.0, std::min(1.0, p_S_to_S));
        p_E_to_E = std::max(0.0, std::min(1.0, p_E_to_E));
        p_I_to_I = std::max(0.0, std::min(1.0, p_I_to_I));

        double p_S_to_E = 1.0 - p_S_to_S;
        double p_E_to_I = 1.0 - p_E_to_E;
        double p_I_to_R = 1.0 - p_I_to_I;

        mobile_population[S_i] = mobile_pop_initial[S_i] * p_S_to_S;
        mobile_population[E_i] = mobile_pop_initial[E_i] * p_E_to_E + mobile_pop_initial[S_i] * p_S_to_E;
        mobile_population[I_i] = mobile_pop_initial[I_i] * p_I_to_I + mobile_pop_initial[E_i] * p_E_to_I;
        mobile_population[R_i] = mobile_pop_initial[R_i] + mobile_pop_initial[I_i] * p_I_to_R;
    }
    mobile_population = mobile_population.cwiseMax(0.0);
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

} // namespace examples
} // namespace mio

#endif // MIO_EXAMPLES_STATE_ESTIMATORS_H
