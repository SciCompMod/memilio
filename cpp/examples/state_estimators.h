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

// create model with explicit commuter compartments

enum class InfectionStateExplicit
{
    S_NonCommuter = 0,
    E_NonCommuter,
    I_NonCommuter,
    R_NonCommuter,
    S_Commuter,
    E_Commuter,
    I_Commuter,
    R_Commuter,

    Count
};

using ParametersExplicit = mio::oseir::Parameters<ScalarType>;

// clang‑format off
using FlowsExplicit =
    mio::TypeList<mio::Flow<InfectionStateExplicit::S_NonCommuter, InfectionStateExplicit::E_NonCommuter>,
                  mio::Flow<InfectionStateExplicit::E_NonCommuter, InfectionStateExplicit::I_NonCommuter>,
                  mio::Flow<InfectionStateExplicit::I_NonCommuter, InfectionStateExplicit::R_NonCommuter>,
                  mio::Flow<InfectionStateExplicit::S_Commuter, InfectionStateExplicit::E_Commuter>,
                  mio::Flow<InfectionStateExplicit::E_Commuter, InfectionStateExplicit::I_Commuter>,
                  mio::Flow<InfectionStateExplicit::I_Commuter, InfectionStateExplicit::R_Commuter>>;
// clang‑format on

template <typename FP = ScalarType>
class ModelExplicit
    : public mio::FlowModel<FP, InfectionStateExplicit, mio::Populations<FP, mio::AgeGroup, InfectionStateExplicit>,
                            ParametersExplicit, FlowsExplicit>
{
    using Base = mio::FlowModel<FP, InfectionStateExplicit, mio::Populations<FP, mio::AgeGroup, InfectionStateExplicit>,
                                ParametersExplicit, FlowsExplicit>;

public:
    using typename Base::ParameterSet;
    using typename Base::Populations;

    ModelExplicit(int num_agegroups)
        : Base(Populations({mio::AgeGroup(num_agegroups), InfectionStateExplicit::Count}),
               ParametersExplicit(num_agegroups))
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const auto& params = this->parameters;
        mio::AgeGroup g(0);

        const size_t SN = this->populations.get_flat_index({g, InfectionStateExplicit::S_NonCommuter});
        const size_t EN = this->populations.get_flat_index({g, InfectionStateExplicit::E_NonCommuter});
        const size_t IN = this->populations.get_flat_index({g, InfectionStateExplicit::I_NonCommuter});
        const size_t RN = this->populations.get_flat_index({g, InfectionStateExplicit::R_NonCommuter});
        const size_t SC = this->populations.get_flat_index({g, InfectionStateExplicit::S_Commuter});
        const size_t EC = this->populations.get_flat_index({g, InfectionStateExplicit::E_Commuter});
        const size_t IC = this->populations.get_flat_index({g, InfectionStateExplicit::I_Commuter});
        const size_t RC = this->populations.get_flat_index({g, InfectionStateExplicit::R_Commuter});

        FP N    = pop[SN] + pop[EN] + pop[IN] + pop[RN] + pop[SC] + pop[EC] + pop[IC] + pop[RC];
        FP invN = (N < mio::Limits<FP>::zero_tolerance()) ? FP(0) : FP(1) / N;

        FP beta    = params.template get<mio::oseir::TransmissionProbabilityOnContact<FP>>()[g];
        FP Te      = params.template get<mio::oseir::TimeExposed<FP>>()[g];
        FP Ti      = params.template get<mio::oseir::TimeInfected<FP>>()[g];
        auto contM = params.template get<mio::oseir::ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(t);

        FP lambda = contM(g.get(), g.get()) * beta * (pop[IN] + pop[IC]) * invN;

        // S→E flows
        flows[Base::template get_flat_flow_index<InfectionStateExplicit::S_NonCommuter,
                                                 InfectionStateExplicit::E_NonCommuter>(g)] = lambda * y[SN];
        flows[Base::template get_flat_flow_index<InfectionStateExplicit::S_Commuter,
                                                 InfectionStateExplicit::E_Commuter>(g)]    = lambda * y[SC];

        // E→I flows
        flows[Base::template get_flat_flow_index<InfectionStateExplicit::E_NonCommuter,
                                                 InfectionStateExplicit::I_NonCommuter>(g)] = (FP(1) / Te) * y[EN];
        flows[Base::template get_flat_flow_index<InfectionStateExplicit::E_Commuter,
                                                 InfectionStateExplicit::I_Commuter>(g)]    = (FP(1) / Te) * y[EC];

        // I→R flows
        flows[Base::template get_flat_flow_index<InfectionStateExplicit::I_NonCommuter,
                                                 InfectionStateExplicit::R_NonCommuter>(g)] = (FP(1) / Ti) * y[IN];
        flows[Base::template get_flat_flow_index<InfectionStateExplicit::I_Commuter,
                                                 InfectionStateExplicit::R_Commuter>(g)]    = (FP(1) / Ti) * y[IC];
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

    for (mio::AgeGroup group(0); group < mio::AgeGroup(num_groups); ++group) {
        size_t flow_base_idx = static_cast<size_t>(group) * 3;
        size_t S_i           = model.populations.get_flat_index({group, InfState::Susceptible});
        size_t E_i           = model.populations.get_flat_index({group, InfState::Exposed});
        size_t I_i           = model.populations.get_flat_index({group, InfState::Infected});
        size_t R_i           = model.populations.get_flat_index({group, InfState::Recovered});

        // Calculate compartment-specific scaling factors.
        double scale_S = (total[S_i] > 1e-10) ? (mobile_population[S_i] / total[S_i]) : 0.0;
        double scale_E = (total[E_i] > 1e-10) ? (mobile_population[E_i] / total[E_i]) : 0.0;
        double scale_I = (total[I_i] > 1e-10) ? (mobile_population[I_i] / total[I_i]) : 0.0;

        scale_S = std::max(0.0, std::min(1.0, scale_S));
        scale_E = std::max(0.0, std::min(1.0, scale_E));
        scale_I = std::max(0.0, std::min(1.0, scale_I));

        // Estimate flows for the mobile population using compartment-specific scaling
        double flow_SE_mobile = diff[flow_base_idx] * scale_S;
        double flow_EI_mobile = diff[flow_base_idx + 1] * scale_E;
        double flow_IR_mobile = diff[flow_base_idx + 2] * scale_I;

        // Apply the estimated flows to mobile population
        mobile_population[S_i] -= flow_SE_mobile;
        mobile_population[E_i] += flow_SE_mobile - flow_EI_mobile;
        mobile_population[I_i] += flow_EI_mobile - flow_IR_mobile;
        mobile_population[R_i] += flow_IR_mobile;
    }
    // Ensure no negative populations due to numerical issues
    mobile_population = mobile_population.cwiseMax(0.0);
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
