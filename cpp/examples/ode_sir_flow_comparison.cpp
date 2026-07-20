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

#include "memilio/compartments/compartmental_model.h"
#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "ode_sir/infection_state.h"
#include "ode_sir/model.h"
#include "ode_sir/parameters.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

namespace
{

using FP = ScalarType;
using SirState = mio::osir::InfectionState;
using SirPop   = mio::Populations<FP, mio::AgeGroup, SirState>;
using SirParams = mio::osir::Parameters<FP>;
using SirFlows  = mio::TypeList<mio::Flow<SirState::Susceptible, SirState::Infected>,
                                mio::Flow<SirState::Infected, SirState::Recovered>>;

enum class AuxiliarySirState
{
    Susceptible,
    Infected,
    Recovered,
    CumulativeInfections,
    CumulativeRecoveries,
    Count
};

using AuxiliarySirPop = mio::Populations<FP, mio::AgeGroup, AuxiliarySirState>;

struct FlowValues {
    FP infections = 0.0;
    FP recoveries = 0.0;
};

FP total_population(Eigen::Ref<const Eigen::VectorX<FP>> y, size_t s, size_t i, size_t r)
{
    return y[s] + y[i] + y[r];
}

FP infection_coefficient(const SirParams& params, FP t, mio::AgeGroup group)
{
    const auto& contact_matrix = params.template get<mio::osir::ContactPatterns<FP>>().get_cont_freq_mat();
    const auto contact_rate = contact_matrix.get_matrix_at(mio::SimulationTime<FP>(t))(
        static_cast<Eigen::Index>((size_t)group), static_cast<Eigen::Index>((size_t)group));
    return contact_rate * params.template get<mio::osir::TransmissionProbabilityOnContact<FP>>()[group];
}

FlowValues compute_sir_rates(const SirParams& params, Eigen::Ref<const Eigen::VectorX<FP>> pop,
                             Eigen::Ref<const Eigen::VectorX<FP>> y, FP t, size_t s, size_t i, size_t r)
{
    const FP n = total_population(pop, s, i, r);
    if (n < mio::Limits<FP>::zero_tolerance()) {
        return {};
    }

    const auto group = mio::AgeGroup(0);
    const FP infection_rate =
        infection_coefficient(params, t, group) * y[s] * pop[i] / n;
    const FP recovery_rate = y[i] / params.template get<mio::osir::TimeInfected<FP>>()[group];
    return {infection_rate, recovery_rate};
}

class FlowSirModel : public mio::FlowModel<FP, SirState, SirPop, SirParams, SirFlows>
{
public:
    using Base = mio::FlowModel<FP, SirState, SirPop, SirParams, SirFlows>;

    FlowSirModel(const SirPop& initial_populations, const SirParams& model_parameters)
        : Base(initial_populations, model_parameters)
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const auto group = mio::AgeGroup(0);
        const size_t s   = this->populations.get_flat_index({group, SirState::Susceptible});
        const size_t i   = this->populations.get_flat_index({group, SirState::Infected});
        const size_t r   = this->populations.get_flat_index({group, SirState::Recovered});

        const auto rates = compute_sir_rates(this->parameters, pop, y, t, s, i, r);
        flows[this->template get_flat_flow_index<SirState::Susceptible, SirState::Infected>(group)] =
            rates.infections;
        flows[this->template get_flat_flow_index<SirState::Infected, SirState::Recovered>(group)] =
            rates.recoveries;
    }
};

class AuxiliarySirModel
    : public mio::CompartmentalModel<FP, AuxiliarySirState, AuxiliarySirPop, SirParams>
{
public:
    using Base = mio::CompartmentalModel<FP, AuxiliarySirState, AuxiliarySirPop, SirParams>;

    AuxiliarySirModel(const AuxiliarySirPop& initial_populations, const SirParams& model_parameters)
        : Base(initial_populations, model_parameters)
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                         Eigen::Ref<Eigen::VectorX<FP>> dydt) const override
    {
        const auto group = mio::AgeGroup(0);
        const size_t s = this->populations.get_flat_index({group, AuxiliarySirState::Susceptible});
        const size_t i = this->populations.get_flat_index({group, AuxiliarySirState::Infected});
        const size_t r = this->populations.get_flat_index({group, AuxiliarySirState::Recovered});
        const size_t cumulative_infections =
            this->populations.get_flat_index({group, AuxiliarySirState::CumulativeInfections});
        const size_t cumulative_recoveries =
            this->populations.get_flat_index({group, AuxiliarySirState::CumulativeRecoveries});

        const auto rates = compute_sir_rates(this->parameters, pop, y, t, s, i, r);
        dydt[s] -= rates.infections;
        dydt[i] += rates.infections - rates.recoveries;
        dydt[r] += rates.recoveries;
        dydt[cumulative_infections] += rates.infections;
        dydt[cumulative_recoveries] += rates.recoveries;
    }
};

std::unique_ptr<mio::OdeIntegratorCore<FP>> make_integrator()
{
    return std::make_unique<mio::EulerIntegratorCore<FP>>();
}

SirParams make_parameters()
{
    SirParams params(mio::AgeGroup(1));
    params.set<mio::osir::TimeInfected<FP>>(4.0);
    params.set<mio::osir::TransmissionProbabilityOnContact<FP>>(0.08);
    params.get<mio::osir::ContactPatterns<FP>>().get_cont_freq_mat()[0].get_baseline().setConstant(8.0);
    return params;
}

SirPop make_sir_populations()
{
    SirPop populations({mio::AgeGroup(1), SirState::Count});
    const auto group = mio::AgeGroup(0);
    populations[{group, SirState::Susceptible}] = 9990.0;
    populations[{group, SirState::Infected}]    = 10.0;
    populations[{group, SirState::Recovered}]   = 0.0;
    return populations;
}

AuxiliarySirPop make_auxiliary_populations()
{
    AuxiliarySirPop populations({mio::AgeGroup(1), AuxiliarySirState::Count}, 0.0);
    const auto group = mio::AgeGroup(0);
    populations[{group, AuxiliarySirState::Susceptible}] = 9990.0;
    populations[{group, AuxiliarySirState::Infected}]    = 10.0;
    populations[{group, AuxiliarySirState::Recovered}]   = 0.0;
    return populations;
}

mio::osir::Model<FP> make_classical_model()
{
    mio::osir::Model<FP> model(make_sir_populations(), make_parameters());
    model.check_constraints();
    return model;
}

FlowSirModel make_flow_model()
{
    FlowSirModel model(make_sir_populations(), make_parameters());
    model.check_constraints();
    return model;
}

AuxiliarySirModel make_auxiliary_model()
{
    AuxiliarySirModel model(make_auxiliary_populations(), make_parameters());
    model.check_constraints();
    return model;
}

FlowValues get_last_flow_values_from_flow_simulation(const mio::TimeSeries<FP>& flows)
{
    const auto& last = flows.get_last_value();
    return {last[0], last[1]};
}

FlowValues get_last_flow_values_from_auxiliary_simulation(const mio::TimeSeries<FP>& result)
{
    const auto& last = result.get_last_value();
    return {last[3], last[4]};
}

FlowValues reconstruct_flows_from_sir_balances(const mio::TimeSeries<FP>& result)
{
    const auto& initial = result.get_value(0);
    const auto& last    = result.get_last_value();
    return {initial[0] - last[0], last[2] - initial[2]};
}

Eigen::Vector2<FP> compute_rates_from_classical_state(const mio::osir::Model<FP>& model,
                                                      Eigen::Ref<const Eigen::VectorX<FP>> y, FP t)
{
    const auto group = mio::AgeGroup(0);
    const size_t s   = model.populations.get_flat_index({group, SirState::Susceptible});
    const size_t i   = model.populations.get_flat_index({group, SirState::Infected});
    const size_t r   = model.populations.get_flat_index({group, SirState::Recovered});

    const auto rates = compute_sir_rates(model.parameters, y, y, t, s, i, r);
    return Eigen::Vector2<FP>(rates.infections, rates.recoveries);
}

mio::TimeSeries<FP> reconstruct_flows_a_posteriori_trapezoid(const mio::osir::Model<FP>& model,
                                                             const mio::TimeSeries<FP>& result)
{
    mio::TimeSeries<FP> flows(result.get_time(0), Eigen::Vector2<FP>::Zero());
    for (Eigen::Index k = 1; k < result.get_num_time_points(); ++k) {
        const FP t_prev = result.get_time(k - 1);
        const FP t      = result.get_time(k);
        const FP dt     = t - t_prev;

        const auto rates_prev = compute_rates_from_classical_state(model, result.get_value(k - 1), t_prev);
        const auto rates      = compute_rates_from_classical_state(model, result.get_value(k), t);
        const auto next_value = (flows.get_last_value() + 0.5 * dt * (rates_prev + rates)).eval();
        flows.add_time_point(t, next_value);
    }
    return flows;
}

void print_flow_summary(const std::string& label, const FlowValues& values)
{
    std::cout << std::left << std::setw(36) << label << std::right << std::fixed << std::setprecision(6)
              << std::setw(16) << values.infections << std::setw(16) << values.recoveries << '\n';
}

} // namespace

int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    const FP t0   = 0.0;
    const FP tmax = 20.0;
    const FP dt   = 0.1;

    auto classical_model = make_classical_model();
    auto flow_model      = make_flow_model();
    auto auxiliary_model = make_auxiliary_model();

    auto classical_result = mio::simulate<FP>(t0, tmax, dt, classical_model, make_integrator());
    auto flow_result      = mio::simulate_flows<FP>(t0, tmax, dt, flow_model, make_integrator());
    auto auxiliary_result = mio::simulate<FP>(t0, tmax, dt, auxiliary_model, make_integrator());

    const auto direct_flows = get_last_flow_values_from_flow_simulation(flow_result[1]);
    const auto auxiliary_flows = get_last_flow_values_from_auxiliary_simulation(auxiliary_result);
    const auto balance_flows = reconstruct_flows_from_sir_balances(classical_result);
    const auto aposteriori_flows_ts =
        reconstruct_flows_a_posteriori_trapezoid(classical_model, classical_result);
    const auto aposteriori_flows = get_last_flow_values_from_flow_simulation(aposteriori_flows_ts);

    std::cout << "SIR cumulative flows without mobility, t = " << t0 << " ... " << tmax << ", dt = " << dt
              << "\n\n";
    std::cout << std::left << std::setw(36) << "Method" << std::right << std::setw(16) << "S->I"
              << std::setw(16) << "I->R" << '\n';
    std::cout << std::string(68, '-') << '\n';
    print_flow_summary("FlowSimulation", direct_flows);
    print_flow_summary("Auxiliary compartments", auxiliary_flows);
    print_flow_summary("A-posteriori trapezoid", aposteriori_flows);
    print_flow_summary("SIR balance identities", balance_flows);

    std::cout << "\nThe balance identities use S(0)-S(t) and R(t)-R(0). They work for this simple SIR chain, "
                 "but do not generalize to compartments with multiple incoming or outgoing transitions.\n";

    return 0;
}
