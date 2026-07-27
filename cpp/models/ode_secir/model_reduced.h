/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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

#ifndef MIO_MODELS_OSECIR_MODEL_REDUCED_H
#define MIO_MODELS_OSECIR_MODEL_REDUCED_H

#include "memilio/compartments/flow_model.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/math/smoother.h"
#include "ode_secir/parameters.h"

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace mio
{
namespace osecir
{
namespace reduced
{

enum class InfectionState
{
    Susceptible,
    Exposed,
    InfectedNoSymptoms,
    InfectedSymptoms,
    InfectedSevere,
    InfectedCritical,
    Recovered,
    Dead,
    Count
};

enum class Transition
{
    Infection,
    ExposedToNoSymptoms,
    NoSymptomsToSymptoms,
    NoSymptomsToRecovered,
    SymptomsToSevere,
    SymptomsToRecovered,
    SevereToCritical,
    SevereToRecovered,
    CriticalToDead,
    CriticalToRecovered
};

inline constexpr std::array<Transition, 10> transitions{
    Transition::Infection,
    Transition::ExposedToNoSymptoms,
    Transition::NoSymptomsToSymptoms,
    Transition::NoSymptomsToRecovered,
    Transition::SymptomsToSevere,
    Transition::SymptomsToRecovered,
    Transition::SevereToCritical,
    Transition::SevereToRecovered,
    Transition::CriticalToDead,
    Transition::CriticalToRecovered,
};

using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::Exposed>,
                       Flow<InfectionState::Exposed, InfectionState::InfectedNoSymptoms>,
                       Flow<InfectionState::InfectedNoSymptoms, InfectionState::InfectedSymptoms>,
                       Flow<InfectionState::InfectedNoSymptoms, InfectionState::Recovered>,
                       Flow<InfectionState::InfectedSymptoms, InfectionState::InfectedSevere>,
                       Flow<InfectionState::InfectedSymptoms, InfectionState::Recovered>,
                       Flow<InfectionState::InfectedSevere, InfectionState::InfectedCritical>,
                       Flow<InfectionState::InfectedSevere, InfectionState::Recovered>,
                       Flow<InfectionState::InfectedCritical, InfectionState::Dead>,
                       Flow<InfectionState::InfectedCritical, InfectionState::Recovered>>;

template <typename FP>
class Model : public FlowModel<FP, InfectionState, Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>
{
    using PopulationType = mio::Populations<FP, AgeGroup, InfectionState>;
    using ParameterType  = Parameters<FP>;
    using Base           = FlowModel<FP, InfectionState, PopulationType, ParameterType, Flows>;

public:
    Model(const PopulationType& initial_populations, const ParameterType& model_parameters)
        : Base(initial_populations, model_parameters)
    {
    }

    explicit Model(int num_age_groups)
        : Model(PopulationType({AgeGroup(num_age_groups), InfectionState::Count}),
                ParameterType(AgeGroup(num_age_groups)))
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const auto num_groups            = this->parameters.get_num_groups();
        const FP test_and_trace_required = compute_test_and_trace_required(pop);

        for (AgeGroup group = 0; group < num_groups; ++group) {
            const auto exposed     = state_index(group, InfectionState::Exposed);
            const auto no_symptoms = state_index(group, InfectionState::InfectedNoSymptoms);
            const auto symptoms    = state_index(group, InfectionState::InfectedSymptoms);
            const auto severe      = state_index(group, InfectionState::InfectedSevere);
            const auto critical    = state_index(group, InfectionState::InfectedCritical);

            const FP exposed_rate          = 1.0 / this->parameters.template get<TimeExposed<FP>>()[group];
            const FP no_symptoms_rate      = 1.0 / this->parameters.template get<TimeInfectedNoSymptoms<FP>>()[group];
            const FP symptoms_rate         = 1.0 / this->parameters.template get<TimeInfectedSymptoms<FP>>()[group];
            const FP severe_rate           = 1.0 / this->parameters.template get<TimeInfectedSevere<FP>>()[group];
            const FP critical_rate         = 1.0 / this->parameters.template get<TimeInfectedCritical<FP>>()[group];
            const FP recovered_no_symptoms = this->parameters.template get<RecoveredPerInfectedNoSymptoms<FP>>()[group];
            const FP severe_symptoms       = this->parameters.template get<SeverePerInfectedSymptoms<FP>>()[group];
            const FP critical_severe       = this->parameters.template get<CriticalPerSevere<FP>>()[group];
            const FP deaths_critical       = this->parameters.template get<DeathsPerCritical<FP>>()[group];

            flows[static_cast<Eigen::Index>(get_transition_index(Transition::Infection, group))] =
                compute_infection_rate(group, pop, y, t, test_and_trace_required);
            flows[static_cast<Eigen::Index>(get_transition_index(Transition::ExposedToNoSymptoms, group))] =
                exposed_rate * y[exposed];
            flows[static_cast<Eigen::Index>(get_transition_index(Transition::NoSymptomsToSymptoms, group))] =
                (1.0 - recovered_no_symptoms) * no_symptoms_rate * y[no_symptoms];
            flows[static_cast<Eigen::Index>(get_transition_index(Transition::NoSymptomsToRecovered, group))] =
                recovered_no_symptoms * no_symptoms_rate * y[no_symptoms];
            flows[static_cast<Eigen::Index>(get_transition_index(Transition::SymptomsToSevere, group))] =
                severe_symptoms * symptoms_rate * y[symptoms];
            flows[static_cast<Eigen::Index>(get_transition_index(Transition::SymptomsToRecovered, group))] =
                (1.0 - severe_symptoms) * symptoms_rate * y[symptoms];
            flows[static_cast<Eigen::Index>(get_transition_index(Transition::SevereToCritical, group))] =
                critical_severe * severe_rate * y[severe];
            flows[static_cast<Eigen::Index>(get_transition_index(Transition::SevereToRecovered, group))] =
                (1.0 - critical_severe) * severe_rate * y[severe];
            flows[static_cast<Eigen::Index>(get_transition_index(Transition::CriticalToDead, group))] =
                deaths_critical * critical_rate * y[critical];
            flows[static_cast<Eigen::Index>(get_transition_index(Transition::CriticalToRecovered, group))] =
                (1.0 - deaths_critical) * critical_rate * y[critical];
        }
    }

    FP get_transition_rate(Transition transition, AgeGroup group, Eigen::Ref<const Eigen::VectorX<FP>> pop,
                           Eigen::Ref<const Eigen::VectorX<FP>> y, FP t) const
    {
        const auto value = [&](InfectionState state) {
            return y[state_index(group, state)];
        };

        switch (transition) {
        case Transition::Infection:
            return compute_infection_rate(group, pop, y, t, compute_test_and_trace_required(pop));
        case Transition::ExposedToNoSymptoms:
            return value(InfectionState::Exposed) / this->parameters.template get<TimeExposed<FP>>()[group];
        case Transition::NoSymptomsToSymptoms:
            return (1.0 - this->parameters.template get<RecoveredPerInfectedNoSymptoms<FP>>()[group]) *
                   value(InfectionState::InfectedNoSymptoms) /
                   this->parameters.template get<TimeInfectedNoSymptoms<FP>>()[group];
        case Transition::NoSymptomsToRecovered:
            return this->parameters.template get<RecoveredPerInfectedNoSymptoms<FP>>()[group] *
                   value(InfectionState::InfectedNoSymptoms) /
                   this->parameters.template get<TimeInfectedNoSymptoms<FP>>()[group];
        case Transition::SymptomsToSevere:
            return this->parameters.template get<SeverePerInfectedSymptoms<FP>>()[group] *
                   value(InfectionState::InfectedSymptoms) /
                   this->parameters.template get<TimeInfectedSymptoms<FP>>()[group];
        case Transition::SymptomsToRecovered:
            return (1.0 - this->parameters.template get<SeverePerInfectedSymptoms<FP>>()[group]) *
                   value(InfectionState::InfectedSymptoms) /
                   this->parameters.template get<TimeInfectedSymptoms<FP>>()[group];
        case Transition::SevereToCritical:
            return this->parameters.template get<CriticalPerSevere<FP>>()[group] *
                   value(InfectionState::InfectedSevere) /
                   this->parameters.template get<TimeInfectedSevere<FP>>()[group];
        case Transition::SevereToRecovered:
            return (1.0 - this->parameters.template get<CriticalPerSevere<FP>>()[group]) *
                   value(InfectionState::InfectedSevere) /
                   this->parameters.template get<TimeInfectedSevere<FP>>()[group];
        case Transition::CriticalToDead:
            return this->parameters.template get<DeathsPerCritical<FP>>()[group] *
                   value(InfectionState::InfectedCritical) /
                   this->parameters.template get<TimeInfectedCritical<FP>>()[group];
        case Transition::CriticalToRecovered:
            return (1.0 - this->parameters.template get<DeathsPerCritical<FP>>()[group]) *
                   value(InfectionState::InfectedCritical) /
                   this->parameters.template get<TimeInfectedCritical<FP>>()[group];
        }
        throw std::runtime_error("Unhandled reduced SECIR transition.");
    }

    size_t get_transition_index(Transition transition, AgeGroup group) const
    {
        switch (transition) {
        case Transition::Infection:
            return this->template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>({group});
        case Transition::ExposedToNoSymptoms:
            return this->template get_flat_flow_index<InfectionState::Exposed, InfectionState::InfectedNoSymptoms>(
                {group});
        case Transition::NoSymptomsToSymptoms:
            return this
                ->template get_flat_flow_index<InfectionState::InfectedNoSymptoms, InfectionState::InfectedSymptoms>(
                    {group});
        case Transition::NoSymptomsToRecovered:
            return this->template get_flat_flow_index<InfectionState::InfectedNoSymptoms, InfectionState::Recovered>(
                {group});
        case Transition::SymptomsToSevere:
            return this->template get_flat_flow_index<InfectionState::InfectedSymptoms, InfectionState::InfectedSevere>(
                {group});
        case Transition::SymptomsToRecovered:
            return this->template get_flat_flow_index<InfectionState::InfectedSymptoms, InfectionState::Recovered>(
                {group});
        case Transition::SevereToCritical:
            return this->template get_flat_flow_index<InfectionState::InfectedSevere, InfectionState::InfectedCritical>(
                {group});
        case Transition::SevereToRecovered:
            return this->template get_flat_flow_index<InfectionState::InfectedSevere, InfectionState::Recovered>(
                {group});
        case Transition::CriticalToDead:
            return this->template get_flat_flow_index<InfectionState::InfectedCritical, InfectionState::Dead>({group});
        case Transition::CriticalToRecovered:
            return this->template get_flat_flow_index<InfectionState::InfectedCritical, InfectionState::Recovered>(
                {group});
        }
        throw std::runtime_error("Unhandled reduced SECIR transition.");
    }

private:
    size_t state_index(AgeGroup group, InfectionState state) const
    {
        return this->populations.get_flat_index({group, state});
    }

    FP compute_test_and_trace_required(Eigen::Ref<const Eigen::VectorX<FP>> pop) const
    {
        FP required = 0.0;
        for (AgeGroup group = 0; group < this->parameters.get_num_groups(); ++group) {
            required += (1.0 - this->parameters.template get<RecoveredPerInfectedNoSymptoms<FP>>()[group]) *
                        pop[state_index(group, InfectionState::InfectedNoSymptoms)] /
                        this->parameters.template get<TimeInfectedNoSymptoms<FP>>()[group];
        }
        return required;
    }

    FP compute_infection_rate(AgeGroup group, Eigen::Ref<const Eigen::VectorX<FP>> pop,
                              Eigen::Ref<const Eigen::VectorX<FP>> y, FP t, FP test_and_trace_required) const
    {
        const auto susceptible = state_index(group, InfectionState::Susceptible);
        const FP seasonality   = 1.0 + this->parameters.template get<Seasonality<FP>>() *
                                         std::sin(std::numbers::pi_v<FP> *
                                                  ((this->parameters.template get<StartDay<FP>>() + t) / 182.5 + 0.5));
        const ContactMatrixGroup<FP>& contact_matrix = this->parameters.template get<ContactPatterns<FP>>();

        FP rate = 0.0;
        for (AgeGroup source_group = 0; source_group < this->parameters.get_num_groups(); ++source_group) {
            const FP source_population = pop[state_index(source_group, InfectionState::Susceptible)] +
                                         pop[state_index(source_group, InfectionState::Exposed)] +
                                         pop[state_index(source_group, InfectionState::InfectedNoSymptoms)] +
                                         pop[state_index(source_group, InfectionState::InfectedSymptoms)] +
                                         pop[state_index(source_group, InfectionState::InfectedSevere)] +
                                         pop[state_index(source_group, InfectionState::InfectedCritical)] +
                                         pop[state_index(source_group, InfectionState::Recovered)];
            const FP inverse_population =
                source_population < Limits<FP>::zero_tolerance() ? FP(0.0) : FP(1.0 / source_population);
            const FP symptomatic_risk = smoother_cosine<FP>(
                test_and_trace_required, this->parameters.template get<TestAndTraceCapacity<FP>>(),
                this->parameters.template get<TestAndTraceCapacity<FP>>() *
                    this->parameters.template get<TestAndTraceCapacityMaxRisk<FP>>(),
                this->parameters.template get<RiskOfInfectionFromSymptomatic<FP>>()[source_group],
                this->parameters.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[source_group]);
            const FP contact_rate = seasonality * contact_matrix.get_matrix_at(SimulationTime<FP>(t))(
                                                      static_cast<Eigen::Index>((size_t)group),
                                                      static_cast<Eigen::Index>((size_t)source_group));
            rate += y[susceptible] * contact_rate * inverse_population *
                    this->parameters.template get<TransmissionProbabilityOnContact<FP>>()[group] *
                    (this->parameters.template get<RelativeTransmissionNoSymptoms<FP>>()[source_group] *
                         pop[state_index(source_group, InfectionState::InfectedNoSymptoms)] +
                     symptomatic_risk * pop[state_index(source_group, InfectionState::InfectedSymptoms)]);
        }
        return rate;
    }
};

} // namespace reduced
} // namespace osecir
} // namespace mio

#endif
