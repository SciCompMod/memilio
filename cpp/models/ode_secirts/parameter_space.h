/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Henrik Zunker, Wadim Koslow, Daniel Abele, Martin J. Kühn
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
#ifndef MIO_ODE_SECIRTS_PARAMETER_SPACE_H
#define MIO_ODE_SECIRTS_PARAMETER_SPACE_H

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/memory.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/math/math_utils.h"
#include "ode_secirts/model.h"

#include <assert.h>
#include <string>
#include <vector>
#include <random>
#include <memory>

namespace mio
{
namespace osecirts
{
/**
 * Draws a sample from the specified distributions for all parameters
 * related to the demographics, e.g., population.
 * @tparam FP Floating point type, e.g., double.
 * @param[inout] model Model including contact patterns for all age groups
 */
template <typename FP>
void draw_sample_demographics(Model<FP>& model)
{
    model.parameters.template get<ICUCapacity<FP>>().draw_sample();
    model.parameters.template get<TestAndTraceCapacity<FP>>().draw_sample();

    const static std::vector<InfectionState> naive_immunity_states = {
        InfectionState::SusceptibleNaive,
        InfectionState::ExposedNaive,
        InfectionState::InfectedNoSymptomsNaive,
        InfectionState::InfectedNoSymptomsNaiveConfirmed,
        InfectionState::InfectedSymptomsNaive,
        InfectionState::InfectedSymptomsNaiveConfirmed,
        InfectionState::InfectedSevereNaive,
        InfectionState::InfectedCriticalNaive,
        InfectionState::DeadNaive,
    };

    const static std::vector<InfectionState> partial_immunity_states = {
        InfectionState::SusceptiblePartialImmunity,        InfectionState::ExposedPartialImmunity,
        InfectionState::InfectedNoSymptomsPartialImmunity, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,
        InfectionState::InfectedSymptomsPartialImmunity,   InfectionState::InfectedSymptomsPartialImmunityConfirmed,
        InfectionState::InfectedSeverePartialImmunity,     InfectionState::InfectedCriticalPartialImmunity,
        InfectionState::TemporaryImmunePartialImmunity,    InfectionState::DeadPartialImmunity,
    };

    const static std::vector<InfectionState> improved_immunity_states = {
        InfectionState::SusceptibleImprovedImmunity,        InfectionState::ExposedImprovedImmunity,
        InfectionState::InfectedNoSymptomsImprovedImmunity, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed,
        InfectionState::InfectedSymptomsImprovedImmunity,   InfectionState::InfectedSymptomsImprovedImmunityConfirmed,
        InfectionState::InfectedSevereImprovedImmunity,     InfectionState::InfectedCriticalImprovedImmunity,
        InfectionState::TemporaryImmuneImprovedImmunity,    InfectionState::DeadImprovedImmunity,
    };

    // helper function to calculate the total population of a layer for a given age group
    auto calculate_layer_total = [&model](const std::vector<InfectionState>& states, AgeGroup ageGroup) {
        return std::accumulate(states.begin(), states.end(), FP(0.0),
                               [&model, &ageGroup](FP sum, const InfectionState& state) {
                                   return evaluate_intermediate<FP>(sum + model.populations[{ageGroup, state}]);
                               });
    };

    // helper function to adjust the susceptible population of a layer for a given age group
    auto adjust_susceptible_population = [&model](AgeGroup i, FP diff, InfectionState susceptibleState) {
        model.populations[{i, susceptibleState}] += diff;
        if (model.populations[{i, susceptibleState}] < 0) {
            mio::log_warning("Negative population in State " + std::to_string(static_cast<size_t>(susceptibleState)) +
                             " for age group " + std::to_string(static_cast<size_t>(i)) + ". Setting to 0.");
            model.populations[{i, susceptibleState}] = 0;
        }
    };

    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {

        const FP group_naive_total    = calculate_layer_total(naive_immunity_states, i);
        const FP group_partial_total  = calculate_layer_total(partial_immunity_states, i);
        const FP group_improved_total = calculate_layer_total(improved_immunity_states, i);

        //sample initial compartments (with exceptions)
        for (auto inf_state = Index<InfectionState>(0); inf_state < InfectionState::Count; ++inf_state) {
            if (inf_state != InfectionState::DeadNaive && //not sampled, fixed from data
                inf_state != InfectionState::DeadPartialImmunity && //not sampled, fixed from data
                inf_state != InfectionState::DeadImprovedImmunity) { //not sampled, fixed from data
                model.populations[{i, inf_state}].draw_sample();
            }
        }
        const FP diff_naive    = group_naive_total - calculate_layer_total(naive_immunity_states, i);
        const FP diff_partial  = group_partial_total - calculate_layer_total(partial_immunity_states, i);
        const FP diff_improved = group_improved_total - calculate_layer_total(improved_immunity_states, i);

        adjust_susceptible_population(i, diff_naive, InfectionState::SusceptibleNaive);
        adjust_susceptible_population(i, diff_partial, InfectionState::SusceptiblePartialImmunity);
        adjust_susceptible_population(i, diff_improved, InfectionState::SusceptibleImprovedImmunity);
    }
}

/**
 * Draws a sample from the specified distributions for all parameters
 * related to the infection.
 *
 * @tparam FP Floating point type, e.g., double.
 * @param[inout] model Model including contact patterns for all age groups.
 */
template <typename FP>
void draw_sample_infection(Model<FP>& model)
{
    model.parameters.template get<Seasonality<FP>>().draw_sample();

    //not age dependent
    model.parameters.template get<TimeExposed<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<TimeInfectedNoSymptoms<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<RelativeTransmissionNoSymptoms<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<RiskOfInfectionFromSymptomatic<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<TimeTemporaryImmunityPI<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<TimeTemporaryImmunityII<FP>>()[AgeGroup(0)].draw_sample();

    model.parameters.template get<ReducExposedPartialImmunity<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<ReducExposedImprovedImmunity<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<ReducTimeInfectedMild<FP>>()[AgeGroup(0)].draw_sample();

    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        //not age dependent
        model.parameters.template get<TimeExposed<FP>>()[i] =
            model.parameters.template get<TimeExposed<FP>>()[AgeGroup(0)];
        model.parameters.template get<TimeInfectedNoSymptoms<FP>>()[i] =
            model.parameters.template get<TimeInfectedNoSymptoms<FP>>()[AgeGroup(0)];
        model.parameters.template get<RelativeTransmissionNoSymptoms<FP>>()[i] =
            model.parameters.template get<RelativeTransmissionNoSymptoms<FP>>()[AgeGroup(0)];
        model.parameters.template get<RiskOfInfectionFromSymptomatic<FP>>()[i] =
            model.parameters.template get<RiskOfInfectionFromSymptomatic<FP>>()[AgeGroup(0)];
        model.parameters.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[i] =
            model.parameters.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[AgeGroup(0)];

        model.parameters.template get<ReducExposedPartialImmunity<FP>>()[i] =
            model.parameters.template get<ReducExposedPartialImmunity<FP>>()[AgeGroup(0)];
        model.parameters.template get<ReducExposedImprovedImmunity<FP>>()[i] =
            model.parameters.template get<ReducExposedImprovedImmunity<FP>>()[AgeGroup(0)];
        model.parameters.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[i] =
            model.parameters.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[AgeGroup(0)];
        model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[i] =
            model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[AgeGroup(0)];
        model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i] =
            model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[AgeGroup(0)];
        model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i] =
            model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[AgeGroup(0)];
        model.parameters.template get<ReducTimeInfectedMild<FP>>()[i] =
            model.parameters.template get<ReducTimeInfectedMild<FP>>()[AgeGroup(0)];

        //age dependent
        model.parameters.template get<TimeInfectedSymptoms<FP>>()[i].draw_sample();
        model.parameters.template get<TimeInfectedSevere<FP>>()[i].draw_sample();
        model.parameters.template get<TimeInfectedCritical<FP>>()[i].draw_sample();

        model.parameters.template get<TransmissionProbabilityOnContact<FP>>()[i].draw_sample();
        model.parameters.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i].draw_sample();
        model.parameters.template get<DeathsPerCritical<FP>>()[i].draw_sample();
        model.parameters.template get<SeverePerInfectedSymptoms<FP>>()[i].draw_sample();
        model.parameters.template get<CriticalPerSevere<FP>>()[i].draw_sample();
    }
}

/**
 * Draws a sample from model parameter distributions and stores sample values
 * as parameters values (cf. UncertainValue and Parameters classes).
 *
 * @tparam FP Floating point type, e.g., double.
 * @param[inout] model Model including contact patterns for all age groups.
 */
template <typename FP>
void draw_sample(Model<FP>& model)
{
    draw_sample_infection<FP>(model);
    draw_sample_demographics<FP>(model);
    model.parameters.template get<ContactPatterns<FP>>().draw_sample();
    model.apply_constraints();
}

/**
 * Draws samples for each model node in a graph.
 * Some parameters are shared between nodes and are only sampled once.
 *
 * @tparam FP Floating point type, e.g., double.
 * @param graph Graph to be sampled.
 * @return Graph with nodes and edges from the input graph sampled.
 */
template <typename FP>
Graph<Model<FP>, MobilityParameters<FP>> draw_sample(Graph<Model<FP>, MobilityParameters<FP>>& graph)
{
    Graph<Model<FP>, MobilityParameters<FP>> sampled_graph;

    //sample global parameters
    auto& shared_params_model = graph.nodes()[0].property;
    draw_sample_infection<FP>(shared_params_model);
    auto& shared_contacts = shared_params_model.parameters.template get<ContactPatterns<FP>>();
    shared_contacts.draw_sample_dampings();
    auto& shared_dynamic_npis = shared_params_model.parameters.template get<DynamicNPIsInfectedSymptoms<FP>>();
    shared_dynamic_npis.draw_sample();

    for (auto& params_node : graph.nodes()) {
        auto& node_model = params_node.property;

        //sample local parameters
        draw_sample_demographics<FP>(params_node.property);

        //copy global parameters
        //save demographic parameters so they aren't overwritten
        auto local_icu_capacity = node_model.parameters.template get<ICUCapacity<FP>>();
        auto local_tnt_capacity = node_model.parameters.template get<TestAndTraceCapacity<FP>>();
        auto local_holidays     = node_model.parameters.template get<ContactPatterns<FP>>().get_school_holidays();
        auto local_daily_v1     = node_model.parameters.template get<DailyPartialVaccinations<FP>>();
        auto local_daily_v2     = node_model.parameters.template get<DailyFullVaccinations<FP>>();
        auto local_daily_v3     = node_model.parameters.template get<DailyBoosterVaccinations<FP>>();
        node_model.parameters   = shared_params_model.parameters;
        node_model.parameters.template get<ICUCapacity<FP>>()                           = local_icu_capacity;
        node_model.parameters.template get<TestAndTraceCapacity<FP>>()                  = local_tnt_capacity;
        node_model.parameters.template get<ContactPatterns<FP>>().get_school_holidays() = local_holidays;
        node_model.parameters.template get<DailyPartialVaccinations<FP>>()              = local_daily_v1;
        node_model.parameters.template get<DailyFullVaccinations<FP>>()                 = local_daily_v2;
        node_model.parameters.template get<DailyBoosterVaccinations<FP>>()              = local_daily_v3;

        node_model.parameters.template get<ContactPatterns<FP>>().make_matrix();
        node_model.apply_constraints();

        sampled_graph.add_node(params_node.id, node_model);
    }

    for (auto& edge : graph.edges()) {
        auto edge_params = edge.property;
        //no dynamic NPIs
        //TODO: add switch to optionally enable dynamic NPIs to edges
        sampled_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge_params);
    }

    return sampled_graph;
}

} // namespace osecirts
} // namespace mio

#endif // MIO_ODE_SECIRTS_PARAMETER_SPACE_H
