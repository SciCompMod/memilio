/* 
* Copyright (C) 2020-2024 MEmilio
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
#include "ode_secirvvs/parameter_space.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/model.h"

#include <numeric>
#include <algorithm>
#include <iterator>

namespace mio
{
namespace osecirvvs
{

void draw_sample_demographics(Model& model)
{
    model.parameters.get<ICUCapacity>().draw_sample();
    model.parameters.get<TestAndTraceCapacity>().draw_sample();

    std::vector<InfectionState> naive_states = {
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

    std::vector<InfectionState> partial_states = {
        InfectionState::SusceptiblePartialImmunity,        InfectionState::ExposedPartialImmunity,
        InfectionState::InfectedNoSymptomsPartialImmunity, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,
        InfectionState::InfectedSymptomsPartialImmunity,   InfectionState::InfectedSymptomsPartialImmunityConfirmed,
        InfectionState::InfectedSeverePartialImmunity,     InfectionState::InfectedCriticalPartialImmunity,
        InfectionState::TemporaryImmunPartialImmunity,     InfectionState::DeadPartialImmunity,
    };

    std::vector<InfectionState> improved_states = {
        InfectionState::SusceptibleImprovedImmunity,        InfectionState::ExposedImprovedImmunity,
        InfectionState::InfectedNoSymptomsImprovedImmunity, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed,
        InfectionState::InfectedSymptomsImprovedImmunity,   InfectionState::InfectedSymptomsImprovedImmunityConfirmed,
        InfectionState::InfectedSevereImprovedImmunity,     InfectionState::InfectedCriticalImprovedImmunity,
        InfectionState::TemporaryImmunImprovedImmunity,     InfectionState::DeadImprovedImmunity,
    };

    // helper function to calculate the total population of a layer for a given age group
    auto calculate_layer_total = [&model](const std::vector<InfectionState>& states, AgeGroup ageGroup) {
        return std::accumulate(states.begin(), states.end(), 0.0,
                               [&model, &ageGroup](double sum, const InfectionState& state) {
                                   return sum + model.populations[{ageGroup, state}];
                               });
    };

    // helper function to adjust the susceptible population of a layer for a given age group
    auto adjust_susceptible_population = [&model](AgeGroup i, double diff, InfectionState susceptibleState) {
        model.populations[{i, susceptibleState}] += diff;
        if (model.populations[{i, susceptibleState}] < 0) {
            mio::log_warning("Negative population in State " + std::to_string(static_cast<size_t>(susceptibleState)) +
                             " for age group " + std::to_string(static_cast<size_t>(i)) + ". Setting to 0.");
            model.populations[{i, susceptibleState}] = 0;
        }
    };

    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {

        double group_naive_total    = calculate_layer_total(naive_states, i);
        double group_partial_total  = calculate_layer_total(partial_states, i);
        double group_improved_total = calculate_layer_total(improved_states, i);

        //sample initial compartments (with exceptions)
        for (auto inf_state = Index<InfectionState>(0); inf_state < InfectionState::Count; ++inf_state) {
            if (inf_state != InfectionState::DeadNaive && //not sampled, fixed from data
                inf_state != InfectionState::DeadPartialImmunity && //not sampled, fixed from data
                inf_state != InfectionState::DeadImprovedImmunity) { //not sampled, fixed from data
                model.populations[{i, inf_state}].draw_sample();
            }
        }
        double diff_naive    = group_naive_total - calculate_layer_total(naive_states, i);
        double diff_partial  = group_partial_total - calculate_layer_total(partial_states, i);
        double diff_improved = group_improved_total - calculate_layer_total(improved_states, i);

        adjust_susceptible_population(i, diff_naive, InfectionState::SusceptibleNaive);
        adjust_susceptible_population(i, diff_partial, InfectionState::SusceptiblePartialImmunity);
        adjust_susceptible_population(i, diff_improved, InfectionState::SusceptibleImprovedImmunity);
    }
}

void draw_sample_infection(Model& model)
{
    model.parameters.get<Seasonality>().draw_sample();

    //not age dependent
    model.parameters.get<TimeExposed>()[AgeGroup(0)].draw_sample();
    model.parameters.get<TimeInfectedNoSymptoms>()[AgeGroup(0)].draw_sample();
    model.parameters.get<RelativeTransmissionNoSymptoms>()[AgeGroup(0)].draw_sample();
    model.parameters.get<RiskOfInfectionFromSymptomatic>()[AgeGroup(0)].draw_sample();
    model.parameters.get<MaxRiskOfInfectionFromSymptomatic>()[AgeGroup(0)].draw_sample();
    model.parameters.get<TimeTemporaryImmunityPI>()[AgeGroup(0)].draw_sample();
    model.parameters.get<TimeTemporaryImmunityII>()[AgeGroup(0)].draw_sample();

    model.parameters.get<ReducExposedPartialImmunity>()[AgeGroup(0)].draw_sample();
    model.parameters.get<ReducExposedImprovedImmunity>()[AgeGroup(0)].draw_sample();
    model.parameters.get<ReducInfectedSymptomsPartialImmunity>()[AgeGroup(0)].draw_sample();
    model.parameters.get<ReducInfectedSymptomsImprovedImmunity>()[AgeGroup(0)].draw_sample();
    model.parameters.get<ReducInfectedSevereCriticalDeadPartialImmunity>()[AgeGroup(0)].draw_sample();
    model.parameters.get<ReducInfectedSevereCriticalDeadImprovedImmunity>()[AgeGroup(0)].draw_sample();
    model.parameters.get<ReducTimeInfectedMild>()[AgeGroup(0)].draw_sample();

    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        //not age dependent
        model.parameters.get<TimeExposed>()[i]            = model.parameters.get<TimeExposed>()[AgeGroup(0)];
        model.parameters.get<TimeInfectedNoSymptoms>()[i] = model.parameters.get<TimeInfectedNoSymptoms>()[AgeGroup(0)];
        model.parameters.get<RelativeTransmissionNoSymptoms>()[i] =
            model.parameters.get<RelativeTransmissionNoSymptoms>()[AgeGroup(0)];
        model.parameters.get<RiskOfInfectionFromSymptomatic>()[i] =
            model.parameters.get<RiskOfInfectionFromSymptomatic>()[AgeGroup(0)];
        model.parameters.get<MaxRiskOfInfectionFromSymptomatic>()[i] =
            model.parameters.get<MaxRiskOfInfectionFromSymptomatic>()[AgeGroup(0)];

        model.parameters.get<ReducExposedPartialImmunity>()[i] =
            model.parameters.get<ReducExposedPartialImmunity>()[AgeGroup(0)];
        model.parameters.get<ReducExposedImprovedImmunity>()[i] =
            model.parameters.get<ReducExposedImprovedImmunity>()[AgeGroup(0)];
        model.parameters.get<ReducInfectedSymptomsPartialImmunity>()[i] =
            model.parameters.get<ReducInfectedSymptomsPartialImmunity>()[AgeGroup(0)];
        model.parameters.get<ReducInfectedSymptomsImprovedImmunity>()[i] =
            model.parameters.get<ReducInfectedSymptomsImprovedImmunity>()[AgeGroup(0)];
        model.parameters.get<ReducInfectedSevereCriticalDeadPartialImmunity>()[i] =
            model.parameters.get<ReducInfectedSevereCriticalDeadPartialImmunity>()[AgeGroup(0)];
        model.parameters.get<ReducInfectedSevereCriticalDeadImprovedImmunity>()[i] =
            model.parameters.get<ReducInfectedSevereCriticalDeadImprovedImmunity>()[AgeGroup(0)];
        model.parameters.get<ReducTimeInfectedMild>()[i] = model.parameters.get<ReducTimeInfectedMild>()[AgeGroup(0)];

        //age dependent
        model.parameters.get<TimeInfectedSymptoms>()[i].draw_sample();
        model.parameters.get<TimeInfectedSevere>()[i].draw_sample();
        model.parameters.get<TimeInfectedCritical>()[i].draw_sample();

        model.parameters.get<TransmissionProbabilityOnContact>()[i].draw_sample();
        model.parameters.get<RecoveredPerInfectedNoSymptoms>()[i].draw_sample();
        model.parameters.get<DeathsPerCritical>()[i].draw_sample();
        model.parameters.get<SeverePerInfectedSymptoms>()[i].draw_sample();
        model.parameters.get<CriticalPerSevere>()[i].draw_sample();
    }
}

void draw_sample(Model& model)
{
    draw_sample_infection(model);
    draw_sample_demographics(model);
    model.parameters.get<ContactPatterns>().draw_sample();
    model.apply_constraints();
}

Graph<Model, MigrationParameters> draw_sample(Graph<Model, MigrationParameters>& graph)
{
    Graph<Model, MigrationParameters> sampled_graph;

    //sample global parameters
    auto& shared_params_model = graph.nodes()[0].property;
    draw_sample_infection(shared_params_model);
    auto& shared_contacts = shared_params_model.parameters.template get<ContactPatterns>();
    shared_contacts.draw_sample_dampings();
    auto& shared_dynamic_npis = shared_params_model.parameters.template get<DynamicNPIsInfectedSymptoms>();
    shared_dynamic_npis.draw_sample();

    for (auto& params_node : graph.nodes()) {
        auto& node_model = params_node.property;

        //sample local parameters
        draw_sample_demographics(params_node.property);

        //copy global parameters
        //save demographic parameters so they aren't overwritten
        auto local_icu_capacity = node_model.parameters.template get<ICUCapacity>();
        auto local_tnt_capacity = node_model.parameters.template get<TestAndTraceCapacity>();
        auto local_holidays     = node_model.parameters.template get<ContactPatterns>().get_school_holidays();
        auto local_daily_v1     = node_model.parameters.template get<DailyPartialVaccination>();
        auto local_daily_v2     = node_model.parameters.template get<DailyFullVaccination>();
        auto local_daily_v3     = node_model.parameters.template get<DailyBoosterVaccination>();
        node_model.parameters   = shared_params_model.parameters;
        node_model.parameters.template get<ICUCapacity>()                           = local_icu_capacity;
        node_model.parameters.template get<TestAndTraceCapacity>()                  = local_tnt_capacity;
        node_model.parameters.template get<ContactPatterns>().get_school_holidays() = local_holidays;
        node_model.parameters.template get<DailyPartialVaccination>()               = local_daily_v1;
        node_model.parameters.template get<DailyFullVaccination>()                  = local_daily_v2;
        node_model.parameters.template get<DailyBoosterVaccination>()               = local_daily_v3;

        node_model.parameters.template get<ContactPatterns>().make_matrix();
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

} // namespace osecirvvs
} // namespace mio
