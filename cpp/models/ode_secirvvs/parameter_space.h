/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kühn
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
#ifndef ODESECIRVVS_PARAMETER_SPACE_H
#define ODESECIRVVS_PARAMETER_SPACE_H

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/logging.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/infection_state.h"

#include <assert.h>

namespace mio
{
namespace osecirvvs
{
/**
     * draws a sample from the specified distributions for all parameters related to the demographics, e.g. population.
     * @tparam FP floating point type, e.g., double
     * @param[inout] model Model including contact patterns for alle age groups
     */
template <typename FP = double>
void draw_sample_demographics(Model<FP>& model)
{
    model.parameters.template get<ICUCapacity<FP>>().draw_sample();
    model.parameters.template get<TestAndTraceCapacity<FP>>().draw_sample();

    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        double group_total = model.populations.get_group_total(i);

        //sample initial compartments (with exceptions)
        for (auto inf_state = Index<InfectionState>(0); inf_state < InfectionState::Count; ++inf_state) {
            if (inf_state != InfectionState::SusceptibleNaive && //not sampled, fixed after sampling everything else
                inf_state != InfectionState::DeadNaive && //not sampled, fixed from data
                inf_state != InfectionState::DeadPartialImmunity && //not sampled, fixed from data
                inf_state != InfectionState::DeadImprovedImmunity) { //not sampled, fixed from data
                model.populations[{i, inf_state}].draw_sample();
            }
        }

        //set susceptibles so the total number stays the same as before sampling.
        //if the new total without susceptibles is already bigger than the previous total
        //subtract the overflow from SusceptibleImprovedImmunity, susceptibles will then be approximately zero.
        model.populations[{i, InfectionState::SusceptibleNaive}] = 0;
        double diff                                              = model.populations.get_group_total(i) - group_total;
        if (diff > 0) {
            model.populations[{i, InfectionState::SusceptibleImprovedImmunity}] -= diff;
            if (model.populations[{i, InfectionState::SusceptibleImprovedImmunity}] < 0.0) {
                log_error("Negative Compartment after sampling.");
            }
            assert(std::abs(group_total - model.populations.get_group_total(i)) < 1e-10 && "Sanity check.");
        }
        model.populations.template set_difference_from_group_total<AgeGroup>({i, InfectionState::SusceptibleNaive},
                                                                             group_total);
    }
}

/**
     * draws a sample from the specified distributions for all parameters related to the infection.
     * @tparam FP floating point type, e.g., double
     * @param[inout] model Model including contact patterns for alle age groups
     */
template <typename FP = double>
void draw_sample_infection(Model<FP>& model)
{
    model.parameters.template get<Seasonality<FP>>().draw_sample();

    //not age dependent
    model.parameters.template get<TimeExposed<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<TimeInfectedNoSymptoms<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<RelativeTransmissionNoSymptoms<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<RiskOfInfectionFromSymptomatic<FP>>()[AgeGroup(0)].draw_sample();
    model.parameters.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[AgeGroup(0)].draw_sample();

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

/** Draws a sample from Model parameter distributions and stores sample values
    * as Parameters parameter values (cf. UncertainValue and Parameters classes)
    * @tparam FP floating point type, e.g., double
    * @param[inout] model Model including contact patterns for alle age groups
    */
template <typename FP = double>
void draw_sample(Model<FP>& model)
{
    draw_sample_infection(model);
    draw_sample_demographics(model);
    model.parameters.template get<ContactPatterns<FP>>().draw_sample();
    model.apply_constraints();
}

/**
    * Draws samples for each model node in a graph.
    * Some parameters are shared between nodes and only sampled once.
    * @tparam FP floating point type, e.g., double
    * @param graph Graph to be sampled.
    * @param variant_high If true, use high value for infectiousness of variant.
    * @return Graph with nodes and edges from the input graph sampled.
    */
template <typename FP = double>
Graph<Model<FP>, MobilityParameters<FP>> draw_sample(Graph<Model<FP>, MobilityParameters<FP>>& graph, bool variant_high)
{
    Graph<Model<FP>, MobilityParameters<FP>> sampled_graph;

    //sample global parameters
    auto& shared_params_model = graph.nodes()[0].property;
    draw_sample_infection(shared_params_model);
    auto& shared_contacts = shared_params_model.parameters.template get<ContactPatterns<FP>>();
    shared_contacts.draw_sample_dampings();
    auto& shared_dynamic_npis = shared_params_model.parameters.template get<DynamicNPIsInfectedSymptoms<FP>>();
    shared_dynamic_npis.draw_sample();
    auto& shared_dynamic_npis_delay = shared_params_model.parameters.template get<DynamicNPIsImplementationDelay<FP>>();
    shared_dynamic_npis_delay.draw_sample();

    double delta_fac;
    if (variant_high) {
        delta_fac = 1.6;
    }
    else {
        delta_fac = 1.4;
    }

    //infectiousness of virus variants is not sampled independently but depend on base infectiousness
    for (auto i = AgeGroup(0); i < shared_params_model.parameters.get_num_groups(); ++i) {
        shared_params_model.parameters.template get<InfectiousnessNewVariant<FP>>()[i] = delta_fac;
    }

    for (auto& params_node : graph.nodes()) {
        auto& node_model = params_node.property;

        //sample local parameters
        draw_sample_demographics(params_node.property);

        //copy global parameters
        //save demographic parameters so they aren't overwritten
        auto local_icu_capacity = node_model.parameters.template get<ICUCapacity<FP>>();
        auto local_tnt_capacity = node_model.parameters.template get<TestAndTraceCapacity<FP>>();
        auto local_holidays     = node_model.parameters.template get<ContactPatterns<FP>>().get_school_holidays();
        auto local_daily_v1     = node_model.parameters.template get<DailyFirstVaccination<FP>>();
        auto local_daily_v2     = node_model.parameters.template get<DailyFullVaccination<FP>>();
        node_model.parameters   = shared_params_model.parameters;
        node_model.parameters.template get<ICUCapacity<FP>>()                           = local_icu_capacity;
        node_model.parameters.template get<TestAndTraceCapacity<FP>>()                  = local_tnt_capacity;
        node_model.parameters.template get<ContactPatterns<FP>>().get_school_holidays() = local_holidays;
        node_model.parameters.template get<DailyFirstVaccination<FP>>()                 = local_daily_v1;
        node_model.parameters.template get<DailyFullVaccination<FP>>()                  = local_daily_v2;

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
} // namespace osecirvvs
} // namespace mio

#endif // ODESECIRVVS_PARAMETER_SPACE_H
