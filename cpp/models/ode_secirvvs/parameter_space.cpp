/* 
* Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
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
#include "ode_secirvvs/parameter_space.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/model.h"

namespace mio
{
namespace osecirvvs
{

    void draw_sample_demographics(Model& model)
    {
        model.parameters.get<ICUCapacity>().draw_sample();
        model.parameters.get<TestAndTraceCapacity>().draw_sample();

        for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
            double group_total = model.populations.get_group_total(i);

            //sample initial compartments (with exceptions)
            for (auto inf_state = Index<InfectionState>(0); inf_state < InfectionState::Count; ++inf_state) {
                if (inf_state != InfectionState::SusceptibleNaive && //not sampled, fixed after sampling everything else
                    inf_state != InfectionState::Dead && //not sampled, fixed from data
                    inf_state != InfectionState::TotalInfections) { //not sampled, only for record keeping
                    model.populations[{i, inf_state}].draw_sample();
                }
            }

            //set susceptibles so the total number stays the same as before sampling.
            //if the new total without susceptibles is already bigger than the previous total
            //subtract the overflow from recovered, susceptibles will then be approximately zero.
            model.populations[{i, InfectionState::SusceptibleNaive}] = 0;
            double diff = model.populations.get_group_total(i) - group_total;
            if (diff > 0) {
                model.populations[{i, InfectionState::Recovered}] -= diff;
                if (model.populations[{i, InfectionState::Recovered}] < 0.0) {
                    log_error("Negative Compartment after sampling.");
                }
                assert(std::abs(group_total - model.populations.get_group_total(i)) < 1e-10 && "Sanity check.");
            }
            model.populations.set_difference_from_group_total<AgeGroup>({i, InfectionState::SusceptibleNaive},
                                                                        group_total);
        }
    }

    void draw_sample_infection(Model& model)
    {
        model.parameters.get<Seasonality>().draw_sample();

        //not age dependent
        model.parameters.get<IncubationTime>()[AgeGroup(0)].draw_sample();
        model.parameters.get<SerialInterval>()[AgeGroup(0)].draw_sample();
        model.parameters.get<TimeInfectedSymptoms>()[AgeGroup(0)].draw_sample();
        model.parameters.get<HospitalizedToICUTime>()[AgeGroup(0)].draw_sample();
        model.parameters.get<RelativeCarrierInfectability>()[AgeGroup(0)].draw_sample();
        model.parameters.get<RiskOfInfectionFromSympomatic>()[AgeGroup(0)].draw_sample();
        model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[AgeGroup(0)].draw_sample();

        model.parameters.get<ExposedFactorPartialImmunity>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ExposedFactorImprovedImmunity>()[AgeGroup(0)].draw_sample();
        model.parameters.get<InfectedFactorPartialImmunity>()[AgeGroup(0)].draw_sample();
        model.parameters.get<InfectedFactorImprovedImmunity>()[AgeGroup(0)].draw_sample();
        model.parameters.get<HospitalizedFactorPartialImmunity>()[AgeGroup(0)].draw_sample();
        model.parameters.get<HospitalizedFactorImprovedImmunity>()[AgeGroup(0)].draw_sample();
        model.parameters.get<InfectiousTimeFactorImmune>()[AgeGroup(0)].draw_sample();

        for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
            //not age dependent
            model.parameters.get<IncubationTime>()[i]     = model.parameters.get<IncubationTime>()[AgeGroup(0)];
            model.parameters.get<SerialInterval>()[i]     = model.parameters.get<SerialInterval>()[AgeGroup(0)];
            model.parameters.get<TimeInfectedSymptoms>()[i] = model.parameters.get<TimeInfectedSymptoms>()[AgeGroup(0)];
            model.parameters.get<HospitalizedToICUTime>()[i] =
                model.parameters.get<HospitalizedToICUTime>()[AgeGroup(0)];
            model.parameters.get<RelativeCarrierInfectability>()[i] =
                model.parameters.get<RelativeCarrierInfectability>()[AgeGroup(0)];
            model.parameters.get<RiskOfInfectionFromSympomatic>()[i] =
                model.parameters.get<RiskOfInfectionFromSympomatic>()[AgeGroup(0)];
            model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[i] =
                model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[AgeGroup(0)];

            model.parameters.get<ExposedFactorPartialImmunity>()[i]     = model.parameters.get<ExposedFactorPartialImmunity>()[AgeGroup(0)];
            model.parameters.get<ExposedFactorImprovedImmunity>()[i]   = model.parameters.get<ExposedFactorImprovedImmunity>()[AgeGroup(0)];
            model.parameters.get<InfectedFactorPartialImmunity>()[i]     = model.parameters.get<InfectedFactorPartialImmunity>()[AgeGroup(0)];
            model.parameters.get<InfectedFactorImprovedImmunity>()[i]   = model.parameters.get<InfectedFactorImprovedImmunity>()[AgeGroup(0)];
            model.parameters.get<HospitalizedFactorPartialImmunity>()[i]    = model.parameters.get<HospitalizedFactorPartialImmunity>()[AgeGroup(0)];
            model.parameters.get<HospitalizedFactorImprovedImmunity>()[i]  = model.parameters.get<HospitalizedFactorImprovedImmunity>()[AgeGroup(0)];
            model.parameters.get<InfectiousTimeFactorImmune>()[i] = model.parameters.get<InfectiousTimeFactorImmune>()[AgeGroup(0)];

            //age dependent
            model.parameters.get<HospitalizedToHomeTime>()[i].draw_sample(); // here: home=recovered
            model.parameters.get<ICUToDeathTime>()[i].draw_sample();
            model.parameters.get<ICUToHomeTime>()[i].draw_sample();

            model.parameters.get<InfectionProbabilityFromContact>()[i].draw_sample();
            model.parameters.get<AsymptoticCasesPerInfectious>()[i].draw_sample();
            model.parameters.get<DeathsPerICU>()[i].draw_sample();
            model.parameters.get<HospitalizedCasesPerInfectious>()[i].draw_sample();
            model.parameters.get<ICUCasesPerHospitalized>()[i].draw_sample();
        }
    }

    void draw_sample(Model& model)
    {
        draw_sample_infection(model);
        draw_sample_demographics(model);
        model.parameters.get<ContactPatterns>().draw_sample();
        model.apply_constraints();
    }

    Graph<Model, MigrationParameters> draw_sample(Graph<Model, MigrationParameters>& graph, bool variant_high)
    {
        Graph<Model, MigrationParameters> sampled_graph;

        //sample global parameters
        auto& shared_params_model = graph.nodes()[0].property;
        draw_sample_infection(shared_params_model);
        auto& shared_contacts = shared_params_model.parameters.template get<ContactPatterns>();
        shared_contacts.draw_sample_dampings();
        auto& shared_dynamic_npis = shared_params_model.parameters.template get<DynamicNPIsInfected>();
        shared_dynamic_npis.draw_sample();

        double delta_fac;
        if (variant_high) {
            delta_fac = 1.6;
        }
        else {
            delta_fac = 1.4;
        }

        //infectiousness of virus variants is not sampled independently but depend on base infectiousness
        for (auto i = AgeGroup(0); i < shared_params_model.parameters.get_num_groups(); ++i) {
            shared_params_model.parameters.template get<BaseInfectiousnessB117>()[i] =
                shared_params_model.parameters.template get<InfectionProbabilityFromContact>()[i];
            shared_params_model.parameters.template get<BaseInfectiousnessB161>()[i] =
                shared_params_model.parameters.template get<InfectionProbabilityFromContact>()[i] * delta_fac;
        }

        for (auto& params_node : graph.nodes()) {
            auto& node_model = params_node.property;

            //sample local parameters
            draw_sample_demographics(params_node.property);

            //copy global parameters
            //save demographic parameters so they aren't overwritten
            auto local_icu_capacity = node_model.parameters.template get<ICUCapacity>();
            auto local_tnt_capacity = node_model.parameters.template get<TestAndTraceCapacity>();
            auto local_holidays     = node_model.parameters.template get<ContactPatterns>().get_school_holidays();
            auto local_daily_v1     = node_model.parameters.template get<DailyFirstVaccination>();
            auto local_daily_v2     = node_model.parameters.template get<DailyFullVaccination>();
            node_model.parameters   = shared_params_model.parameters;
            node_model.parameters.template get<ICUCapacity>()                           = local_icu_capacity;
            node_model.parameters.template get<TestAndTraceCapacity>()                  = local_tnt_capacity;
            node_model.parameters.template get<ContactPatterns>().get_school_holidays() = local_holidays;
            node_model.parameters.template get<DailyFirstVaccination>()                 = local_daily_v1;
            node_model.parameters.template get<DailyFullVaccination>()                  = local_daily_v2;

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
