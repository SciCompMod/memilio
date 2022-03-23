/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "secir_vaccine/parameter_space.h"
#include "memilio/utils/parameter_distributions.h"
#include "secir_vaccine/secir.h"

namespace mio
{
namespace vaccinated
{

    void draw_sample_demographics(SecirModel& model)
    {
        model.parameters.get<ICUCapacity>().draw_sample();
        model.parameters.get<TestAndTraceCapacity>().draw_sample();

        for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
            double group_total = model.populations.get_group_total(i);

            model.populations[{i, InfectionState::Exposed}].draw_sample();
            model.populations[{i, InfectionState::Carrier}].draw_sample();
            model.populations[{i, InfectionState::Infected}].draw_sample();
            model.populations[{i, InfectionState::Hospitalized}].draw_sample();
            model.populations[{i, InfectionState::ICU}].draw_sample();
            model.populations[{i, InfectionState::Recovered}].draw_sample();

            // no sampling for dead and total numbers
            // [...]

            model.populations[{i, InfectionStateV::Susceptible}] = 0;
            double group_total_dummy                             = model.populations.get_group_total(i);
            if (group_total_dummy < group_total) {
                model.populations.set_difference_from_group_total<AgeGroup>({i, epi::InfectionStateV::Susceptible},
                                                                            group_total);
            }
            else {
                double diff = group_total_dummy - group_total;
                model.populations[{i, InfectionStateV::Recovered}] =
                    model.populations[{i, InfectionStateV::Recovered}] - diff;
                assert(std::abs(group_total - model.populations.get_group_total(i)) < 1e-10);
            }

            model.populations.set_difference_from_group_total<AgeGroup>({i, InfectionStateV::Susceptible},
                                                                        model.populations.get_group_total(i));
        }
    }

    void draw_sample_infection(SecirModel& model)
    {
        model.parameters.get<Seasonality>().draw_sample();

        //not age dependent
        model.parameters.get<IncubationTime>()[AgeGroup(0)].draw_sample();
        model.parameters.get<SerialInterval>()[AgeGroup(0)].draw_sample();
        model.parameters.get<InfectiousTimeMild>()[AgeGroup(0)].draw_sample();
        model.parameters.get<HospitalizedToICUTime>()[AgeGroup(0)].draw_sample();
        model.parameters.get<RelativeCarrierInfectability>()[AgeGroup(0)].draw_sample();
        model.parameters.get<RiskOfInfectionFromSympomatic>()[AgeGroup(0)].draw_sample();
        model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[AgeGroup(0)].draw_sample();

        model.parameters.get<ReducVaccExp>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducImmuneExp>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducExpInf>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducImmuneExpInf>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducInfHosp>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducImmuneInfHosp>()[AgeGroup(0)].draw_sample();
        model.parameters.get<ReducTime>()[AgeGroup(0)].draw_sample();

        for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
            //not age dependent
            model.parameters.get<IncubationTime>()[i]     = model.parameters.get<IncubationTime>()[AgeGroup(0)];
            model.parameters.get<SerialInterval>()[i]     = model.parameters.get<SerialInterval>()[AgeGroup(0)];
            model.parameters.get<InfectiousTimeMild>()[i] = model.parameters.get<InfectiousTimeMild>()[AgeGroup(0)];
            model.parameters.get<HospitalizedToICUTime>()[i] =
                model.parameters.get<HospitalizedToICUTime>()[AgeGroup(0)];
            model.parameters.get<RelativeCarrierInfectability>()[i] =
                model.parameters.get<RelativeCarrierInfectability>()[AgeGroup(0)];
            model.parameters.get<RiskOfInfectionFromSympomatic>()[i] =
                model.parameters.get<RiskOfInfectionFromSympomatic>()[AgeGroup(0)];
            model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[i] =
                model.parameters.get<MaxRiskOfInfectionFromSympomatic>()[AgeGroup(0)];

            model.parameters.get<ReducVaccExp>()[i]       = model.parameters.get<ReducVaccExp>()[AgeGroup(0)];
            model.parameters.get<ReducImmuneExp>()[i]     = model.parameters.get<ReducImmuneExp>()[AgeGroup(0)];
            model.parameters.get<ReducExpInf>()[i]        = model.parameters.get<ReducExpInf>()[AgeGroup(0)];
            model.parameters.get<ReducImmuneExpInf>()[i]  = model.parameters.get<ReducImmuneExpInf>()[AgeGroup(0)];
            model.parameters.get<ReducInfHosp>()[i]       = model.parameters.get<ReducInfHosp>()[AgeGroup(0)];
            model.parameters.get<ReducImmuneInfHosp>()[i] = model.parameters.get<ReducImmuneInfHosp>()[AgeGroup(0)];
            model.parameters.get<ReducTime>()[i]          = model.parameters.get<ReducTime>()[AgeGroup(0)];

            //age dependent
            model.parameters.get<HospitalizedToHomeTime>()[i].draw_sample(); // here: home=recovered
            model.parameters.get<HomeToHospitalizedTime>()[i].draw_sample(); // here: home=infectious
            model.parameters.get<InfectiousTimeAsymptomatic>()[i].draw_sample();
            model.parameters.get<ICUToDeathTime>()[i].draw_sample();
            model.parameters.get<ICUToHomeTime>()[i].draw_sample();

            model.parameters.get<InfectionProbabilityFromContact>()[i].draw_sample();
            model.parameters.get<AsymptoticCasesPerInfectious>()[i].draw_sample();
            model.parameters.get<DeathsPerHospitalized>()[i].draw_sample();
            model.parameters.get<HospitalizedCasesPerInfectious>()[i].draw_sample();
            model.parameters.get<ICUCasesPerHospitalized>()[i].draw_sample();
        }
    }

    void draw_sample(SecirModel& model)
    {
        draw_sample_infection(model);
        draw_sample_demographics(model);
        model.parameters.get<ContactPatterns>().draw_sample();
        model.apply_constraints();
    }

    Graph<SecirModel, MigrationParameters> draw_sample(Graph<SecirModel, MigrationParameters>& graph)
    {
        Graph<SecirModel, MigrationParameters> sampled_graph;

        //sample global parameters
        auto& shared_params_model = graph.nodes()[0].property;
        draw_sample_infection(shared_params_model);
        auto& shared_contacts = shared_params_model.parameters.template get<ContactPatterns>();
        shared_contacts.draw_sample_dampings();
        auto& shared_dynamic_npis = shared_params_model.parameters.template get<DynamicNPIsInfected>();
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
            node_model.parameters   = shared_params_model.parameters;
            node_model.parameters.template get<ICUCapacity>()                           = local_icu_capacity;
            node_model.parameters.template get<TestAndTraceCapacity>()                  = local_tnt_capacity;
            node_model.parameters.template get<ContactPatterns>().get_school_holidays() = local_holidays;

            node_model.parameters.template get<ContactPatterns>().make_matrix();
            node_model.apply_constraints();

            sampled_graph.add_node(params_node.id, node_model);
        }

        for (auto& edge : graph.edges()) {
            auto edge_params = edge.property;
            apply_dampings(edge_params.get_coefficients(), shared_contacts.get_dampings(), [&edge_params](auto& v) {
                return make_migration_damping_vector(edge_params.get_coefficients().get_shape(), v);
            });
            edge_params.set_dynamic_npis_infected(shared_dynamic_npis);
            sampled_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge_params);
        }

        return sampled_graph;
    }

} // namespace vaccinated
} // namespace mio
