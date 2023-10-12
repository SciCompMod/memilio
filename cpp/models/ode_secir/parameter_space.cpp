/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "ode_secir/parameter_space.h"
#include "memilio/utils/parameter_distributions.h"
#include "ode_secir/model.h"

namespace mio
{
namespace osecir
{

void set_params_distributions_normal(Model& model, double t0, double tmax, double dev_rel)
{
    auto set_distribution = [dev_rel](UncertainValue& v, double min_val = 0.001) {
        v.set_distribution(ParameterDistributionNormal(
            //add add limitsfornonsense big values. Also mscv has a problem with a few doublesso this fixes it
            std::min(std::max(min_val, (1 - dev_rel * 2.6) * v), 0.1 * std::numeric_limits<double>::max()),
            std::min(std::max(min_val, (1 + dev_rel * 2.6) * v), 0.5 * std::numeric_limits<double>::max()),
            std::min(std::max(min_val, double(v)), 0.3 * std::numeric_limits<double>::max()),
            std::min(std::max(min_val, dev_rel * v), std::numeric_limits<double>::max())));
    };

    set_distribution(model.parameters.get<Seasonality>(), 0.0);
    set_distribution(model.parameters.get<ICUCapacity>());
    set_distribution(model.parameters.get<TestAndTraceCapacity>());

    // populations
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        for (auto j = Index<InfectionState>(0); j < Index<InfectionState>(InfectionState::Count); j++) {

            // don't touch S and D
            if (j == InfectionState::Susceptible || j == InfectionState::Dead) {
                continue;
            }

            // variably sized groups
            set_distribution(model.populations[{i, j}], 0.0);
        }
    }

    // times
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {

        set_distribution(model.parameters.get<IncubationTime>()[i]);
        set_distribution(model.parameters.get<SerialInterval>()[i]);
        set_distribution(model.parameters.get<TimeInfectedSymptoms>()[i]);
        set_distribution(model.parameters.get<TimeInfectedSevere>()[i]);
        set_distribution(model.parameters.get<TimeInfectedCritical>()[i]);

        set_distribution(model.parameters.get<TransmissionProbabilityOnContact>()[i]);
        set_distribution(model.parameters.get<RelativeTransmissionNoSymptoms>()[i]);
        set_distribution(model.parameters.get<RecoveredPerInfectedNoSymptoms>()[i]);
        set_distribution(model.parameters.get<RiskOfInfectionFromSymptomatic>()[i]);
        set_distribution(model.parameters.get<MaxRiskOfInfectionFromSymptomatic>()[i]);
        set_distribution(model.parameters.get<DeathsPerCritical>()[i]);
        set_distribution(model.parameters.get<SeverePerInfectedSymptoms>()[i]);
        set_distribution(model.parameters.get<CriticalPerSevere>()[i]);
    }

    // dampings
    auto matrices = std::vector<size_t>();
    for (size_t i = 0; i < model.parameters.get<ContactPatterns>().get_cont_freq_mat().get_num_matrices(); ++i) {
        matrices.push_back(i);
    }
    auto groups = Eigen::VectorXd::Constant(Eigen::Index(model.parameters.get_num_groups().get()), 1.0);
    model.parameters.get<ContactPatterns>().get_dampings().emplace_back(
        mio::UncertainValue(0.5), mio::DampingLevel(0), mio::DampingType(0), mio::SimulationTime(t0 + (tmax - t0) / 2),
        matrices, groups);
    set_distribution(model.parameters.get<ContactPatterns>().get_dampings()[0].get_value(), 0.0);
}

void draw_sample_demographics(Model& model)
{
    model.parameters.get<ICUCapacity>().draw_sample();
    model.parameters.get<TestAndTraceCapacity>().draw_sample();

    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        double group_total = model.populations.get_group_total(i);

        model.populations[{i, InfectionState::Exposed}].draw_sample();
        model.populations[{i, InfectionState::InfectedNoSymptoms}].draw_sample();
        model.populations[{i, InfectionState::InfectedSymptoms}].draw_sample();
        model.populations[{i, InfectionState::InfectedSevere}].draw_sample();
        model.populations[{i, InfectionState::InfectedCritical}].draw_sample();
        model.populations[{i, InfectionState::Recovered}].draw_sample();

        // no sampling for dead and total numbers
        // [...]

        model.populations.set_difference_from_group_total<AgeGroup>({i, InfectionState::Susceptible}, group_total);
        model.populations.set_difference_from_group_total<AgeGroup>({i, InfectionState::Susceptible},
                                                                    model.populations.get_group_total(i));
    }
}

void draw_sample_infection(Model& model)
{
    model.parameters.get<Seasonality>().draw_sample();

    //not age dependent
    model.parameters.get<IncubationTime>()[AgeGroup(0)].draw_sample();
    model.parameters.get<SerialInterval>()[AgeGroup(0)].draw_sample();
    model.parameters.get<TimeInfectedSymptoms>()[AgeGroup(0)].draw_sample();
    model.parameters.get<RelativeTransmissionNoSymptoms>()[AgeGroup(0)].draw_sample();
    model.parameters.get<RiskOfInfectionFromSymptomatic>()[AgeGroup(0)].draw_sample();
    model.parameters.get<MaxRiskOfInfectionFromSymptomatic>()[AgeGroup(0)].draw_sample();

    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        //not age dependent
        model.parameters.get<IncubationTime>()[i] = model.parameters.get<IncubationTime>()[AgeGroup(0)];
        model.parameters.get<SerialInterval>()[i] = model.parameters.get<SerialInterval>()[AgeGroup(0)];
        model.parameters.get<RelativeTransmissionNoSymptoms>()[i] =
            model.parameters.get<RelativeTransmissionNoSymptoms>()[AgeGroup(0)];
        model.parameters.get<RiskOfInfectionFromSymptomatic>()[i] =
            model.parameters.get<RiskOfInfectionFromSymptomatic>()[AgeGroup(0)];
        model.parameters.get<MaxRiskOfInfectionFromSymptomatic>()[i] =
            model.parameters.get<MaxRiskOfInfectionFromSymptomatic>()[AgeGroup(0)];

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

} // namespace osecir
} // namespace mio
