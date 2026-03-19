/* 
* Copyright (C) 2020-2026 MEmilio
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
#ifndef ODEMENG_PARAMETER_SPACE_H
#define ODEMENG_PARAMETER_SPACE_H

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/memory.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"
#include "ode_meningitis/model.h"

#include <assert.h>
#include <limits>
#include <string>
#include <vector>
#include <random>
#include <memory>

namespace mio
{
namespace omeng
{

/**
 * @brief Sets all meningitis model parameters normally distributed,
 *  using the current value and a given standard deviation.
 * @param[inout] model Model including contact patterns
 * @param[in] t0 start time
 * @param[in] tmax end time
 * @param[in] dev_rel maximum relative deviation from particular value(s) given in params
 */
template <typename FP>
void set_params_distributions_normal(Model<FP>& model, FP t0, FP tmax, FP dev_rel)
{
    using std::max;
    using std::min;
    auto set_distribution = [dev_rel](UncertainValue<FP>& v, FP min_val = 0.001) {
        auto lower_bound = min<FP>(max<FP>(min_val, (1 - dev_rel * 2.6) * v), 0.1 * std::numeric_limits<FP>::max());
        auto upper_bound = min<FP>(max<FP>(min_val, (1 + dev_rel * 2.6) * v), 0.5 * std::numeric_limits<FP>::max());

        if (mio::floating_point_equal<FP>(lower_bound, upper_bound, mio::Limits<FP>::zero_tolerance())) {
            mio::log_debug("Bounded ParameterDistribution has standard deviation close to zero. Therefore constant "
                           "distribution is used.");
            v.set_distribution(ParameterDistributionConstant(
                min<FP>(max<FP>(min_val, v.value()), 0.3 * std::numeric_limits<FP>::max())));
        }
        else {
            v.set_distribution(ParameterDistributionNormal(
                lower_bound, upper_bound, min<FP>(max<FP>(min_val, v.value()), 0.3 * std::numeric_limits<FP>::max()),
                min<FP>(max<FP>(min_val, dev_rel * v), std::numeric_limits<FP>::max())));
        }
    };

    set_distribution(model.parameters.template get<DynamicNPIsImplementationDelay<FP>>());

    // populations – skip Incoming, Dead and DeadNatural
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        for (auto j = Index<InfectionState>(0); j < Index<InfectionState>(InfectionState::Count); j++) {
            if (j == InfectionState::Incoming || j == InfectionState::Dead ||
                j == InfectionState::DeadNatural) {
                continue;
            }
            set_distribution(model.populations[{i, j}], 0.0);
        }
    }

    // per-group infection parameters
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        set_distribution(model.parameters.template get<RateCarrierToInfected<FP>>()[i]);
        set_distribution(model.parameters.template get<RateCarrierToRecovered<FP>>()[i]);
        set_distribution(model.parameters.template get<RateInfectedToDead<FP>>()[i]);
        set_distribution(model.parameters.template get<RateInfectedToRecovered<FP>>()[i]);
        set_distribution(model.parameters.template get<RateNaturalDeath<FP>>()[i]);
        set_distribution(model.parameters.template get<RateImmunityLoss<FP>>()[i]);
        set_distribution(model.parameters.template get<ProbabilityImmunityLossSusLow<FP>>()[i], 0.0);
        set_distribution(model.parameters.template get<TransmissionProbabilityOnContact<FP>>()[i]);
        set_distribution(model.parameters.template get<RiskOfInfectionFromFromCarrier<FP>>()[i]);
        set_distribution(model.parameters.template get<RiskOfInfectionFromFromInfected<FP>>()[i]);
        set_distribution(model.parameters.template get<ModificationRate<FP>>()[i], 0.0);
        set_distribution(model.parameters.template get<IncomeFractionSusLow<FP>>()[i], 0.0);
        set_distribution(model.parameters.template get<IncomeRate<FP>>()[i]);
    }

    // dampings
    auto matrices = std::vector<size_t>();
    for (size_t i = 0; i < model.parameters.template get<ContactPatterns<FP>>().get_cont_freq_mat().get_num_matrices();
         ++i) {
        matrices.push_back(i);
    }
    auto groups = Eigen::VectorX<FP>::Constant(Eigen::Index(model.parameters.get_num_groups().get()), 1.0);
    model.parameters.template get<ContactPatterns<FP>>().get_dampings().emplace_back(
        mio::UncertainValue<FP>(0.5), mio::DampingLevel(0), mio::DampingType(0),
        mio::SimulationTime<FP>(t0 + (tmax - t0) / 2), matrices, groups);
    set_distribution(model.parameters.template get<ContactPatterns<FP>>().get_dampings()[0].get_value(), 0.0);
}

/**
 * @brief Draws a sample from distributions for all demographic parameters (population).
 * @param[inout] model Meningitis model
 */
template <typename FP>
void draw_sample_demographics(Model<FP>& model)
{
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        FP group_total = model.populations.get_group_total(i);

        model.populations[{i, InfectionState::SusceptibleLow}].draw_sample();
        model.populations[{i, InfectionState::Carrier}].draw_sample();
        model.populations[{i, InfectionState::Infected}].draw_sample();
        model.populations[{i, InfectionState::Recovered}].draw_sample();

        // SusceptibleHigh is set as remainder – Dead/DeadNatural/Incoming are not sampled
        model.populations.template set_difference_from_group_total<AgeGroup>({i, InfectionState::SusceptibleHigh},
                                                                              group_total);
        model.populations.template set_difference_from_group_total<AgeGroup>({i, InfectionState::SusceptibleHigh},
                                                                              model.populations.get_group_total(i));
    }
}

/**
 * @brief Draws a sample from distributions for all infection-related parameters.
 * @param[inout] model Meningitis model
 */
template <typename FP>
void draw_sample_infection(Model<FP>& model)
{
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); i++) {
        model.parameters.template get<RateCarrierToInfected<FP>>()[i].draw_sample();
        model.parameters.template get<RateCarrierToRecovered<FP>>()[i].draw_sample();
        model.parameters.template get<RateInfectedToDead<FP>>()[i].draw_sample();
        model.parameters.template get<RateInfectedToRecovered<FP>>()[i].draw_sample();
        model.parameters.template get<RateNaturalDeath<FP>>()[i].draw_sample();
        model.parameters.template get<RateImmunityLoss<FP>>()[i].draw_sample();
        model.parameters.template get<ProbabilityImmunityLossSusLow<FP>>()[i].draw_sample();
        model.parameters.template get<TransmissionProbabilityOnContact<FP>>()[i].draw_sample();
        model.parameters.template get<RiskOfInfectionFromFromCarrier<FP>>()[i].draw_sample();
        model.parameters.template get<RiskOfInfectionFromFromInfected<FP>>()[i].draw_sample();
        model.parameters.template get<ModificationRate<FP>>()[i].draw_sample();
        model.parameters.template get<IncomeFractionSusLow<FP>>()[i].draw_sample();
        model.parameters.template get<IncomeRate<FP>>()[i].draw_sample();
    }
}

/**
 * @brief Draws a sample from all Model parameter distributions.
 * @param[inout] model Meningitis model
 */
template <typename FP>
void draw_sample(Model<FP>& model)
{
    draw_sample_infection(model);
    draw_sample_demographics(model);
    model.parameters.template get<ContactPatterns<FP>>().draw_sample();
    model.apply_constraints();
}

/**
 * @brief Draws samples for a graph of meningitis models.
 * @param[inout] graph Graph of meningitis models with mobility parameters
 * @return New graph with sampled parameter values
 */
template <typename FP>
Graph<Model<FP>, MobilityParameters<FP>> draw_sample(Graph<Model<FP>, MobilityParameters<FP>>& graph)
{
    Graph<Model<FP>, MobilityParameters<FP>> sampled_graph;

    // sample global parameters from first node
    auto& shared_params_model = graph.nodes()[0].property;
    draw_sample_infection(shared_params_model);
    auto& shared_contacts = shared_params_model.parameters.template get<ContactPatterns<FP>>();
    shared_contacts.draw_sample_dampings();
    auto& shared_dynamic_npis = shared_params_model.parameters.template get<DynamicNPIsInfectedSymptoms<FP>>();
    shared_dynamic_npis.draw_sample();
    auto& shared_dynamic_npis_delay = shared_params_model.parameters.template get<DynamicNPIsImplementationDelay<FP>>();
    shared_dynamic_npis_delay.draw_sample();

    for (auto& params_node : graph.nodes()) {
        auto& node_model = params_node.property;

        // sample local demographics
        draw_sample_demographics(params_node.property);

        // copy global parameters but keep node-local income rate and contact holidays
        auto local_income_rate = node_model.parameters.template get<IncomeRate<FP>>();
        auto local_holidays    = node_model.parameters.template get<ContactPatterns<FP>>().get_school_holidays();
        node_model.parameters  = shared_params_model.parameters;
        node_model.parameters.template get<IncomeRate<FP>>()                             = local_income_rate;
        node_model.parameters.template get<ContactPatterns<FP>>().get_school_holidays()  = local_holidays;

        node_model.parameters.template get<ContactPatterns<FP>>().make_matrix();
        node_model.apply_constraints();

        sampled_graph.add_node(params_node.id, node_model);
    }

    for (auto& edge : graph.edges()) {
        auto edge_params = edge.property;
        apply_dampings(edge_params.get_coefficients(), shared_contacts.get_dampings(), [&edge_params](auto& v) {
            return make_mobility_damping_vector(edge_params.get_coefficients().get_shape(), v);
        });
        edge_params.set_dynamic_npis_infected(shared_dynamic_npis);
        sampled_graph.add_edge(edge.start_node_idx, edge.end_node_idx, edge_params);
    }

    return sampled_graph;
}

} // namespace omeng
} // namespace mio

#endif // ODEMENG_PARAMETER_SPACE_H
