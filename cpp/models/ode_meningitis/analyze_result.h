/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Martin J. Kuehn
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
#ifndef ODEMENG_ANALYZE_RESULT_H
#define ODEMENG_ANALYZE_RESULT_H

#include "ode_meningitis/infection_state.h"
#include "ode_meningitis/parameters.h"

namespace mio
{

namespace omeng
{

/**
 * @brief computes the p percentile of the parameters for each node.
 * @param ensemble_params vector of multiple simulation runs, each a vector of Models per node
 * @param p percentile value in open interval (0, 1)
 * @return p percentile of the parameters over all runs
 */
template <typename FP, class Model>
std::vector<Model> ensemble_params_percentile(const std::vector<std::vector<Model>>& ensemble_params, FP p)
{
    assert(p > 0.0 && p < 1.0 && "Invalid percentile value.");

    auto num_runs   = ensemble_params.size();
    auto num_nodes  = ensemble_params[0].size();
    auto num_groups = (size_t)ensemble_params[0][0].parameters.get_num_groups();

    std::vector<Model> percentile(num_nodes, Model((int)num_groups));
    auto param_percentil = [&ensemble_params, p, num_runs, &percentile](auto n, auto get_param) mutable {
        std::vector<FP> single_element(num_runs);
        for (size_t run = 0; run < num_runs; run++) {
            auto const& params  = ensemble_params[run][n];
            single_element[run] = get_param(params);
        }
        std::sort(single_element.begin(), single_element.end());
        auto& new_params = get_param(percentile[n]);
        new_params       = single_element[static_cast<size_t>(num_runs * p)];
    };

    for (size_t node = 0; node < num_nodes; node++) {
        for (auto i = AgeGroup(0); i < AgeGroup(num_groups); i++) {
            // Population
            for (size_t compart = 0; compart < (size_t)InfectionState::Count; ++compart) {
                param_percentil(
                    node, [ compart, i ](auto&& model) -> auto& {
                        return model.populations[{i, (InfectionState)compart}];
                    });
            }
            // per group infection parameters
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<RateCarrierToInfected<FP>>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<RateCarrierToRecovered<FP>>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<RateInfectedToDead<FP>>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<RateInfectedToRecovered<FP>>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.template get<RateNaturalDeath<FP>>()[i]; });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.template get<RateImmunityLoss<FP>>()[i]; });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<ProbabilityImmunityLossSusLow<FP>>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<TransmissionProbabilityOnContact<FP>>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<RiskOfInfectionFromFromCarrier<FP>>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<RiskOfInfectionFromFromInfected<FP>>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.template get<ModificationRate<FP>>()[i]; });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<IncomeFractionSusLow<FP>>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.template get<IncomeRate<FP>>()[i]; });
        }
        // group independent params
        param_percentil(
            node,
            [](auto&& model) -> auto& { return model.parameters.template get<DynamicNPIsImplementationDelay<FP>>(); });
    }
    return percentile;
}

} // namespace omeng
} // namespace mio

#endif // ODEMENG_ANALYZE_RESULT_H
