/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#ifndef EPI_SECIR_ANALYZE_RESULT_H
#define EPI_SECIR_ANALYZE_RESULT_H

#include "secir/secir.h"
#include "memilio/data/analyze_result.h"

#include <functional>
#include <vector>

namespace mio
{

/**
 * @brief computes the p percentile of the parameters for each node.
 * @param ensemble_result graph of multiple simulation runs
 * @param p percentile value in open interval (0, 1)
 * @return p percentile of the parameters over all runs
 */
template <class Model>
std::vector<Model> ensemble_params_percentile(const std::vector<std::vector<Model>>& ensemble_params, double p)
{
    assert(p > 0.0 && p < 1.0 && "Invalid percentile value.");

    auto num_runs   = ensemble_params.size();
    auto num_nodes  = ensemble_params[0].size();
    auto num_groups = (size_t)ensemble_params[0][0].parameters.get_num_groups();

    std::vector<double> single_element_ensemble(num_runs);

    // lamda function that calculates the percentile of a single paramter
    std::vector<Model> percentile(num_nodes, Model((int)num_groups));
    auto param_percentil = [&ensemble_params, p, num_runs, &percentile](auto n, auto get_param) mutable {
        std::vector<double> single_element(num_runs);
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
            //Population
            for (size_t compart = 0; compart < (size_t)InfectionState::Count; ++compart) {
                param_percentil(
                    node, [ compart, i ](auto&& model) -> auto& {
                        return model.populations[{i, (InfectionState)compart}];
                    });
            }
            // times
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.template get<IncubationTime>()[i]; });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.template get<SerialInterval>()[i]; });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.template get<TimeInfectedSymptoms>()[i]; });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.template get<HospitalizedToICUTime>()[i]; });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<HospitalizedToHomeTime>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<HomeToHospitalizedTime>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.template get<ICUToDeathTime>()[i]; });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.template get<ICUToHomeTime>()[i]; });
            //probs
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<RelativeCarrierInfectability>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<RiskOfInfectionFromSymptomatic>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<MaxRiskOfInfectionFromSymptomatic>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<AsymptomaticCasesPerInfectious>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<HospitalizedCasesPerInfectious>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.template get<ICUCasesPerHospitalized>()[i];
                });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.template get<mio::DeathsPerICU>()[i]; });
        }
        // group independent params
        param_percentil(
            node, [](auto&& model) -> auto& { return model.parameters.template get<mio::Seasonality>(); });
        param_percentil(
            node, [](auto&& model) -> auto& { return model.parameters.template get<mio::TestAndTraceCapacity>(); });

        for (size_t run = 0; run < num_runs; run++) {

            auto const& params = ensemble_params[run][node];
            single_element_ensemble[run] =
                params.parameters.template get<mio::ICUCapacity>() * params.populations.get_total();
        }
        std::sort(single_element_ensemble.begin(), single_element_ensemble.end());
        percentile[node].parameters.template set<mio::ICUCapacity>(
            single_element_ensemble[static_cast<size_t>(num_runs * p)]);
    }
    return percentile;
}

} // namespace mio

#endif //EPI_SECIR_ANALYZE_RESULT_H
