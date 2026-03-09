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
#ifndef MIO_ODE_SECIRTS_ANALYZE_RESULT_H
#define MIO_ODE_SECIRTS_ANALYZE_RESULT_H

#include "ode_secirts/model.h"
#include "memilio/data/analyze_result.h"
#include "ode_secirts/parameters.h"

#include <functional>
#include <vector>

namespace mio
{
namespace osecirts
{
/**
    * @brief computes the p percentile of the parameters for each node.
    * @param ensemble_result graph of multiple simulation runs
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
    auto num_days   = ensemble_params[0][0]
                        .parameters.template get<DailyPartialVaccinations<FP>>()
                        .template size<mio::SimulationDay>();

    std::vector<FP> single_element_ensemble(num_runs);

    // lambda function that calculates the percentile of a single parameter
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
        percentile[node].parameters.template get<DailyPartialVaccinations<FP>>().resize(num_days);
        percentile[node].parameters.template get<DailyFullVaccinations<FP>>().resize(num_days);
        percentile[node].parameters.template get<DailyBoosterVaccinations<FP>>().resize(num_days);

        for (auto i = AgeGroup(0); i < AgeGroup(num_groups); i++) {
            //Population
            for (auto compart = Index<InfectionState>(0); compart < InfectionState::Count; ++compart) {
                param_percentil(node, [compart, i](auto&& model) -> auto& {
                    return model.populations[{i, compart}];
                });
            }
            // times
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<TimeExposed<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<TimeInfectedNoSymptoms<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<TimeInfectedSymptoms<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<TimeInfectedSevere<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<TimeInfectedCritical<FP>>()[i];
            });
            //probs
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<TransmissionProbabilityOnContact<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<RelativeTransmissionNoSymptoms<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<RiskOfInfectionFromSymptomatic<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<SeverePerInfectedSymptoms<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<CriticalPerSevere<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<DeathsPerCritical<FP>>()[i];
            });
            //vaccinations
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<ReducExposedPartialImmunity<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<ReducExposedImprovedImmunity<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<ReducTimeInfectedMild<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<DaysUntilEffectivePartialVaccination<FP>>()[i];
            });
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<DaysUntilEffectiveImprovedVaccination<FP>>()[i];
            });

            for (auto day = SimulationDay(0); day < num_days; ++day) {
                param_percentil(node, [i, day](auto&& model) -> auto& {
                    return model.parameters.template get<DailyPartialVaccinations<FP>>()[{i, day}];
                });
                param_percentil(node, [i, day](auto&& model) -> auto& {
                    return model.parameters.template get<DailyFullVaccinations<FP>>()[{i, day}];
                });
                param_percentil(node, [i, day](auto&& model) -> auto& {
                    return model.parameters.template get<DailyBoosterVaccinations<FP>>()[{i, day}];
                });
            }
            //virus variants
            param_percentil(node, [i](auto&& model) -> auto& {
                return model.parameters.template get<InfectiousnessNewVariant<FP>>()[i];
            });
        }
        // group independent params
        param_percentil(node, [](auto&& model) -> auto& {
            return model.parameters.template get<Seasonality<FP>>();
        });
        param_percentil(node, [](auto&& model) -> auto& {
            return model.parameters.template get<TestAndTraceCapacity<FP>>();
        });
        param_percentil(node, [](auto&& model) -> auto& {
            return model.parameters.template get<ICUCapacity<FP>>();
        });
        param_percentil(node, [](auto&& model) -> auto& {
            return model.parameters.template get<DynamicNPIsImplementationDelay<FP>>();
        });

        for (size_t run = 0; run < num_runs; run++) {

            auto const& params = ensemble_params[run][node];
            single_element_ensemble[run] =
                params.parameters.template get<ICUCapacity<FP>>() * params.populations.get_total();
        }
        std::sort(single_element_ensemble.begin(), single_element_ensemble.end());
        percentile[node].parameters.template set<ICUCapacity<FP>>(
            single_element_ensemble[static_cast<size_t>(num_runs * p)]);
    }
    return percentile;
}

} // namespace osecirts
} // namespace mio

#endif //MIO_ODE_SECIRTS_ANALYZE_RESULT_H
