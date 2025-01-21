/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Khoa Nguyen
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
#ifndef ABM_ANALYZE_RESULT_H
#define ABM_ANALYZE_RESULT_H

#include "abm/parameters.h"

#include <vector>

namespace mio
{
namespace abm
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

    auto num_runs  = ensemble_params.size();
    auto num_nodes = ensemble_params[0].size();
    std::vector<double> single_element_ensemble(num_runs);
    auto num_groups = (size_t)ensemble_params[0][0].parameters.get_num_groups();

    // Lambda function that calculates the percentile of a single parameter
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
        for (auto age_group = AgeGroup(0); age_group < AgeGroup(num_groups); age_group++) {
            for (auto virus_variant : enum_members<VirusVariant>()) {
                // Global infection parameters
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<IncubationPeriod>()[{virus_variant, age_group}].params.m();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<IncubationPeriod>()[{virus_variant, age_group}].params.s();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedNoSymptomsToSymptoms>()[{virus_variant, age_group}]
                            .params.m();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedNoSymptomsToSymptoms>()[{virus_variant, age_group}]
                            .params.s();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedNoSymptomsToRecovered>()[{virus_variant, age_group}]
                            .params.m();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedNoSymptomsToRecovered>()[{virus_variant, age_group}]
                            .params.s();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedSymptomsToSevere>()[{virus_variant, age_group}]
                            .params.m();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedSymptomsToSevere>()[{virus_variant, age_group}]
                            .params.s();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedSymptomsToRecovered>()[{virus_variant, age_group}]
                            .params.m();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedSymptomsToRecovered>()[{virus_variant, age_group}]
                            .params.s();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedSevereToCritical>()[{virus_variant, age_group}]
                            .params.m();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedSevereToCritical>()[{virus_variant, age_group}]
                            .params.s();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedSevereToRecovered>()[{virus_variant, age_group}]
                            .params.m();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedSevereToRecovered>()[{virus_variant, age_group}]
                            .params.s();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedCriticalToDead>()[{virus_variant, age_group}]
                            .params.m();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedCriticalToDead>()[{virus_variant, age_group}]
                            .params.s();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedCriticalToRecovered>()[{virus_variant, age_group}]
                            .params.m();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<TimeInfectedCriticalToRecovered>()[{virus_variant, age_group}]
                            .params.s();
                    return result;
                });

                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    return model.parameters.template get<SymptomsPerInfectedNoSymptoms>()[{virus_variant, age_group}];
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    return model.parameters.template get<SeverePerInfectedSymptoms>()[{virus_variant, age_group}];
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    return model.parameters.template get<CriticalPerInfectedSevere>()[{virus_variant, age_group}];
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    return model.parameters.template get<DeathsPerInfectedCritical>()[{virus_variant, age_group}];
                });

                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    return model.parameters.template get<DetectInfection>()[{virus_variant, age_group}];
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<ViralLoadDistributions>()[{virus_variant, age_group}]
                            .viral_load_incline.params.a();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<ViralLoadDistributions>()[{virus_variant, age_group}]
                            .viral_load_incline.params.b();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<ViralLoadDistributions>()[{virus_variant, age_group}]
                            .viral_load_decline.params.a();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<ViralLoadDistributions>()[{virus_variant, age_group}]
                            .viral_load_decline.params.b();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<ViralLoadDistributions>()[{virus_variant, age_group}]
                            .viral_load_peak.params.a();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<ViralLoadDistributions>()[{virus_variant, age_group}]
                            .viral_load_peak.params.b();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<InfectivityDistributions>()[{virus_variant, age_group}]
                            .infectivity_alpha.params.a();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<InfectivityDistributions>()[{virus_variant, age_group}]
                            .infectivity_alpha.params.b();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<InfectivityDistributions>()[{virus_variant, age_group}]
                            .infectivity_beta.params.a();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<InfectivityDistributions>()[{virus_variant, age_group}]
                            .infectivity_beta.params.b();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<VirusShedFactor>()[{virus_variant, age_group}].params.a();
                    return result;
                });
                param_percentil(node, [age_group, virus_variant](auto&& model) -> auto& {
                    static auto result =
                        model.parameters.template get<VirusShedFactor>()[{virus_variant, age_group}].params.b();
                    return result;
                });
                param_percentil(node, [virus_variant](auto&& model) -> auto& {
                    return model.parameters.template get<AerosolTransmissionRates>()[{virus_variant}];
                });
            }
            param_percentil(node, [age_group](auto&& model) -> auto& {
                return model.parameters.template get<BasicShoppingRate>()[{age_group}];
            });
            param_percentil(node, [age_group](auto&& model) -> auto& {
                static auto result = model.parameters.template get<GotoWorkTimeMinimum>()[{age_group}].hours();
                return result;
            });
            param_percentil(node, [age_group](auto&& model) -> auto& {
                static auto result = model.parameters.template get<GotoWorkTimeMaximum>()[{age_group}].hours();
                return result;
            });
            param_percentil(node, [age_group](auto&& model) -> auto& {
                static auto result = model.parameters.template get<GotoSchoolTimeMinimum>()[{age_group}].hours();
                return result;
            });
            param_percentil(node, [age_group](auto&& model) -> auto& {
                static auto result = model.parameters.template get<GotoSchoolTimeMaximum>()[{age_group}].hours();
                return result;
            });
        }
        param_percentil(node, [](auto&& model) -> auto& {
            return model.parameters.template get<MaskProtection>()[MaskType::Community];
        });
        param_percentil(node, [](auto&& model) -> auto& {
            return model.parameters.template get<MaskProtection>()[MaskType::FFP2];
        });
        param_percentil(node, [](auto&& model) -> auto& {
            return model.parameters.template get<MaskProtection>()[MaskType::Surgical];
        });
        param_percentil(node, [](auto&& model) -> auto& {
            static auto result = model.parameters.template get<LockdownDate>().days();
            return result;
        });
    }

    return percentile;
}

} // namespace abm
} // namespace mio

#endif
