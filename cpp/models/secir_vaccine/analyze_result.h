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

#include "memilio/utils/time_series.h"
#include "secir_vaccine/secir.h"
#include "memilio/mobility/mobility.h"

#include <functional>
#include <vector>

namespace mio
{
namespace vaccinated
{

    /**
 * interpolate time series with evenly spaced, integer time points.
 * time points [t0, t1, t2, ..., tmax] interpolated as [floor(t0), floor(t0) + 1,...,ceil(tmax)].
 * values at new time points are linearly interpolated from their immediate neighbors from the old time points.
 * @param simulation_result time series to interpolate
 * @return interpolated time series
 */
    TimeSeries<double> interpolate_simulation_result(const TimeSeries<double>& simulation_result);

    /**
 * helper template, type returned by overload interpolate_simulation_result(T t)
 */
    template <class T>
    using InterpolateResultT = std::decay_t<decltype(interpolate_simulation_result(std::declval<T>()))>;

    /**
 * @brief Interpolates results of all runs with evenly spaced, integer time points.
 * @see interpolate_simulation_result
 * @param ensemble_result result of multiple simulations (single TimeSeries or Graph)
 * @return interpolated time series, one (or as many as nodes in the graph) per result in the ensemble
 */
    template <class T>
    std::vector<InterpolateResultT<T>> interpolate_ensemble_results(const std::vector<T>& ensemble_results)
    {
        std::vector<InterpolateResultT<T>> interpolated;
        interpolated.reserve(ensemble_results.size());
        std::transform(ensemble_results.begin(), ensemble_results.end(), std::back_inserter(interpolated),
                       [](auto& run) {
                           return interpolate_simulation_result(run);
                       });
        return interpolated;
    }

    std::vector<std::vector<TimeSeries<double>>>
    sum_nodes(const std::vector<std::vector<TimeSeries<double>>>& ensemble_result);

    /**
 * @brief computes mean of each compartment, node, and time point over all runs
 * input must be uniform as returned by interpolated_ensemble_result:
 * same number of nodes, same time points and elements.
 * @see interpolated_ensemble_result
 * @param ensemble_results uniform results of multiple simulation runs
 * @return mean of the results over all runs
 */
    std::vector<TimeSeries<double>> ensemble_mean(const std::vector<std::vector<TimeSeries<double>>>& ensemble_results);

    /**
 * @brief computes the p percentile of the result for each compartment, node, and time point.
 * Produces for each compartment the value that that is bigger than approximately
 * a p-th share of the values of this compartment over all runs.
 * input must be uniform as returned by interpolated_ensemble_result:
 * same number of nodes, same time points and elements.
 * @see interpolated_ensemble_result
 * @param ensemble_result uniform results of multiple simulation runs
 * @param p percentile value in open interval (0, 1)
 * @return p percentile of the results over all runs
 */
    std::vector<TimeSeries<double>>
    ensemble_percentile(const std::vector<std::vector<TimeSeries<double>>>& ensemble_result, double p);
    /**
 * interpolate time series with evenly spaced, integer time points for each node.
 * @see interpolate_simulation_result
 * @param graph_result graph of simulations whose results will be interpolated
 * @return one interpolated time series per node
 */
    template <class Simulation>
    std::vector<TimeSeries<double>>
    interpolate_simulation_result(const Graph<SimulationNode<Simulation>, MigrationEdge>& graph_result)
    {
        std::vector<TimeSeries<double>> interpolated;
        interpolated.reserve(graph_result.nodes().size());
        std::transform(graph_result.nodes().begin(), graph_result.nodes().end(), std::back_inserter(interpolated),
                       [](auto& n) {
                           return interpolate_simulation_result(n.property.get_result());
                       });
        return interpolated;
    }

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
                    node, [i](auto&& model) -> auto& {
                        return model.parameters.template get<InfectiousTimeMild>()[i];
                    });
                param_percentil(
                    node, [i](auto&& model) -> auto& {
                        return model.parameters.template get<HospitalizedToICUTime>()[i];
                    });
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
                        return model.parameters.template get<RiskOfInfectionFromSympomatic>()[i];
                    });
                param_percentil(
                    node, [i](auto&& model) -> auto& {
                        return model.parameters.template get<MaxRiskOfInfectionFromSympomatic>()[i];
                    });
                param_percentil(
                    node, [i](auto&& model) -> auto& {
                        return model.parameters.template get<AsymptoticCasesPerInfectious>()[i];
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
                    node, [i](auto&& model) -> auto& {
                        return model.parameters.template get<mio::DeathsPerHospitalized>()[i];
                    });
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

    /**
 * Compute the distance between two SECIR simulation results.
 * The distance is the 2-norm of the element-wise difference of the two results.
 * The two results (e.g. output of interpolate_simulation_result) must have the same dimensions and number of time points.
 * @param result1 first result.
 * @param result2 second result.
 * @return Computed distance between result1 and result2.
 */
    double result_distance_2norm(const std::vector<mio::TimeSeries<double>>& result1,
                                 const std::vector<mio::TimeSeries<double>>& result2);

    /**
 * Compute the distance between two SECIR simulation results in one compartment.
 * The distance is the 2-norm of the element-wise difference of the two results in the specified compartment.
 * The two results (e.g. output of interpolate_simulation_result) must have the same dimensions and number of time points.
 * @param result1 first result.
 * @param result2 second result.
 * @param compartment the compartment to compare.
 * @return Computed distance between result1 and result2.
 */
    double result_distance_2norm(const std::vector<mio::TimeSeries<double>>& result1,
                                 const std::vector<mio::TimeSeries<double>>& result2, InfectionState compartment);

} // namespace vaccinated
} // namespace mio

#endif //EPI_SECIR_ANALYZE_RESULT_H
