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
#include "memilio/mobility/mobility.h"
#include "models/abm/time.h"

#include <functional>
#include <vector>

namespace mio
{

/**
 * interpolate time series with evenly spaced, integer time points.
 * time points [t0, t1, t2, ..., tmax] interpolated as [floor(t0), floor(t0) + 1,...,ceil(tmax)].
 * values at new time points are linearly interpolated from their immediate neighbors from the old time points.
 * @param simulation_result time series to interpolate
 * @param interpolations_times vector of time points that represent the evaluation of interpolation
 * @return interpolated time series
 */
TimeSeries<double> interpolate_simulation_result(const TimeSeries<double>& simulation_result, const std::vector<TimePoint> interpolation_times);

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
    std::transform(ensemble_results.begin(), ensemble_results.end(), std::back_inserter(interpolated), [](auto& run) {
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
std::vector<TimeSeries<double>> ensemble_percentile(const std::vector<std::vector<TimeSeries<double>>>& ensemble_result,
                                                    double p);
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

} // namespace mio

#endif //EPI_SECIR_ANALYZE_RESULT_H
