/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Wadim Koslow, Daniel Abele, David Kerkmann, Sascha Korf
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
#ifndef MEMILIO_DATA_ANALYZE_RESULT_H
#define MEMILIO_DATA_ANALYZE_RESULT_H

#include "memilio/utils/time_series.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/math/interpolation.h"

#include <functional>
#include <vector>

namespace mio
{

/**
 * @brief interpolate time series with evenly spaced, integer time points that represent whole days.
 * time points [t0, t1, t2, ..., tmax] interpolated as days [ceil(t0), floor(t0) + 1,...,floor(tmax)].
 * tolerances in the first and last time point (t0 and t_max) are accounted for.
 * values at new time points are linearly interpolated from their immediate neighbors from the old time points.
 * @see interpolate_simulation_result
 * @param simulation_result time series to interpolate
 * @param abs_tol  absolute tolerance given for doubles t0 and tmax to account for small deviations from whole days.
 * @return interpolated time series
 */
template <typename FP = double>
TimeSeries<FP> interpolate_simulation_result(const TimeSeries<FP>& simulation_result,
                                                 const FP abs_tol = 1e-14)
{
    using std::floor;
    using std::ceil;
    const auto t0    = simulation_result.get_time(0);
    const auto t_max = simulation_result.get_last_time();
    // add another day if the first time point is equal to day_0 up to absolute tolerance tol
    const auto day0 = (t0 - abs_tol < ceil(t0) - 1) ? floor(t0) : ceil(t0);
    // add another day if the last time point is equal to day_max up to absolute tolerance tol
    const auto day_max = (t_max + abs_tol > floor(t_max) + 1) ? ceil(t_max) : floor(t_max);

    // create interpolation_times vector with all days between day0 and day_max
    std::vector<FP> tps(static_cast<int>(day_max) - static_cast<int>(day0) + 1);
    std::iota(tps.begin(), tps.end(), day0);

    return interpolate_simulation_result<FP>(simulation_result, tps);
}

/**
 * @brief interpolate time series with freely chosen time points that lie in between the time points of the given time series up to a given tolerance.
 * values at new time points are linearly interpolated from their immediate neighbors from the old time points.
 * @param simulation_result time series to interpolate
 * @param interpolations_times std::vector of time points at which simulation results are interpolated.
 * @return interpolated time series at given interpolation points
 */
template <typename FP = double>
TimeSeries<FP> interpolate_simulation_result(const TimeSeries<FP>& simulation_result,
                                                 const std::vector<FP>& interpolation_times)
{
    assert(simulation_result.get_num_time_points() > 0 && "TimeSeries must not be empty.");

    assert(std::is_sorted(interpolation_times.begin(), interpolation_times.end()) &&
           "Time points for interpolation have to be sorted in non-descending order.");

    if (interpolation_times.size() >= 2) {
        assert((interpolation_times[1] > simulation_result.get_time(0) &&
                interpolation_times.rbegin()[1] <= simulation_result.get_last_time()) &&
               "All but the first and the last time point of interpolation have lie between simulation times (strictly "
               "for lower boundary).");
    }

    TimeSeries<FP> interpolated(simulation_result.get_num_elements());

    if (interpolation_times.size() == 0) {
        return interpolated;
    }

    size_t interp_idx = 0;
    // add first time point of interpolation times in case it is smaller than the first time point of simulation_result
    // this is used for the case that it equals the first time point of simulation up to tolerance
    // this is necessary even if the tolerance is 0 due to the way the comparison in the loop is implemented (< and >=)
    if (simulation_result.get_time(0) >= interpolation_times[0]) {
        interpolated.add_time_point(interpolation_times[0], simulation_result[0]);
        ++interp_idx;
    }

    //interpolate between pair of time points that lie on either side of each interpolation point
    for (Eigen::Index sim_idx = 0;
         sim_idx < simulation_result.get_num_time_points() - 1 && interp_idx < interpolation_times.size();) {
        //only go to next pair of time points if no time point is added.
        //otherwise check the same time points again
        //in case there is more than one interpolation point between the two time points
        if (simulation_result.get_time(sim_idx) < interpolation_times[interp_idx] &&
            simulation_result.get_time(sim_idx + 1) >= interpolation_times[interp_idx]) {
            
            std::cout << "debug 1\n";
            interpolated.add_time_point(
                interpolation_times[interp_idx],
                linear_interpolation(interpolation_times[interp_idx], simulation_result.get_time(sim_idx),
                                     simulation_result.get_time(sim_idx + 1), simulation_result[sim_idx],
                                     simulation_result[sim_idx + 1]));
            std::cout << "debug 1\n";
            ++interp_idx;
        }
        else {
            ++sim_idx;
        }
    }

    // add last time point of interpolation times in case it is larger than the last time point of simulation_result
    // this is used for the case that it equals the last time point of simulation up to tolerance
    if (interp_idx < interpolation_times.size() &&
        simulation_result.get_last_time() < interpolation_times[interp_idx]) {
        interpolated.add_time_point(interpolation_times[interp_idx], simulation_result.get_last_value());
    }

    return interpolated;
}

/**
 * helper template, type returned by overload interpolate_simulation_result(T t)
 */
template <class T>
using InterpolateResultT = std::decay_t<decltype(interpolate_simulation_result(std::declval<T>()))>;

/**
 * @brief Interpolates results of all runs with evenly spaced, integer time points that represent whole days.
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
template <class Simulation, typename FP = double>
std::vector<TimeSeries<FP>>
interpolate_simulation_result(const Graph<SimulationNode<Simulation>, MobilityEdge<FP>>& graph_result)
{
    std::vector<TimeSeries<FP>> interpolated;
    interpolated.reserve(graph_result.nodes().size());
    std::transform(graph_result.nodes().begin(), graph_result.nodes().end(), std::back_inserter(interpolated),
                   [](auto& n) {
                       return interpolate_simulation_result(n.property.get_result());
                   });
    return interpolated;
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
 * Compute the distance between two compartment model simulation results in one compartment.
 * The distance is the 2-norm of the element-wise difference of the two results in the specified compartment.
 * The two results (e.g. output of interpolate_simulation_result) must have the same dimensions and number of time points.
 * @tparam InfectionState enum type that defines the compartments of the model that produced the results.
 * @param result1 first result.
 * @param result2 second result.
 * @param compartment the compartment to compare.
 * @return Computed distance between result1 and result2.
 */
template <class InfectionState>
double result_distance_2norm(const std::vector<mio::TimeSeries<double>>& result1,
                             const std::vector<mio::TimeSeries<double>>& result2, InfectionState compartment)
{
    assert(result1.size() == result2.size());
    assert(result1.size() > 0);
    assert(result1[0].get_num_time_points() > 0);
    assert(result1[0].get_num_elements() > 0);

    auto num_compartments = Eigen::Index(InfectionState::Count);
    auto num_age_groups   = result1[0].get_num_elements() / num_compartments;

    auto norm_sqr = 0.0;
    for (auto iter_node1 = result1.begin(), iter_node2 = result2.begin(); iter_node1 < result1.end();
         ++iter_node1, ++iter_node2) {
        for (Eigen::Index time_idx = 0; time_idx < iter_node1->get_num_time_points(); ++time_idx) {
            auto v1 = (*iter_node1)[time_idx];
            auto v2 = (*iter_node2)[time_idx];
            for (Eigen::Index age_idx = 0; age_idx < num_age_groups; ++age_idx) {
                auto d1 = v1[age_idx * num_compartments + Eigen::Index(compartment)];
                auto d2 = v2[age_idx * num_compartments + Eigen::Index(compartment)];
                norm_sqr += (d1 - d2) * (d1 - d2);
            }
        }
    }
    return std::sqrt(norm_sqr);
}

} // namespace mio

#endif //MEMILIO_DATA_ANALYZE_RESULT_H
