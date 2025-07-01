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
#include "memilio/io/io.h"

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
template <typename FP>
TimeSeries<FP> interpolate_simulation_result(const TimeSeries<FP>& simulation_result, const FP abs_tol = 1e-14);

/**
 * @brief interpolate time series with freely chosen time points that lie in between the time points of the given time series up to a given tolerance.
 * values at new time points are linearly interpolated from their immediate neighbors from the old time points.
 * @param simulation_result time series to interpolate
 * @param interpolations_times std::vector of time points at which simulation results are interpolated.
 * @return interpolated time series at given interpolation points
 */
template <typename FP>
TimeSeries<FP> interpolate_simulation_result(const TimeSeries<FP>& simulation_result,

                                             const std::vector<FP>& interpolation_times);

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

template <typename FP>
std::vector<std::vector<TimeSeries<FP>>> sum_nodes(const std::vector<std::vector<TimeSeries<FP>>>& ensemble_result);

/**
 * @brief computes mean of each compartment, node, and time point over all runs
 * input must be uniform as returned by interpolated_ensemble_result:
 * same number of nodes, same time points and elements.
 * @see interpolated_ensemble_result
 * @param ensemble_results uniform results of multiple simulation runs
 * @return mean of the results over all runs
 */
template <typename FP>
std::vector<TimeSeries<FP>> ensemble_mean(const std::vector<std::vector<TimeSeries<FP>>>& ensemble_results);

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
template <typename FP>
std::vector<TimeSeries<FP>> ensemble_percentile(const std::vector<std::vector<TimeSeries<FP>>>& ensemble_result, FP p);
/**
 * interpolate time series with evenly spaced, integer time points for each node.
 * @see interpolate_simulation_result
 * @param graph_result graph of simulations whose results will be interpolated
 * @return one interpolated time series per node
 */
template <typename FP, class Simulation>
std::vector<TimeSeries<FP>>
interpolate_simulation_result(const Graph<SimulationNode<FP, Simulation>, MobilityEdge<FP>>& graph_result)
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
template <typename FP>
FP result_distance_2norm(const std::vector<mio::TimeSeries<FP>>& result1,
                         const std::vector<mio::TimeSeries<FP>>& result2);

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
template <typename FP, class InfectionState>
FP result_distance_2norm(const std::vector<mio::TimeSeries<FP>>& result1,
                         const std::vector<mio::TimeSeries<FP>>& result2, InfectionState compartment)
{
    using std::sqrt;

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
    return sqrt(norm_sqr);
}

template <typename FP>
TimeSeries<FP> interpolate_simulation_result(const TimeSeries<FP>& simulation_result, const FP abs_tol)
{
    const auto t0    = simulation_result.get_time(0);
    const auto t_max = simulation_result.get_last_time();
    // add another day if the first time point is equal to day_0 up to absolute tolerance tol
    const auto day0 = (t0 - abs_tol < std::ceil(t0) - 1) ? std::floor(t0) : std::ceil(t0);
    // add another day if the last time point is equal to day_max up to absolute tolerance tol
    const auto day_max = (t_max + abs_tol > std::floor(t_max) + 1) ? std::ceil(t_max) : std::floor(t_max);

    // create interpolation_times vector with all days between day0 and day_max
    std::vector<FP> tps(static_cast<int>(day_max) - static_cast<int>(day0) + 1);
    std::iota(tps.begin(), tps.end(), day0);

    return interpolate_simulation_result<FP>(simulation_result, tps);
}

template <typename FP>
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
            interpolated.add_time_point(
                interpolation_times[interp_idx],
                linear_interpolation<FP>(interpolation_times[interp_idx], simulation_result.get_time(sim_idx),
                                         simulation_result.get_time(sim_idx + 1), simulation_result[sim_idx],
                                         simulation_result[sim_idx + 1]));
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

template <typename FP>
std::vector<std::vector<TimeSeries<FP>>> sum_nodes(const std::vector<std::vector<TimeSeries<FP>>>& ensemble_result)
{
    auto num_runs        = ensemble_result.size();
    auto num_nodes       = ensemble_result[0].size();
    auto num_time_points = ensemble_result[0][0].get_num_time_points();
    auto num_elements    = ensemble_result[0][0].get_num_elements();

    std::vector<std::vector<TimeSeries<FP>>> sum_result(
        num_runs, std::vector<TimeSeries<FP>>(1, TimeSeries<FP>::zero(num_time_points, num_elements)));

    for (size_t run = 0; run < num_runs; run++) {
        for (Eigen::Index time = 0; time < num_time_points; time++) {
            sum_result[run][0].get_time(time) = ensemble_result[run][0].get_time(time);
            for (size_t node = 0; node < num_nodes; node++) {
                sum_result[run][0][time] += ensemble_result[run][node][time];
            }
        }
    }
    return sum_result;
}

template <typename FP>
std::vector<TimeSeries<FP>> ensemble_mean(const std::vector<std::vector<TimeSeries<FP>>>& ensemble_result)
{
    auto num_runs        = ensemble_result.size();
    auto num_nodes       = ensemble_result[0].size();
    auto num_time_points = ensemble_result[0][0].get_num_time_points();
    auto num_elements    = ensemble_result[0][0].get_num_elements();

    std::vector<TimeSeries<FP>> mean(num_nodes, TimeSeries<FP>::zero(num_time_points, num_elements));

    for (size_t run = 0; run < num_runs; run++) {
        assert(ensemble_result[run].size() == num_nodes && "ensemble results not uniform.");
        for (size_t node = 0; node < num_nodes; node++) {
            assert(ensemble_result[run][node].get_num_time_points() == num_time_points &&
                   "ensemble results not uniform.");
            for (Eigen::Index time = 0; time < num_time_points; time++) {
                assert(ensemble_result[run][node].get_num_elements() == num_elements &&
                       "ensemble results not uniform.");
                mean[node].get_time(time) = ensemble_result[run][node].get_time(time);
                mean[node][time] += ensemble_result[run][node][time] / num_runs;
            }
        }
    }

    return mean;
}

template <typename FP>
std::vector<TimeSeries<FP>> ensemble_percentile(const std::vector<std::vector<TimeSeries<FP>>>& ensemble_result, FP p)
{
    assert(p > 0.0 && p < 1.0 && "Invalid percentile value.");

    auto num_runs        = ensemble_result.size();
    auto num_nodes       = ensemble_result[0].size();
    auto num_time_points = ensemble_result[0][0].get_num_time_points();
    auto num_elements    = ensemble_result[0][0].get_num_elements();

    std::vector<TimeSeries<FP>> percentile(num_nodes, TimeSeries<FP>::zero(num_time_points, num_elements));

    std::vector<FP> single_element_ensemble(num_runs); //reused for each element
    for (size_t node = 0; node < num_nodes; node++) {
        for (Eigen::Index time = 0; time < num_time_points; time++) {
            percentile[node].get_time(time) = ensemble_result[0][node].get_time(time);
            for (Eigen::Index elem = 0; elem < num_elements; elem++) {
                std::transform(ensemble_result.begin(), ensemble_result.end(), single_element_ensemble.begin(),
                               [=](auto& run) {
                                   return run[node][time][elem];
                               });
                std::sort(single_element_ensemble.begin(), single_element_ensemble.end());
                percentile[node][time][elem] = single_element_ensemble[static_cast<size_t>(num_runs * p)];
            }
        }
    }
    return percentile;
}

template <typename FP>
FP result_distance_2norm(const std::vector<mio::TimeSeries<FP>>& result1,
                         const std::vector<mio::TimeSeries<FP>>& result2)
{
    using std::sqrt;
    assert(result1.size() == result2.size());
    assert(result1.size() > 0);
    assert(result1[0].get_num_time_points() > 0);
    assert(result1[0].get_num_elements() > 0);

    auto norm_sqr = 0.0;
    for (auto iter_node1 = result1.begin(), iter_node2 = result2.begin(); iter_node1 < result1.end();
         ++iter_node1, ++iter_node2) {
        for (Eigen::Index time_idx = 0; time_idx < iter_node1->get_num_time_points(); ++time_idx) {
            auto v1 = (*iter_node1)[time_idx];
            auto v2 = (*iter_node2)[time_idx];
            norm_sqr += ((v1 - v2).array() * (v1 - v2).array()).sum();
        }
    }
    return sqrt(norm_sqr);
}

/**
 * @brief This function merges two TimeSeries by copying their time points and values to a new TimeSeries in the correct order.
 * If both TimeSeries have values for the same time point, their values are either added or only one value is taken.
 * @param[in] ts1 First TimeSeries.
 * @param[in] ts2 Second TimeSeries.
 * @param[in] add_values Boolean specifying whether the values should be added if both TimeSeries contain the same time point. If false, the value of just the first TimeSeries is taken.
 * @tparam FP A floating point type.
 * @return A TimeSeries containing all time points and values from both input TimeSeries.
 */
template <class FP>
IOResult<TimeSeries<FP>> merge_time_series(const TimeSeries<FP>& ts1, const TimeSeries<FP>& ts2,
                                           bool add_values = false)
{
    TimeSeries<FP> merged_ts(ts1.get_num_elements());
    if (ts1.get_num_elements() != ts2.get_num_elements()) {
        log_error("TimeSeries have a different number of elements.");
        return failure(mio::StatusCode::InvalidValue);
    }
    else {
        Eigen::Index t1_iterator = 0;
        Eigen::Index t2_iterator = 0;
        bool t1_finished         = false;
        bool t2_finished         = false;
        while (!t1_finished || !t2_finished) {
            if (!t1_finished) {
                if (ts1.get_time(t1_iterator) < ts2.get_time(t2_iterator) ||
                    t2_finished) { // Current time point of first TimeSeries is smaller than current time point of second TimeSeries or second TimeSeries has already been copied entirely
                    merged_ts.add_time_point(ts1.get_time(t1_iterator), ts1.get_value(t1_iterator));
                    t1_iterator += 1;
                }
                else if (!t2_finished && ts1.get_time(t1_iterator) ==
                                             ts2.get_time(t2_iterator)) { // Both TimeSeries have the current time point
                    if (add_values) {
                        merged_ts.add_time_point(ts1.get_time(t1_iterator),
                                                 ts1.get_value(t1_iterator) + ts2.get_value(t2_iterator));
                    }
                    else {
                        merged_ts.add_time_point(ts1.get_time(t1_iterator), ts1.get_value(t1_iterator));
                        log_warning("Both TimeSeries have values for t={}. The value of the first TimeSeries is used",
                                    ts1.get_time(t1_iterator));
                    }
                    t1_iterator += 1;
                    t2_iterator += 1;
                    if (t2_iterator >=
                        ts2.get_num_time_points()) { // Check if all values of second TimeSeries have been copied
                        t2_finished = true;
                        t2_iterator = ts2.get_num_time_points() - 1;
                    }
                }
                if (t1_iterator >=
                    ts1.get_num_time_points()) { // Check if all values of first TimeSeries have been copied
                    t1_finished = true;
                    t1_iterator = ts1.get_num_time_points() - 1;
                }
            }
            if (!t2_finished) {
                if (ts2.get_time(t2_iterator) < ts1.get_time(t1_iterator) ||
                    t1_finished) { // Current time point of second TimeSeries is smaller than current time point of first TimeSeries or first TimeSeries has already been copied entirely
                    merged_ts.add_time_point(ts2.get_time(t2_iterator), ts2.get_value(t2_iterator));
                    t2_iterator += 1;
                }
                else if (!t1_finished && ts2.get_time(t2_iterator) ==
                                             ts1.get_time(t1_iterator)) { // Both TimeSeries have the current time point
                    if (add_values) {
                        merged_ts.add_time_point(ts1.get_time(t1_iterator),
                                                 ts1.get_value(t1_iterator) + ts2.get_value(t2_iterator));
                    }
                    else {
                        merged_ts.add_time_point(ts1.get_time(t1_iterator), ts1.get_value(t1_iterator));
                        log_warning("Both TimeSeries have values for t={}. The value of the first TimeSeries is used",
                                    ts1.get_time(t1_iterator));
                    }
                    t1_iterator += 1;
                    t2_iterator += 1;
                    if (t1_iterator >=
                        ts1.get_num_time_points()) { // Check if all values of first TimeSeries have been copied
                        t1_finished = true;
                        t1_iterator = ts1.get_num_time_points() - 1;
                    }
                }
                if (t2_iterator >=
                    ts2.get_num_time_points()) { // Check if all values of second TimeSeries have been copied
                    t2_finished = true;
                    t2_iterator = ts2.get_num_time_points() - 1;
                }
            }
        }
    }
    return success(merged_ts);
}

} // namespace mio

#endif //MEMILIO_DATA_ANALYZE_RESULT_H
