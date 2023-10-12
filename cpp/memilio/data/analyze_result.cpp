/* 
* Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
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
#include "memilio/data/analyze_result.h"
#include "memilio/math/interpolation.h"

#include <algorithm>
#include <cassert>

namespace mio
{
TimeSeries<double> interpolate_simulation_result(const TimeSeries<double>& simulation_result, const double abs_tol)
{
    const auto t0    = simulation_result.get_time(0);
    const auto t_max = simulation_result.get_last_time();
    // add another day if the first time point is equal to day_0 up to absolute tolerance tol
    const auto day0 = (t0 - abs_tol < std::ceil(t0) - 1) ? std::floor(t0) : std::ceil(t0);
    // add another day if the last time point is equal to day_max up to absolute tolerance tol
    const auto day_max = (t_max + abs_tol > std::floor(t_max) + 1) ? std::ceil(t_max) : std::floor(t_max);

    // create interpolation_times vector with all days between day0 and day_max
    std::vector<double> tps(static_cast<int>(day_max) - static_cast<int>(day0) + 1);
    std::iota(tps.begin(), tps.end(), day0);

    return interpolate_simulation_result(simulation_result, tps);
}

TimeSeries<double> interpolate_simulation_result(const TimeSeries<double>& simulation_result,
                                                 const std::vector<double>& interpolation_times)
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

    TimeSeries<double> interpolated(simulation_result.get_num_elements());

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
                linear_interpolation(interpolation_times[interp_idx], simulation_result.get_time(sim_idx),
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

std::vector<std::vector<TimeSeries<double>>>
sum_nodes(const std::vector<std::vector<TimeSeries<double>>>& ensemble_result)
{
    auto num_runs        = ensemble_result.size();
    auto num_nodes       = ensemble_result[0].size();
    auto num_time_points = ensemble_result[0][0].get_num_time_points();
    auto num_elements    = ensemble_result[0][0].get_num_elements();

    std::vector<std::vector<TimeSeries<double>>> sum_result(
        num_runs, std::vector<TimeSeries<double>>(1, TimeSeries<double>::zero(num_time_points, num_elements)));

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

std::vector<TimeSeries<double>> ensemble_mean(const std::vector<std::vector<TimeSeries<double>>>& ensemble_result)
{
    auto num_runs        = ensemble_result.size();
    auto num_nodes       = ensemble_result[0].size();
    auto num_time_points = ensemble_result[0][0].get_num_time_points();
    auto num_elements    = ensemble_result[0][0].get_num_elements();

    std::vector<TimeSeries<double>> mean(num_nodes, TimeSeries<double>::zero(num_time_points, num_elements));

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

std::vector<TimeSeries<double>> ensemble_percentile(const std::vector<std::vector<TimeSeries<double>>>& ensemble_result,
                                                    double p)
{
    assert(p > 0.0 && p < 1.0 && "Invalid percentile value.");

    auto num_runs        = ensemble_result.size();
    auto num_nodes       = ensemble_result[0].size();
    auto num_time_points = ensemble_result[0][0].get_num_time_points();
    auto num_elements    = ensemble_result[0][0].get_num_elements();

    std::vector<TimeSeries<double>> percentile(num_nodes, TimeSeries<double>::zero(num_time_points, num_elements));

    std::vector<double> single_element_ensemble(num_runs); //reused for each element
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

double result_distance_2norm(const std::vector<mio::TimeSeries<double>>& result1,
                             const std::vector<mio::TimeSeries<double>>& result2)
{
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
    return std::sqrt(norm_sqr);
}

} // namespace mio
