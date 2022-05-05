/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, David Kerkmann, Sascha Korf
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
#include "memilio/utils/analyze_result.h"
#include "memilio/math/interpolation.h"
#include <algorithm>
#include <cassert>

namespace mio
{

TimeSeries<double> interpolate_simulation_result_days(const TimeSeries<double>& simulation_result, const double abs_tol)
{
    const auto t0 = simulation_result.get_time(0);
    const auto t_max = simulation_result.get_last_time();
    // add another day if the first time point is equal to day_0 up to absolute tolerance tol
    const auto day0 = (t0 - abs_tol < std::ceil(t0) - 1) ? std::floor(t0) : std::ceil(t0);
    // add another day if the last time point is equal to day_max up to absolute tolerance tol
    const auto day_max = (t_max + abs_tol > std::floor(t_max) + 1) ? std::ceil(t_max) : std::floor(t_max);
    
    std::vector<double> tps(day_max - day0 + 1);
    std::iota(tps.begin(), tps.end(), day0);
    
    return interpolate_simulation_result(simulation_result, tps, abs_tol);
}


TimeSeries<double> interpolate_simulation_result(const TimeSeries<double>& simulation_result, const std::vector<double>& interpolation_times, const double abs_tol)
{
    assert(simulation_result.get_num_time_points() > 0 && "TimeSeries must not be empty.");
    
    assert((std::min(interpolation_times) + abs_tol >= simulation_result.get_time(0) && std::max(interpolation_times) - abs_tol <= simulation_result.get_last_time()) && "Interpolation times have lie between simulation times (up to given tolerance).");
    unused(abs_tol);
    
    if ((int)interpolation_times.size() == 0) {
        std::cout << "Warning: Vector of interpolation times is empty. Returning empty TimeSeries.";
    }
    
    TimeSeries<double> interpolated(simulation_result.get_num_elements());

    int pointer_interp{};
    // add first time point of interpolation times in case it equals the first time point of simulation (up to tolerance)
    // this is necessary even if the tolerance is 0 due to the way the comparison in the loop is implemented (< and >=)
    if (simulation_result.get_time(0) >= interpolation_times[0]) {
        interpolated.add_time_point(interpolation_times[0], simulation_result[0]);
        ++pointer_interp;
    }
    
    //interpolate between pair of time points that lie on either side of each interpolation point
    for (int pointer_sim = 0; pointer_sim < simulation_result.get_num_time_points() - 1 && pointer_interp < (int)interpolation_times.size();) {
        //only go to next pair of time points if no time point is added.
        //otherwise check the same time points again
        //in case there is more than one interpolation point between the two time points
        if (simulation_result.get_time(pointer_sim) < interpolation_times[pointer_interp] && simulation_result.get_time(pointer_sim + 1) >= interpolation_times[pointer_interp] )
        {
            interpolated.add_time_point(interpolation_times[pointer_interp],
                                        linear_interpolation(interpolation_times[pointer_interp], simulation_result.get_time(pointer_sim), simulation_result.get_time(pointer_sim + 1), simulation_result[pointer_sim], simulation_result[pointer_sim + 1]));
            ++pointer_interp;
        }
        else {
            ++pointer_sim;
        }
    }
    
    // add last time point of interpolation times in case it equals the first time point of simulation (up to tolerance)
    if (simulation_result.get_last_time() < interpolation_times[pointer_interp]) {
        interpolated.add_time_point(interpolation_times[pointer_interp], simulation_result.get_last_value());
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


} // namespace mio
