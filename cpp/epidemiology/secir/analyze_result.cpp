#include "epidemiology/secir/analyze_result.h"

#include <algorithm>
#include <cassert>

namespace epi
{

/**
 * TODO: extrapolate first and last point
 */
TimeSeries<double> interpolate_simulation_result(const TimeSeries<double>& simulation_result)
{
    assert(simulation_result.get_num_time_points() > 0 && "TimeSeries must not be empty.");

    const auto t0      = simulation_result.get_time(0);
    const auto tmax    = simulation_result.get_last_time();
    const auto day0    = static_cast<int>(floor(t0));
    const auto day_max = static_cast<int>(ceil(tmax));

    auto day = day0;
    TimeSeries<double> interpolated(simulation_result.get_num_elements());
    interpolated.reserve(day_max - day0 + 1);
    interpolated.add_time_point(day, simulation_result.get_value(0));
    day++;

    //interpolate between pair of time points that lie on either side of each integer day
    for (int i = 0; i < simulation_result.get_num_time_points() - 1;) {
        //only go to next pair of time points if no time point is added.
        //otherwise check the same time points again
        //in case there is more than one day between the two time points
        if (simulation_result.get_time(i) < day && simulation_result.get_time(i + 1) >= day) {
            auto weight = (day - simulation_result.get_time(i)) / (simulation_result.get_time(i + 1) - simulation_result.get_time(i));
            interpolated.add_time_point(day, simulation_result[i] + (simulation_result[i + 1] - simulation_result[i]) * weight);
            ++day;
        }
        else {
            ++i;
        }
    }

    if (day_max > tmax) {
        interpolated.add_time_point(day, simulation_result.get_last_value());
    }

    return interpolated;
}

std::vector<TimeSeries<double>>
interpolate_simulation_result(const Graph<ModelNode<SecirSimulation>, MigrationEdge>& graph_result)
{
    std::vector<TimeSeries<double>> interpolated;
    interpolated.reserve(graph_result.nodes().size());
    std::transform(graph_result.nodes().begin(), graph_result.nodes().end(), std::back_inserter(interpolated),
                   [](auto& n) {
                       return interpolate_simulation_result(n.get_result());
                   });
    return interpolated;
}

std::vector<TimeSeries<double>> ensemble_mean(const std::vector<std::vector<TimeSeries<double>>>& ensemble_result)
{
    auto num_runs        = ensemble_result.size();
    auto num_nodes       = ensemble_result[0].size();
    auto num_time_points = ensemble_result[0][0].get_num_time_points();
    auto num_elements    = ensemble_result[0][0].get_num_elements();

    std::vector<TimeSeries<double>> mean(num_nodes, TimeSeries<double>::zero(num_time_points, num_elements));

    for (int run = 0; run < num_runs; run++) {
        assert(ensemble_result[run].size() == num_nodes && "ensemble results not uniform.");
        for (int node = 0; node < num_nodes; node++) {
            assert(ensemble_result[run][node].get_num_time_points() == num_time_points &&
                   "ensemble results not uniform.");
            for (int time = 0; time < num_time_points; time++) {
                assert(ensemble_result[run][node].get_num_elements() == num_elements &&
                       "ensemble results not uniform.");
                mean[node].get_time(time) = time;
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

} // namespace epi