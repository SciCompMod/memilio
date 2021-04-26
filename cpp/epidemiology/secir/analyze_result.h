#ifndef EPI_SECIR_ANALYZE_RESULT_H
#define EPI_SECIR_ANALYZE_RESULT_H

#include "epidemiology/utils/time_series.h"
#include "epidemiology/secir/secir.h"
#include "epidemiology/migration/migration.h"

#include <functional>
#include <vector>

namespace epi
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
interpolate_simulation_result(const Graph<ModelNode<Simulation>, MigrationEdge>& graph_result)
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
template <class Simulation>
std::vector<Simulation> ensemble_params_percentile(const std::vector<std::vector<Simulation>>& ensemble_params,
                                                   double p)
{

    assert(p > 0.0 && p < 1.0 && "Invalid percentile value.");

    auto num_runs   = ensemble_params.size();
    auto num_nodes  = ensemble_params[0].size();
    auto num_groups = ensemble_params[0][0].parameters.get_num_groups();

    std::vector<double> single_element_ensemble(num_runs);

    // lamda function that calculates the percentile of a single paramter
    std::vector<Simulation> percentile(num_nodes, Simulation(num_groups));
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
        for (size_t i = 0; i < num_groups; i++) {
            //Population
            for (size_t compart = 0; compart < (size_t)InfectionState::Count; ++compart) {
                param_percentil(
                    node, [ compart, i ](auto&& model) -> auto& {
                        return model.populations[{epi::AgeGroup(i), (InfectionState)compart}];
                    });
            }
            // times
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.times[i].get_incubation(); });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.times[i].get_serialinterval(); });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.times[i].get_infectious_mild(); });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.times[i].get_hospitalized_to_icu(); });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.times[i].get_hospitalized_to_home(); });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.times[i].get_home_to_hospitalized(); });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.times[i].get_icu_to_dead(); });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.times[i].get_icu_to_home(); });
            //probs
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.probabilities[i].get_carrier_infectability();
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.probabilities[i].get_risk_from_symptomatic();
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.probabilities[i].get_test_and_trace_max_risk_from_symptomatic();
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.probabilities[i].get_asymp_per_infectious();
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.probabilities[i].get_hospitalized_per_infectious();
                });
            param_percentil(
                node, [i](auto&& model) -> auto& {
                    return model.parameters.probabilities[i].get_icu_per_hospitalized();
                });
            param_percentil(
                node, [i](auto&& model) -> auto& { return model.parameters.probabilities[i].get_dead_per_icu(); });
        }
        // group independent params
        param_percentil(
            node, [](auto&& model) -> auto& { return model.parameters.get_seasonality(); });
        param_percentil(
            node, [](auto&& model) -> auto& { return model.parameters.get_test_and_trace_capacity(); });

        for (size_t run = 0; run < num_runs; run++) {

            auto const& params           = ensemble_params[run][node];
            single_element_ensemble[run] = params.parameters.get_icu_capacity() * params.populations.get_total();
        }
        std::sort(single_element_ensemble.begin(), single_element_ensemble.end());
        percentile[node].parameters.set_icu_capacity(single_element_ensemble[static_cast<size_t>(num_runs * p)]);
    }
    return percentile;
}

} // namespace epi

#endif //EPI_SECIR_ANALYZE_RESULT_H
