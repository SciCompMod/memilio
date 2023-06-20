/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kuehn
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
#ifndef MEMILIO_IO_RESULT_IO_H
#define MEMILIO_IO_RESULT_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_HDF5

#include "memilio/math/eigen_util.h"
#include "memilio/utils/time_series.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/mobility_io.h"
#include "memilio/io/io.h"
#include "boost/filesystem.hpp"

namespace fs = boost::filesystem;
namespace mio
{

/**
 * @brief Save the results of a graph simulation run.
 * @param result Simulation results per node of the graph.
 * @param ids Identifiers for each node of the graph. 
 * @param num_groups Number of groups in the results.
 * @param filename Name of file
 * @return Any io errors that occur during writing of the files. 
 */
IOResult<void> save_result(const std::vector<TimeSeries<double>>& result, const std::vector<int>& ids, int num_groups,
                           const std::string& filename);

class SimulationResult
{
public:
    /**
     * @brief Standard constructor of SimulationResult
     * @param num_groups Number of groups or subpopulations in the simulation.
     * @param num_infectionstates Number of infection states in the considered simulation.
     */
    SimulationResult(int num_groups, int num_infectionstates)
        : m_groups(num_groups * num_infectionstates)
        , m_totals(num_infectionstates)
    {
    }

    /**
     * @brief Constructor of SimulationResult storing time, groups, and total sums of all groups
     * @param groups Simulation results of individual groups.
     * @param total Simulation results as the sum over all groups.
     */
    SimulationResult(const TimeSeries<double>& groups, const TimeSeries<double>& totals)
        : m_groups(groups)
        , m_totals(totals)
    {
    }

    /**
     * @brief Simulation results of individual groups.
     */
    const TimeSeries<double>& get_groups() const
    {
        return m_groups;
    }

    /**
     * @brief Simulation results of the sum over all groups.
     */
    const TimeSeries<double>& get_totals() const
    {
        return m_totals;
    }

private:
    TimeSeries<double> m_groups;
    TimeSeries<double> m_totals;
};

/**
 * @brief Read simulation result from h5 file.
 * @param filename name of the file to be read.
 */
IOResult<std::vector<SimulationResult>> read_result(const std::string& filename);

/**
 * Save the results and the parameters of a single graph simulation run.
 * Creates a new subdirectory for each run according to run_idx.
 * @param result Simulation results per node of the graph.
 * @param params Parameters used for the simulation run.
 * @param ids Identifiers for each node of the graph. 
 * @param result_dir Top level directory for all results of the parameter study.
 * @param run_idx Index of the run; used in directory name.
 * @return Any io errors that occur during writing of the files.
 */
template <class Model>
IOResult<void> save_result_with_params(const std::vector<TimeSeries<double>>& result, const std::vector<Model>& params,
                                       const std::vector<int>& county_ids, const fs::path& result_dir, size_t run_idx)
{
    auto result_dir_run = result_dir / ("run" + std::to_string(run_idx));
    BOOST_OUTCOME_TRY(create_directory(result_dir_run.string()));
    BOOST_OUTCOME_TRY(save_result(result, county_ids, (int)(size_t)params[0].parameters.get_num_groups(),
                                  (result_dir_run / "Result.h5").string()));
    BOOST_OUTCOME_TRY(write_graph(create_graph_without_edges<Model, MigrationParameters>(params, county_ids),
                                  result_dir_run.string(), IOF_OmitDistributions));
    return success();
}

/**
 * @brief Save the results of a parameter study.
 * Stores different percentiles (p5, p25, p50, p75, p90) and sums of the results and parameters. 
 * @param ensemble_results Result of each simulation run.
 * @param ensemble_params Parameters used for each simulation run.
 * @param county_ids Ids of the county nodes.
 * @param result_dir Top level directory for all results of the parameter study.
 * @param save_single_runs [Default: true] Defines if single run results are written to the disk.
 * @param save_single_runs [Default: true] Defines if percentiles are written to the disk.
 * @return Any io errors that occur during writing of the files.
 */
template <class Model>
IOResult<void> save_results(const std::vector<std::vector<TimeSeries<double>>>& ensemble_results,
                            const std::vector<std::vector<Model>>& ensemble_params, const std::vector<int>& county_ids,
                            const fs::path& result_dir, bool save_single_runs = true, bool save_percentiles = true)
{
    //save results and sum of results over nodes
    auto ensemble_result_sum = sum_nodes(ensemble_results);
    auto num_groups          = (int)(size_t)ensemble_params[0][0].parameters.get_num_groups();
    if (save_single_runs) {
        for (size_t i = 0; i < ensemble_result_sum.size(); ++i) {
            BOOST_OUTCOME_TRY(save_result(ensemble_result_sum[i], {0}, num_groups,
                                          (result_dir / ("results_run" + std::to_string(i) + "_sum.h5")).string()));
            BOOST_OUTCOME_TRY(save_result(ensemble_results[i], county_ids, num_groups,
                                          (result_dir / ("results_run" + std::to_string(i) + ".h5")).string()));
        }
    }

    if (save_percentiles) {
        // make directories for percentiles
        auto result_dir_p05 = result_dir / "p05";
        auto result_dir_p25 = result_dir / "p25";
        auto result_dir_p50 = result_dir / "p50";
        auto result_dir_p75 = result_dir / "p75";
        auto result_dir_p95 = result_dir / "p95";
        BOOST_OUTCOME_TRY(create_directory(result_dir_p05.string()));
        BOOST_OUTCOME_TRY(create_directory(result_dir_p25.string()));
        BOOST_OUTCOME_TRY(create_directory(result_dir_p50.string()));
        BOOST_OUTCOME_TRY(create_directory(result_dir_p75.string()));
        BOOST_OUTCOME_TRY(create_directory(result_dir_p95.string()));

        // save percentiles of results, summed over nodes
        {
            auto ensemble_results_sum_p05 = ensemble_percentile(ensemble_result_sum, 0.05);
            auto ensemble_results_sum_p25 = ensemble_percentile(ensemble_result_sum, 0.25);
            auto ensemble_results_sum_p50 = ensemble_percentile(ensemble_result_sum, 0.50);
            auto ensemble_results_sum_p75 = ensemble_percentile(ensemble_result_sum, 0.75);
            auto ensemble_results_sum_p95 = ensemble_percentile(ensemble_result_sum, 0.95);

            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_sum_p05, {0}, num_groups, (result_dir_p05 / "Results_sum.h5").string()));
            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_sum_p25, {0}, num_groups, (result_dir_p25 / "Results_sum.h5").string()));
            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_sum_p50, {0}, num_groups, (result_dir_p50 / "Results_sum.h5").string()));
            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_sum_p75, {0}, num_groups, (result_dir_p75 / "Results_sum.h5").string()));
            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_sum_p95, {0}, num_groups, (result_dir_p95 / "Results_sum.h5").string()));
        }

        // save percentiles of results
        {
            auto ensemble_results_p05 = ensemble_percentile(ensemble_results, 0.05);
            auto ensemble_results_p25 = ensemble_percentile(ensemble_results, 0.25);
            auto ensemble_results_p50 = ensemble_percentile(ensemble_results, 0.50);
            auto ensemble_results_p75 = ensemble_percentile(ensemble_results, 0.75);
            auto ensemble_results_p95 = ensemble_percentile(ensemble_results, 0.95);

            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_p05, county_ids, num_groups, (result_dir_p05 / "Results.h5").string()));
            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_p25, county_ids, num_groups, (result_dir_p25 / "Results.h5").string()));
            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_p50, county_ids, num_groups, (result_dir_p50 / "Results.h5").string()));
            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_p75, county_ids, num_groups, (result_dir_p75 / "Results.h5").string()));
            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_p95, county_ids, num_groups, (result_dir_p95 / "Results.h5").string()));
        }

        // save percentiles of parameters
        {
            auto ensemble_params_p05 = ensemble_params_percentile(ensemble_params, 0.05);
            auto ensemble_params_p25 = ensemble_params_percentile(ensemble_params, 0.25);
            auto ensemble_params_p50 = ensemble_params_percentile(ensemble_params, 0.50);
            auto ensemble_params_p75 = ensemble_params_percentile(ensemble_params, 0.75);
            auto ensemble_params_p95 = ensemble_params_percentile(ensemble_params, 0.95);

            auto make_graph = [&county_ids](auto&& params) {
                return create_graph_without_edges<Model, MigrationParameters>(params, county_ids);
            };
            BOOST_OUTCOME_TRY(
                write_graph(make_graph(ensemble_params_p05), result_dir_p05.string(), IOF_OmitDistributions));
            BOOST_OUTCOME_TRY(
                write_graph(make_graph(ensemble_params_p25), result_dir_p25.string(), IOF_OmitDistributions));
            BOOST_OUTCOME_TRY(
                write_graph(make_graph(ensemble_params_p50), result_dir_p50.string(), IOF_OmitDistributions));
            BOOST_OUTCOME_TRY(
                write_graph(make_graph(ensemble_params_p75), result_dir_p75.string(), IOF_OmitDistributions));
            BOOST_OUTCOME_TRY(
                write_graph(make_graph(ensemble_params_p95), result_dir_p95.string(), IOF_OmitDistributions));
        }
    }
    return success();
}

} // namespace mio

#endif // MEMILIO_HAS_HDF5

#endif // MEMILIO_IO_RESULT_IO_H
