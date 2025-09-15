/* 
* Copyright (C) 2020-2025 MEmilio
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
template <typename FP>
IOResult<void> save_result(const std::vector<TimeSeries<FP>>& results, const std::vector<int>& ids, int num_groups,
                           const std::string& filename)
{
    int region_idx = 0;
    H5File file{H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(file.id, StatusCode::FileNotFound, filename);
    for (auto& result : results) {
        auto h5group_name = "/" + std::to_string(ids[region_idx]);
        H5Group region_h5group{H5Gcreate(file.id, h5group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
        MEMILIO_H5_CHECK(region_h5group.id, StatusCode::UnknownError,
                         "Group could not be created (" + h5group_name + ")");

        const int num_timepoints      = static_cast<int>(result.get_num_time_points());
        const int num_infectionstates = (int)result.get_num_elements() / num_groups;

        hsize_t dims_t[] = {static_cast<hsize_t>(num_timepoints)};
        H5DataSpace dspace_t{H5Screate_simple(1, dims_t, NULL)};
        MEMILIO_H5_CHECK(dspace_t.id, StatusCode::UnknownError, "Time DataSpace could not be created.");
        H5DataSet dset_t{H5Dcreate(region_h5group.id, "Time", H5T_NATIVE_DOUBLE, dspace_t.id, H5P_DEFAULT, H5P_DEFAULT,
                                   H5P_DEFAULT)};
        MEMILIO_H5_CHECK(dset_t.id, StatusCode::UnknownError, "Time DataSet could not be created (Time).");
        auto values_t = std::vector<FP>(result.get_times().begin(), result.get_times().end());
        MEMILIO_H5_CHECK(H5Dwrite(dset_t.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values_t.data()),
                         StatusCode::UnknownError, "Time data could not be written.");

        auto total = Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(num_timepoints,
                                                                                              num_infectionstates)
                         .eval();

        for (int group_idx = 0; group_idx <= num_groups; ++group_idx) {
            auto group = Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(num_timepoints,
                                                                                                  num_infectionstates)
                             .eval();
            if (group_idx < num_groups) {
                for (Eigen::Index t_idx = 0; t_idx < result.get_num_time_points(); ++t_idx) {
                    auto v           = result[t_idx].transpose().eval();
                    auto group_slice = mio::slice(v, {group_idx * num_infectionstates, num_infectionstates});
                    mio::slice(group, {t_idx, 1}, {0, num_infectionstates}) = group_slice;
                    mio::slice(total, {t_idx, 1}, {0, num_infectionstates}) += group_slice;
                }
            }

            hsize_t dims_values[] = {static_cast<hsize_t>(num_timepoints), static_cast<hsize_t>(num_infectionstates)};
            H5DataSpace dspace_values{H5Screate_simple(2, dims_values, NULL)};
            MEMILIO_H5_CHECK(dspace_values.id, StatusCode::UnknownError, "Values DataSpace could not be created.");
            auto dset_name = group_idx == num_groups ? std::string("Total") : "Group" + std::to_string(group_idx + 1);
            H5DataSet dset_values{H5Dcreate(region_h5group.id, dset_name.c_str(), H5T_NATIVE_DOUBLE, dspace_values.id,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
            MEMILIO_H5_CHECK(dset_values.id, StatusCode::UnknownError, "Values DataSet could not be created.");

            MEMILIO_H5_CHECK(H5Dwrite(dset_values.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                      group_idx == num_groups ? total.data() : group.data()),
                             StatusCode::UnknownError, "Values data could not be written.");
        }
        region_idx++;
    }
    return success();
}

template <typename FP>
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
    SimulationResult(const TimeSeries<FP>& groups, const TimeSeries<FP>& totals)
        : m_groups(groups)
        , m_totals(totals)
    {
    }

    /**
     * @brief Simulation results of individual groups.
     */
    const TimeSeries<FP>& get_groups() const
    {
        return m_groups;
    }

    /**
     * @brief Simulation results of the sum over all groups.
     */
    const TimeSeries<FP>& get_totals() const
    {
        return m_totals;
    }

private:
    TimeSeries<FP> m_groups;
    TimeSeries<FP> m_totals;
};

// Forward declaration of store_group_name
herr_t store_group_name(hid_t /*id*/, const char* name, const H5L_info_t* /*linfo*/, void* opdata);
/**
 * @brief Read simulation result from h5 file.
 * @param filename name of the file to be read.
 */
template <typename FP>
IOResult<std::vector<SimulationResult<FP>>> read_result(const std::string& filename)
{
    std::vector<SimulationResult<FP>> results;

    H5File file{H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(file.id, StatusCode::FileNotFound, filename);

    std::vector<std::string> h5group_names;
    MEMILIO_H5_CHECK(H5Literate(file.id, H5_INDEX_NAME, H5_ITER_INC, NULL, &store_group_name, &h5group_names),
                     StatusCode::UnknownError, "Group names could not be read.");

    for (auto& h5group_name : h5group_names) {
        H5Group region_h5group{H5Gopen(file.id, h5group_name.c_str(), H5P_DEFAULT)};
        MEMILIO_H5_CHECK(region_h5group.id, StatusCode::UnknownError,
                         "Group could not be opened (" + h5group_name + ")");

        std::vector<std::string> h5dset_names;
        MEMILIO_H5_CHECK(
            H5Literate(region_h5group.id, H5_INDEX_NAME, H5_ITER_INC, NULL, &store_group_name, &h5dset_names),
            StatusCode::UnknownError, "Dataset names could not be read.");

        std::string h5_key = std::any_of(h5dset_names.begin(), h5dset_names.end(),
                                         [](const std::string& str) {
                                             return str.find("Group") == 0;
                                         })
                                 ? "Group"
                                 : "End";

        auto num_groups = (Eigen::Index)std::count_if(h5dset_names.begin(), h5dset_names.end(), [&h5_key](auto&& str) {
            return str.find(h5_key) != std::string::npos;
        });

        H5DataSet dataset_t{H5Dopen(region_h5group.id, "Time", H5P_DEFAULT)};
        MEMILIO_H5_CHECK(dataset_t.id, StatusCode::UnknownError, "Time DataSet could not be read.");

        // dataset dimensions
        H5DataSpace dataspace_t{H5Dget_space(dataset_t.id)};
        MEMILIO_H5_CHECK(dataspace_t.id, StatusCode::UnknownError, "Time DataSpace could not be read.");
        const auto n_dims_t = 1;
        hsize_t dims_t[n_dims_t];
        H5Sget_simple_extent_dims(dataspace_t.id, dims_t, NULL);

        auto num_timepoints = Eigen::Index(dims_t[0]);
        auto time           = std::vector<FP>(num_timepoints);
        MEMILIO_H5_CHECK(H5Dread(dataset_t.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time.data()),
                         StatusCode::UnknownError, "Time data could not be read.");

        auto dataset_name_total("/" + h5group_name + "/Total");
        H5DataSet dataset_total{H5Dopen(file.id, dataset_name_total.c_str(), H5P_DEFAULT)};
        MEMILIO_H5_CHECK(dataset_total.id, StatusCode::UnknownError, "Totals DataSet could not be read.");

        //read data space dimensions
        H5DataSpace dataspace_total{H5Dget_space(dataset_total.id)};
        MEMILIO_H5_CHECK(dataspace_total.id, StatusCode::UnknownError, "Totals DataSpace could not be read.");
        const auto n_dims_total = 2;
        hsize_t dims_total[n_dims_total];
        H5Sget_simple_extent_dims(dataspace_total.id, dims_total, NULL);
        if (num_timepoints != Eigen::Index(dims_total[0])) {
            return failure(StatusCode::InvalidFileFormat, "Number of time points does not match.");
        }
        auto num_infectionstates = Eigen::Index(dims_total[1]);

        auto total_values =
            Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(num_timepoints, num_infectionstates);
        MEMILIO_H5_CHECK(
            H5Dread(dataset_total.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, total_values.data()),
            StatusCode::UnknownError, "Totals data could not be read");

        auto totals = TimeSeries<FP>(num_infectionstates);
        totals.reserve(num_timepoints);
        for (auto t_idx = 0; t_idx < num_timepoints; ++t_idx) {
            totals.add_time_point(time[t_idx], slice(total_values, {t_idx, 1}, {0, num_infectionstates}).transpose());
        }

        auto groups = TimeSeries<FP>(num_infectionstates * num_groups);
        groups.reserve(num_timepoints);
        for (Eigen::Index t_idx = 0; t_idx < num_timepoints; ++t_idx) {
            groups.add_time_point(time[t_idx]);
        }

        std::vector<int> h5_key_indices;
        // Extract group indices from h5dset_names
        for (const auto& name : h5dset_names) {
            if (name.find(h5_key) == 0) {
                h5_key_indices.push_back(std::stoi(name.substr(h5_key.size())));
            }
        }

        for (auto h5_key_indx = size_t(0); h5_key_indx < h5_key_indices.size(); h5_key_indx++) {
            auto group_name = "/" + h5group_name + "/" + h5_key + std::to_string(h5_key_indices[h5_key_indx]);
            H5DataSet dataset_values{H5Dopen(file.id, group_name.c_str(), H5P_DEFAULT)};
            MEMILIO_H5_CHECK(dataset_values.id, StatusCode::UnknownError, "Values DataSet could not be read.");

            //read data space dimensions
            H5DataSpace dataspace_values{H5Dget_space(dataset_values.id)};
            MEMILIO_H5_CHECK(dataspace_values.id, StatusCode::UnknownError, "Values DataSpace could not be read.");
            const auto n_dims_values = 2;
            hsize_t dims_values[n_dims_values];
            H5Sget_simple_extent_dims(dataspace_values.id, dims_values, NULL);
            if (num_timepoints != Eigen::Index(dims_values[0])) {
                return failure(StatusCode::InvalidFileFormat, "Number of time points does not match.");
            }
            if (num_infectionstates != Eigen::Index(dims_values[1])) {
                return failure(StatusCode::InvalidFileFormat, "Number of infection states does not match.");
            }

            auto group_values =
                Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(num_timepoints, num_infectionstates);
            MEMILIO_H5_CHECK(
                H5Dread(dataset_values.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_values.data()),
                StatusCode::UnknownError, "Values data could not be read");

            for (Eigen::Index idx_t = 0; idx_t < num_timepoints; idx_t++) {
                for (Eigen::Index idx_c = 0; idx_c < num_infectionstates; idx_c++) {
                    groups[idx_t][num_infectionstates * h5_key_indx + idx_c] = group_values(idx_t, idx_c);
                }
            }
        }

        results.push_back(SimulationResult<FP>(groups, totals));
    }
    return success(results);
}

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
template <typename FP, class Model>
IOResult<void> save_result_with_params(const std::vector<TimeSeries<FP>>& result, const std::vector<Model>& params,
                                       const std::vector<int>& county_ids, const fs::path& result_dir, size_t run_idx)
{
    auto result_dir_run = result_dir / ("run" + std::to_string(run_idx));
    BOOST_OUTCOME_TRY(create_directory(result_dir_run.string()));
    BOOST_OUTCOME_TRY(save_result<FP>(result, county_ids, (int)(size_t)params[0].parameters.get_num_groups(),
                                      (result_dir_run / "Result.h5").string()));
    BOOST_OUTCOME_TRY(write_graph<FP>(create_graph_without_edges<Model, MobilityParameters<FP>>(params, county_ids),
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
template <typename FP, class Model>
IOResult<void> save_results(const std::vector<std::vector<TimeSeries<FP>>>& ensemble_results,
                            const std::vector<std::vector<Model>>& ensemble_params, const std::vector<int>& county_ids,
                            const fs::path& result_dir, bool save_single_runs = true, bool save_percentiles = true)
{
    //save results and sum of results over nodes
    auto ensemble_result_sum = sum_nodes(ensemble_results);
    auto num_groups          = (int)(size_t)ensemble_params[0][0].parameters.get_num_groups();
    if (save_single_runs) {
        for (size_t i = 0; i < ensemble_result_sum.size(); ++i) {
            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_result_sum[i], {0}, num_groups,
                                              (result_dir / ("results_run" + std::to_string(i) + "_sum.h5")).string()));
            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_results[i], county_ids, num_groups,
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
            auto ensemble_results_sum_p05 = ensemble_percentile<FP>(ensemble_result_sum, 0.05);
            auto ensemble_results_sum_p25 = ensemble_percentile<FP>(ensemble_result_sum, 0.25);
            auto ensemble_results_sum_p50 = ensemble_percentile<FP>(ensemble_result_sum, 0.50);
            auto ensemble_results_sum_p75 = ensemble_percentile<FP>(ensemble_result_sum, 0.75);
            auto ensemble_results_sum_p95 = ensemble_percentile<FP>(ensemble_result_sum, 0.95);

            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_results_sum_p05, {0}, num_groups,
                                              (result_dir_p05 / "Results_sum.h5").string()));
            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_results_sum_p25, {0}, num_groups,
                                              (result_dir_p25 / "Results_sum.h5").string()));
            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_results_sum_p50, {0}, num_groups,
                                              (result_dir_p50 / "Results_sum.h5").string()));
            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_results_sum_p75, {0}, num_groups,
                                              (result_dir_p75 / "Results_sum.h5").string()));
            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_results_sum_p95, {0}, num_groups,
                                              (result_dir_p95 / "Results_sum.h5").string()));
        }

        // save percentiles of results
        {
            auto ensemble_results_p05 = ensemble_percentile<FP>(ensemble_results, 0.05);
            auto ensemble_results_p25 = ensemble_percentile<FP>(ensemble_results, 0.25);
            auto ensemble_results_p50 = ensemble_percentile<FP>(ensemble_results, 0.50);
            auto ensemble_results_p75 = ensemble_percentile<FP>(ensemble_results, 0.75);
            auto ensemble_results_p95 = ensemble_percentile<FP>(ensemble_results, 0.95);

            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_results_p05, county_ids, num_groups,
                                              (result_dir_p05 / "Results.h5").string()));
            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_results_p25, county_ids, num_groups,
                                              (result_dir_p25 / "Results.h5").string()));
            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_results_p50, county_ids, num_groups,
                                              (result_dir_p50 / "Results.h5").string()));
            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_results_p75, county_ids, num_groups,
                                              (result_dir_p75 / "Results.h5").string()));
            BOOST_OUTCOME_TRY(save_result<FP>(ensemble_results_p95, county_ids, num_groups,
                                              (result_dir_p95 / "Results.h5").string()));
        }

        // save percentiles of parameters
        {
            auto ensemble_params_p05 = ensemble_params_percentile<FP>(ensemble_params, 0.05);
            auto ensemble_params_p25 = ensemble_params_percentile<FP>(ensemble_params, 0.25);
            auto ensemble_params_p50 = ensemble_params_percentile<FP>(ensemble_params, 0.50);
            auto ensemble_params_p75 = ensemble_params_percentile<FP>(ensemble_params, 0.75);
            auto ensemble_params_p95 = ensemble_params_percentile<FP>(ensemble_params, 0.95);

            auto make_graph = [&county_ids](auto&& params) {
                return create_graph_without_edges<Model, MobilityParameters<FP>>(params, county_ids);
            };
            BOOST_OUTCOME_TRY(
                write_graph<FP>(make_graph(ensemble_params_p05), result_dir_p05.string(), IOF_OmitDistributions));
            BOOST_OUTCOME_TRY(
                write_graph<FP>(make_graph(ensemble_params_p25), result_dir_p25.string(), IOF_OmitDistributions));
            BOOST_OUTCOME_TRY(
                write_graph<FP>(make_graph(ensemble_params_p50), result_dir_p50.string(), IOF_OmitDistributions));
            BOOST_OUTCOME_TRY(
                write_graph<FP>(make_graph(ensemble_params_p75), result_dir_p75.string(), IOF_OmitDistributions));
            BOOST_OUTCOME_TRY(
                write_graph<FP>(make_graph(ensemble_params_p95), result_dir_p95.string(), IOF_OmitDistributions));
        }
    }
    return success();
}

} // namespace mio

#endif // MEMILIO_HAS_HDF5

#endif // MEMILIO_IO_RESULT_IO_H
