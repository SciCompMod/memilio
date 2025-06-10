/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Wadim Koslow, Henrik Zunker, Martin J. Kuehn
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
#include "memilio/io/mobility_io.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"

#ifdef MEMILIO_HAS_HDF5
#include "memilio/io/hdf5_cpp.h"
#endif // MEMILIO_HAS_HDF5

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace mio
{

std::vector<std::string> split(const std::string& s, char delimitor)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimitor)) {
        tokens.push_back(token);
    }
    return tokens;
}

IOResult<int> count_lines(const std::string& filename)
{
    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, filename);
    }

    int count = 0;
    std::string line;
    while (getline(file, line)) {
        count++;
    }
    return success(count);
}

IOResult<Eigen::MatrixXd> read_mobility_formatted(const std::string& filename)
{
    BOOST_OUTCOME_TRY(auto&& num_lines, count_lines(filename));

    if (num_lines == 0) {
        return success(Eigen::MatrixXd(0, 0));
    }

    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, filename);
    }

    Eigen::MatrixXd txt_matrix(num_lines - 1, 3);
    std::vector<int> ids;

    try {
        std::string tp;
        int linenumber = 0;
        while (getline(file, tp)) {
            if (linenumber > 0) {
                auto line = split(tp, '\t');
                if (line.size() < 5) {
                    return failure(StatusCode::InvalidFileFormat,
                                   filename + ":" + std::to_string(linenumber) + ": Not enough entries in line.");
                }
                ids.push_back(std::stoi(line[2]));
                txt_matrix(linenumber - 1, 0) = std::stoi(line[2]);
                txt_matrix(linenumber - 1, 1) = std::stoi(line[3]);
                txt_matrix(linenumber - 1, 2) = std::stod(line[4]);
            }
            linenumber++;
        }
    }
    catch (std::runtime_error& ex) {
        return failure(StatusCode::InvalidFileFormat, filename + ": " + ex.what());
    }

    sort(ids.begin(), ids.end());
    std::vector<int>::iterator iter = std::unique(ids.begin(), ids.end());
    ids.resize(std::distance(ids.begin(), iter));

    Eigen::MatrixXd mobility = Eigen::MatrixXd::Zero(ids.size(), ids.size());

    for (int k = 0; k < num_lines - 1; k++) {
        int row_ind = 0;
        int col_ind = 0;
        while (txt_matrix(k, 0) != ids[row_ind]) {
            row_ind++;
        }
        while (txt_matrix(k, 1) != ids[col_ind]) {
            col_ind++;
        }
        mobility(row_ind, col_ind) = txt_matrix(k, 2);
    }

    return success(mobility);
}

IOResult<Eigen::MatrixXd> read_mobility_plain(const std::string& filename)
{
    BOOST_OUTCOME_TRY(auto&& num_lines, count_lines(filename));

    if (num_lines == 0) {
        return success(Eigen::MatrixXd(0, 0));
    }

    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, filename);
    }

    Eigen::MatrixXd mobility(num_lines, num_lines);

    try {
        std::string tp;
        int linenumber = 0;
        while (getline(file, tp)) {
            auto line = split(tp, ' ');
            if (line.size() != size_t(num_lines)) {
                return failure(StatusCode::InvalidFileFormat, filename + ": Not a square matrix.");
            }
            Eigen::Index i = static_cast<Eigen::Index>(linenumber);
            for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(line.size()); j++) {
                mobility(i, j) = std::stod(line[j]);
            }
            linenumber++;
        }
    }
    catch (std::runtime_error& ex) {
        return failure(StatusCode::InvalidFileFormat, filename + ": " + ex.what());
    }

    return success(mobility);
}

#ifdef MEMILIO_HAS_HDF5
IOResult<void> save_edges(const std::vector<std::vector<TimeSeries<double>>>& ensemble_edges,
                          const std::vector<std::pair<int, int>>& pairs_edges, const fs::path& result_dir,
                          bool save_single_runs, bool save_percentiles)
{
    //save results and sum of results over nodes
    auto ensemble_edges_sum = sum_nodes(ensemble_edges);
    if (save_single_runs) {
        for (size_t i = 0; i < ensemble_edges_sum.size(); ++i) {
            BOOST_OUTCOME_TRY(save_edges(ensemble_edges[i], pairs_edges,
                                         (result_dir / ("Edges_run" + std::to_string(i) + ".h5")).string()));
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

        // save percentiles of Edges
        {
            auto ensemble_edges_p05 = ensemble_percentile(ensemble_edges, 0.05);
            auto ensemble_edges_p25 = ensemble_percentile(ensemble_edges, 0.25);
            auto ensemble_edges_p50 = ensemble_percentile(ensemble_edges, 0.50);
            auto ensemble_edges_p75 = ensemble_percentile(ensemble_edges, 0.75);
            auto ensemble_edges_p95 = ensemble_percentile(ensemble_edges, 0.95);

            BOOST_OUTCOME_TRY(save_edges(ensemble_edges_p05, pairs_edges, (result_dir_p05 / "Edges.h5").string()));
            BOOST_OUTCOME_TRY(save_edges(ensemble_edges_p25, pairs_edges, (result_dir_p25 / "Edges.h5").string()));
            BOOST_OUTCOME_TRY(save_edges(ensemble_edges_p50, pairs_edges, (result_dir_p50 / "Edges.h5").string()));
            BOOST_OUTCOME_TRY(save_edges(ensemble_edges_p75, pairs_edges, (result_dir_p75 / "Edges.h5").string()));
            BOOST_OUTCOME_TRY(save_edges(ensemble_edges_p95, pairs_edges, (result_dir_p95 / "Edges.h5").string()));
        }
    }
    return success();
}

IOResult<void> save_edges(const std::vector<TimeSeries<double>>& results, const std::vector<std::pair<int, int>>& ids,
                          const std::string& filename)
{
    const auto num_edges = results.size();
    size_t edge_indx     = 0;
    // H5Fcreate creates a new HDF5 file.
    // H5F_ACC_TRUNC: If the file already exists, H5Fcreate fails. If the file does not exist, it is created and opened with read-write access.
    // H5P_DEFAULT: default data transfer properties are used.
    H5File file{H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(file.id, StatusCode::FileNotFound, "Failed to open the HDF5 file: " + filename);

    while (edge_indx < num_edges) {
        const auto& result = results[edge_indx];
        auto h5group_name  = "/" + std::to_string(ids[edge_indx].first);
        H5Group start_node_h5group{H5Gcreate(file.id, h5group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
        MEMILIO_H5_CHECK(start_node_h5group.id, StatusCode::UnknownError,
                         "Group could not be created (" + h5group_name + ") in the file: " + filename);

        const auto num_timepoints = result.get_num_time_points();
        if (num_timepoints > 0) {
            const auto num_elements = result.get_value(0).size();

            hsize_t dims_t[] = {static_cast<hsize_t>(num_timepoints)};
            H5DataSpace dspace_t{H5Screate_simple(1, dims_t, NULL)};
            MEMILIO_H5_CHECK(dspace_t.id, StatusCode::UnknownError,
                             "Failed to create the DataSpace for 'Time' in group " + h5group_name +
                                 " in the file: " + filename);

            H5DataSet dset_t{H5Dcreate(start_node_h5group.id, "Time", H5T_NATIVE_DOUBLE, dspace_t.id, H5P_DEFAULT,
                                       H5P_DEFAULT, H5P_DEFAULT)};
            MEMILIO_H5_CHECK(dset_t.id, StatusCode::UnknownError,
                             "Failed to create the 'Time' DataSet in group " + h5group_name +
                                 " in the file: " + filename);

            auto values_t = std::vector<double>(result.get_times().begin(), result.get_times().end());
            MEMILIO_H5_CHECK(H5Dwrite(dset_t.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values_t.data()),
                             StatusCode::UnknownError,
                             "Failed to write 'Time' data in group " + h5group_name + " in the file: " + filename);

            int start_id = ids[edge_indx].first;
            auto total   = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(num_timepoints,
                                                                                                        num_elements)
                             .eval();
            while (edge_indx < num_edges && ids[edge_indx].first == start_id) {
                const auto& result_edge = results[edge_indx];
                auto edge_result        = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(
                                       num_timepoints, num_elements)
                                       .eval();
                for (Eigen::Index t_idx = 0; t_idx < result_edge.get_num_time_points(); ++t_idx) {
                    auto v                 = result_edge.get_value(t_idx).transpose().eval();
                    edge_result.row(t_idx) = v;
                    total.row(t_idx) += v;
                }

                hsize_t dims_values[] = {static_cast<hsize_t>(num_timepoints), static_cast<hsize_t>(num_elements)};
                H5DataSpace dspace_values{H5Screate_simple(2, dims_values, NULL)};
                MEMILIO_H5_CHECK(dspace_values.id, StatusCode::UnknownError,
                                 "Failed to create the DataSpace for End" + std::to_string(ids[edge_indx].second) +
                                     " in group " + h5group_name + " in the file: " + filename);

                // End is the target node of the edge
                auto dset_name = "End" + std::to_string(ids[edge_indx].second);
                H5DataSet dset_values{H5Dcreate(start_node_h5group.id, dset_name.c_str(), H5T_NATIVE_DOUBLE,
                                                dspace_values.id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
                MEMILIO_H5_CHECK(dset_values.id, StatusCode::UnknownError,
                                 "Failed to create the DataSet for " + dset_name + " in group " + h5group_name +
                                     " in the file: " + filename);

                MEMILIO_H5_CHECK(
                    H5Dwrite(dset_values.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, edge_result.data()),
                    StatusCode::UnknownError,
                    "Failed to write data for " + dset_name + " in group " + h5group_name +
                        " in the file: " + filename);

                // In the final iteration, we also save the total values
                if (edge_indx == num_edges - 1 || ids[edge_indx + 1].first != start_id) {
                    dset_name = "Total";
                    H5DataSet dset_total{H5Dcreate(start_node_h5group.id, dset_name.c_str(), H5T_NATIVE_DOUBLE,
                                                   dspace_values.id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
                    MEMILIO_H5_CHECK(dset_total.id, StatusCode::UnknownError,
                                     "Failed to create the Total DataSet in group " + h5group_name +
                                         " in the file: " + filename);

                    MEMILIO_H5_CHECK(
                        H5Dwrite(dset_total.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, total.data()),
                        StatusCode::UnknownError,
                        "Failed to write Total data in group " + h5group_name + " in the file: " + filename);
                }
                edge_indx++;
            }
        }
        else {
            log_error("No time points in the TimeSeries for Edge combination " + std::to_string(ids[edge_indx].first) +
                      " -> " + std::to_string(ids[edge_indx].second));
            return failure(mio::StatusCode::InvalidValue);
        }
    }
    return success();
}
#endif // MEMILIO_HAS_HDF5

} // namespace mio
