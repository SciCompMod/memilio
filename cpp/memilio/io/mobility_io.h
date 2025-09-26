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
#ifndef MEMILIO_IO_MOBILITY_IO_H
#define MEMILIO_IO_MOBILITY_IO_H

#include "memilio/io/json_serializer.h"
#include "memilio/mobility/graph.h"
#include "memilio/data/analyze_result.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
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

/**
 * @brief Splits string into a Vector of strings according to delimiter
 * @param s string which is splitted
 * @param delimiter sign at which to split string
 */
std::vector<std::string> split(const std::string& s, char delimiter);

/**
 * @brief Counts lines of txt file
 * @param filename name of file which is counted
 */
IOResult<int> count_lines(const std::string& filename);

/**
 * @brief Reads formatted mobility or contact data which is given in columns
 *          from_str	to_str	from_rs	    to_rs	count_abs
 *        and separated by tabs. Writes it into a NxN Eigen Matrix,
 *        where N is the number of regions
 * @param filename name of file to be read
 */
template <typename FP>
IOResult<Eigen::MatrixX<FP>> read_mobility_formatted(const std::string& filename)
{
    BOOST_OUTCOME_TRY(auto&& num_lines, count_lines(filename));

    if (num_lines == 0) {
        return success(Eigen::MatrixX<FP>(0, 0));
    }

    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, filename);
    }

    Eigen::MatrixX<FP> txt_matrix(num_lines - 1, 3);
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
                txt_matrix(linenumber - 1, 0) = static_cast<FP>(std::stoi(line[2]));
                txt_matrix(linenumber - 1, 1) = static_cast<FP>(std::stoi(line[3]));
                txt_matrix(linenumber - 1, 2) = static_cast<FP>(std::stod(line[4]));
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

    Eigen::MatrixX<FP> mobility = Eigen::MatrixX<FP>::Zero(ids.size(), ids.size());

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

/**
 * @brief Reads txt mobility data or contact which is given by values only
 *        and separated by spaces. Writes it into a NxN Eigen
 *        Matrix, where N is the number of regions
 * @param filename name of file to be read
 */
template <typename FP>
IOResult<Eigen::MatrixX<FP>> read_mobility_plain(const std::string& filename)
{
    BOOST_OUTCOME_TRY(auto&& num_lines, count_lines(filename));

    if (num_lines == 0) {
        return success(Eigen::MatrixX<FP>(0, 0));
    }

    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, filename);
    }

    Eigen::MatrixX<FP> mobility(num_lines, num_lines);

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
                mobility(i, j) = static_cast<FP>(std::stod(line[j]));
            }
            linenumber++;
        }
    }
    catch (std::runtime_error& ex) {
        return failure(StatusCode::InvalidFileFormat, filename + ": " + ex.what());
    }

    return success(mobility);
}

#ifdef MEMILIO_HAS_JSONCPP

/**
 * @brief creates json files for each node in a simulation graph.
 * Creates two files per node: one contains the models and its parameters, one contains the outgoing edges.
 * @param graph Graph which should be written
 * @param directory directory where files should be stored
 * @param ioflags flags that set the behavior of serialization; see mio::IOFlags
 */
template <typename FP, class Model>
IOResult<void> write_graph(const Graph<Model, MobilityParameters<FP>>& graph, const std::string& directory,
                           int ioflags = IOF_None)
{
    assert(graph.nodes().size() > 0 && "Graph Nodes are empty");

    std::string abs_path;
    BOOST_OUTCOME_TRY(auto&& created, create_directory(directory, abs_path));

    if (created) {
        log_info("Results are stored in {:s}/results.", mio::get_current_dir_name());
    }
    else {
        log_info("Results are stored in {:s}/results. Files from previous "
                 "graph will be "
                 "overwritten",
                 mio::get_current_dir_name());
    }

    //write two files per node
    //one file that contains outgoing edges from the node
    //one file for the model (parameters and population)
    for (auto inode = size_t(0); inode < graph.nodes().size(); ++inode) {
        //node
        const auto node = graph.nodes()[inode];
        BOOST_OUTCOME_TRY(auto&& js_node_model, serialize_json(node.property, ioflags));
        Json::Value js_node(Json::objectValue);
        js_node["NodeId"]  = node.id;
        js_node["Model"]   = js_node_model;
        auto node_filename = path_join(abs_path, "GraphNode" + std::to_string(inode) + ".json");
        BOOST_OUTCOME_TRY(write_json(node_filename, js_node));

        //list of edges
        auto out_edges = graph.out_edges(inode);
        if (out_edges.size()) {
            Json::Value js_edges(Json::arrayValue);
            for (auto& e : graph.out_edges(inode)) {
                BOOST_OUTCOME_TRY(auto&& js_edge_params, serialize_json(e.property, ioflags));
                Json::Value js_edge{Json::objectValue};
                js_edge["StartNodeIndex"] = Json::UInt64(e.start_node_idx);
                js_edge["EndNodeIndex"]   = Json::UInt64(e.end_node_idx);
                js_edge["Parameters"]     = js_edge_params;
                js_edges.append(std::move(js_edge));
            }
            auto edge_filename = path_join(abs_path, "GraphEdges_node" + std::to_string(inode) + ".json");
            BOOST_OUTCOME_TRY(write_json(edge_filename, js_edges));
        }
    }

    return success();
}

/**
 * @brief reads graph json files and returns a simulation graph.
 * See write_graph() for information of expected files.
 * @tparam the type of the simulation model.
 * @param directory directory from where graph should be read.
 * @param ioflags flags that set the behavior of serialization; see mio::IOFlags
 * @param read_edges boolean value that decides whether the edges of the graph should also be read in.
 */
template <typename FP, class Model>
IOResult<Graph<Model, MobilityParameters<FP>>> read_graph(const std::string& directory, int ioflags = IOF_None,
                                                          bool read_edges = true)
{
    std::string abs_path;
    if (!file_exists(directory, abs_path)) {
        log_error("Directory {} does not exist.", directory);
        return failure(StatusCode::FileNotFound, directory);
    }

    auto graph = Graph<Model, MobilityParameters<FP>>{};

    //read nodes, as many as files are available
    for (auto inode = 0;; ++inode) {
        auto node_filename = path_join(abs_path, "GraphNode" + std::to_string(inode) + ".json");
        if (!file_exists(node_filename, node_filename)) {
            break;
        }
        BOOST_OUTCOME_TRY(auto&& js_node, read_json(node_filename));
        if (!js_node["NodeId"].isInt()) {
            log_error("NodeId field must be an integer.");
            return failure(StatusCode::InvalidType, node_filename + ", NodeId must be an integer.");
        }
        auto node_id = js_node["NodeId"].asInt();
        BOOST_OUTCOME_TRY(auto&& model, deserialize_json(js_node["Model"], Tag<Model>{}, ioflags));
        graph.add_node(node_id, model);
    }

    //read edges; nodes must already be available for that)
    if (read_edges) {
        for (auto inode = size_t(0); inode < graph.nodes().size(); ++inode) {
            //list of edges
            auto edge_filename = path_join(abs_path, "GraphEdges_node" + std::to_string(inode) + ".json");
            BOOST_OUTCOME_TRY(auto&& js_edges, read_json(edge_filename));

            for (auto& e : js_edges) {
                auto start_node_idx  = inode;
                auto js_end_node_idx = e["EndNodeIndex"];
                if (!js_end_node_idx.isUInt64()) {
                    log_error("EndNodeIndex must be an integer.");
                    return failure(StatusCode::InvalidType, edge_filename + ", EndNodeIndex must be an integer.");
                }
                auto end_node_idx = js_end_node_idx.asUInt64();
                if (end_node_idx >= graph.nodes().size()) {
                    log_error("EndNodeIndex not in range of number of graph nodes.");
                    return failure(StatusCode::OutOfRange,
                                   edge_filename + ", EndNodeIndex not in range of number of graph nodes.");
                }
                BOOST_OUTCOME_TRY(auto&& parameters,
                                  deserialize_json(e["Parameters"], Tag<MobilityParameters<FP>>{}, ioflags));
                graph.add_edge(start_node_idx, end_node_idx, parameters);
            }
        }
    }

    return success(graph);
}

#endif //MEMILIO_HAS_JSONCPP
#ifdef MEMILIO_HAS_HDF5
/**
 * @brief Save the results of the edges for a single graph simulation run.
 * @param[in] results Simulation results per edge of the graph.
 * @param[in] ids Identifiers for the start and end node of the edges.
 * @param[in] filename Name of the file where the results will be saved.
 * @return Any io errors that occur during writing of the files.
 */
template <typename FP>
IOResult<void> save_edges(const std::vector<TimeSeries<FP>>& results, const std::vector<std::pair<int, int>>& ids,
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

            auto values_t = std::vector<FP>(result.get_times().begin(), result.get_times().end());
            MEMILIO_H5_CHECK(H5Dwrite(dset_t.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values_t.data()),
                             StatusCode::UnknownError,
                             "Failed to write 'Time' data in group " + h5group_name + " in the file: " + filename);

            int start_id = ids[edge_indx].first;
            auto total =
                Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(num_timepoints, num_elements)
                    .eval();
            while (edge_indx < num_edges && ids[edge_indx].first == start_id) {
                const auto& result_edge = results[edge_indx];
                auto edge_result        = Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(
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

/**
 * @brief Saves the results of a simulation for each edge in the graph.
 * @param[in] ensemble_edges Simulation results for each run for each edge.
 * @param[in] pairs_edges Identifiers for the start and end node of the edges.
 * @param[in] result_dir Top level directory for all results of the parameter study.
 * @param[in] save_single_runs [Default: true] Defines if single run results are written.
 * @param[in] save_percentiles [Default: true] Defines if percentiles are written.
 * @return Any io errors that occur during writing of the files.
 */
template <typename FP>
IOResult<void> save_edges(const std::vector<std::vector<TimeSeries<FP>>>& ensemble_edges,
                          const std::vector<std::pair<int, int>>& pairs_edges, const fs::path& result_dir,
                          bool save_single_runs = true, bool save_percentiles = true)
{
    //save results and sum of results over nodes
    auto ensemble_edges_sum = sum_nodes(ensemble_edges);
    if (save_single_runs) {
        for (size_t i = 0; i < ensemble_edges_sum.size(); ++i) {
            BOOST_OUTCOME_TRY(save_edges<FP>(ensemble_edges[i], pairs_edges,
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
            auto ensemble_edges_p05 = ensemble_percentile<FP>(ensemble_edges, 0.05);
            auto ensemble_edges_p25 = ensemble_percentile<FP>(ensemble_edges, 0.25);
            auto ensemble_edges_p50 = ensemble_percentile<FP>(ensemble_edges, 0.50);
            auto ensemble_edges_p75 = ensemble_percentile<FP>(ensemble_edges, 0.75);
            auto ensemble_edges_p95 = ensemble_percentile<FP>(ensemble_edges, 0.95);

            BOOST_OUTCOME_TRY(save_edges<FP>(ensemble_edges_p05, pairs_edges, (result_dir_p05 / "Edges.h5").string()));
            BOOST_OUTCOME_TRY(save_edges<FP>(ensemble_edges_p25, pairs_edges, (result_dir_p25 / "Edges.h5").string()));
            BOOST_OUTCOME_TRY(save_edges<FP>(ensemble_edges_p50, pairs_edges, (result_dir_p50 / "Edges.h5").string()));
            BOOST_OUTCOME_TRY(save_edges<FP>(ensemble_edges_p75, pairs_edges, (result_dir_p75 / "Edges.h5").string()));
            BOOST_OUTCOME_TRY(save_edges<FP>(ensemble_edges_p95, pairs_edges, (result_dir_p95 / "Edges.h5").string()));
        }
    }
    return success();
}

#endif //MEMILIO_HAS_HDF5

} // namespace mio

#endif // MEMILIO_IO_MOBILITY_IO_H
