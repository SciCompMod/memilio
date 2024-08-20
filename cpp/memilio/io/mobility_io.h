/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Wadim Koslow, Martin J. Kuehn
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
#ifndef MIO_MOBILITY_IO_H
#define MIO_MOBILITY_IO_H

#include "memilio/io/json_serializer.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"

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
IOResult<Eigen::MatrixXd> read_mobility_formatted(const std::string& filename);

/**
 * @brief Reads txt mobility data or contact which is given by values only
 *        and separated by spaces. Writes it into a NxN Eigen 
 *        Matrix, where N is the number of regions
 * @param filename name of file to be read
 */
IOResult<Eigen::MatrixXd> read_mobility_plain(const std::string& filename);

/**
 * @brief Reads txt file containing the duration of stay in each county.
          Writes it into a Eigen vector of size N, where N is the number of local entites.
 * @param filename name of file to be read
 * @return IOResult<Eigen::MatrixXd> containing the duration of stay in each local entity
 */
IOResult<Eigen::MatrixXd> read_duration_stay(const std::string& filename);

/**
 * @brief For each edge we have the path from the start node to the end node. This functions reads the file and returns the path for each edge.
 * 
 * @param filename Filename of the file containing the paths.
 * @return IOResult<std::vector<std::vector<std::vector<int>>>> contains the paths for each edge. 
 */
IOResult<std::vector<std::vector<std::vector<int>>>> read_path_mobility(const std::string& filename);

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

} // namespace mio

#endif // MIO_MOBILITY_IO_H
