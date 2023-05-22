/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#ifndef GRAPH_H
#define GRAPH_H

#include <functional>
#include "memilio/utils/stl_util.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/date.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/geography/regions.h"
#include <iostream>

#include "boost/filesystem.hpp"

//is used to provide some paths as function arguments
namespace fs = boost::filesystem;

namespace mio
{

struct OutEdgeBase {
    OutEdgeBase(size_t start)
        : start_node_idx(start)
    {
    }
    size_t start_node_idx;
};
struct InEdgeBase {
    InEdgeBase(size_t end)
        : end_node_idx(end)
    {
    }
    size_t end_node_idx;
};
struct EdgeBase : public OutEdgeBase, public InEdgeBase {
    EdgeBase(size_t start, size_t end)
        : OutEdgeBase(start)
        , InEdgeBase(end)
    {
    }
};

/**
 * @brief represents a node of the graph
 * @param id unique id for the node (e.g. regional keys)
 */
template <class NodePropertyT>
struct Node {
    template <class... Args>
    Node(int node_id, Args&&... args)
        : id{node_id}
        , property(std::forward<Args>(args)...)
    {
    }
    int id;
    NodePropertyT property;
};

/**
 * @brief represents an edge of the graph
 */
template <class EdgePropertyT>
struct Edge : public EdgeBase {
    template <class... Args>
    Edge(size_t start, size_t end, Args&&... args)
        : EdgeBase{start, end}
        , property(std::forward<Args>(args)...)
    {
    }

    EdgePropertyT property;
};

/**
 * @brief comparison operator if node property type is equality comparable
 */
template <class T>
std::enable_if_t<has_eq_op<T>::value, bool> operator==(const Node<T>& n1, const Node<T>& n2)
{
    return n1.id == n2.id && n1.property == n2.property;
}
template <class T>
std::enable_if_t<has_eq_op<T>::value, bool> operator!=(const Node<T>& n1, const Node<T>& n2)
{
    return !(n1 == n2);
}

/**
 * @brief comparison operator if edge property type is equality comparable
 */
template <class T>
std::enable_if_t<has_eq_op<T>::value, bool> operator==(const Edge<T>& e1, const Edge<T>& e2)
{
    return e1.start_node_idx == e2.start_node_idx && e1.end_node_idx == e2.end_node_idx && e1.property == e2.property;
}
template <class T>
std::enable_if_t<has_eq_op<T>::value, bool> operator!=(const Edge<T>& e1, const Edge<T>& e2)
{
    return !(e1 == e2);
}

/**
 * @brief out stream operator for edges if edge property type has stream operator defined
 */
template <class T>
std::enable_if_t<has_ostream_op<T>::value, std::ostream&> operator<<(std::ostream& os, const Edge<T>& e)
{
    os << e.start_node_idx << " > " << e.end_node_idx << " : " << e.property;
    return os;
}

/**
 * @brief out stream operator for edges if edge property type does not have stream operator defined
 */
template <class T>
std::enable_if_t<!has_ostream_op<T>::value, std::ostream&> operator<<(std::ostream& os, const Edge<T>& e)
{
    os << e.start_node_idx << " > " << e.end_node_idx;
    return os;
}

/**
 * @brief generic graph structure
 */
template <class NodePropertyT, class EdgePropertyT>
class Graph
    //ensure correct std::is_copy_constructible; it's not correct by default because the nodes and edges are stored in std::vector
    : not_copyable_if_t<
          !conjunction<std::is_copy_constructible<NodePropertyT>, std::is_copy_constructible<EdgePropertyT>>::value>
{
public:
    using NodeProperty = NodePropertyT;
    using EdgeProperty = EdgePropertyT;

    /**
     * @brief add a node to the graph. property of the node is constructed from arguments.
     */
    template <class... Args>
    Node<NodePropertyT>& add_node(int id, Args&&... args)
    {
        m_nodes.emplace_back(id, std::forward<Args>(args)...);
        return m_nodes.back();
    }

    /**
     * @brief add an edge to the graph. property of the edge is constructed from arguments.
     */
    template <class... Args>
    Edge<EdgePropertyT>& add_edge(size_t start_node_idx, size_t end_node_idx, Args&&... args)
    {
        assert(m_nodes.size() > start_node_idx && m_nodes.size() > end_node_idx);
        return *insert_sorted_replace(m_edges,
                                      Edge<EdgePropertyT>(start_node_idx, end_node_idx, std::forward<Args>(args)...),
                                      [](auto&& e1, auto&& e2) {
                                          return e1.start_node_idx == e2.start_node_idx
                                                     ? e1.end_node_idx < e2.end_node_idx
                                                     : e1.start_node_idx < e2.start_node_idx;
                                      });
    }

    /**
     * @brief range of nodes
     */
    auto nodes()
    {
        return make_range(begin(m_nodes), end(m_nodes));
    }

    /**
     * @brief range of nodes
     */
    auto nodes() const
    {
        return make_range(begin(m_nodes), end(m_nodes));
    };

    /**
     * @brief range of edges
     */
    auto edges()
    {
        return make_range(begin(m_edges), end(m_edges));
    }

    /**
     * @brief range of edges
     */
    auto edges() const
    {
        return make_range(begin(m_edges), end(m_edges));
    }

    /**
     * @brief range of edges going out from a specific node
     */
    auto out_edges(size_t node_idx)
    {
        return out_edges(begin(m_edges), end(m_edges), node_idx);
    }

    /**
     * @brief range of edges going out from a specific node
     */
    auto out_edges(size_t node_idx) const
    {
        return out_edges(begin(m_edges), end(m_edges), node_idx);
    }

private:
    template <typename Iter>
    static auto out_edges(Iter b, Iter e, size_t idx)
    {
        return make_range(std::equal_range(b, e, OutEdgeBase(idx), [](auto&& e1, auto&& e2) {
            return e1.start_node_idx < e2.start_node_idx;
        }));
    }

private:
    std::vector<Node<NodePropertyT>> m_nodes;
    std::vector<Edge<EdgePropertyT>> m_edges;
}; // namespace mio

/**
 * @brief Sets the graph nodes for counties or districts.
 * Reads the node ids which could refer to districts or counties and the epidemiological
 * data from json files and creates one node for each id. Every node contains a model.
 * @param[in] params Model Parameters that are used for every node.
 * @param[in] start_date Start date for which the data should be read.
 * @param[in] end_data End date for which the data should be read.
 * @param[in] data_dir Directory that contains the data files.
 * @param[in] population_data_path Path to json file containing the population data.
 * @param[in] is_node_for_county Specifies whether the node ids should be county ids (true) or district ids (false).
 * @param[in, out] params_graph Graph whose nodes are set by the function.
 * @param[in] read_func Function that reads input data for german counties and sets Model compartments.
 * @param[in] node_func Function that returns the county ids.
 * @param[in] scaling_factor_inf Factor of confirmed cases to account for undetected cases in each county.
 * @param[in] scaling_factor_icu Factor of ICU cases to account for underreporting.
 * @param[in] tnt_capacity_factor Factor for test and trace capacity.
 * @param[in] num_days Number of days to be simulated; required to load data for vaccinations during the simulation.
 * @param[in] export_time_series If true, reads data for each day of simulation and writes it in the same directory as the input files.
 */
template <class TestAndTrace, class ContactPattern, class Model, class MigrationParams, class Parameters,
          class ReadFunction, class NodeIdFunction>
IOResult<void>
set_nodes(const Parameters& params, Date start_date, Date end_date, const fs::path& data_dir,
          const std::string& population_data_path, bool is_node_for_county, Graph<Model, MigrationParams>& params_graph,
          ReadFunction&& read_func, NodeIdFunction&& node_func, const std::vector<double>& scaling_factor_inf,
          double scaling_factor_icu, double tnt_capacity_factor, int num_days = 0, bool export_time_series = false)
{
    BOOST_OUTCOME_TRY(node_ids, node_func(population_data_path, is_node_for_county));
    std::vector<Model> nodes(node_ids.size(), Model(int(size_t(params.get_num_groups()))));
    for (auto& node : nodes) {
        node.parameters = params;
    }

    BOOST_OUTCOME_TRY(read_func(nodes, start_date, node_ids, scaling_factor_inf, scaling_factor_icu, data_dir.string(),
                                num_days, export_time_series));

    for (size_t node_idx = 0; node_idx < nodes.size(); ++node_idx) {

        auto tnt_capacity = nodes[node_idx].populations.get_total() * tnt_capacity_factor;

        //local parameters
        auto& tnt_value = nodes[node_idx].parameters.template get<TestAndTrace>();
        tnt_value       = UncertainValue(0.5 * (1.2 * tnt_capacity + 0.8 * tnt_capacity));
        tnt_value.set_distribution(mio::ParameterDistributionUniform(0.8 * tnt_capacity, 1.2 * tnt_capacity));

        auto id = 0;
        if (is_node_for_county) {
            id = int(regions::CountyId(node_ids[node_idx]));
        }
        else {
            id = int(regions::DistrictId(node_ids[node_idx]));
        }
        //holiday periods
        auto holiday_periods = regions::get_holidays(regions::get_state_id(id), start_date, end_date);
        auto& contacts       = nodes[node_idx].parameters.template get<ContactPattern>();
        contacts.get_school_holidays() =
            std::vector<std::pair<mio::SimulationTime, mio::SimulationTime>>(holiday_periods.size());
        std::transform(
            holiday_periods.begin(), holiday_periods.end(), contacts.get_school_holidays().begin(), [=](auto& period) {
                return std::make_pair(mio::SimulationTime(mio::get_offset_in_days(period.first, start_date)),
                                      mio::SimulationTime(mio::get_offset_in_days(period.second, start_date)));
            });

        //uncertainty in populations
        for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
            for (auto j = Index<typename Model::Compartments>(0); j < Model::Compartments::Count; ++j) {
                auto& compartment_value = nodes[node_idx].populations[{i, j}];
                compartment_value =
                    UncertainValue(0.5 * (1.1 * double(compartment_value) + 0.9 * double(compartment_value)));
                compartment_value.set_distribution(mio::ParameterDistributionUniform(0.9 * double(compartment_value),
                                                                                     1.1 * double(compartment_value)));
            }
        }

        params_graph.add_node(node_ids[node_idx], nodes[node_idx]);
    }
    return success();
}

/**
 * @brief Sets the graph edges.
 * Reads the commuting matrices from txt files and sets the graph edges with that.
 * @param[in] data_dir Directory that contains the data files.
 * @param[in, out] params_graph Graph whose nodes are set by the function.
 * @param[in] migrating_compartments Compartments that commute.
 * @param[in] contact_locations_size Number of contact locations.
 * @param[in] read_func Function that reads commuting matrices.
 * @param[in] commuting_weights Vector with a commuting weight for every AgeGroup.
 */
template <class ContactLocation, class Model, class MigrationParams, class MigrationCoefficientGroup,
          class InfectionState, class ReadFunction>
IOResult<void> set_edges(const fs::path& data_dir, Graph<Model, MigrationParams>& params_graph,
                         std::initializer_list<InfectionState>& migrating_compartments, size_t contact_locations_size,
                         ReadFunction&& read_func,
                         std::vector<ScalarType> commuting_weights = std::vector<ScalarType>{})
{
    // mobility between nodes
    BOOST_OUTCOME_TRY(mobility_data_commuter,
                      read_func((data_dir / "mobility" / "commuter_migration_scaled.txt").string()));
    BOOST_OUTCOME_TRY(mobility_data_twitter, read_func((data_dir / "mobility" / "twitter_scaled_1252.txt").string()));
    if (mobility_data_commuter.rows() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data_commuter.cols() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data_twitter.rows() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data_twitter.cols() != Eigen::Index(params_graph.nodes().size())) {
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Mobility matrices do not have the correct size. You may need to run "
                            "transformMobilitydata.py from pycode memilio epidata package.");
    }

    for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
            auto& populations = params_graph.nodes()[county_idx_i].property.populations;
            // mobility coefficients have the same number of components as the contact matrices.
            // so that the same NPIs/dampings can be used for both (e.g. more home office => fewer commuters)
            auto mobility_coeffs = MigrationCoefficientGroup(contact_locations_size, populations.numel());
            auto num_age_groups  = (size_t)params_graph.nodes()[county_idx_i].property.parameters.get_num_groups();
            commuting_weights =
                (commuting_weights.size() == 0 ? std::vector<ScalarType>(num_age_groups, 1.0) : commuting_weights);
            //commuters
            auto working_population = 0.0;
            for (auto age = AgeGroup(0); age < populations.template size<mio::AgeGroup>(); ++age) {
                working_population += populations.get_group_total(age) * commuting_weights[size_t(age)];
            }
            auto commuter_coeff_ij = mobility_data_commuter(county_idx_i, county_idx_j) /
                                     working_population; //data is absolute numbers, we need relative
            for (auto age = AgeGroup(0); age < populations.template size<mio::AgeGroup>(); ++age) {
                for (auto compartment : migrating_compartments) {
                    auto coeff_index = populations.get_flat_index({age, compartment});
                    mobility_coeffs[size_t(ContactLocation::Work)].get_baseline()[coeff_index] =
                        commuter_coeff_ij * commuting_weights[size_t(age)];
                }
            }
            //others
            auto total_population = populations.get_total();
            auto twitter_coeff    = mobility_data_twitter(county_idx_i, county_idx_j) /
                                 total_population; //data is absolute numbers, we need relative
            for (auto age = AgeGroup(0); age < populations.template size<mio::AgeGroup>(); ++age) {
                for (auto compartment : migrating_compartments) {
                    auto coeff_idx = populations.get_flat_index({age, compartment});
                    mobility_coeffs[size_t(ContactLocation::Other)].get_baseline()[coeff_idx] = twitter_coeff;
                }
            }

            //only add edges with mobility above thresholds for performance
            //thresholds are chosen empirically so that more than 99% of mobility is covered, approx. 1/3 of the edges
            if (commuter_coeff_ij > 4e-5 || twitter_coeff > 1e-5) {
                params_graph.add_edge(county_idx_i, county_idx_j, std::move(mobility_coeffs));
            }
        }
    }

    return success();
}

/**
 * Create an unconnected graph.
 * Can be used to save space on disk when writing parameters if the edges are not required.
 * @param node_properties Vector of node properties of all nodes, e.g., parameters in each model node.
 * @param node_ids Indices for the nodes.
 * @return Graph with nodes only having no edges.
 */
template <class NodePropertyT, class EdgePropertyT>
auto create_graph_without_edges(const std::vector<NodePropertyT>& node_properties, const std::vector<int>& node_ids)
{
    // create a graph without edges for writing to file
    auto graph = mio::Graph<NodePropertyT, EdgePropertyT>();
    for (auto i = size_t(0); i < node_ids.size(); ++i) {
        graph.add_node(node_ids[i], node_properties[i]);
    }
    return graph;
}

template <class T>
std::enable_if_t<!has_ostream_op<T>::value, void> print_graph_object(std::ostream& os, size_t idx, const T&)
{
    os << idx;
}

template <class T>
std::enable_if_t<has_ostream_op<T>::value, void> print_graph_object(std::ostream& os, size_t idx, const T& o)
{
    os << idx << " [" << o << "]";
}

template <class Graph>
void print_graph(std::ostream& os, const Graph& g)
{
    auto nodes = g.nodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        os << "NODE ";
        print_graph_object(os, i, nodes[i]);
        os << '\n';
    }

    auto edges = g.edges();
    for (size_t i = 0; i < edges.size(); ++i) {
        auto& e = edges[i];
        os << "EDGE ";
        print_graph_object(os, i, e.property);
        os << " FROM NODE ";
        print_graph_object(os, e.start_node_idx, nodes[e.start_node_idx]);
        os << " TO ";
        print_graph_object(os, e.end_node_idx, nodes[e.end_node_idx]);
        os << '\n';
    }
}

} // namespace mio

#endif //GRAPH_H
