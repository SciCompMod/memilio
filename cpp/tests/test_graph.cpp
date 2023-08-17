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
#include "memilio/mobility/graph.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/date.h"
#include "memilio/epidemiology/damping.h"
#include "models/ode_secir/parameters_io.h"
#include "models/ode_secir/parameters.h"
#include "models/ode_secir/infection_state.h"
#include "models/ode_secir/model.h"
#include "models/ode_secirvvs/parameters_io.h"
#include "models/ode_secirvvs/parameter_space.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "models/ode_secirvvs/model.h"
#include "memilio/io/io.h"
#include "matchers.h"
#include "memilio/utils/stl_util.h"
#include "gmock/gmock-matchers.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <type_traits>
#include <string>

namespace fs = boost::filesystem;

enum class MockContactLocation
{
    Work,
    Other,
    Count,
};

mio::IOResult<std::vector<int>> mock_node_function(const std::string& path, int node_id)
{
    mio::unused(path);
    mio::unused(node_id);
    std::vector<int> id = {1001, 1002};
    return mio::success(id);
}

template <class Model>
mio::IOResult<void> mock_read_function(std::vector<Model>& model, mio::Date date, const std::vector<int>& node_ids,
                                       const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                       const std::string& dir, int num_days, bool export_time_series)
{
    mio::unused(model);
    mio::unused(date);
    mio::unused(node_ids);
    mio::unused(scaling_factor_inf);
    mio::unused(scaling_factor_icu);
    mio::unused(dir);
    mio::unused(num_days);
    mio::unused(export_time_series);

    return mio::success();
}

template <class Model>
mio::IOResult<void> mock_read_function_county(std::vector<Model>& model, mio::Date date,
                                              const std::vector<int>& node_ids,
                                              const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                                              const std::string& dir, int num_days, bool export_time_series)
{
    mio::unused(model);
    mio::unused(date);
    mio::unused(node_ids);
    mio::unused(scaling_factor_inf);
    mio::unused(scaling_factor_icu);
    mio::unused(dir);
    mio::unused(num_days);
    mio::unused(export_time_series);

    return mio::success();
}

mio::IOResult<Eigen::MatrixXd> mock_read_mobility(const std::string& filename)
{
    mio::unused(filename);
    Eigen::MatrixXd migration(2, 2);
    migration(0, 0) = 0;
    migration(1, 0) = 2;
    migration(0, 1) = 2;
    migration(1, 1) = 0;

    return mio::success(migration);
}

TEST(TestGraph, creation)
{
    mio::Graph<int, int> g;
    g.add_node(0, 6);
    g.add_node(1, 4);
    g.add_node(2, 8);
    g.add_edge(0, 1, 1);
    g.add_edge(2, 1, 3);
    g.add_edge(1, 2, 2);

    std::vector<mio::Node<int>> n = {{0, 6}, {1, 4}, {2, 8}};
    EXPECT_THAT(g.nodes(), testing::ElementsAreArray(n));

    std::vector<mio::Edge<int>> v = {{0, 1, 1}, {1, 2, 2}, {2, 1, 3}};
    EXPECT_THAT(g.edges(), testing::ElementsAreArray(v));
}

TEST(TestGraph, duplicate_edge)
{
    mio::Graph<int, int> g;
    g.add_node(6);
    g.add_node(4);
    g.add_node(8);
    g.add_edge(0, 1, 1);
    g.add_edge(2, 1, 3);
    g.add_edge(1, 2, 2);
    g.add_edge(1, 2, 3);

    EXPECT_EQ(g.edges()[1], (mio::Edge<int>{1, 2, 3}));
}

TEST(TestGraph, graph_without_edges)
{
    struct MockModel {
    };
    struct MockMobility {
    };
    std::vector<MockModel> models = {MockModel(), MockModel()};
    std::vector<int> ids          = {1, 2};

    auto g = mio::create_graph_without_edges<MockModel, MockMobility>(models, ids);

    EXPECT_EQ(g.edges().size(), 0);
    EXPECT_EQ(g.nodes().size(), 2);
    EXPECT_EQ(g.nodes()[0].id, 1);
    EXPECT_EQ(g.nodes()[1].id, 2);
    // Node property is the model itself
    auto model_type_true = std::is_same<decltype(g.nodes()[0].property), MockModel>::value;
    EXPECT_TRUE(model_type_true);
    auto model_type_false = std::is_same<decltype(g.nodes()[0].property), MockMobility>::value;
    EXPECT_FALSE(model_type_false);
}

TEST(TestGraph, set_nodes_secir)
{

    mio::osecir::Parameters params(1);
    mio::Graph<mio::osecir::Model, mio::MigrationParameters> params_graph;
    const auto& read_function_nodes = mock_read_function<mio::osecir::Model>;
    const auto& node_id_function    = mock_node_function;

    const fs::path& dir = " ";

    auto result = mio::set_nodes<mio::osecir::TestAndTraceCapacity, mio::osecir::ContactPatterns, mio::osecir::Model,
                                 mio::MigrationParameters, mio::osecir::Parameters, decltype(read_function_nodes),
                                 decltype(node_id_function)>(
        params, mio::Date(2020, 5, 10), mio::Date(2020, 5, 11), dir, " ", 0, params_graph, read_function_nodes,
        node_id_function, std::vector<double>(size_t(1), 1.0), 1.0, 0.01);

    EXPECT_EQ(params_graph.nodes().size(), 2);
    EXPECT_EQ(params_graph.nodes()[0].id, 1001);
    EXPECT_EQ(params_graph.nodes()[1].id, 1002);
    auto model_type_true1 = std::is_same<decltype(params_graph.nodes()[0].property), mio::osecir::Model>::value;
    EXPECT_TRUE(model_type_true1);
    auto model_type_true2 = std::is_same<decltype(params_graph.nodes()[1].property), mio::osecir::Model>::value;
    EXPECT_TRUE(model_type_true2);
}

TEST(TestGraph, set_nodes_secirvvs)
{

    mio::osecirvvs::Parameters params(1);
    mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters> params_graph;
    const auto& read_function_nodes = mock_read_function<mio::osecirvvs::Model>;
    const auto& node_id_function    = mock_node_function;

    const fs::path& dir = " ";

    auto result = mio::set_nodes<mio::osecirvvs::TestAndTraceCapacity, mio::osecirvvs::ContactPatterns,
                                 mio::osecirvvs::Model, mio::MigrationParameters, mio::osecirvvs::Parameters,
                                 decltype(read_function_nodes), decltype(node_id_function)>(
        params, mio::Date(2020, 5, 10), mio::Date(2020, 5, 11), dir, " ", 0, params_graph, read_function_nodes,
        node_id_function, std::vector<double>(size_t(1), 1.0), 1.0, 0.01);

    EXPECT_EQ(params_graph.nodes().size(), 2);
    EXPECT_EQ(params_graph.nodes()[0].id, 1001);
    EXPECT_EQ(params_graph.nodes()[1].id, 1002);
    auto model_type_true1 = std::is_same<decltype(params_graph.nodes()[0].property), mio::osecirvvs::Model>::value;
    EXPECT_TRUE(model_type_true1);
    auto model_type_true2 = std::is_same<decltype(params_graph.nodes()[1].property), mio::osecirvvs::Model>::value;
    EXPECT_TRUE(model_type_true2);
}

TEST(TestGraph, set_edges)
{
    mio::osecir::Model model(6);
    model.populations[{mio::AgeGroup(3), mio::osecir::InfectionState::Exposed}] = 1;
    mio::Graph<mio::osecir::Model, mio::MigrationParameters> params_graph;
    const fs::path& dir         = " ";
    auto migrating_compartments = {mio::osecir::InfectionState::Susceptible, mio::osecir::InfectionState::Exposed,
                                   mio::osecir::InfectionState::InfectedNoSymptoms,
                                   mio::osecir::InfectionState::InfectedSymptoms,
                                   mio::osecir::InfectionState::Recovered};

    params_graph.add_node(1001, model);
    params_graph.add_node(1002, model);

    const auto& read_function_edges = mock_read_mobility;

    auto result =
        mio::set_edges<MockContactLocation, mio::osecir::Model, mio::MigrationParameters,
                       mio::MigrationCoefficientGroup, mio::osecir::InfectionState, decltype(read_function_edges)>(
            dir, params_graph, migrating_compartments, size_t(2), read_function_edges,
            std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.});

    auto e_work = (Eigen::ArrayXd(6 * Eigen::Index(mio::osecir::InfectionState::Count)) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 0, 0.66, 0.66,
                   0.66, 0, 0.66, 0, 0, 0, 0.66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
                      .finished();
    auto e_other = (Eigen::ArrayXd(6 * Eigen::Index(mio::osecir::InfectionState::Count)) << 2, 2, 2, 0, 2, 0, 0, 0, 2,
                    0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 0, 2, 2,
                    2, 0, 2, 0, 0, 0, 2, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 0)
                       .finished();

    EXPECT_EQ(params_graph.edges().size(), 2);

    ASSERT_THAT(print_wrap(params_graph.edges()[0]
                               .property.get_coefficients()[size_t(MockContactLocation::Work)]
                               .get_baseline()
                               .array()
                               .cast<double>()),
                MatrixNear(print_wrap(e_work.array().cast<double>()), 1e-5, 1e-5));

    ASSERT_THAT(print_wrap(params_graph.edges()[0]
                               .property.get_coefficients()[size_t(MockContactLocation::Other)]
                               .get_baseline()
                               .array()
                               .cast<double>()),
                MatrixNear(print_wrap(e_other.array().cast<double>()), 1e-5, 1e-5));

    ASSERT_THAT(print_wrap(params_graph.edges()[1]
                               .property.get_coefficients()[size_t(MockContactLocation::Work)]
                               .get_baseline()
                               .array()
                               .cast<double>()),
                MatrixNear(print_wrap(e_work.array().cast<double>()), 1e-5, 1e-5));

    ASSERT_THAT(print_wrap(params_graph.edges()[1]
                               .property.get_coefficients()[size_t(MockContactLocation::Other)]
                               .get_baseline()
                               .array()
                               .cast<double>()),
                MatrixNear(print_wrap(e_other.array().cast<double>()), 1e-5, 1e-5));
}

TEST(TestGraph, ot_edges)
{
    mio::Graph<int, int> g;
    g.add_node(0);
    g.add_node(1);
    g.add_node(2);
    g.add_node(3);
    g.add_edge(0, 1, 0);
    g.add_edge(1, 2, 1);
    g.add_edge(0, 2, 2);
    g.add_edge(3, 0, 3);

    std::vector<mio::Edge<int>> v0 = {{0, 1, 0}, {0, 2, 2}};
    EXPECT_THAT(g.out_edges(0), testing::ElementsAreArray(v0));

    std::vector<mio::Edge<int>> v1 = {{1, 2, 1}};
    EXPECT_THAT(g.out_edges(1), testing::ElementsAreArray(v1));
}

namespace
{

struct MoveOnly {
    MoveOnly();
    MoveOnly(const MoveOnly&)            = delete;
    MoveOnly& operator=(const MoveOnly&) = delete;
    MoveOnly(MoveOnly&&)                 = default;
    MoveOnly& operator=(MoveOnly&&)      = default;
};
using MoveOnlyGraph = mio::Graph<MoveOnly, MoveOnly>;

template <class G>
using add_node_expr_t = decltype(std::declval<G>().add_node(int()));
template <class G>
using add_edge_expr_t = decltype(std::declval<G>().add_edge(size_t(), size_t()));

} // namespace

static_assert(std::is_constructible<MoveOnlyGraph>::value, "Graph should support move-only node and edge properties.");
static_assert(std::is_move_constructible<MoveOnlyGraph>::value && std::is_move_assignable<MoveOnlyGraph>::value,
              "Graph should support move-only node and edge properties.");
static_assert(mio::is_expression_valid<add_node_expr_t, MoveOnlyGraph>::value,
              "Graph should support move-only node and edge properties.");
static_assert(mio::is_expression_valid<add_edge_expr_t, MoveOnlyGraph>::value,
              "Graph should support move-only node and edge properties.");
