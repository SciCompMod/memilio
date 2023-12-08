
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_detailed.h"
#include "memilio/compartments/simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "ode_secir/model.h"
#include "ode_secir/parameters.h"

namespace mio
{
mio::IOResult<Eigen::MatrixXd> read_duration_stay(const std::string& filename)
{
    mio::unused(filename);
    return mio::success(Eigen::MatrixXd(0, 0));
}

TEST(TestMobilityDetailed, AddNodeTest)
{
    mio::osecir::Model model(1);
    GraphDetailed<mio::SimulationNode<mio::Simulation<mio::osecir::Model>>, int> sim_graph;

    auto& node = sim_graph.add_node(1, 0.5, model, model, 0.0, 0.1);
    EXPECT_EQ(node.stay_duration, 0.5);
}

TEST(TestMobilityDetailed, AddNodeWithTimeIntegrationTest)
{
    mio::osecir::Model model(1);
    GraphDetailed<mio::osecir::Model, int> model_graph;
    auto& node = model_graph.add_node(1, 0.5, model, model);
    EXPECT_EQ(node.stay_duration, 0.5);
}

TEST(TestMobilityDetailed, AddEdgeTest)
{
    mio::osecir::Model model(1);
    GraphDetailed<mio::osecir::Model, int> model_graph;
    model_graph.add_node(0, 0.5, model, model);
    model_graph.add_node(1, 0.5, model, model);
    auto& edge = model_graph.add_edge(0, 1, 0.5, {1, 2});
    EXPECT_EQ(edge.traveltime, 0.5);
}

} // namespace mio