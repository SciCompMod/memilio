
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"
#include "memilio/io/result_io.h"

#include <chrono>

mio::IOResult<void> run(const fs::path& data_dir, double t0, double tmax, double dt)
{
    // global parameters
    const int num_age_groups = 1;
    mio::oseir::Parameters<double> params(num_age_groups);
    mio::Populations<double, mio::AgeGroup, mio::oseir::InfectionState> population(
        {mio::AgeGroup(num_age_groups), mio::oseir::InfectionState::Count});
    params.set<mio::oseir::TransmissionProbabilityOnContact<>>(1.);

    // set transition times
    params.set<mio::oseir::TimeExposed<>>(3.);
    params.set<mio::oseir::TimeInfected<>>(5.);

    // set contact matrix
    params.get<mio::oseir::ContactPatterns<>>().get_cont_freq_mat()[0].get_baseline().setConstant(7.95);

    population[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 10000;

    // graph of counties with populations and local parameters
    // and mobility between counties
    mio::Graph<mio::SimulationNode<mio::Simulation<ScalarType, mio::oseir::Model<>>>, mio::MobilityEdge<>> params_graph;

    std::vector<mio::oseir::Model<double>> nodes(2, mio::oseir::Model(int(size_t(params.get_num_groups()))));
    for (auto& node : nodes) {
        node.parameters  = params;
        node.populations = population;
    }
    nodes[0].populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 9990;
    nodes[0].populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = 10;

    for (size_t node_idx = 0; node_idx < nodes.size(); ++node_idx) {
        params_graph.add_node(node_idx, nodes[node_idx]);
    }

    BOOST_OUTCOME_TRY(auto&& mobility_data_commuter,
                      mio::read_mobility_plain((data_dir / "mobility" / "commuter_migration_test.txt").string()));
    if (mobility_data_commuter.rows() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data_commuter.cols() != Eigen::Index(params_graph.nodes().size())) {
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Mobility matrices do not have the correct size. You may need to run "
                            "transformMobilitydata.py from pycode memilio epidata package.");
    }

    for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
            auto& populations = params_graph.nodes()[county_idx_i].property.get_simulation().get_model().populations;

            auto commuter_coeff_ij = mobility_data_commuter(county_idx_i, county_idx_j) / populations.get_total();
            params_graph.add_edge(
                county_idx_i, county_idx_j,
                Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, commuter_coeff_ij));
        }
    }

    for (auto& node : params_graph.nodes()) {
        node.property.get_simulation().set_integrator(std::make_shared<mio::EulerIntegratorCore<ScalarType>>());
    }

    auto sim = mio::make_mobility_sim(t0, dt, std::move(params_graph));

    printf("Start Simulation\n");
    auto t1 = std::chrono::high_resolution_clock::now();
    sim.advance(tmax);
    auto t2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> ms_double = t2 - t1;

    printf("Runtime: %f\n", ms_double.count());

    auto result_graph = std::move(sim).get_graph();
    auto result       = mio::interpolate_simulation_result(result_graph);

    std::vector<int> county_ids(result_graph.nodes().size());
    std::transform(result_graph.nodes().begin(), result_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    auto save_result_status = save_result(result, county_ids, 1, "graph_result.h5");
    result_graph.nodes()[0].property.get_result().print_table();
    result_graph.nodes()[1].property.get_result().print_table();
    // for (auto&& node : result_graph.nodes()) {
    //     node.property.get_result().print_table();
    // }

    return mio::success();
}

int main()
{
    const auto t0   = 0.;
    const auto tmax = 1.;
    const auto dt   = 0.5; //time step of mobility, daily mobility every second step

    const std::string& data_dir = "";

    auto result = run(data_dir, t0, tmax, dt);

    return 0;
}
