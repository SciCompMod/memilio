/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele
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
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"
#include "memilio/io/result_io.h"

int main()
{
    const auto t0   = 0.;
    const auto tmax = 15.;
    const auto dt   = 0.5; //time step of mobility, daily mobility every second step

    mio::oseir::Model<> model(1);

    // set population
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 10000;

    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<>>(1.);

    // set transition times
    model.parameters.set<mio::oseir::TimeExposed<>>(3.);
    model.parameters.set<mio::oseir::TimeInfected<>>(5.);

    // set contact matrix
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::oseir::ContactPatterns<>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(2.7);

    //two mostly identical groups
    auto model_group1 = model;
    auto model_group2 = model;

    //some contact restrictions in group 1
    // mio::ContactMatrixGroup& contact_matrix1 =
    //     model_group1.parameters.get<mio::oseir::ContactPatterns<>>().get_cont_freq_mat();
    // contact_matrix1[0].add_damping(0.5, mio::SimulationTime(5));

    //infection starts in group 1
    model_group1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 9990;
    model_group1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = 10;

    mio::Graph<mio::SimulationNode<mio::Simulation<ScalarType, mio::oseir::Model<>>>, mio::MobilityEdge<>> g;
    g.add_node(1001, model_group1, t0);
    g.add_node(1002, model_group2, t0);
    for (auto& node : g.nodes()) {
        node.property.get_simulation().set_integrator(std::make_shared<mio::EulerIntegratorCore<ScalarType>>());
        node.property.get_simulation().get_dt() = dt;
    }
    g.add_edge(0, 1, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, 0.05));
    g.add_edge(1, 0, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, 0.01));

    auto sim = mio::make_mobility_sim(t0, dt, std::move(g));

    sim.advance(tmax);

    auto result_graph = std::move(sim).get_graph();
    auto result       = mio::interpolate_simulation_result(result_graph);

    std::vector<int> county_ids(result_graph.nodes().size());
    std::transform(result_graph.nodes().begin(), result_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    auto save_result_status = save_result(result, county_ids, 1, "graph_result.h5");

    for (auto&& node : result_graph.nodes()) {
        node.property.get_result().print_table();
    }

    return 0;
}
