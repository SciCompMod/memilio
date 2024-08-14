/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Henrik Zunker
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
#include "memilio/mobility/metapopulation_mobility_detailed.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/data/analyze_result.h"

int main()
{
    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 0.5; //time step of migration, daily migration every second step

    mio::oseir::Model<> model(1);

    // set population
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 10000;

    // set transition times
    model.parameters.set<mio::oseir::TimeExposed<>>(1);
    model.parameters.set<mio::oseir::TimeInfected<>>(1);

    // set contact matrix
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::oseir::ContactPatterns<>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(2.7);

    //two mostly identical groups
    auto model_group1 = model;
    auto model_group2 = model;

    //some contact restrictions in group 1
    mio::ContactMatrixGroup& contact_matrix1 =
        model_group1.parameters.get<mio::oseir::ContactPatterns<>>().get_cont_freq_mat();
    contact_matrix1[0].add_damping(0.5, mio::SimulationTime(5));

    //infection starts in group 1
    model_group1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 9990;
    model_group1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = 10;

    auto sim1          = mio::FlowSimulation<ScalarType, mio::oseir::Model<>>(model_group1, t0, dt);
    auto sim2          = mio::FlowSimulation<ScalarType, mio::oseir::Model<>>(model_group2, t0, dt);
    double stay_time_1 = 0.3;
    double stay_time_2 = 0.3;
    // mio::ExtendedGraph<mio::Simulation<ScalarType, mio::oseir::Model<>>> g;
    mio::ExtendedGraph<mio::FlowSimulation<ScalarType, mio::oseir::Model<double>>> g;
    g.add_node(1, sim1, sim2, stay_time_1);
    g.add_node(2, sim1, sim2, stay_time_2);

    double traveltime      = 0.1;
    std::vector<int> path1 = {0, 1};
    std::vector<int> path2 = {1, 0};
    g.add_edge(0, 1, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, 0.01), traveltime, path1);
    g.add_edge(1, 0, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, 0.01), traveltime, path2);

    auto sim = mio::make_migration_sim(t0, dt, std::move(g));

    sim.advance(tmax);
    // results node 1
    std::cout << "Results node 1" << std::endl;
    auto interpolated_sim1 =
        mio::interpolate_simulation_result(sim.get_graph().nodes()[0].property.base_sim.get_result());
    interpolated_sim1.print_table({"S", "E", "I", "R"});

    // results node 1 mobility_sim
    std::cout << "Mobility results node 1" << std::endl;
    auto interpolated_sim1_mobility = sim.get_graph().nodes()[0].property.mobility_sim.get_result();
    interpolated_sim1_mobility.print_table({"S", "E", "I", "R"});

    // results node 2
    std::cout << "Results node 2" << std::endl;
    auto interpolated_sim2 = sim.get_graph().nodes()[1].property.base_sim.get_result();
    interpolated_sim2.print_table({"S", "E", "I", "R"});

    // results node 2 mobility_sim
    std::cout << "Mobility results node 2" << std::endl;
    auto interpolated_sim2_mobility = sim.get_graph().nodes()[1].property.mobility_sim.get_result();
    interpolated_sim2_mobility.print_table({"S", "E", "I", "R"});

    return 0;
}
