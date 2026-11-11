/* 
* Copyright (C) 2020-2026 MEmilio
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
#include "memilio/data/analyze_result.h"
#include "ode_secir/model.h"
#include "memilio/compartments/feedback_simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/utils/logging.h"
#include <iostream>

// alias for the type of the simulation with feedback
using FeedbackSim = mio::FeedbackSimulation<ScalarType, mio::Simulation<ScalarType, mio::osecir::Model<ScalarType>>,
                                            mio::osecir::ContactPatterns<ScalarType>>;

// helper function to initialize the model with population and parameters
void initialize_model(mio::osecir::Model<ScalarType>& model, ScalarType cont_freq)
{
    model.parameters.set<mio::osecir::StartDay<ScalarType>>(60);
    model.parameters.set<mio::osecir::Seasonality<ScalarType>>(0.2);

    // Mean stay times per compartment
    model.parameters.get<mio::osecir::TimeExposed<ScalarType>>()            = 3.2;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>() = 2.0;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()   = 5.8;
    model.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()     = 9.5;
    model.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()   = 7.1;

    // Set transmission and isolation parameters
    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()  = 0.05;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()    = 0.7;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()    = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()    = 0.25;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<ScalarType>>() = 0.45;
    model.parameters.get<mio::osecir::TestAndTraceCapacity<ScalarType>>()              = 35;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()         = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()                 = 0.25;
    model.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()                 = 0.3;

    // contact matrix
    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, cont_freq));
}

// helper function to initialize the feedback mechanism parameters for a simulation
void initialize_feedback(FeedbackSim& feedback_simulation)
{
    // nominal ICU capacity
    feedback_simulation.get_parameters().template get<mio::NominalICUCapacity<ScalarType>>() = 10;

    // ICU occupancy in the past for memory kernel
    auto& icu_occupancy     = feedback_simulation.get_parameters().template get<mio::ICUOccupancyHistory<ScalarType>>();
    Eigen::VectorXd icu_day = Eigen::VectorXd::Constant(1, 1);
    const auto cutoff       = static_cast<int>(feedback_simulation.get_parameters().template get<mio::GammaCutOff>());
    for (int t = -cutoff; t <= 0; ++t) {
        icu_occupancy.add_time_point(t, icu_day);
    }

    // bounds for contact reduction measures
    feedback_simulation.get_parameters().template get<mio::ContactReductionMin<ScalarType>>() = {0.1};
    feedback_simulation.get_parameters().template get<mio::ContactReductionMax<ScalarType>>() = {0.8};

    // Set blending factors. The global blending factor is implicitly defined as 1 - local - regional.
    feedback_simulation.get_parameters().template get<mio::BlendingFactorLocal<ScalarType>>()    = 0.5;
    feedback_simulation.get_parameters().template get<mio::BlendingFactorRegional<ScalarType>>() = 0.3;
}

// helper function to create the graph with nodes and edges
mio::Graph<mio::SimulationNode<ScalarType, FeedbackSim>, mio::MobilityEdge<ScalarType>>
create_graph(int num_nodes, int total_population, ScalarType cont_freq)
{
    // Create a graph for the metapopulation simulation
    mio::Graph<mio::SimulationNode<ScalarType, FeedbackSim>, mio::MobilityEdge<ScalarType>> g;

    // Create models and add nodes to the graph
    for (int i = 0; i < num_nodes; ++i) {
        mio::osecir::Model<ScalarType> model(1);
        initialize_model(model, cont_freq);

        // Set initial populations (infection starts in the first node)
        if (i == 0) {
            model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}] = total_population * 0.1;
            model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] =
                total_population * 0.1;
            model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] =
                total_population * 0.05;
            model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}] =
                total_population * 0.02;
            model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}] =
                total_population * 0.01;
            model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}] = 0;
            model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                        total_population);
        }
        else {
            model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] = total_population;
        }
        // The function apply_constraints() ensures that all parameters are within their defined bounds.
        // Note that negative values are set to zero instead of stopping the simulation.
        model.apply_constraints();

        // Determine the index for the ICU state (InfectedCritical) for the feedback mechanism
        auto icu_index = std::vector<size_t>{
            model.populations.get_flat_index({mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical})};

        // Create feedback simulation
        auto feedback_sim = FeedbackSim(mio::Simulation<ScalarType, mio::osecir::Model<ScalarType>>(model), icu_index);
        initialize_feedback(feedback_sim);

        // Node-ID-Logic: 1001-1005, 2001-2005, ...
        const int region_id = i / 5;
        const int local_id  = i % 5;
        const int node_id   = (region_id + 1) * 1000 + (local_id + 1);
        g.add_node(node_id, std::move(feedback_sim));
    }

    // Define complete graph, i.e. each node is connected to every other node
    std::vector<std::vector<size_t>> mobile_compartments(2);
    for (size_t i = 0; i < g.nodes().size(); ++i) {
        for (size_t j = 0; j < g.nodes().size(); ++j) {
            if (i != j) {
                g.add_edge(i, j, Eigen::VectorXd::Constant((size_t)mio::osecir::InfectionState::Count, 0.01));
            }
        }
    }

    return g;
}

int main()
{
    // This example demonstrates the implementation of a feedback mechanism for a ODE SECIR model in a graph.
    // It shows how the perceived risk dynamically impacts contact reduction measures in different regions (nodes).
    mio::set_log_level(mio::LogLevel::err);

    const auto t0              = 0.;
    const auto tmax            = 10.;
    const auto dt              = 0.5;
    const int total_population = 1000;
    const ScalarType cont_freq = 2.7;
    const int num_nodes        = 10;

    // Create the graph
    auto g = create_graph(num_nodes, total_population, cont_freq);

    // Create and run the simulation
    using Graph = decltype(g);
    auto sim    = mio::FeedbackGraphSimulation<ScalarType, Graph>(std::move(g), t0, dt);
    sim.advance(tmax);

    // The output shows the compartments sizes for a node without any initial infections.
    auto& node         = sim.get_graph().nodes()[1];
    auto& results_node = node.property.get_simulation().get_result();
    // interpolate results
    auto interpolated_results_node = mio::interpolate_simulation_result(results_node);

    // print result with print_table
    std::cout << "Node ID: " << node.id << "\n";
    std::vector<std::string> cols = {"S", "E", "C", "C_confirmed", "I", "I_confirmed", "H", "U", "R", "D"};
    interpolated_results_node.print_table(cols);

    return 0;
}
