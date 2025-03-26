/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/config.h"
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/flow_simulation.h"

#include <iostream>

// Type definitions for the simulation and graph
using SimType   = mio::FlowSimulation<ScalarType, mio::oseir::Model<ScalarType>>;
using GraphType = mio::Graph<mio::SimulationNode<SimType>, mio::MobilityEdge<ScalarType>>;

// Custom mobility return function to adjust mobile population based on flow differences
void flow_based_mobility_returns(Eigen::Ref<Eigen::VectorXd> mobile_population, const SimType& sim,
                                 Eigen::Ref<const Eigen::VectorXd> total, ScalarType t, ScalarType dt)
{
    mio::unused(t, dt);
    // Retrieve flow data from the simulation
    auto flows = sim.get_flows();

    // Check if there are at least two flow data points; if not, use the standard calculation
    if (flows.get_num_time_points() < 2) {
        mio::log_warning("Flow data not available for mobility return calculation.");
        return;
    }

    // Get the last two flow values and compute the difference
    auto last_flow        = flows.get_last_value();
    auto second_last_flow = flows.get_value(flows.get_num_time_points() - 2);
    Eigen::VectorXd diff  = last_flow - second_last_flow;

    // Scale flows by the ratio of the mobile population to the total population
    double scale_factor         = mobile_population.sum() / total.sum();
    Eigen::VectorXd diff_mobile = diff * scale_factor;

    // Access the model to determine compartment indices
    const auto& model       = sim.get_model();
    const size_t num_groups = static_cast<size_t>(model.parameters.get_num_groups());
    using InfState          = mio::oseir::InfectionState;

    // Apply the flows to each age group for the SEIR compartments
    for (mio::AgeGroup group(0); group < mio::AgeGroup(num_groups); ++group) {
        size_t flow_base_idx = static_cast<size_t>(group) * 3; // Three flows per age group: S->E, E->I, I->R
        size_t S_idx         = model.populations.get_flat_index({group, InfState::Susceptible});
        size_t E_idx         = model.populations.get_flat_index({group, InfState::Exposed});
        size_t I_idx         = model.populations.get_flat_index({group, InfState::Infected});
        size_t R_idx         = model.populations.get_flat_index({group, InfState::Recovered});

        if (flow_base_idx + 2 < static_cast<size_t>(diff_mobile.size())) {
            // Flow S -> E
            mobile_population[S_idx] -= diff_mobile[flow_base_idx];
            mobile_population[E_idx] += diff_mobile[flow_base_idx];
            // Flow E -> I
            mobile_population[E_idx] -= diff_mobile[flow_base_idx + 1];
            mobile_population[I_idx] += diff_mobile[flow_base_idx + 1];
            // Flow I -> R
            mobile_population[I_idx] -= diff_mobile[flow_base_idx + 2];
            mobile_population[R_idx] += diff_mobile[flow_base_idx + 2];
        }
    }
}

/**
 * Integrates the dynamics of a mobile population during their stay at a destination node.
 */
void integrate_mobile_population_euler(Eigen::Ref<Eigen::VectorXd> mobile_population, const SimType& sim,
                                       Eigen::Ref<const Eigen::VectorXd> total, ScalarType t, ScalarType dt)
{
    auto y0 = mobile_population.eval();
    auto y1 = mobile_population;
    mio::EulerIntegratorCore<ScalarType>().step(
        [&](auto&& y, auto&& t_, auto&& dydt) {
            sim.get_model().get_derivatives(total, y, t_, dydt);
        },
        y0, t, dt, y1);
}

// Function to set up the base SEIR model with initial populations and parameters
mio::oseir::Model<ScalarType> setupBaseModel(const ScalarType total_population, const size_t num_groups)
{
    mio::oseir::Model<ScalarType> model(num_groups);

    // Set initial population distribution
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}];

    // Set model parameters
    model.parameters.set<mio::oseir::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<ScalarType>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.1);

    // Configure baseline contact patterns
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(2.7);

    return model;
}

// Function to set up node-specific models using modifications to the base model
void setupNodeModels(mio::oseir::Model<ScalarType>& model_group1, mio::oseir::Model<ScalarType>& model_group2,
                     mio::oseir::Model<ScalarType>& model_group3, const mio::oseir::Model<ScalarType>& base_model)
{
    model_group1 = base_model;
    // mio::ContactMatrixGroup& contact_matrix_m1 = model_group1.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
    // contact_matrix_m1[0].add_damping(0.7, mio::SimulationTime(15.));
    model_group1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 9700;
    model_group1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = 100;

    model_group2 = base_model;

    model_group3                                                                          = base_model;
    model_group3.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 9795;
    model_group3.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = 5;
}

// Function to set up indices for tracking the SEIR compartments (Susceptible, Exposed, Infected, Recovered)
std::vector<std::vector<size_t>> setupIndices(const mio::oseir::Model<ScalarType>& model, const size_t num_groups)
{
    std::vector<std::vector<size_t>> indices(4);
    for (auto& vec : indices) {
        vec.reserve(num_groups);
    }
    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(num_groups); ++group) {
        indices[0].emplace_back(model.populations.get_flat_index({group, mio::oseir::InfectionState::Susceptible}));
        indices[1].emplace_back(model.populations.get_flat_index({group, mio::oseir::InfectionState::Exposed}));
        indices[2].emplace_back(model.populations.get_flat_index({group, mio::oseir::InfectionState::Infected}));
        indices[3].emplace_back(model.populations.get_flat_index({group, mio::oseir::InfectionState::Recovered}));
    }
    return indices;
}

// Function to create and configure the graph with nodes and mobility edges
GraphType createGraph(const mio::oseir::Model<ScalarType>& model_group1,
                      const mio::oseir::Model<ScalarType>& model_group2,
                      const mio::oseir::Model<ScalarType>& model_group3,
                      const std::vector<std::vector<size_t>>& indices_save_edges, const ScalarType t0,
                      bool costum_mobility)
{
    GraphType graph;
    // Add nodes to the graph
    graph.add_node(0, model_group1, t0);
    graph.add_node(1, model_group2, t0);
    graph.add_node(2, model_group3, t0);

    // Define mobility rates for different connections
    const ScalarType mobility_rate_1_2 = 0.1;
    const ScalarType mobility_rate_1_3 = 0.1;
    const ScalarType mobility_rate_2_3 = 0.1;

    // Add edges representing population movement between nodes
    graph.add_edge(0, 1, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, mobility_rate_1_2),
                   indices_save_edges);
    graph.add_edge(1, 0, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, mobility_rate_1_2),
                   indices_save_edges);

    graph.add_edge(0, 2, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, mobility_rate_1_3),
                   indices_save_edges);
    graph.add_edge(2, 0, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, mobility_rate_1_3),
                   indices_save_edges);

    graph.add_edge(1, 2, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, mobility_rate_2_3),
                   indices_save_edges);
    graph.add_edge(2, 1, Eigen::VectorXd::Constant((size_t)mio::oseir::InfectionState::Count, mobility_rate_2_3),
                   indices_save_edges);

    // Set the custom mobility return function for each edge
    if (costum_mobility) {
        for (auto& edge : graph.edges()) {
            edge.property.set_custom_mobility_return<SimType>(integrate_mobile_population_euler);
        }
    }
    return graph;
}

// Function to process simulation results and print flow differences and summary statistics
template <typename GraphSim>
void processSimulationResults(const GraphSim& sim)
{
    std::vector<double> time_points;
    std::vector<std::vector<double>> all_edges_values(4); // For S, E, I, and R

    // Collect time points and combine data from all edges
    for (size_t edge_idx = 0; edge_idx < sim.get_graph().edges().size(); ++edge_idx) {
        auto& edge         = sim.get_graph().edges()[edge_idx];
        auto& edge_results = edge.property.get_mobility_results();

        // Initialize time points and data arrays on the first edge
        if (edge_idx == 0) {
            for (size_t i = 0; i < static_cast<size_t>(edge_results.get_num_time_points()); ++i) {
                time_points.push_back(edge_results.get_time(i));
            }
            for (auto& vec : all_edges_values) {
                vec.resize(time_points.size(), 0.0);
            }
        }
        // Sum the compartment values from this edge into the combined arrays
        for (size_t i = 0; i < static_cast<size_t>(edge_results.get_num_time_points()); ++i) {
            all_edges_values[0][i] += edge_results[i][0]; // Susceptible
            all_edges_values[1][i] += edge_results[i][1]; // Exposed
            all_edges_values[2][i] += edge_results[i][2]; // Infected
            all_edges_values[3][i] += edge_results[i][3]; // Recovered
        }
    }

    // Calculate differences between consecutive time steps and print the results
    std::vector<std::vector<double>> all_edges_diff(4);
    size_t num_periods = time_points.size() / 2;
    for (auto& vec : all_edges_diff) {
        vec.resize(num_periods, 0.0);
    }

    std::cout << "Time Period      S Change         E Change         I Change         R Change" << std::endl;
    for (size_t i = 1; i < time_points.size(); i += 2) {
        double diff_S = all_edges_values[0][i] - all_edges_values[0][i - 1];
        double diff_E = all_edges_values[1][i] - all_edges_values[1][i - 1];
        double diff_I = all_edges_values[2][i] - all_edges_values[2][i - 1];
        double diff_R = all_edges_values[3][i] - all_edges_values[3][i - 1];

        size_t idx             = i / 2;
        all_edges_diff[0][idx] = diff_S;
        all_edges_diff[1][idx] = diff_E;
        all_edges_diff[2][idx] = diff_I;
        all_edges_diff[3][idx] = diff_R;

        std::string time_period = std::to_string(time_points[i - 1]) + " -> " + std::to_string(time_points[i]);
        printf("%-15s %15.5f %15.5f %15.5f %15.5f\n", time_period.c_str(), diff_S, diff_E, diff_I, diff_R);
    }

    double sum_diff_S = 0.0, sum_diff_E = 0.0, sum_diff_I = 0.0, sum_diff_R = 0.0;
    for (size_t i = 0; i < all_edges_diff[0].size(); ++i) {
        sum_diff_S += all_edges_diff[0][i];
        sum_diff_E += all_edges_diff[1][i];
        sum_diff_I += all_edges_diff[2][i];
        sum_diff_R += all_edges_diff[3][i];
    }

    std::cout << "\n--- Summary of Flow Changes Over Entire Simulation ---" << std::endl;
    printf("Net S Change:     %15.5f\n", sum_diff_S);
    printf("Net E Change:     %15.5f\n", sum_diff_E);
    printf("Net I Change:     %15.5f\n", sum_diff_I);
    printf("Net R Change:     %15.5f\n", sum_diff_R);
}

int main()
{
    // --- Simulation setup ---
    mio::set_log_level(mio::LogLevel::warn);
    const ScalarType t0        = 0.0;
    const ScalarType tmax      = 40.0;
    const ScalarType dt        = 0.5; // Time step for mobility simulation
    const bool costum_mobility = false;

    // --- Base model configuration ---
    const size_t num_groups           = 1;
    const ScalarType total_population = 10000;
    auto base_model                   = setupBaseModel(total_population, num_groups);

    // --- Set up node-specific models ---
    mio::oseir::Model<ScalarType> model_group1(num_groups), model_group2(num_groups), model_group3(num_groups);
    setupNodeModels(model_group1, model_group2, model_group3, base_model);

    // --- Set up indices for tracking compartments on edges ---
    auto indices_save_edges = setupIndices(base_model, num_groups);

    // --- Create graph with nodes and mobility edges ---
    auto graph = createGraph(model_group1, model_group2, model_group3, indices_save_edges, t0, costum_mobility);

    // --- Run simulation ---
    auto sim = mio::make_mobility_sim(t0, dt, std::move(graph));
    sim.advance(tmax);

    // --- Process and display simulation results ---
    processSimulationResults(sim);

    return 0;
}
