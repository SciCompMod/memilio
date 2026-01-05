#pragma once

#include <iostream>
#include <vector>
#include <cstddef>
#include <cassert>

#include "constraints.h"
#include "infection_state_utils.h"

#include "../optimization_settings/secir_optimization.h"

#include <vector>

template <typename FP, class Sim>
void update_path_constraint(
    const SecirOptimization& settings,
    mio::GraphSimulation<FP, mio::Graph<mio::SimulationNode<FP, Sim>, mio::MobilityEdge<FP>>, FP, FP,
                         void (*)(FP, FP, mio::MobilityEdge<FP>&, mio::SimulationNode<FP, Sim>&,
                                  mio::SimulationNode<FP, Sim>&),
                         void (*)(FP, FP, mio::SimulationNode<FP, Sim>&)>& graph_simulation,
    std::vector<FP>& path_constraint_values)
{
    assert(path_constraint_values.size() == settings.num_path_constraints());
    using std::max;

    const auto& graph = graph_simulation.get_graph();

    for (size_t i = 0; i < settings.num_path_constraints(); i++) {
        const Constraint& constraint  = settings.path_constraints()[i];
        auto states                   = query_infection_states(constraint.name());
        std::optional<size_t> node_id = constraint.node_id();

        FP value = 0.0;
        for (size_t node_idx = 0; node_idx < graph.nodes().size(); node_idx++) {
            const auto& node            = graph.nodes()[node_idx].property;
            const auto& node_simulation = node.get_simulation();
            const auto& node_model      = node_simulation.get_model();
            size_t graph_node_id        = graph.nodes()[node_idx].id;
            size_t num_age_groups       = node_model.parameters.get_num_groups();

            const mio::TimeSeries<FP>& simulation_result      = node_simulation.get_result();
            Eigen::Ref<const Eigen::VectorX<FP>> latest_state = simulation_result.get_last_value();

            if (!node_id.has_value() || node_id.value() == graph_node_id) {
                for (size_t age_group = 0; age_group < num_age_groups; age_group++) {
                    for (auto state : states) {
                        size_t index = age_group * num_infection_states() + static_cast<size_t>(state);
                        value += latest_state[index];
                    }
                }
            }
        }
        path_constraint_values[i] = max<FP>(path_constraint_values[i], value);
    }
}

template <typename FP, class Sim>
void update_terminal_constraint(
    const SecirOptimization& settings,
    mio::GraphSimulation<FP, mio::Graph<mio::SimulationNode<FP, Sim>, mio::MobilityEdge<FP>>, FP, FP,
                         void (*)(FP, FP, mio::MobilityEdge<FP>&, mio::SimulationNode<FP, Sim>&,
                                  mio::SimulationNode<FP, Sim>&),
                         void (*)(FP, FP, mio::SimulationNode<FP, Sim>&)>& graph_sim_mobility,
    std::vector<FP>& terminal_constraint_values)
{
    assert(terminal_constraint_values.size() == settings.num_terminal_constraints());
    using std::max;

    const auto& graph = graph_sim_mobility.get_graph();

    for (size_t i = 0; i < settings.num_terminal_constraints(); i++) {
        const Constraint& constraint  = settings.terminal_constraints()[i];
        auto states                   = query_infection_states(constraint.name());
        std::optional<size_t> node_id = constraint.node_id();

        FP value = 0.0;
        for (size_t node_idx = 0; node_idx < graph.nodes().size(); node_idx++) {
            const auto& node            = graph.nodes()[node_idx].property;
            const auto& node_simulation = node.get_simulation();
            const auto& node_model      = node_simulation.get_model();
            size_t graph_node_id        = graph.nodes()[node_idx].id;
            size_t num_age_groups       = node_model.parameters.get_num_groups();

            const mio::TimeSeries<FP>& simulation_result      = node_simulation.get_result();
            Eigen::Ref<const Eigen::VectorX<FP>> latest_state = simulation_result.get_last_value();

            if (!node_id.has_value() || node_id.value() == graph_node_id) {
                for (size_t age_group = 0; age_group < num_age_groups; age_group++) {
                    for (auto state : states) {
                        size_t index = age_group * num_infection_states() + static_cast<size_t>(state);
                        value += latest_state[index];
                    }
                }
            }
        }
        terminal_constraint_values[i] = value;
    }
}

template <typename FP, class Sim>
FP get_tests_used(const SecirOptimization& settings,
                  mio::GraphSimulation<FP, mio::Graph<mio::SimulationNode<FP, Sim>, mio::MobilityEdge<FP>>, FP, FP,
                                       void (*)(FP, FP, mio::MobilityEdge<FP>&, mio::SimulationNode<FP, Sim>&,
                                                mio::SimulationNode<FP, Sim>&),
                                       void (*)(FP, FP, mio::SimulationNode<FP, Sim>&)>& graph_sim_mobility)
{
    const auto& graph = graph_sim_mobility.get_graph();
    FP total_tests    = 0.0;

    for (const auto& node_entry : graph.nodes()) {
        const auto& node_simulation = node_entry.property.get_simulation();
        const auto& node_model      = node_simulation.get_model();
        total_tests += node_model.commuter_tests_used;
    }

    return total_tests;
}
