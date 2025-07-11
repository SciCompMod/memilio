#pragma once

#include <vector>

#include "../optimization_settings/secirvvs_optimization.h"
#include "models/ode_secirvvs/model.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"

#include "../constraints/infection_state_utils.h"

template <typename FP>
FP objective_function(
    mio::Graph<mio::SimulationNode<FP, mio::Simulation<FP, mio::osecirvvs::Model<FP>>>, mio::MobilityEdge<FP>>
        graph_model,
    const SecirvvsOptimization& settings, const FP* ptr_parameters, size_t n)
{
    // ------------------------------------------------------------------ //
    // Evaluate the objective function of the model.                      //
    // Step 1. Define dampings based on 'const FP* parameters'.           //
    // Step 2. Evaluate the objective function based on                   //
    //         the parameters and the infection states in the simulation. //
    // ------------------------------------------------------------------ //
    assert(n == settings.num_control_parameters() * settings.num_control_intervals());
    const size_t num_graph_nodes = graph_model.nodes().size();

    FP objective = 0.0;

    std::vector<std::vector<FP>> parameters(
        num_graph_nodes, std::vector<FP>(settings.num_control_parameters() * settings.num_control_intervals()));

    for (size_t control_interval = 0; control_interval < settings.num_control_intervals(); control_interval++) {
        for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
            size_t idx = control_index + control_interval * settings.num_control_parameters();
            for (size_t node_index = 0; node_index < num_graph_nodes; node_index++) {
                size_t global_idx =
                    idx + node_index * settings.num_control_intervals() * settings.num_control_parameters();
                parameters[node_index][idx] = settings.activation_function()(ptr_parameters[global_idx]);
            }
        }
    }

    set_control_dampings<FP>(settings, graph_model, parameters);

    for (auto& node : graph_model.nodes()) {
        node.property.get_simulation().set_integrator(make_integrator<FP>(settings.integrator_type(), settings.dt()));
    }

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());
    auto graph_sim_mobility    = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), graph_model);

    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {

        graph_sim_mobility.advance(time_steps[interval + 1]);

        size_t control_interval = interval / settings.pc_resolution();

        auto param_at = [&](const std::string& name, size_t nodex_index) {
            size_t control_index = static_cast<size_t>(string_to_control(name));
            return parameters[nodex_index][control_index + control_interval * settings.num_control_parameters()];
        };

        auto cost = [&](const std::string& name) {
            size_t control_index = static_cast<size_t>(string_to_control(name));
            return settings.control_parameters()[control_index].cost();
        };

        FP interval_cost = 0.0;
        for (size_t node_index = 0; node_index < num_graph_nodes; node_index++) {
            interval_cost += cost("SchoolClosure") * param_at("SchoolClosure", node_index);
            interval_cost += cost("HomeOffice") * param_at("HomeOffice", node_index);
            interval_cost += cost("PhysicalDistancingSchool") * param_at("PhysicalDistancingSchool", node_index);
            interval_cost += cost("PhysicalDistancingWork") * param_at("PhysicalDistancingWork", node_index);
            interval_cost += cost("PhysicalDistancingOther") * param_at("PhysicalDistancingOther", node_index);
        }

        objective += interval_cost / settings.num_intervals();
    }

    return objective;
}
