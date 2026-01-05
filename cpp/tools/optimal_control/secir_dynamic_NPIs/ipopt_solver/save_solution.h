
#pragma once

#include "../optimization_settings/optimization_settings.h"
#include "models/ode_secir/model.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"

#include "../constraints/infection_state_utils.h"

#include <vector>
#include <iostream>

template <typename FP>
void save_solution(mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>> graph_model,
                   const SecirOptimization& settings, size_t n, const FP* ptr_parameters, const FP* z_L, const FP* z_U,
                   size_t m, const FP* ptr_constraints, const FP* lambda, FP obj_value)
{
    std::vector<FP> dynamic_NPI_strengths(settings.num_control_parameters());
    for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
        dynamic_NPI_strengths[control_index] = ptr_parameters[control_index];
    }

    std::cout << "\nPath Constraints:\n";
    for (size_t i = 0; i < settings.num_path_constraints(); ++i) {
        std::cout << settings.path_constraints()[i].name() << ": " << ptr_constraints[i] << std::endl;
    }
    std::cout << "\nTerminal Constraints:\n";
    for (size_t i = 0; i < settings.num_terminal_constraints(); ++i) {
        std::cout << settings.terminal_constraints()[i].name() << ": "
                  << ptr_constraints[i + settings.num_path_constraints()] << std::endl;
    }

    // Print optimal control parameters
    std::cout << "\nOptimal Control Parameters:\n";
    size_t index = 0; // index into dynamic_NPI_strengths
    for (auto& cp : settings.control_parameters()) {
        std::cout << cp.name() << ", " << "Threshold: " << cp.threshold() << "\n";

        // Print each damping with its corresponding name
        const auto& damping_names = cp.damping_names();
        const auto& dampings      = cp.dampings();

        for (size_t d = 0; d < dampings.size(); ++d) {
            FP strength = std::max(dynamic_NPI_strengths[index++], static_cast<FP>(0));
            double tol  = 1e-7;
            if (strength < tol) {
                strength = 0.0;
            }
            std::cout << "    Damping: " << damping_names[d] << ", Strength: " << strength << "\n";
        }
    }

    // Set the control parameters for each node simulation
    for (auto& node : graph_model.nodes()) {
        set_dynamic_NPIs<FP>(settings, node.property.get_simulation().get_model(), dynamic_NPI_strengths);
    }

    // Set the integrator for each node simulation
    for (auto& node : graph_model.nodes()) {
        node.property.get_simulation().set_integrator_core(std::make_unique<mio::EulerIntegratorCore<FP>>());
    }

    // Create graph simulation (move graph_model since parameter was passed by value)
    auto graph_simulation = mio::make_mobility_sim<FP>(settings.t0(), 0.01, graph_model);

    graph_simulation.advance(settings.tmax());

    std::vector<mio::TimeSeries<double>> results = mio::interpolate_simulation_result(graph_simulation.get_graph());

    std::vector<int> county_ids(graph_simulation.get_graph().nodes().size());
    for (int i = 0; i < county_ids.size(); i++) {
        county_ids[i] = graph_simulation.get_graph().nodes()[i].id;
    }

    auto res = mio::save_result(results, county_ids, 1, "result_optimal_controls.h5");
}
